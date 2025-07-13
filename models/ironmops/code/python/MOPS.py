import numpy as np
from petsc4py import PETSc
import tmm4py as TMM
from numba import njit
from numba.typed import List

from INSOLATION import insolation
from CARCHEM import car_ini, co2_surfforcing
from OXYGEN import o2_surfforcing, o2saturation

from MOPS_PARAMS import *

class DotDict(dict):
    """dot.notation access to dictionary attributes"""

    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, d):
        return self.__dict__.update(d)

bgc_ktotal = 60

class MOPS:

    p = DotDict()

    def __init__(self, efctxId, prefix):
      self.efctxId = efctxId
      self.prefix = prefix
      self.ef = DotDict()
      
    def iniExternalForcingFn(self, tc, Iter, state, *args, **kwargs): 
      p=MOPS.p
      myId = PETSc.COMM_WORLD.getRank()
      OptDB=PETSc.Options()
      ef=self.ef

#     First set data common to all instances of MOPS
      if self.efctxId==1:
        p.myId = myId
        p.profileConfig=TMM.getProfileConfig()
        p.timeConfig=TMM.getTimeConfig()
        p.bgc_zu=np.zeros(bgc_ktotal+1) # geometry needed by all MOPS instances

      lSize=p.profileConfig['lSize']
      lNumProfiles=p.profileConfig['lNumProfiles']
      lProfileLength=p.profileConfig['lProfileLength']
      lStartIndices=p.profileConfig['lStartIndices']
      maxSteps=p.timeConfig['maxSteps']
      Iter0=p.timeConfig['Iter0']
      
      if self.efctxId==1:
        p.useSeparateBiogeochemTimeStepping=OptDB.hasName('separate_biogeochem_time_stepping')
        if p.useSeparateBiogeochemTimeStepping:
          PETSc.Sys.Print("Biogeochem model will be time-stepped independently")
          p.numBiogeochemStepsPerOceanStep=OptDB.getInt('num_biogeochem_steps_per_ocean_step')
          PETSc.Sys.Print("Number of biogeochem model time steps per ocean time step = %d" % p.numBiogeochemStepsPerOceanStep)

        p.useCarbon = OptDB.hasName('use_carbon')
        p.useIMPRO = OptDB.hasName('use_impro')
        
        try:
          p.nzeuph=OptDB.getInt('nzeuph')
        except:
          PETSc.Error("Must indicate number of euphotic zone layers with the -nzeuph option")
        PETSc.Sys.Print("Number of euphotic zone layers is %d" % p.nzeuph)

        try:
          p.DeltaT=OptDB.getReal('biogeochem_deltat')
        except:
          PETSc.Error("Must indicate biogeochemical time step in seconds with the -biogeochem_deltat option")
        PETSc.Sys.Print("Ocean time step for BGC length is  %12.7f seconds" % p.DeltaT)

        p.daysPerYear=OptDB.getReal('days_per_year',360.)
        PETSc.Sys.Print("Number of days per year is %12.7f" % p.daysPerYear)
        p.secondsPerYear = 86400.0*p.daysPerYear

#/* Need this for atmospheric exchange, river runoff, ... */
        p.localdA=TMM.readProfileScalarData('dA.bin')
        p.totalA=TMM.globalSum(p.localdA)

#/* Grid arrays */
        p.localdz=TMM.loadVecIntoArray("dz.petsc")
        f = open("drF.bin", 'rb')
        p.nzmax = np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0]
        drF = np.fromfile(f, dtype=np.dtype('>f8'), count=p.nzmax)
        f.close()
        PETSc.Sys.Print("Number of vertical layers is %d" % p.nzmax)

        p.periodicBiogeochemForcing=OptDB.hasName('periodic_biogeochem_forcing')

        if p.periodicBiogeochemForcing:
          PETSc.Sys.Print("Periodic biogeochemical forcing specified")
          p.biogeochemTimer=TMM.PeriodicTimer()
          p.biogeochemTimer.create(prefix="periodic_biogeochem_")
  
#/*   Read T and S */
        if p.periodicBiogeochemForcing:
          p.localTs=np.zeros(lSize)
          p.localSs=np.zeros(lSize)
          p.localFes=np.zeros(lSize)
          p.Tsp=TMM.PeriodicVec()
          p.Tsp.create(typ='tracer', buf=p.localTs)
          p.Ssp=TMM.PeriodicVec()
          p.Ssp.create(typ='tracer', buf=p.localSs)
          p.Fesp=TMM.PeriodicVec()
          p.Fesp.create(typ='tracer', buf=p.localFes)
        else:
          p.localTs=TMM.loadVecIntoArray("Ts.petsc")
          p.localSs=TMM.loadVecIntoArray("Ss.petsc")
          p.localFes=TMM.loadVecIntoArray("Fes.petsc")
          PETSc.Sys.Print("Done reading T/S/Fe")

        if p.useCarbon:
          p.useVirtualFlux=OptDB.hasName("use_virtual_flux")
          if p.periodicBiogeochemForcing:
            p.localEmP=np.zeros(lNumProfiles)
            p.localEmPp=TMM.PeriodicArray()
            p.localEmPp.create(lNumProfiles, buf=p.localEmP)
          else:
            p.localEmP=TMM.readProfileScalarData("EmP.bin")
          if p.useVirtualFlux:
            p.surfVolFrac=TMM.loadVecIntoArray("surface_volume_fraction.petsc")

# /* fraction of global river runoff in each box, divided by the box volume (a 3D field) */
        p.useRunoff = OptDB.hasName('use_runoff')
        p.writeRunoff = OptDB.hasName('write_runoff')
        if p.useRunoff:
          p.localrunoffvol=TMM.loadVecIntoArray("runoff_volume_annual.petsc")
          PETSc.Sys.Print("Runoff will be supplied via rivers")
        else:
          PETSc.Sys.Print("Runoff will be distributed over total ocean area of %g" % p.totalA)
          p.localrunoffvol=TMM.loadVecIntoArray("dz.petsc")
# 	/* IK: loading dz.petsc is just a dummy for now; runoff will be divided by dz(1) in BGC_MODEL.F */ 

# /* Forcing fields */
#         p.INS = INSOLATION()
        p.localswrad=np.zeros(lNumProfiles)
        p.localtau=np.zeros(lNumProfiles)
        p.locallatitude=TMM.readProfileScalarData("latitude.bin")
        if p.periodicBiogeochemForcing:
          p.localfice=np.zeros(lNumProfiles)
          p.localficep=TMM.PeriodicArray()
          p.localficep.create(lNumProfiles,buf=p.localfice)
        else:
          p.localfice=TMM.readProfileScalarData("fice.bin")

        if p.periodicBiogeochemForcing:
          p.localwind=np.zeros(lNumProfiles)
          p.localwindp=TMM.PeriodicArray()
          p.localwindp.create(lNumProfiles,buf=p.localwind)
        else:
          p.localwind=TMM.readProfileScalarData("wind.bin")

        if p.periodicBiogeochemForcing:
          p.localatmosp=np.zeros(lNumProfiles)
          p.localatmospp=TMM.PeriodicArray()
          p.localatmospp.create(lNumProfiles,buf=p.localatmosp)
        else:
          p.localatmosp=TMM.readProfileScalarData("atmosp.bin")

#SPK This is done in mops_biogeochem_ini.F
# /* Initialize biogeochem model */
# C Biogeochemical time step length in days
        p.bgc_timesteps = p.numBiogeochemStepsPerOceanStep
        p.bgc_dt=p.DeltaT/(86400.0*p.numBiogeochemStepsPerOceanStep)

# C Depth of layers for biogeochemistry

        p.bgc_kmax = p.nzmax
        p.bgc_keuph = p.nzeuph

        p.bgc_zu[1:p.bgc_kmax+1]=np.cumsum(drF)
#         for k in range(1,p.bgc_kmax+1):
#           bgc_zu[k]=bgc_zu[k-1]+drF[k-1]
# do k=2,bgc_kmax+1
#   bgc_zu(k)=bgc_zu(k-1)+drFloc(k-1)
# enddo        

#   } /* end common data */

#     Now set data for this MOPS instance
      ef.config=state.getConfig()
      numTracers=ef.config['numTracers']

      for itr in range(numTracers):
        state.qef[itr][:]=0.0

      OptDBState=PETSc.Options(self.prefix)
# // Some defaults
      ef.numDiag=12
      ef.appendDiagnostics = False

      if p.useCarbon:
# // Defaults
        ef.fixedAtmosCO2 = True
        ef.Focean=0.0
        ef.Foceanint = 0.0

        ef.useAtmModel = OptDBState.hasName('use_atm_model')

        if (ef.useAtmModel):
          PETSc.Sys.Print("Using interactive atmospheric model")
# /* overwrite default value */
          if OptDBState.hasName('pco2atm_ini'):
            pCO2atm_ini = OptDBState.getReal('pco2atm_ini', 280.0)
          else:
             if OptDBState.hasName("pco2atm_ini_file"):
               pCO2atmIniFile = OptDBState.getString("pco2atm_ini_file")
               pCO2atm_ini = np.fromfile(pCO2atmIniFile, dtype=np.dtype('>r8'), count=1)
          ef.pCO2atm = pCO2atm_ini
          PETSc.Sys.Print("Using initial atmospheric pCO2 of %g ppm" % ef.pCO2atm)
            
          ef.atmWriteTimer = TMM.StepTimer()
          ef.atmWriteTimer.create(prefix="atm_write_",startTimeStep=Iter0)

          ef.atmAppendOutput = OptDBState.hasName('atm_append')
          if ef.atmAppendOutput:
            PETSc.Sys.Print("Atmospheric model output will be appended")
          else:
            PETSc.Sys.Print("Atmospheric model output will overwrite existing file(s)")

         # /* Output times */
          ef.atmOutTimeFile = OptDBState.getString('atm_time_file',"atm_output_time.txt")
          PETSc.Sys.Print("Atmospheric model output times will be written to %s" % atmOutTimeFile)

          if not ef.atmAppendOutput:
            if p.myId==0:
              ef.atmfptime = open(atmOutTimeFile, 'wt')
              ef.atmfptime.write("{:d}   {:.5f}\n".format(Iter0,time0))
            PETSc.Sys.Print("Writing atmospheric output at time %10.5f, step %d" % (tc,Iter))
            TMM.writeBinaryArray("pCO2atm_output.bin",np.atleast_1d(ef.pCO2atm),False)
            TMM.writeBinaryArray("Foceanint_output.bin",np.atleast_1d(ef.Focean),False)
          else:
            if p.myId==0:
              ef.atmfptime = open(atmOutTimeFile, 'at')
            PETSc.Sys.Print("Atmospheric model output will be appended. Initial condition will NOT be written\n")

          ef.atmModelDeltaT = p.DeltaT/p.secondsPerYear #/* time step in years */

        else:  #/* not using atm model */
          PETSc.Sys.Print("Using prescribed atmospheric pCO2")
#            /* prescribed atmospheric CO2 */
          if OptDBState.hasName('pco2atm_history'):
            pCO2atmFiles = OptDBState.getString("pco2atm_history").split(",")
            if len(pCO2atmFiles) != 2:
              PETSc.Error("Insufficient number of file names specified for atmospheric pCO2 history")
            ef.fixedAtmosCO2 = False
            PETSc.Sys.Print("Reading time-dependent atmospheric pCO2 history")
#                /* read time data */
            f = open(pCO2atmFiles[0], 'rb')
            ef.numpCO2atm_hist = np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0]
            PETSc.Sys.Print("Number of points in atmospheric history file is %d" % ef.numpCO2atm_hist)
            ef.TpCO2atm_hist = np.fromfile(f, dtype=np.dtype('>f8'), count=ef.numpCO2atm_hist)
            f.close()
#                /* read atmospheric pCO2 data */
            ef.pCO2atm_hist = np.fromfile(pCO2atmFiles[1], dtype=np.dtype('>f8'), count=ef.numpCO2atm_hist)
            ef.pCO2atm = ef.pCO2atm_hist[0];
          else:
            ef.pCO2atm = OptDBState.getReal('pco2atm', 280.0)
            PETSc.Sys.Print("Using fixed atmospheric pCO2 of %g ppm" % ef.pCO2atm)
  
        ef.localph = np.zeros(lNumProfiles)
        ef.localco2airseaflux=np.zeros(lNumProfiles)

        # finished useCarbon

# // Defaults
      ef.localFburial = 0.0;
      ef.Fburial=0.0;
      PETSc.Sys.Print("Using Burial-Runoff model")

# /* Define the interval over which to integrate global burial */
      try:
        ef.burialSumSteps = OptDBState.getInt("burial_sum_steps")
      except:
        PETSc.Error("Must indicate burial integration interval with the -burial_sum_steps option")
        
      if ((maxSteps % ef.burialSumSteps)!=0):
        PETSc.Error("maxSteps not divisible by burialSumSteps!")

      PETSc.Sys.Print("Runoff will be integrated over every %d time steps" % ef.burialSumSteps)

      if p.writeRunoff:
# /* set the name of the runoff time file */
        runoffOutTimeFile = OptDBState.getString('runoff_time_file',"runoff_output_time.txt")
        PETSc.Sys.Print("Runoff output times will be written to %s\n" % runoffOutTimeFile)

# /* set inititial runoff: overwrite default value with value from command line*/
      runoff_ini = 0.0
      try:
        runoff_ini = OptDBState.getReal('runoff_ini', 0.0)
# /* set inititial runoff: overwrite default value with value from file*/
      except:
        if OptDBState.hasName("runoff_ini_file"):
          runoffIniFile = OptDBState.getString("runoff_ini_file")
          runoff_ini = np.fromfile(runoffIniFile, dtype=np.dtype('>f8'), count=1)

      ef.GRunoff = runoff_ini
      PETSc.Sys.Print("Using initial runoff of %g Gmol P/d" % ef.GRunoff)

      ef.runoffIniOutFile =  OptDBState.getString('pickup_runoff_out',"pickup_runoff.bin")
      PETSc.Sys.Print("Final runoff output will be written to %s" % ef.runoffIniOutFile)

      if p.writeRunoff:
# /* if run is continued, always append runoff and output times */
        if (Iter0>0):
          PETSc.Sys.Print("Runoff output will be appended")
          if p.myId==0:
            ef.runofffptime = open(runoffOutTimeFile, 'at')
          PETSc.Sys.Print("Initial runoff output will not be written\n");CHKERRQ(ierr);
        else:  
          PETSc.Sys.Print("Runoff output will overwrite existing file(s)\n");CHKERRQ(ierr);
          if p.myId==0:
            ef.runofffptime = open(runoffOutTimeFile, 'wt')
            ef.runofffptime.write("{:d}   {:.5f}\n".format(Iter0,time0))
          TMM.writeBinaryArray("Grunoff_output.bin",np.atleast_1d(ef.GRunoff),False)
          PETSc.Sys.Print("Writing runoff output at time %10.5f, step %d\n", tc,Iter);CHKERRQ(ierr);  

# /* Initialize biogeochem model */
# /* Read and overwrite default parameter values here */
#SPKSPK Need a way to transfer bgcparams to parameter variables that can be used during initialization
  
      ef.readBGCParams = OptDBState.hasName('bgc_params_file')
      if (ef.readBGCParams):
        bgcParamsFile = OptDBState.getString('bgc_params_file')
        try:
          numBGCParams = OptDBState.getInt('num_bgc_params')
        except:
          PETSc.Error("Must indicate number of BGC parameters to read with the -num_bgc_params option")
          
        PETSc.Sys.Print("Reading %d parameters from file\n" % numBGCParams)

        if OptDBState.hasName('ascii_params'):
          f = open(bgcParamsFile, 'rt')
          ef.bgcparams = np.fromfile(f, count=numBGCParams)
          f.close()
        else:
          ef.bgcparams = np.fromfile(bgcParamsFile, dtype=np.dtype('>f8'), count=numBGCParams)

#SPK These are currently global but need to be moved into ef after parameters have been optionally read
# C REMINERALISATION LENGTH SCALES
      ef.wdet=np.zeros(p.bgc_kmax)
      for k in range(p.bgc_kmax):
        ef.wdet[k] = detwb + (p.bgc_zu[k]+p.bgc_zu[k+1])*0.5*detwa

      if p.useIMPRO:
# C Use implicit profiles for particle sinking to overcome numerical diffusion,
# C as explained in Kriest and Oschlies, 2011, Ocean Model., 39, 275-283.
# C Only do this for deeper layers, as many other processes (physical and biological)
# C beside sinking and remineralization might play a role in the euphotic zone.
        anafac0 = (1.0+detwa*p.bgc_dt)**(detlambda/detwa)-1.0
        for k in range(p.bgc_keuph,p.bgc_kmax):
          anafacz = (p.bgc_zu(k)/p.bgc_zu[k+1])**(detlambda/detwa)
          ef.wdet[k] = ((p.bgc_zu[k+1]-p.bgc_zu[k])*anafac0* \
            anafacz/(1.0-anafacz))/p.bgc_dt
#       enddo

      if p.useCarbon:
        ef.fcaco3=np.zeros(p.bgc_kmax)

        for k in range(p.bgc_kmax):
          ef.fcaco3[k] =np.exp(0.0-p.bgc_zu[k]/length_caco3) # from OCMIP2: this is the flux fraction through the top of each box

      myTime = p.DeltaT*Iter #/* Iter should start at 0 */

#     Common data only needs to be updated once
      if (self.efctxId==1):
        if p.periodicBiogeochemForcing:
          p.Tsp.interp(tc, p.biogeochemTimer, "Ts_")
          p.Ssp.interp(tc, p.biogeochemTimer, "Ss_")
          p.Fesp.interp(tc, p.biogeochemTimer, "Fes_")
          insolation(myTime,p.locallatitude,p.daysPerYear,p.localswrad,p.localtau)
          p.localficep.interp(tc, p.biogeochemTimer, "fice_")
          p.localwindp.interp(tc, p.biogeochemTimer, "wind_")
          p.localatmospp.interp(tc, p.biogeochemTimer, "atmosp_")
          if p.useCarbon:
            p.localEmPp.interp(tc, p.biogeochemTimer, "EmP_")
        else:
          insolation(myTime,p.locallatitude,p.daysPerYear,p.localswrad,p.localtau)
# finished common initialization

      for ip in range(lNumProfiles):
        nzloc=lProfileLength[ip]
        kl=lStartIndices[ip]

        if p.useCarbon:

# ! Global mean silicate for surface layer
          ssil=ocmip_silfac

# ! Surface total alkalinity follows the OCMIP protocol:
# ! ta = 2310*s/sbar where sbar is the model's annual mean surface salinity      
          salk = ocmip_alkfac*p.localSs[kl]
              
# ! Initialize carbonate system.
# ! The constants and initial coefficients will be written into a common block.
          sdic = state.c[idic][kl]
          spho = state.c[ipo4][kl]

          ef.localph[ip]=car_ini(p.localTs[kl],p.localSs[kl],sdic,spho,salk,ssil)

      ef.calcDiagnostics = OptDBState.hasName('calc_diagnostics')
      if ef.calcDiagnostics:
#     /*Data for diagnostics */
        ef.diagTimer = TMM.StepTimer()
        ef.diagTimer.create(prefix="diag_",startTimeStep=Iter0)
        PETSc.Sys.Print("Diagnostics will be computed starting at (and including) time step: %d" % ef.diagTimer.startTimeStep)
        PETSc.Sys.Print("Diagnostics will be computed over %d time steps" % ef.diagTimer.numTimeSteps)

        try:
          ef.diagOutFile = OptDBState.getString("diag_files").split(",")
        except:
          PETSc.Error("Must indicate file name(s) for writing diagnostics with the -diag_files option")
        if len(ef.diagOutFile) != ef.numDiag:
          PETSc.Error("Insufficient number of time average file names specified")

        ef.diagMode = "w"
        ef.localfbgc=[np.zeros(lSize) for _ in range(ef.numDiag)]
        ef.localfbgcavg=[np.zeros(lSize) for _ in range(ef.numDiag)]
        
        if p.useCarbon:
          ef.localco2airseafluxdiagavg=np.zeros(lNumProfiles)

          try:
            ef.co2airseafluxFile = OptDBState.getString('co2airseaflux_file')
          except:
            PETSc.Error("Must indicate file name for writing co2airsea flux with the -co2airseaflux_file option")
      
    def calcExternalForcingFn(self, tc, Iter, iLoop, state, *args, **kwargs):

      p=MOPS.p
      ef=self.ef

      lSize=p.profileConfig['lSize']
      lNumProfiles=p.profileConfig['lNumProfiles']
      lProfileLength=p.profileConfig['lProfileLength']
      lStartIndices=p.profileConfig['lStartIndices']

      Iter0=p.timeConfig['Iter0']

      numTracers=ef.config['numTracers']
      
#       if p.useCarbon:
#         localco2airseaflux = 0.0
#       else:
# #       Dummy values to pass to BGC_MODEL for when not using carbon      
# #         sph=None
# #         emp=None
#         pco2atm=None
#         dicgave=None
#         alkgave=None
      
      localFocean=0.0
      localburial = 0.0
             
      myTime = p.DeltaT*Iter # /* Iter should start at 0 */

#     Common data only needs to be updated once
      if self.efctxId==1:
        if p.periodicBiogeochemForcing:
          p.Tsp.interp(tc, p.biogeochemTimer, "Ts_")
          p.Ssp.interp(tc, p.biogeochemTimer, "Ss_")
          p.Fesp.interp(tc, p.biogeochemTimer, "Fes_")
          insolation(myTime,p.locallatitude,p.daysPerYear,p.localswrad,p.localtau)
          p.localficep.interp(tc, p.biogeochemTimer, "fice_")
          p.localwindp.interp(tc, p.biogeochemTimer, "wind_")
          p.localatmospp.interp(tc, p.biogeochemTimer, "atmosp_")
          if p.useCarbon:
            p.localEmPp.interp(tc, p.biogeochemTimer, "EmP_")
        else:
          insolation(myTime,p.locallatitude,p.daysPerYear,p.localswrad,p.localtau)

      if p.useCarbon:
        if not ef.useAtmModel:
# /* Interpolate atmospheric pCO2   */
          if not ef.fixedAtmosCO2:
            if (tc>=ef.TpCO2atm_hist[0]):
              ef.pCO2atm = np.interp(tc, ef.TpCO2atm_hist, ef.pCO2atm_hist, left=ef.pCO2atm_hist[0])
            elif (tc>ef.TpCO2atm_hist[-1]):
              PETSc.Error("Error: time out of bounds")
            else:  
              PETSc.Sys.Print("Warning: time < %10.5f. Assuming pCO2atm=%g" % (ef.TpCO2atm_hist[0],ef.pCO2atm))

        if (p.useVirtualFlux): #{ /* use the global surface mean value to calculate E-P contribution */
          DICemp = TMM.dotProd(p.surfVolFrac,state.c[idic]) # /* volume weighted mean surface DIC */									              
          ALKemp = TMM.dotProd(p.surfVolFrac,state.c[ialk]) # /* volume weighted mean surface ALK */									                  
      
      wdet = ef.wdet
      
      fcaco3 = None
      if p.useCarbon:
        fcaco3 = ef.fcaco3

      fbgc = None
      doDiagnostics = False
      if ef.calcDiagnostics:
        if (Iter0+iLoop>=ef.diagTimer.startTimeStep): #{ /* start time averaging (note: startTimeStep is ABSOLUTE time step) */
          doDiagnostics=True
      if doDiagnostics:
        fbgc=List(ef.localfbgc)

      # Pack physics and geometry into lists

      bgc_geomi=List([lProfileLength,lStartIndices])
      bgc_forcing2d=List([p.localswrad,p.localtau,p.localfice,p.localwind,p.localatmosp])
      bgc_forcing3d=List([p.localdz,p.localTs,p.localSs,p.localFes,p.localrunoffvol])
      tracer = List(state.c)
      bgc_globalrunoff = ef.GRunoff
      bgc_flags=List([p.useCarbon,p.useRunoff,p.useVirtualFlux])

      pco2atm=None
      dicgave=None
      alkgave=None
      if p.useCarbon:
        bgc_carb=List([ef.localph,p.localEmP,ef.localco2airseaflux])
        pco2atm=ef.pCO2atm
        if p.useVirtualFlux:
          dicgave=DICemp
          alkgave=ALKemp
      else:
        bgc_carb=None

      localtotFburial = BGC_MODEL(bgc_flags,bgc_globalrunoff,p.localdA,p.bgc_dt,p.bgc_timesteps,tracer,bgc_geomi,p.bgc_keuph,wdet,fcaco3, \
                                  bgc_forcing2d,bgc_forcing3d,bgc_carb,pco2atm,dicgave,alkgave,fbgc)

# # /* burial in sediment is already integrated within BGC_MODEL for this processor */
      # accumulate
      ef.localFburial = ef.localFburial  + localtotFburial
       
      if p.useCarbon:
        if (ef.useAtmModel):
          ef.Focean=TMM.dotProd(localco2airseaflux,p.localdA)*(1.0/p.DeltaT)*(12.0/1.e18)*secondsPerYear #/* PgC/y */
          # time step atmosphere
          ef.pCO2atm = ef.pCO2atm + ef.atmModelDeltaT*(-ef.Focean)/ppmToPgC
          # reset values
          ef.Foceanint = ef.Foceanint + ef.atmModelDeltaT*ef.Focean #/* calculate the time integrated flux */

# /* sum burial in sediment over all processors, and scale by time step etc.*/
# /* do this only once every burialSumSteps , and then take this value for next year's runoff */

      if ((iLoop % ef.burialSumSteps)==0):
        ef.Fburial=TMM.globalScalarSum(ef.localFburial)

        if p.useRunoff:
          ef.GRunoff = ef.Fburial/(1.e12*ef.burialSumSteps)*(86400.0/p.DeltaT) # /* This is Gmol P/day. 
#         Note: localrunoff is scaled with 1e12. Note: GRunoff will be scaled with bgc_dt.*/
        else:
          ef.GRunoff = ef.Fburial/(p.totalA*ef.burialSumSteps)*(86400.0/p.DeltaT) # /* This is mmol P/m2/day. 
#         Note: this will later be divided by depth of first layer. Note: GRunoff will be scaled with bgc_dt.*/ 

        ef.localFburial = 0.0

    def writeExternalForcingFn(self, tc, Iter, iLoop, state, *args, **kwargs):

      p=MOPS.p
      ef=self.ef

      lSize=p.profileConfig['lSize']
      lNumProfiles=p.profileConfig['lNumProfiles']
      lProfileLength=p.profileConfig['lProfileLength']
      lStartIndices=p.profileConfig['lStartIndices']

      Iter0=p.timeConfig['Iter0']

      numTracers=ef.config['numTracers']
      

      if p.useCarbon:
        if (ef.useAtmModel):
#       /* write instantaneous atmos model state */
          if (Iter0+iLoop>=(ef.atmWriteTimer.startTimeStep)): # { /* note: startTimeStep is ABSOLUTE time step */
            if ((ef.atmWriteTimer.count)<=(ef.atmWriteTimer.numTimeSteps)):
              ef.atmWriteTimer.incr()
            if ((ef.atmWriteTimer.count)==(ef.atmWriteTimer.numTimeSteps)): # { /* time to write out */
              PETSc.Sys.Print("Writing atmospheric model output at time %10.5f, step %d" % (tc, Iter0+iLoop))
              if p.myId==0:
                ef.atmfptime.write("{:d}   {:.5f}\n".format(Iter0+iLoop,tc))
              TMM.writeBinaryArray("pCO2atm_output.bin",np.atleast_1d(ef.pCO2atm),True)
              TMM.writeBinaryArray("Foceanint_output.bin",np.atleast_1d(ef.Focean),True)
              ef.Foceanint = 0.0
              ef.atmWriteTimer.update(Iter0+iLoop)

      if p.writeRunoff:
        if ((iLoop % ef.burialSumSteps)==0): #/*  time to write out */
          PETSc.Sys.Print("Writing runoff output at time %10.5f, step %d" % (tc, Iter0+iLoop))
          if p.myId==0:
            ef.runofffptime.write("{:d}   {:.5f}\n".format(Iter0+iLoop,tc))
          TMM.writeBinaryArray("Grunoff_output.bin",np.atleast_1d(ef.GRunoff),True)

      if (ef.calcDiagnostics):
        if (Iter0+iLoop>=(ef.diagTimer.startTimeStep)):  #/* start time averaging (note: startTimeStep is ABSOLUTE time step) */  
  
          if (ef.diagTimer.count<=ef.diagTimer.numTimeSteps): #/* still within same averaging block so accumulate */
            for i in range(ef.numDiag):
              ef.localfbgcavg[i][:] = ef.localfbgcavg[i][:] + ef.localfbgc[i][:]

            if p.useCarbon:
              ef.localco2airseafluxdiagavg[:]=ef.localco2airseaflux[:]+ef.localco2airseafluxdiagavg[:]

            ef.diagTimer.incr()
          
          if ((ef.diagTimer.count)==(ef.diagTimer.numTimeSteps)): #/* time to write averages to file */
            PETSc.Sys.Print("Writing diagnostics time average at time %10.5f, step %d" % (tc, Iter0+iLoop))
            for i in range(ef.numDiag):
              ef.localfbgcavg[i][:] = ef.localfbgcavg[i][:]/ef.diagTimer.count
              TMM.writeArrayToVec(ef.diagOutFile[i], ef.localfbgcavg[i], ef.diagMode)
              ef.localfbgcavg[i][:] = 0.0
            ef.diagMode = "a"
            if p.useCarbon:
               ef.localco2airseafluxdiagavg[:]=ef.localco2airseafluxdiagavg[:]/ef.diagTimer.count
               TMM.writeProfileScalarData(ef.co2airseafluxFile, ef.localco2airseafluxdiagavg, 1, ef.appendDiagnostics)
#       /*      reset diagnostic arrays */
               ef.localco2airseafluxdiagavg[:]=0.0

            ef.appendDiagnostics=True
            ef.diagTimer.update(Iter0+iLoop)
      

    def finalizeExternalForcingFn(self, tc, Iter, state, *args, **kwargs):

      p=MOPS.p
      ef=self.ef

      lSize=p.profileConfig['lSize']
      lNumProfiles=p.profileConfig['lNumProfiles']
      lProfileLength=p.profileConfig['lProfileLength']
      lStartIndices=p.profileConfig['lStartIndices']

      Iter0=p.timeConfig['Iter0']

      numTracers=ef.config['numTracers']

#       /* write final pickup */
      if p.useCarbon:
        if (ef.useAtmModel):
#       /* write instantaneous atmos model state */
          TMM.writeBinaryArray("pickup_pCO2atm.bin",np.atleast_1d(ef.pCO2atm),False)
      
      TMM.writeBinaryArray(ef.runoffIniOutFile,np.atleast_1d(ef.GRunoff),False)

#     Delete common data
      if (ef.efctxId==1):
        if (p.periodicBiogeochemForcing):
          p.Tsp.destroy()
          p.Ssp.destroy()
          p.Fesp.destroy()
          p.localficep.destroy()
          p.localwindp.destroy()
          p.localatmospp.destroy()
          if p.useCarbon:
           p.localEmPp.destroy()

      if (ef.useAtmModel):
        ef.atmfptime.close()

      if p.writeRunoff:
        ef.runofffptime.close()

    def reInitializeExternalForcingFn(self, tc, Iter, iLoop, state, *args, **kwargs):
        pass

@njit(cache=False)
def BGC_MODEL(bgc_flags,bgc_globalrunoff,localdA,bgc_dt,bgc_timesteps,tracer,bgc_geomi,bgc_keuph,wdet,fcaco3, \
              bgc_forcing2d,bgc_forcing3d,bgc_carb,pco2atm,dicgave,alkgave,fbgc):
                                                  
  useCarbon,useRunoff,useVirtualFlux=bgc_flags
  lProfileLength,lStartIndices=bgc_geomi
  localswrad,localtau,localfice,localwind,localatmosp=bgc_forcing2d
  localdz,localTs,localSs,localFes,localrunoffvol=bgc_forcing3d
  
  numTracers=len(tracer)
  lNumProfiles=len(lProfileLength)

  if useCarbon:
    localph,localEmP,localco2airseaflux=bgc_carb

  doDiagnostics=(fbgc is not None)
  if doDiagnostics:
    numDiag=len(fbgc)

  localflux_bury=0.0
 
 #  bgc_globalrunoff = GRunoff

  for ip in range(lNumProfiles):
    nzloc=lProfileLength[ip]
    kl=lStartIndices[ip]

    bgc_tracer = [tracer[itr][kl:kl+nzloc] for itr in range(numTracers)]

    if doDiagnostics:
      localfbgc=[fbgc[i][kl:kl+nzloc] for i in range(numDiag)]

    if useCarbon:
      sph=localph[ip]
      emp=localEmP[ip]
      if not useVirtualFlux: # /* use the local surface value to calculate E-P contribution */
        dicgave=bgc_tracer[idic][0]
        alkgave=bgc_tracer[ialk][0]
 
    bgc_kloc=nzloc
    bgc_swr=localswrad[ip]
    bgc_tau=localtau[ip]
    bgc_seaice=localfice[ip]
    bgc_wind=localwind[ip]
    bgc_atmosp=localatmosp[ip]

    bgc_dz=localdz[kl:kl+nzloc]
    bgc_theta=localTs[kl:kl+nzloc]
    bgc_salt=localSs[kl:kl+nzloc]
    bgc_fedep=localFes[kl:kl+nzloc]
    bgc_runoffvol=localrunoffvol[kl:kl+nzloc]
 
  
  # ! Depth of the euphotic zone
    bgc_keuphloc = np.fmin(bgc_kloc,bgc_keuph)

  # ! Reset the diagnostic fluxes per ocean time step.
    if doDiagnostics:
      f1_out=localfbgc[0]
      f2_out=localfbgc[1]
      f3_out=localfbgc[2]
      f4_out=localfbgc[3]
      f5_out=localfbgc[4]
      f6_out=localfbgc[5]
      f7_out=localfbgc[6]
      f8_out=localfbgc[7]
      f9_out=localfbgc[8]
      f10_out=localfbgc[9]
      f11_out=localfbgc[10]
      f12_out=localfbgc[11]

    if doDiagnostics:
      f1_out[:]=0.0 # PP
      f2_out[:]=0.0 # Grazing
      f3_out[:]=0.0 # Sedimentation and burial
      f4_out[:]=0.0 # Oxic remineralisation of POP and DOP
      f5_out[:]=0.0 # River runoff or surface deposition of buried material
      f6_out[:]=0.0 # Nitrogen fixation
      f7_out[:]=0.0 # Anoxic remineralisation of POP and DOP
      f8_out[:]=0.0 # Production of calcite (only for option DCARBON)
      f9_out[:]=0.0 # Dissolution of calcite  (only for option DCARBON)
      f10_out[:]=0.0 # Precipitation and organic scavenging of dFe  (only for option DCARBON)
      f11_out[:]=0.0 # Sedimentation of pFe (layer 1: burial)
      f12_out[:]=0.0 # Sedimentary release of dFe
  
    # ! arrays to store the sms and other fluxes
  #   topo4=localto[ipo4]
  #   tophy=localto[iphy]
  #   tozoo=localto[izoo]
  #   todop=localto[idop]
  #   todet=localto[idet]
  #   tooxy=localto[ioxy]
  #   todin=localto[idin]
  #   # ! For Fe: arrays to store sms
  #   todfe=localto[idfe]
  #   topfe=localto[ipfe]
  #   if useCarbon:
  #     todic=localto[idic]
  #     toalk=localto[ialk]
    topo4=np.zeros(bgc_kloc)
    tophy=np.zeros(bgc_kloc)
    tozoo=np.zeros(bgc_kloc)
    todop=np.zeros(bgc_kloc)
    todet=np.zeros(bgc_kloc)
    tooxy=np.zeros(bgc_kloc)
    todin=np.zeros(bgc_kloc)
    # ! For Fe: arrays to store sms
    todfe=np.zeros(bgc_kloc)
    topfe=np.zeros(bgc_kloc)
    if useCarbon:
      todic=np.zeros(bgc_kloc)
      toalk=np.zeros(bgc_kloc)

    flux=np.zeros(bgc_kloc)
    fdiv=np.zeros(bgc_kloc)
    feflux=np.zeros(bgc_kloc)
    fefdiv=np.zeros(bgc_kloc)
    runoff=np.zeros(bgc_kloc)
    ciz=np.zeros(bgc_kloc)
  
  # ! Things for air-sea gas exchange that are common to all tracers 
  # ! (O2,CO2)
  # ! vgas660 = exchange coefficient normalized to a Sc of 660,
  # ! Conversion: m/day (0.24=[.01 m/cm]*[24 hr/day] converts cm/hr to m/d)

    vgas660=(0.337*bgc_wind**2)*0.24*(1.0-bgc_seaice)

    if useCarbon:
  #     car_coeffs(bgc_theta[0],bgc_salt[0])
      co2airseaflux=0.0
  #   else:
  # #       Dummy return values for when not using carbon      
  #     sph=None
  #     co2airseaflux=None
    
    flux_bury = 0.0
  
  # C INTERNAL TIME LOOP FOR BGC

    for it in range(bgc_timesteps):
  # ! Reset fluxes.
      topo4[0:bgc_kloc]=0.0
      todop[0:bgc_kloc]=0.0
      tooxy[0:bgc_kloc]=0.0
      tophy[0:bgc_kloc]=0.0
      tozoo[0:bgc_kloc]=0.0
      todet[0:bgc_kloc]=0.0
      flux[0:bgc_kloc] =0.0
      fdiv[0:bgc_kloc] =0.0
      todin[0:bgc_kloc]=0.0

      if useCarbon:
        todic[0:bgc_kloc]=0.0
        toalk[0:bgc_kloc]=0.0

      runoff[0:bgc_kloc]=0.0 

  #           ! For Fe       
      todfe[0:bgc_kloc] = 0.0
      topfe[0:bgc_kloc] = 0.0
      feflux[0:bgc_kloc] =0.0
      fefdiv[0:bgc_kloc] =0.0

  # ! AIR-SEA GAS EXCHANGE OF OXYGEN

      o2sat=o2saturation(bgc_theta[0],bgc_salt[0])    
      o2gasex=o2_surfforcing(vgas660,bgc_atmosp,bgc_theta[0],bgc_salt[0],bgc_tracer[ioxy][0],o2sat)

      if useCarbon:
  # ! AIR-SEA GAS EXCHANGE OF CO2

        surf_dic = bgc_tracer[idic][0]
        surf_pho = bgc_tracer[ipo4][0]

  # ! Surface total alkalinity 
        surf_alk = bgc_tracer[ialk][0]

  # ! Surface silicate from the OCMIP protocol
        surf_sil=ocmip_silfac

        sph, co2gasex, co2emp, alkemp = co2_surfforcing(bgc_theta[0],bgc_salt[0],vgas660,bgc_atmosp,surf_dic,surf_pho,surf_alk,surf_sil,pco2atm,emp,dicgave,alkgave,sphin=sph)

  # !     co2gasex and co2emp are in mmolC/(m^2 d)
  # !     co2airseaflux at the end of the internal time stepping loop will be 
  # !     in mmolC/(m^2 ocean_time_step)
        co2airseaflux = co2airseaflux + (co2gasex + co2emp)*bgc_dt

  # !      CALL CAR_CHEMISTRY(...,dicchem,alkchem)

  # ! EUPHOTIC ZONE  AND ITS EXPORT TO DEEPER LAYERS 
  # 
  # ! Net solar radiation at top of every layer.

      parfrac=0.4
      ciz[0]=bgc_swr*(1.0-bgc_seaice)*parfrac

      for k in range(1,bgc_keuphloc):
        attlim=np.fmax(0.0,bgc_tracer[iphy][k-1])
        atten = (ACkw+ACkchl*attlim)*bgc_dz[k-1]
        ciz[k]=ciz[k-1]*np.exp(-atten)

  # ! Biogeochemical fluxes in euphotic zone and their export
  # ! Take care of negative tracer concentrations.

      for k in range(bgc_keuphloc):
        PO4=bgc_tracer[ipo4][k]
        DIN=bgc_tracer[idin][k]
        PHY=bgc_tracer[iphy][k]
        ZOO=bgc_tracer[izoo][k]
        DFE=bgc_tracer[idfe][k]
        attlim=np.fmax(PHY,0.0)

  # ! temperature dependence of phytoplankton growth (Eppley)
  # ! this affects the light-half-saturation constant via acik=acmuphy/alpha
  # ! the extension to dfe means that muphy and ik decline with declining fe, leaving alpha untouched 
        tempscale = np.exp(bgc_theta[k]/TempB)
        felim  = DFE/(kfe+DFE)
        TACmuphy = ACmuphy*tempscale*felim
        TACik = ACik*tempscale*felim       
 
  # ! The light limitation function of phytoplankton.
  # ! This function corresponds to Evans and Garcon, 1997.
  # ! Note that the initial slope of the P-I curve, alpha, is ACMuPhy/ACIk
  # ! flightlim thus gives the light limited growth rate, averaged over day 
  # ! and layer, normalised by max. growth rate
        atten = (ACkw+ACkchl*attlim)*bgc_dz[k] #attenuation at bottom of layer
        glbygd = 2.0*ciz[k]/(TACik*bgc_tau)   # 2 * G_L/G_D of EG97
        flightlim = bgc_tau/atten*(__phi(glbygd)-__phi(glbygd*np.exp(-atten)))

        if (PHY > 0.0):

          limnut = np.fmin(PO4,DIN/rnp) # because I use the same (scaled) half-sat constant for N and P I here determine the limiting nutrient in this simple way

          if (limnut > vsafe):

  #      ! The nutrient limitation of phytoplankton
            fnutlim = limnut/(ACkpo4+limnut)

  #      ! The growth rate of phytoplankton: minimum of light and nutrients; note that felim is 
  #      ! already in TACMuPhy! Thus, here we assume multiplicative limitation

            phygrow0 = TACmuphy*PHY*np.fmin(flightlim,fnutlim)

  #      ! Make sure not to take up more nutrients than available (mmol P/m3/d)
            phygrowFe = np.fmin(DFE,phygrow0*rfep*bgc_dt)/(bgc_dt*rfep) # mmol P/m3/d
            phygrowP  = np.fmin(limnut,phygrow0*bgc_dt)/bgc_dt          # mmol P/m3/d
            phygrow   = np.fmin(phygrowFe,phygrowP)                     # mmol P/m3/d

          else: #limnut < vsafe

            phygrow=0.0

  #               endif !limnut

  #      ! The exudation of phytoplankton
          phyexu = AClambda * PHY

  #      ! Other losses of phytoplankton       
          phyloss = AComni * PHY * PHY

          if (ZOO > 0.0):

  #      ! Grazing of zooplankton, Holling III
            graz0=ACmuzoo*PHY*PHY/(ACkphy*ACkphy+PHY*PHY)*ZOO

  #      ! Make sure not to graze more phytoplankton than available.
            graz = np.fmin(PHY,graz0*bgc_dt)/bgc_dt

          else: #ZOO < 0

            graz=0.0

  #               endif !ZOO

        else: #PHY < 0

          phygrow=0.0
          phyexu =0.0
          phyloss=0.0
          graz   =0.0

  #             endif !PHY

        if (ZOO > 0.0):

  #      ! Zooplankton exudation
           zooexu = AClambdaz * ZOO

  #      ! Zooplankton mortality 
           zooloss = AComniz * ZOO * ZOO
 
        else: #ZOO < 0

            zooexu = 0.0
            zooloss = 0.0

  #             endif !ZOO

  # ! Relaxation of N:P to Redfield values (mimick cyanobacteria)

        if(PO4 > vsafe):

          ttemp = bgc_theta[k]
          nfixtfac = np.fmax(0.0,tf2*ttemp*ttemp + tf1*ttemp + tf0)/tff
          dinlim = np.fmax(0.0,DIN)
          nfixnfac = np.fmax(0.0, 1.0-dinlim/(PO4*rnp))
          nfixation = nfixtfac*nfixnfac*nfix

        else:

           nfixation = 0.0  

  #             endif  
 
  # ! Photosynthesis stored in this array for diagnostic purposes only.
        if doDiagnostics:
          f1_out[k] = f1_out[k]+phygrow*bgc_dt
          f2_out[k] = f2_out[k]+graz*bgc_dt
          f6_out[k] = f6_out[k]+nfixation*bgc_dt

  # ! Collect all euphotic zone fluxes in these arrays.
        topo4[k]=-phygrow+zooexu
        todop[k]= graztodop*(1.0-ACeff)*graz \
          +graztodop*(phyexu+zooloss) \
          +phyloss
        tooxy[k]= tooxy[k]+(phygrow-zooexu)*ro2ut
        tophy[k]= phygrow-graz-phyexu-phyloss
        tozoo[k]= ACeff*graz-zooexu-zooloss
        todet[k] = (1.0-graztodop)*(1.0-ACeff)*graz \
          + (1.0-graztodop)*(phyexu+zooloss)
        todin[k]=topo4[k]*rnp + nfixation

  # ! For Fe:
        todfe[k]=topo4[k]*rfep + todop[k]*rfep    
        topfe[k]=todet[k]*rfep

  #     !end loop over euphotic zone

  # ! Explicit sinking of detritus in seperate loop. 
      flux_u = 0.0

      for k in range(bgc_kloc-1): #loop over all layers except last one, which accounts for burial (see below)

        DET = np.fmax(bgc_tracer[idet][k]-alimit*alimit,0.0)
        flux_l=wdet[k]*DET
        flux[k]  = flux[k]+flux_u
        fdiv[k] = fdiv[k]+(flux_u-flux_l)/bgc_dz[k]
        flux_u=flux_l          

      flux_l = 0.0

  # ! account for burial in the sediment 
      DET = np.fmax(bgc_tracer[idet][bgc_kloc-1]-alimit*alimit,0.0)
      fDET = wdet[bgc_kloc-1]*DET
      flux_l = np.fmin(1.0,burdige_fac*fDET**burdige_exp)*fDET
      flux_bury = flux_bury + flux_l*bgc_dt

      flux[bgc_kloc-1] = flux[bgc_kloc-1]+flux_u
      fdiv[bgc_kloc-1] = fdiv[bgc_kloc-1]+(flux_u-flux_l)/bgc_dz[bgc_kloc-1]

  # ! Store flux for diagnostic purposes. Flux is at the top of a box; for box # 1 it is the burial
      if doDiagnostics:
        f3_out[0]  = flux_bury

        for k in range(1,bgc_kloc):
          f3_out[k]  = f3_out[k]+flux[k]*bgc_dt

  # ! For Fe: Explicit sinking of pFe in seperate loop 
  # ! Sedimentary release of dFe depends on respiration divided by available oxygen
  # ! In contrast to Dale et al. (2015) and Somes et al. (2021) I don't calculate the dependence from the ratio of Cox and O2, but rather from O2-consumption and O2
  # ! In principle, o2demand/OXY is essentially the height of the water column [m] that would be stripped of oxygen within this time step because of the non-buried organic matter
  # ! This lowers the value for dFe release a bit.
  # ! 2024-04-16: Changed to omitting multiplication by bgc_dt. This is now the oxidative demand over a day. 

      o2demand = (fDET-flux_l)*ro2ut  # the oxidative demand of all organic matter that is NOT buried within a time step [mmol O2/m2]
      OXY = np.fmax(bgc_tracer[ioxy][bgc_kloc-1]-alimit*alimit,0.0) # bottom water oxygen

      flux_u = 0.0

      for k in range(bgc_kloc): #loop over all layers - assumes no boundary parameterisation for burial (in contrast to detritus)

        PFE = np.fmax(bgc_tracer[ipfe][k]-alimit*alimit,0.0)
        flux_l=wdet[k]*PFE
        feflux[k]  = feflux[k]+flux_u
        fefdiv[k]  = fefdiv[k]+(flux_u-flux_l)/bgc_dz[k]
        flux_u=flux_l         

  # ! sedimentary release of dFe: in contrast to Somes et al. (2021) do NOT release more dFE than would be buried
      fesed = np.fmin(flux_l,fesedmax*np.tanh(o2demand/OXY))

  # ! Store flux for diagnostic purposes. Flux is at the top of a box; for box # 1 it is the burial

      if doDiagnostics:
        f11_out[0]  = f11_out[0] + flux_l*bgc_dt

        for k in range(1,bgc_kloc):
          f11_out[k]  = f11_out[k]+feflux[k]*bgc_dt

      flux_l = 0.0

      if useCarbon:
  # ! effect of CaCO3 production on alkalinity and DIC
        caco3_prod = 0.0
        for k in range(bgc_kloc): # IK 2023-03-16: This should be bgc_kloc!
          caco3_prod = caco3_prod + todet[k]*rcp*frac_caco3*bgc_dz[k]
          toalk[k] = toalk[k] - 2.0 * todet[k]*rcp*frac_caco3 # mmol Alk/m3/d
          todic[k] = todic[k] - 1.0 * todet[k]*rcp*frac_caco3 # mmol Alk/m3/d
          if doDiagnostics:
            f9_out[k] = f9_out[k] + todet[k]*rcp*frac_caco3*bgc_dt

  # ! effect of CaCO3 flux divergence in aphotic layers on alkalinity and CaCO3
        for k in range(bgc_kloc-1):
          fdiv_caco3 = caco3_prod*(fcaco3[k]-fcaco3[k+1])/bgc_dz[k]
          toalk[k] = toalk[k] + 2.0 * fdiv_caco3
          todic[k] = todic[k] + 1.0 * fdiv_caco3
          if doDiagnostics:
            f8_out[k] = f8_out[k] + fdiv_caco3*bgc_dt
 
  # ! The bottom box is closed = no exchange with the sediment
        fdiv_caco3 = caco3_prod*fcaco3[bgc_kloc-1]/bgc_dz[bgc_kloc-1]
        toalk[bgc_kloc-1] = toalk[bgc_kloc-1] + 2.0 * fdiv_caco3
        todic[bgc_kloc-1] = todic[bgc_kloc-1] + 1.0 * fdiv_caco3
        if doDiagnostics:
          f8_out[bgc_kloc-1] = f8_out[bgc_kloc-1] + fdiv_caco3*bgc_dt

  # ! PROCESSES AFFECTING THE ENTIRE WATER COLUMN

      for k in range(bgc_kloc):

        DOP = np.fmax(bgc_tracer[idop][k]-alimit*alimit,0.0)
        PHY = np.fmax(bgc_tracer[iphy][k]-alimit*alimit,0.0)
        ZOO = np.fmax(bgc_tracer[izoo][k]-alimit*alimit,0.0)
        DET = np.fmax(bgc_tracer[idet][k]-alimit*alimit,0.0)
        DFE = np.fmax(bgc_tracer[idfe][k]-alimit*alimit,0.0)
        PFE = np.fmax(bgc_tracer[ipfe][k]-alimit*alimit,0.0)

  #     c AEROBIC DECAY
  # 
  #     c In contrast to the older (Kriest&Oschlies, 2013) version, this option:
  #     c (1) does not degrade OM in the absence of O2, i.e. OM can accumulate 
  #     c (2) uses a Michaelis-Menten Kinetic to slow down bacterial remineralisation under low O2
  #     c (2) takes care not to use more O2 per timestep than available
  # 
  #     c Michaelis-Menten limitation for oxic degradation: 

        OXY = np.fmax(bgc_tracer[ioxy][k]-subox,0.0)
        oxymm = OXY*OXY/(OXY*OXY+ACkbaco2*ACkbaco2)

  #     c O2 required for total remineralisation in a time step will then be:

        o2req = oxymm*(dlambda*DOP+detlambda*DET)*ro2ut*bgc_dt      

  #     c restrict remineralisation to amount of vaialable oxygen

        if (o2req > 0.0):
           o2usefrac = np.fmin(OXY,o2req)/o2req
        else:
           o2usefrac = 0.0
  #             endif

        remindop = oxymm*dlambda*DOP*o2usefrac
        remindet = oxymm*detlambda*DET*o2usefrac

  #     c ANAEROBIC DECAY INCL. ANAMMOX ETC.

        if (OXY < 36.0):

          DIN = np.fmax(bgc_tracer[idin][k]-subdin,0.0)
          dinmm = DIN*DIN/(DIN*DIN+ACkbacdin*ACkbacdin)*(1.0-oxymm)

  #       c NO3 required for total remineralisation in a time step will then be:

          dinreq = dinmm*(dlambda*DOP+detlambda*DET)*rhno3ut*bgc_dt

  #       c restrict remineralisation to amount of variable oxygen

          if (dinreq > 0.0):
             dinusefrac = np.fmin(DIN,dinreq)/dinreq
          else:
             dinusefrac = 0.0

  #       c restrict anaerobic processes to regions with low oxygen concentration

          denitdop = dinmm*dlambda*DOP*dinusefrac
          denitdet = dinmm*detlambda*DET*dinusefrac

        else:

          denitdop = 0.0
          denitdet = 0.0

        topo4[k]=topo4[k]+remindop+remindet+denitdop+denitdet
        todop[k]=todop[k]-remindop-denitdop \
          +plambda*PHY \
          +zlambda*ZOO     
        tooxy[k]=tooxy[k]-(remindop+remindet)*ro2ut
        tophy[k]=tophy[k]-plambda*PHY
        tozoo[k]=tozoo[k]-zlambda*ZOO
        todet[k]=todet[k]-remindet-denitdet
        todin[k]=todin[k]+(remindop+remindet)*rnp \
          -(denitdop+denitdet)*rhno3ut 
        if doDiagnostics:
          f4_out[k] = f4_out[k] + (remindop+remindet)*bgc_dt
          f7_out[k] = f7_out[k] + (denitdop+denitdet)*bgc_dt

  #       ! For Fe: release of Fe from within detritus (O2 and NO3 dependent) and additional release
        todfe[k]=todfe[k]+ (
        (remindet+denitdet)*rfep \
          +(plambda*PHY+zlambda*ZOO )*rfep
  #       !     &                 +(remindop+denitdop)*rfep   
          +detlambda*np.fmax(0.0,PFE-DET*rfep) 
        )
        topfe[k]=topfe[k]-(remindet+denitdet)*rfep \
          -detlambda*np.fmax(0.0,PFE-DET*rfep)    

  #       ! For Fe: scavenging by organic matter and precipipation

        o2sat = o2saturation(bgc_theta[k],bgc_salt[k])
        AOU = np.fmax(40.0,o2sat-OXY)
        ligands = np.fmax(fealpha*AOU**0.8+febeta*DOP**0.8,0.5)
        fea = 1.0+kfeleq*(ligands-DFE)
        feprime = ((fea**2.0+4.0*kfeleq*DFE)**0.5-fea) \
          /(2.0*kfeleq)
        feorgads = kfeorg*feprime*DET**0.58
        feprecip = kfepre*feprime**2.0*np.tanh(np.fmax(OXY,0.0)) # not in the paper, but in the code

        todfe[k] = todfe[k] - feorgads - feprecip
        topfe[k] = topfe[k] + feorgads + feprecip

        if doDiagnostics:
          if useCarbon:
            f10_out[k] = f10_out[k] + (feprecip+feorgads)*bgc_dt
          else:
            f8_out[k]  = ligands
            f9_out[k]  = f9_out[k] + feorgads*bgc_dt
            f10_out[k] = f10_out[k] + feprecip*bgc_dt

  #           ENDDO

  # ! RESUPPLY OF BURIED MATTER VIA RIVER RUNOFF OR VIA SURFACE

      if useRunoff:
        for k in range(bgc_kloc):
          runoff[k] = bgc_globalrunoff * bgc_runoffvol[k]
          if doDiagnostics:
            f5_out[k] = f5_out[k] + runoff[k]*bgc_dt
      else:
        runoff[0] = bgc_globalrunoff/bgc_dz[0]
        if doDiagnostics:
          f5_out[0] = f5_out[0] + runoff[0]*bgc_dt

  # ! UPDATE MASS CONCENTRATIONS 

      for k in range(bgc_kloc):
  # ! Update tracer concentrations, by adding the fluxes scaled by 
  # ! time step length.
        bgc_tracer[ipo4][k] = bgc_tracer[ipo4][k] + \
          topo4[k]*bgc_dt + runoff[k]*bgc_dt
        bgc_tracer[idop][k] = bgc_tracer[idop][k] + \
          todop[k]*bgc_dt
        bgc_tracer[ioxy][k] = bgc_tracer[ioxy][k] +  \
          tooxy[k]*bgc_dt
        bgc_tracer[iphy][k]= bgc_tracer[iphy][k]  + \
          tophy[k]*bgc_dt
        bgc_tracer[izoo][k]= bgc_tracer[izoo][k]  + \
          tozoo[k]*bgc_dt
        bgc_tracer[idet][k]= bgc_tracer[idet][k]  + \
          (todet[k]+fdiv[k])*bgc_dt
        bgc_tracer[idin][k]= bgc_tracer[idin][k]  + \
          todin[k]*bgc_dt + rnp*runoff[k]*bgc_dt

        if useCarbon:
          bgc_tracer[idic][k]= bgc_tracer[idic][k]  + (
            rcp*topo4[k]*bgc_dt # this is from organic production and decay   
            + todic[k]*bgc_dt # this is from CaCO3 production and dissolution
            + rcp*runoff[k]*bgc_dt # this is DIC runoff
            )
          bgc_tracer[ialk][k]= bgc_tracer[ialk][k]  - (
            (todin[k] + topo4[k])*bgc_dt) + (  # this is from organic production and decay    
            + toalk[k]*bgc_dt    # this is from CaCO3 production and dissolution
            )
          
  # ! For Fe - including sources from deposition and hydrothermal fluxes
        bgc_tracer[idfe][k]= bgc_tracer[idfe][k]  + \
          (todfe[k]+bgc_fedep[k])*bgc_dt
        bgc_tracer[ipfe][k]= bgc_tracer[ipfe][k]  + \
          (topfe[k]+fefdiv[k])*bgc_dt 

  #       ENDDO 

      bgc_tracer[ioxy][0]= bgc_tracer[ioxy][0] +  \
        o2gasex/bgc_dz[0]*bgc_dt

  # ! For Fe: sedimentary release
      bgc_tracer[idfe][bgc_kloc-1]= bgc_tracer[idfe][bgc_kloc-1] +  \
        fesed/bgc_dz[bgc_kloc-1]*bgc_dt

      if doDiagnostics:
        f12_out[bgc_kloc-1] = f12_out[bgc_kloc-1] + fesed*bgc_dt 

  # ! Air-sea gas exchange of CO2 and E-P
  #ifdef CARBON
      if useCarbon:
        bgc_tracer[idic][0]= bgc_tracer[idic][0] + \
         (co2emp+co2gasex)/bgc_dz[0]*bgc_dt

        bgc_tracer[ialk][0]= bgc_tracer[ialk][0] +  \
         (alkemp)/bgc_dz[0]*bgc_dt

  #endif

  #       ENDDO !internal time loop

    if useCarbon:
      localph[ip]=sph
      localco2airseaflux[ip]=co2airseaflux
    
    localflux_bury=localflux_bury+flux_bury*localdA[ip]
  
  return localflux_bury

@njit
def __phi(u):
# !      phi= u*(0.555588d0+0.004926d0*u)/(1.0d0+0.188721d0*u)
  if (u > 1.0e-6):
    phi= np.log(u+np.sqrt(1.0+u*u))-(np.sqrt(1.0+u*u)-1.0)/u
  else:
    phi=0.0

  return phi