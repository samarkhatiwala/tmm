import numpy as np
from petsc4py import PETSc
import tmm4py as TMM
from numba import njit

class DotDict(dict):
    """dot.notation access to dictionary attributes"""

    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, d):
        return self.__dict__.update(d)


# from PaTh_PARAMS import *

kmax=60
pathtrcmin=5e-12
rhosw=1024.5
avogradroNumber=6.02214076e23

betaPa=2.33e-3 # production rate in ocean in dpm/m^3/y
betaTh=2.52e-2 # production rate in ocean in dpm/m^3/y    
lambdaDecayPa=2.13e-5 # y^-1
lambdaDecayTh=9.22e-6 # y^-1

# Convert units
lambdaDecayPa=lambdaDecayPa/(86400.0*365.0) # y^-1 -> s^-1
lambdaDecayTh=lambdaDecayTh/(86400.0*365.0) # y^-1 -> s^-1

fmoltodpmPa=((1e-15)*avogradroNumber*lambdaDecayPa*60.0) # fmol Pa to dpm
fmoltodpmTh=((1e-15)*avogradroNumber*lambdaDecayTh*60.0) # fmol Th to dpm

betaPa=betaPa/(86400.0*365.0)/fmoltodpmPa # production rate: dpm/m^3/y -> fmol/m^3/sec
betaTh=betaTh/(86400.0*365.0)/fmoltodpmTh # production rate: dpm/m^3/y -> fmol/m^3/sec

# set defaults
PaKref=1.e7
ThKref=1.e7
KPaPOMFac=1.0
KPaCaCO3Fac=1.0/40.0
KPaOpalFac=1.0/6.0
KPaDustFac=1.0/20.0      
KPaLithFac=1.0/20.0
KThPOMFac=1.0
KThCaCO3Fac=1.0
KThOpalFac=1.0/20.0
KThDustFac=1.0/20.0
KThLithFac=1.0/20.0
      
#
KPaPOM=PaKref*KPaPOMFac*np.ones(kmax)
KPaCaCO3=PaKref*KPaCaCO3Fac*np.ones(kmax)
KPaOpal=PaKref*KPaOpalFac*np.ones(kmax)
KPaDust=PaKref*KPaDustFac*np.ones(kmax)
KPaLith=PaKref*KPaLithFac*np.ones(kmax)

KThPOM=ThKref*KThPOMFac*np.ones(kmax)
KThCaCO3=ThKref*KThCaCO3Fac*np.ones(kmax)
KThOpal=ThKref*KThOpalFac*np.ones(kmax)
KThDust=ThKref*KThDustFac*np.ones(kmax)
KThLith=ThKref*KThLithFac*np.ones(kmax)

class PaTh:

    p = DotDict()

    def __init__(self, efctxId, prefix):
      self.efctxId = efctxId
      self.prefix = prefix
      self.ef = DotDict()
      
    def iniExternalForcingFn(self, tc, Iter, state, *args, **kwargs): 
      p=PaTh.p
      myId = PETSc.COMM_WORLD.getRank()
      OptDB=PETSc.Options()
      ef=self.ef

#     First set data common to all instances of PaTh
      if self.efctxId==1:
        p.myId = myId
        p.profileConfig=TMM.getProfileConfig()
        p.timeConfig=TMM.getTimeConfig()

      lSize=p.profileConfig['lSize']
      lNumProfiles=p.profileConfig['lNumProfiles']
      lProfileLength=p.profileConfig['lProfileLength']
      lStartIndices=p.profileConfig['lStartIndices']
      maxSteps=p.timeConfig['maxSteps']
      Iter0=p.timeConfig['Iter0']
      
      if self.efctxId==1:
        try:
          p.DeltaT=OptDB.getReal('biogeochem_deltat')
        except:
          PETSc.Error("Must indicate biogeochemical time step in seconds with the -biogeochem_deltat option")
        PETSc.Sys.Print("Ocean time step for Pa/Th is  %12.7f seconds" % p.DeltaT)

#       Grid data
#       drF.bin is in MOBI's units of cm; we convert to meter here
        f = open("drF.bin", 'rb')
        p.nzmax = np.fromfile(f, dtype=np.dtype('>i4'), count=1)[0]
        drF = np.fromfile(f, dtype=np.dtype('>f8'), count=p.nzmax)/100.
        f.close()
        p.dzmr=1./drF[:] # reciprocal of cell thickness in m^-1

        PETSc.Sys.Print("Number of vertical layers is %d" % p.nzmax)

#       Forcing data
#       Vertical sinking velocity profiles (the astype conversion is for numba)
        with open("wPOM_siddall.bin", 'rb') as f:
          p.wPOM = np.fromfile(f, dtype=np.dtype('>f8'), count=p.nzmax).astype(np.float64)

        with open("wOpal_siddall.bin", 'rb') as f:
          p.wOpal = np.fromfile(f, dtype=np.dtype('>f8'), count=p.nzmax).astype(np.float64)

        with open("wCaCO3_siddall.bin", 'rb') as f:
          p.wCaCO3 = np.fromfile(f, dtype=np.dtype('>f8'), count=p.nzmax).astype(np.float64)

        with open("wLith_siddall.bin", 'rb') as f:
          p.wLith = np.fromfile(f, dtype=np.dtype('>f8'), count=p.nzmax).astype(np.float64)

#       Vertically uniform speed for dust
        p.wDust=1000.0 # m/y
        p.wDust=p.wDust/(86400.0*365.0) # dust sinking speed in m/s

        p.periodicBiogeochemForcing=OptDB.hasName('periodic_biogeochem_forcing')

        if p.periodicBiogeochemForcing:
          PETSc.Sys.Print("Periodic biogeochemical forcing specified")
          p.biogeochemTimer=TMM.PeriodicTimer()
          p.biogeochemTimer.create(prefix="periodic_biogeochem_")

#       Particle fields in kg/m^3
        p.localPaTh_pom=TMM.loadVecIntoArray("PaTh_pom_siddall.petsc")
        p.localPaTh_opal=TMM.loadVecIntoArray("PaTh_opal_siddall.petsc")
        p.localPaTh_caco3=TMM.loadVecIntoArray("PaTh_caco3_siddall.petsc")
        p.localPaTh_lith=TMM.loadVecIntoArray("PaTh_lith_siddall.petsc")
        p.localPaTh_dust=np.zeros(lSize) # this is calculated later from dust deposition

#       dustdep is the surface dust flux in kg/m^2/s
        if p.periodicBiogeochemForcing:
          p.localdustdep=np.zeros(lNumProfiles)
          p.localdustdepp=TMM.PeriodicArray()
          p.localdustdepp.create(lNumProfiles,buf=p.localdustdep)
        else:
          p.localdustdep=TMM.readProfileScalarData("dust_dep.bin")
          for ip in range(lNumProfiles):
            nzloc=lProfileLength[ip]
            kl=lStartIndices[ip]
            p.localPaTh_dust[kl:kl+nzloc]=p.localdustdep[ip]/p.wDust # kg/m^3

#   } /* end common data */

#     Now set data for this PaTh instance
      ef.config=state.getConfig()
      numTracers=ef.config['numTracers']

      for itr in range(numTracers):
        state.qef[itr][:]=0.0

      OptDBState=PETSc.Options(self.prefix)

#     Data for diagnostics
#     Diagnostics are always calculated within protac_thor_driver but we only 
#     accumulate and write them if requested via the -calc_diagnostics option
      ef.num3dDiag=4
      ef.num2dDiag=2
      ef.local3dDiag=[np.zeros(lSize) for _ in range(ef.num3dDiag)]
      ef.local2dDiag=[np.zeros(lNumProfiles) for _ in range(ef.num2dDiag)]
      ef.calcDiagnostics = OptDBState.hasName('calc_diagnostics')
      if ef.calcDiagnostics:
        ef.diagTimer = TMM.StepTimer()
        ef.diagTimer.create(prefix="diag_",startTimeStep=Iter0)
        PETSc.Sys.Print("Diagnostics will be computed starting at (and including) time step: %d" % ef.diagTimer.startTimeStep)
        PETSc.Sys.Print("Diagnostics will be computed over %d time steps\n" % ef.diagTimer.numTimeSteps)

        ef.diag3dOutFile = ["O_protacd.petsc", "O_protacb.petsc", "O_thord.petsc", "O_thorb.petsc"]
        ef.diag2dOutFile = ["O_protac_bottom_flux.bin", "O_thor_bottom_flux.bin"]

        ef.diagMode = "w"
        ef.appendDiagnostics = False
        ef.local3dDiagavg=[np.zeros(lSize) for _ in range(ef.num3dDiag)]
        ef.local2dDiagavg=[np.zeros(lNumProfiles) for _ in range(ef.num2dDiag)]

#     Write out configuration
      PETSc.Sys.Print("betaPa=%g fmol/m^3/sec" % (betaPa))
      PETSc.Sys.Print("betaTh=%g fmol/m^3/sec" % (betaTh))
      PETSc.Sys.Print("lambdaDecayPa=%g 1/sec" % (lambdaDecayPa))
      PETSc.Sys.Print("lambdaDecayTh=%g 1/sec" % (lambdaDecayTh))
      PETSc.Sys.Print("fmoltodpmPa=%g fmol/dpm" % (fmoltodpmPa))
      PETSc.Sys.Print("fmoltodpmTh=%g fmol/dpm" % (fmoltodpmTh))
      PETSc.Sys.Print("Pa scavenging coefficients:")
      PETSc.Sys.Print("KPaPOM=%g" % (KPaPOM[0]))
      PETSc.Sys.Print("KPaCaCO3=%g" % (KPaCaCO3[0]))
      PETSc.Sys.Print("KPaOpal=%g" % (KPaOpal[0]))
      PETSc.Sys.Print("KPaDust=%g" % (KPaDust[0]))
      PETSc.Sys.Print("  Pa lithogenic scavenging coefficient:")
      PETSc.Sys.Print("KPaLith=%g" % (KPaLith[0]))

      PETSc.Sys.Print("Th scavenging coefficients:")
      PETSc.Sys.Print("KThPOM=%g" % (KThPOM[0]))
      PETSc.Sys.Print("KThCaCO3=%g" % (KThCaCO3[0]))
      PETSc.Sys.Print("KThOpal=%g" % (KThOpal[0]))
      PETSc.Sys.Print("KThDust=%g" % (KThDust[0]))
      PETSc.Sys.Print("  Th lithogenic scavenging coefficient:")
      PETSc.Sys.Print("KThLith=%g" % (KThLith[0]))
      PETSc.Sys.Print("Sinking speed profile for Pa-Th:")
      PETSc.Sys.Print("k  wPOM   wCaCO3 wOpal wDust wLith")
      if p.myId==0:
        for k in range(p.nzmax):
          PETSc.Sys.Print("k=%d, %g, %g, %g, %g, %g" % (k,p.wPOM[k],p.wCaCO3[k],p.wOpal[k],p.wDust,p.wLith[k]))
      
      TMM.barrier()

    def calcExternalForcingFn(self, tc, Iter, iLoop, state, *args, **kwargs):

      p=PaTh.p
      ef=self.ef

      lSize=p.profileConfig['lSize']
      lNumProfiles=p.profileConfig['lNumProfiles']
      lProfileLength=p.profileConfig['lProfileLength']
      lStartIndices=p.profileConfig['lStartIndices']

      Iter0=p.timeConfig['Iter0']

      numTracers=ef.config['numTracers']
      
      myTime = p.DeltaT*Iter # /* Iter should start at 0 */

#     Common data only needs to be updated once
      if self.efctxId==1:
        if p.periodicBiogeochemForcing:
          p.localdustdepp.interp(tc, p.biogeochemTimer, "dust_dep_")
          for ip in range(lNumProfiles):
            nzloc=lProfileLength[ip]
            kl=lStartIndices[ip]
            p.localPaTh_dust[kl:kl+nzloc]=p.localdustdep[ip]/p.wDust # kg/m^3

      doDiagnostics = False
      if ef.calcDiagnostics:
        if (Iter0+iLoop>=ef.diagTimer.startTimeStep): #{ /* start time averaging (note: startTimeStep is ABSOLUTE time step) */
          doDiagnostics=True
            
      for ip in range(lNumProfiles):
        nzloc=lProfileLength[ip]
        kl=lStartIndices[ip]

#       Get pointers to data for vertical profile
#       Tracers
        Pa = state.c[0][kl:kl+nzloc]
        Th = state.c[1][kl:kl+nzloc]

#       Forcing
        POM = p.localPaTh_pom[kl:kl+nzloc]
        Opal = p.localPaTh_opal[kl:kl+nzloc]
        CaCO3 = p.localPaTh_caco3[kl:kl+nzloc]
        Lith = p.localPaTh_lith[kl:kl+nzloc]
        Dust = p.localPaTh_dust[kl:kl+nzloc]
#       Diagnostic arrays
        local3dDiag=[ef.local3dDiag[i][kl:kl+nzloc] for i in range(ef.num3dDiag)]

#       Call the main computational function
        FPa_bot, FTh_bot = protac_thor_driver(nzloc, p.dzmr, p.DeltaT, Pa, Th, \
           POM, CaCO3, Opal, Lith, Dust, \
           p.wPOM, p.wCaCO3, p.wOpal, p.wLith, p.wDust, \
           local3dDiag[0], local3dDiag[1], local3dDiag[2], local3dDiag[3])

        ef.local2dDiag[0][ip]=FPa_bot
        ef.local2dDiag[1][ip]=FTh_bot
#     end loop over profiles

      TMM.barrier()
      
    def writeExternalForcingFn(self, tc, Iter, iLoop, state, *args, **kwargs):

      p=PaTh.p
      ef=self.ef

      Iter0=p.timeConfig['Iter0']

      numTracers=ef.config['numTracers']
      
      if (ef.calcDiagnostics):
        if (Iter0+iLoop>=(ef.diagTimer.startTimeStep)):  #/* start time averaging (note: startTimeStep is ABSOLUTE time step) */  
          if (ef.diagTimer.count<=ef.diagTimer.numTimeSteps): #/* still within same averaging block so accumulate */
            for i in range(ef.num3dDiag):
              ef.local3dDiagavg[i][:] = ef.local3dDiagavg[i][:] + ef.local3dDiag[i][:]
            for i in range(ef.num2dDiag):
              ef.local2dDiagavg[i][:]=ef.local2dDiag[i][:]+ef.local2dDiagavg[i][:]
#           Increment the count
            ef.diagTimer.incr()
          
          if ((ef.diagTimer.count)==(ef.diagTimer.numTimeSteps)): # time to write averages to file
            PETSc.Sys.Print("Writing diagnostics time average at time %10.5f, step %d" % (tc, Iter0+iLoop))
            for i in range(ef.num3dDiag):
              ef.local3dDiagavg[i][:] = ef.local3dDiagavg[i][:]/ef.diagTimer.count
              TMM.writeArrayToVec(ef.diag3dOutFile[i], ef.local3dDiagavg[i], ef.diagMode)
#             reset diagnostic arrays
              ef.local3dDiagavg[i][:] = 0.0
            for i in range(ef.num2dDiag):
               ef.local2dDiagavg[i][:]=ef.local2dDiagavg[i][:]/ef.diagTimer.count
               TMM.writeProfileScalarData(ef.diag2dOutFile[i], ef.local2dDiagavg[i], 1, ef.appendDiagnostics)
#              reset diagnostic arrays
               ef.local2dDiagavg[i][:]=0.0

            ef.diagMode = "a"
            ef.appendDiagnostics=True
            ef.diagTimer.update(Iter0+iLoop)
      
    def finalizeExternalForcingFn(self, tc, Iter, state):

      p=PaTh.p
      ef=self.ef

      TMM.barrier()
      
#     Delete common data
      if (self.efctxId==1):
        if (p.periodicBiogeochemForcing):
          p.localdustdepp.destroy()

    def reInitializeExternalForcingFn(self, tc, Iter, iLoop, state):
      pass

@njit(cache=False)
def protac_thor_driver(kmx, dzmr, twodt, Pa, Th, \
      POM, CaCO3, Opal, Lith, Dust, \
      wPOM, wCaCO3, wOpal, wLith, wDust, \
      Pad, Pab, Thd, Thb):

# #     Input POM, CaCO3, Opal, Lith and Dust are in kg/m^3

# Number of internal sub time steps
  ntpath = 2
  dtpath = twodt/ntpath
  
  FPa_bot=0.0
  FTh_bot=0.0

  for it in range(ntpath):
    FPa_in=0.0
    FTh_in=0.0
    for k in range(kmx):
      Pam=np.fmax(Pa[k], pathtrcmin) # total concentration in fmol/m^3		
      SPa_pom=KPaPOM[k]*np.fmax(POM[k],pathtrcmin)/rhosw
      SPa_caco3=KPaCaCO3[k]*np.fmax(CaCO3[k],pathtrcmin)/rhosw
      SPa_opal=KPaOpal[k]*np.fmax(Opal[k],pathtrcmin)/rhosw
      SPa_lith=KPaLith[k]*np.fmax(Lith[k],pathtrcmin)/rhosw
      SPa_dust=KPaDust[k]*np.fmax(Dust[k],pathtrcmin*1.e-3)/rhosw
      SPa=SPa_pom+SPa_caco3+SPa_opal+SPa_lith+SPa_dust
      Pad[k]=(1.0/(1.0+SPa))*Pam # dissolved concentration in fmol/m^3
      Pab[k]=Pam-Pad[k] # particle-associated concentration in fmol/m^3
      Pa_pom=SPa_pom*Pad[k] # particle-associated concentration in fmol/m^3
      Pa_caco3=SPa_caco3*Pad[k] # particle-associated concentration in fmol/m^3
      Pa_opal=SPa_opal*Pad[k] # particle-associated concentration in fmol/m^3
      Pa_lith=SPa_lith*Pad[k] # particle-associated concentration in fmol/m^3
      Pa_dust=SPa_dust*Pad[k] # particle-associated concentration in fmol/m^3
#       if (k > 0): FPa_in = FPa_out
      FPa_out = wPOM[k]*Pa_pom + wCaCO3[k]*Pa_caco3 + \
                wOpal[k]*Pa_opal + wLith[k]*Pa_lith + \
                wDust*Pa_dust
      srcPa = betaPa - lambdaDecayPa*Pam + (FPa_in-FPa_out)*dzmr[k] # source term in fmol/m^3/s

      Thm=np.fmax(Th[k], pathtrcmin) # total concentration in fmol/m^3		
      STh_pom=KThPOM[k]*np.fmax(POM[k],pathtrcmin)/rhosw
      STh_caco3=KThCaCO3[k]*np.fmax(CaCO3[k],pathtrcmin)/rhosw
      STh_opal=KThOpal[k]*np.fmax(Opal[k],pathtrcmin)/rhosw
      STh_lith=KThLith[k]*np.fmax(Lith[k],pathtrcmin)/rhosw
      STh_dust=KThDust[k]*np.fmax(Dust[k],pathtrcmin*1.e-3)/rhosw
      STh=STh_pom+STh_caco3+STh_opal+STh_lith+STh_dust
      Thd[k]=(1.0/(1.0+STh))*Thm # dissolved concentration in fmol/m^3
      Thb[k]=Thm-Thd[k] # particle-associated concentration in fmol/m^3
      Th_pom=STh_pom*Thd[k] # particle-associated concentration in fmol/m^3
      Th_caco3=STh_caco3*Thd[k] # particle-associated concentration in fmol/m^3
      Th_opal=STh_opal*Thd[k] # particle-associated concentration in fmol/m^3
      Th_lith=STh_lith*Thd[k] # particle-associated concentration in fmol/m^3
      Th_dust=STh_dust*Thd[k] # particle-associated concentration in fmol/m^3
#       if (k > 0): FTh_in = FTh_out
      FTh_out = wPOM[k]*Th_pom + wCaCO3[k]*Th_caco3 + \
                wOpal[k]*Th_opal + wLith[k]*Th_lith + \
                wDust*Th_dust
      srcTh = betaTh - lambdaDecayTh*Thm + (FTh_in-FTh_out)*dzmr[k]

#     Flux into next box
      FPa_in = FPa_out
      FTh_in = FTh_out

#     Time step tracers
      Pa[k] = Pa[k] + dtpath*srcPa
      Th[k] = Th[k] + dtpath*srcTh
      
#   Accumulate bottom flux
    FPa_bot = FPa_bot + dtpath*FPa_out
    FTh_bot = FTh_bot + dtpath*FTh_out

# Bottom flux
  FPa_bot = FPa_bot/twodt # fmol/m^2/s
  FTh_bot = FTh_bot/twodt # fmol/m^2/s
      
  return FPa_bot, FTh_bot
