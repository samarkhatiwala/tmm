#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define READ_SWRAD
#include "petscmat.h"
#include "petsc_matvec_utils.h"
#include "tmm_main.h"
#include "tmm_forcing_utils.h"
#include "tmm_profile_utils.h"
#include "tmm_profile_data.h"
#include "tmm_main.h"
#include "MOBI_TMM_OPTIONS.h"
#include "mobi.h"

/* Macros to map tracer names to vectors */
/* v[0]=DIC,v[1]=Alk,v[2]=PO4,v[3]=DOP,v[4]=O2,v[5]=Fe */
#define TR1 v[0]
#define TR2 v[1]
#define TR3 v[2]
#define TR4 v[3]
#define TR5 v[4]
#define TR6 v[5]
#define TR7 v[6]
#define TR8 v[7]
#define TR9 v[8]
#define TR10 v[9]
#define TR11 v[10]
#define TR12 v[11]
#define TR13 v[12]
#define TR14 v[13]
#define TR15 v[14]
#define TR16 v[15]
#define TR17 v[16]
#define TR18 v[17]
#define TR19 v[18]
#define TR20 v[19]
#define TR21 v[20]
#define TR22 v[21]
#define TR23 v[22]
#define TR24 v[23]

#define JTR1 ut[0]
#define JTR2 ut[1]
#define JTR3 ut[2]
#define JTR4 ut[3]
#define JTR5 ut[4]
#define JTR6 ut[5]
#define JTR7 ut[6]
#define JTR8 ut[7]
#define JTR9 ut[8]
#define JTR10 ut[9]
#define JTR11 ut[10]
#define JTR12 ut[11]
#define JTR13 ut[12]
#define JTR14 ut[13]
#define JTR15 ut[14]
#define JTR16 ut[15]
#define JTR17 ut[16]
#define JTR18 ut[17]
#define JTR19 ut[18]
#define JTR20 ut[19]
#define JTR21 ut[20]
#define JTR22 ut[21]
#define JTR23 ut[22]
#define JTR24 ut[23]

Vec Ts,Ss,Fe;
PetscInt *gIndices;
PetscScalar *localTs,*localSs,*localFe;
PetscScalar **localTR, **localJTR;
PetscTruth useEmP = PETSC_FALSE;
PetscScalar *localEmP, EmPglobavg;
PeriodicArray localEmPp;
Vec surfVolFrac;
PetscScalar Sglobavg = 0.0, *TRglobavg;
PetscScalar *localwind,*localaice,*localhice,*localhsno,*localdz,*localatmosp;
PetscScalar *localswrad;
PetscScalar *locallatitude;
PetscScalar *localsgbathy;

PetscScalar daysPerYear, secondsPerYear;

PetscTruth useSeparateBiogeochemTimeStepping = PETSC_FALSE;
PetscInt nzmax;
PetscScalar DeltaT;
PetscScalar *drF, *zt;

PeriodicVec Tsp, Ssp,Fep;
PeriodicArray localwindp,localaicep,localhicep,localhsnop,localatmospp;
#ifdef READ_SWRAD
PeriodicArray localswradp;
#endif

PetscInt numBiogeochemPeriods;
PetscScalar *tdpBiogeochem; /* arrays for periodic forcing */
PetscTruth periodicBiogeochemForcing = PETSC_FALSE;
PetscScalar biogeochemCyclePeriod, biogeochemCycleStep;

/* atmospheric model variables */
PetscScalar *TpCO2atm_hist, *pCO2atm_hist;
PetscInt numpCO2atm_hist = 0;
PetscTruth fixedAtmosCO2 = PETSC_TRUE;
char pCO2atmIniFile[PETSC_MAX_PATH_LEN];  

PetscTruth useAtmModel = PETSC_FALSE;
PetscScalar pCO2atm_ini = 277.0; /* default initial value */
PetscScalar pCO2atm = 277.0; /* default initial value */
PetscScalar *localdA;
PetscScalar ppmToPgC=2.1324;
PetscScalar atmModelDeltaT;
PetscScalar secPerYear=86400.0*360.0;
PetscScalar Focean=0.0;
PetscScalar localFocean=0.0;
PetscScalar Foceanint = 0.0;
PetscInt atmModelUpdateTimeSteps=1;

PetscInt atmWriteSteps;
PetscTruth atmAppendOutput;
FILE *atmfptime;
PetscViewer atmfd;
PetscInt atmfp;
char atmOutTimeFile[PETSC_MAX_PATH_LEN];  

PetscTruth calcDiagnostics = PETSC_FALSE;
PetscInt diagNumTimeSteps, diagStartTimeStep, diagCount;
PetscTruth appendDiagnostics = PETSC_FALSE;
/* Add model specific diagnostic variables below */

#ifdef CARBON
PetscScalar *localco2airseafluxdiag, *localco2airseafluxdiagavg;
#endif

PetscMPIInt myId;
PetscInt debugFlag = 0;

PetscTruth TRUE = PETSC_TRUE, FALSE = PETSC_FALSE;

#if defined (FORSPINUP) || defined (FORJACOBIAN)
PetscScalar relaxTau[50], relaxLambda[50], relaxValue[50];
PetscTruth relaxTracer = PETSC_FALSE;
#endif

/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "iniExternalForcing"
PetscErrorCode iniExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt numTracers, Vec *v, Vec *ut)
{
  PetscErrorCode ierr;
  PetscInt gLow, gHigh, il;
  PetscInt ip, kl, nzloc;
  PetscInt itr;
  PetscViewer fd;
  PetscInt fp;
  PetscTruth flg;
  PetscInt it, m;
  PetscScalar myTime;
  PetscScalar zero = 0.0;
  
  daysPerYear = 360.0;
  secondsPerYear = 86400*daysPerYear;

#if defined (FORSPINUP) || defined (FORJACOBIAN)
  ierr = PetscOptionsHasName(PETSC_NULL,"-relax_tracer",&relaxTracer);CHKERRQ(ierr);
  if (relaxTracer) {  
    PetscInt maxValsToRead, itr;

    maxValsToRead = numTracers;
    ierr = PetscOptionsGetRealArray(PETSC_NULL,"-relax_tau",relaxTau,&maxValsToRead,&flg);
    if (!flg) SETERRQ(1,"Must indicate tracer relaxation tau with the -relax_tau option");
    if (maxValsToRead != numTracers) {
      SETERRQ(1,"Insufficient number of relaxation tau values specified");
    }

    maxValsToRead = numTracers;
    ierr = PetscOptionsGetRealArray(PETSC_NULL,"-relax_value",relaxValue,&maxValsToRead,&flg);
    if (!flg) SETERRQ(1,"Must indicate relaxation values with the -relax_value option");
    if (maxValsToRead != numTracers) {
      SETERRQ(1,"Insufficient number of relaxation values specified");
    }
    
    for (itr=0; itr<numTracers; itr++) {
      if (relaxTau[itr]>0.0) {
        relaxLambda[itr]=1.0/relaxTau[itr];
      } else {
        relaxLambda[itr]=0.0;
        relaxValue[itr]=0.0;
      }      
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d relaxation lambda=%15.11f, relaxation value=%10.8f\n",itr,relaxLambda[itr],relaxValue[itr]);CHKERRQ(ierr);
    }      
    
  }
#endif
    
/*   v[0]=TR1,v[1]=TR2,... */
  m = 0;
# if defined O_carbon
#define dic v[0]
#define sdic ut[0]
  m++;
#  if defined O_carbon_13
#define dic13 v[1]
#define sdic13 ut[1]
  m++;
#  endif
#  if defined O_carbon_14
#define c14 v[2]
#define sc14 ut[2]
  m++;
#  endif
# endif
# if defined O_npzd_o2
#define o2 v[3]
#define so2 ut[3]
  m++;
# endif
# if defined O_npzd_alk
#define alk v[4]
#define salk ut[4]
  m++;
# endif
# if defined O_npzd
#define po4 v[5]
#define spo4 ut[5]
  m++;
#define dop v[6]
#define sdop ut[6]
  m++;
#define phyt v[7]
#define sphyt ut[7]
  m++;
#define zoop v[8]
#define szoop ut[8]
  m++;
#define detr v[9]
#define sdetr ut[9]
  m++;
#  if defined O_npzd_nitrogen
#define no3 v[10]
#define sno3 ut[10]
  m++;
#define don v[11]
#define sdon ut[11]
  m++;
#define diaz v[12]
#define sdiaz ut[12]
  m++;
#   if defined O_npzd_nitrogen_15
#define din15 v[13]
#define sdin15 ut[13]
  m++;
#define don15 v[14]
#define sdon15 ut[14]
  m++;
#define phytn15 v[15]
#define sphytn15 ut[15]
  m++;
#define zoopn15 v[16]
#define szoopn15 ut[16]
  m++;
#define detrn15 v[17]
#define sdetrn15 ut[17]
  m++;
#define diazn15 v[18]
#define sdiazn15 ut[18]
  m++;
#   endif
#  endif
#  if defined O_carbon_13
#define doc13 v[19]
#define sdoc13 ut[19]
  m++;
#define phytc13 v[20]
#define sphytc13 ut[20]
  m++;
#define zoopc13 v[21]
#define szoopc13 ut[21]
  m++;
#define detrc13 v[22]
#define sdetrc13 ut[22]
  m++;
#   if defined O_npzd_nitrogen
#define diazc13 v[23]
#define sdiazc13 ut[23]
  m++;
#   endif
#  endif
# endif
  
  if (m != numTracers) {
    SETERRQ(1,"Error: numTracers does not match expected number of tracers!");    
  }
  
  for (itr=0; itr<numTracers; itr++) {
    ierr = VecSet(ut[itr],zero); CHKERRQ(ierr);
  }

  ierr = VecGetArrays(v,numTracers,&localTR);CHKERRQ(ierr);
  ierr = VecGetArrays(ut,numTracers,&localJTR);CHKERRQ(ierr);

  ierr = PetscOptionsGetTruth(PETSC_NULL,"-separate_biogeochem_time_stepping",&useSeparateBiogeochemTimeStepping,0);CHKERRQ(ierr);
#if defined (FORSPINUP) || defined (FORJACOBIAN)
  if (useSeparateBiogeochemTimeStepping) {
    SETERRQ(1,"Cannot use the -separate_biogeochem_time_stepping option with SPINUP or JACOBIAN ");  
  
  }
#endif
  if (useSeparateBiogeochemTimeStepping) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Biogeochem model will be time-stepped independently\n");CHKERRQ(ierr);
  }  
  
  ierr = PetscOptionsGetInt(PETSC_NULL,"-nzmax",&nzmax,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(1,"Must indicate maximum number of z points with the -nzmax option");  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of vertical layers is %d \n",nzmax);CHKERRQ(ierr);

  ierr = PetscOptionsGetReal(PETSC_NULL,"-biogeochem_deltat",&DeltaT,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(1,"Must indicate biogeochemical time step in seconds with the -biogeochem_deltat option");  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Ocean time step for BGC length is  %12.7f seconds\n",DeltaT);CHKERRQ(ierr);

  ierr = PetscOptionsHasName(PETSC_NULL,"-periodic_biogeochem_forcing",&periodicBiogeochemForcing);CHKERRQ(ierr);

  if (periodicBiogeochemForcing) {    
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic biogeochemical forcing specified\n");CHKERRQ(ierr);

/*  read time data */
/*  IMPORTANT: time units must be the same as that used by the toplevel driver */
    ierr = PetscOptionsGetReal(PETSC_NULL,"-periodic_biogeochem_cycle_period",&biogeochemCyclePeriod,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(1,"Must indicate biogeochemical forcing cycling time with the -periodic_biogeochem_cycle_period option");
    ierr = PetscOptionsGetReal(PETSC_NULL,"-periodic_biogeochem_cycle_step",&biogeochemCycleStep,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(1,"Must indicate biogeochemical forcing cycling step with the -periodic_biogeochem_cycle_step option");
    numBiogeochemPeriods=biogeochemCyclePeriod/biogeochemCycleStep;
/*  array for holding extended time array */
    PetscMalloc((numBiogeochemPeriods+2)*sizeof(PetscScalar), &tdpBiogeochem); 
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Periodic biogeochemical forcing specified at times:\n");CHKERRQ(ierr);            
    for (it=0; it<=numBiogeochemPeriods+1; it++) {
      tdpBiogeochem[it]=(-biogeochemCycleStep/2.0) + it*biogeochemCycleStep;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"tdpBiogeochem=%10.5f\n", tdpBiogeochem[it]);CHKERRQ(ierr);        
    }    
  }
  
/*   Read T and S */
  ierr = VecDuplicate(TR1,&Ts);CHKERRQ(ierr);
  ierr = VecDuplicate(TR1,&Ss);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    Tsp.firstTime = PETSC_TRUE;
    Ssp.firstTime = PETSC_TRUE;
  } else {
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Ts.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoadIntoVector(fd,Ts);CHKERRQ(ierr);  
	ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);    
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Ss.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoadIntoVector(fd,Ss);CHKERRQ(ierr);    
	ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);    
  }  
  ierr = VecGetArray(Ts,&localTs);CHKERRQ(ierr);
  ierr = VecGetArray(Ss,&localSs);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Done reading T/S\n");CHKERRQ(ierr);

/*   Compute global indices for local piece of vectors */
  ierr = VecGetOwnershipRange(Ts,&gLow,&gHigh);CHKERRQ(ierr);
  gHigh = gHigh - 1; /* Note: gHigh is one more than the last local element */
  ierr = PetscMalloc(lSize*sizeof(PetscInt),&gIndices);CHKERRQ(ierr);  
  for (il=0; il<lSize; il++) {
    gIndices[il] = il + gLow;
  }  

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using prescribed atmospheric pCO2\n");CHKERRQ(ierr);  
  ierr = PetscOptionsGetReal(PETSC_NULL,"-pco2_atm",&pCO2atm,&flg);CHKERRQ(ierr); /* overwrite default value */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using fixed atmospheric pCO2 of %g ppm\n",pCO2atm);CHKERRQ(ierr);
      
  ierr = PetscOptionsHasName(PETSC_NULL,"-use_emp",&useEmP);CHKERRQ(ierr);
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localEmP);CHKERRQ(ierr); /* always need this */
  for (ip=0; ip<lNumProfiles; ip++) { /* initialize to zero to be safe */
    localEmP[ip]=0.0;
  }
  if (useEmP) {
    if (periodicBiogeochemForcing) {    
      localEmPp.firstTime = PETSC_TRUE;
      localEmPp.arrayLength = lNumProfiles;
    } else {  
      ierr = readProfileSurfaceScalarData("EmP.bin",localEmP,1);  
    }

	ierr = VecDuplicate(TR1,&surfVolFrac);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"surface_volume_fraction.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoadIntoVector(fd,surfVolFrac);CHKERRQ(ierr);  
	ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);      

	ierr = PetscMalloc(numTracers*sizeof(PetscScalar),&TRglobavg);CHKERRQ(ierr);	
  }

/* Grid arrays */
//  ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localdz);CHKERRQ(ierr);    
//  ierr = VecLoadVecIntoArray(TR1,"dz.petsc",localdz);CHKERRQ(ierr);

  ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localsgbathy);CHKERRQ(ierr);    
  ierr = VecLoadVecIntoArray(TR1,"sgbathy.petsc",localsgbathy);CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"zt.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
  ierr = PetscMalloc(nzmax*sizeof(PetscScalar),&zt);CHKERRQ(ierr); 
  ierr = PetscBinaryRead(fp,zt,nzmax,PETSC_SCALAR);CHKERRQ(ierr);  
  ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"drF.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
  ierr = PetscMalloc(nzmax*sizeof(PetscScalar),&drF);CHKERRQ(ierr); 
  ierr = PetscBinaryRead(fp,drF,nzmax,PETSC_SCALAR);CHKERRQ(ierr);  
  ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);

/* Forcing fields */  
  ierr = VecDuplicate(TR1,&Fe);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    Fep.firstTime = PETSC_TRUE;
  } else {
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Fe.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoadIntoVector(fd,Fe);CHKERRQ(ierr);    
	ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);    
  }  
  ierr = VecGetArray(Fe,&localFe);CHKERRQ(ierr);

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localswrad);CHKERRQ(ierr);  
#ifdef READ_SWRAD
  if (periodicBiogeochemForcing) {    
    localswradp.firstTime = PETSC_TRUE;
    localswradp.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("swrad.bin",localswrad,1);  
  }
#endif	  

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&locallatitude);CHKERRQ(ierr);  
  ierr = readProfileSurfaceScalarData("latitude.bin",locallatitude,1);  
  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localaice);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localaicep.firstTime = PETSC_TRUE;
    localaicep.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("aice.bin",localaice,1);  
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localhice);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localhicep.firstTime = PETSC_TRUE;
    localhicep.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("hice.bin",localhice,1);  
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localhsno);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localhsnop.firstTime = PETSC_TRUE;
    localhsnop.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("hsno.bin",localhsno,1);  
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localwind);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localwindp.firstTime = PETSC_TRUE;
    localwindp.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("wind.bin",localwind,1);  
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localatmosp);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localatmospp.firstTime = PETSC_TRUE;
    localatmospp.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("atmosp.bin",localatmosp,1);  
  }

/* Initialize biogeochem model */
  myTime = DeltaT*Iter; /* Iter should start at 0 */

  if (periodicBiogeochemForcing) {   
	ierr = interpPeriodicVector(tc,&Ts,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Ssp,"Ss_");	
	ierr = interpPeriodicVector(tc,&Fe,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Fep,"Fe_");		
#ifdef READ_SWRAD
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localswradp,"swrad_");
#else
   insolation_(&lNumProfiles,&myTime,&locallatitude[0],&localswrad[0],&localtau[0]);
#endif                                                 
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localaice,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localaicep,"aice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhice,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localhicep,"hice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhsno,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localhsnop,"hsno_");                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemCyclePeriod,numBiogeochemPeriods,
					                              tdpBiogeochem,&localwindp,"wind_");   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemCyclePeriod,numBiogeochemPeriods,
									              tdpBiogeochem,&localatmospp,"atmosp_");					                              
    if (useEmP) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                    tdpBiogeochem,&localEmPp,"EmP_");                                                  
    }
  } else {
#ifndef READ_SWRAD
    insolation_(&lNumProfiles,&myTime,&locallatitude[0],&localswrad[0],&localtau[0]);
#endif    
  }

/* compute global means */
  if (useEmP) {
	ierr = VecDot(surfVolFrac,Ss,&Sglobavg);CHKERRQ(ierr); /* volume weighted mean surface salinity */    
	for (itr=0; itr<numTracers; itr++) {    
      ierr = VecDot(surfVolFrac,v[itr],&TRglobavg[itr]);CHKERRQ(ierr); /* volume weighted mean surface TR */									                    
    }  
  }

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myId);CHKERRQ(ierr);
  if (myId == 0) debugFlag = 1;
  
  mobi_ini_(zt,drF,&DeltaT,&Sglobavg,TRglobavg,&debugFlag);  

  return 0;
}




/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/




#undef __FUNCT__
#define __FUNCT__ "calcExternalForcing"
PetscErrorCode calcExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop, PetscInt numTracers, Vec *v, Vec *ut)
{

  PetscErrorCode ierr;
  PetscInt itr, ip, nzloc, kl;
  PetscScalar myTime;
  PetscScalar relyr, day;
  PetscInt toMobi = 1; 
  PetscInt fromMobi = 2;
  PetscInt itf;
  PetscScalar alpha;
  PetscScalar localco2airseaflux = 0.0;

  myTime = DeltaT*Iter; /* Iter should start at 0 */

  if (periodicBiogeochemForcing) {   
	ierr = interpPeriodicVector(tc,&Ts,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Ssp,"Ss_");	
	ierr = interpPeriodicVector(tc,&Fe,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Fep,"Fe_");	
#ifdef READ_SWRAD
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localswradp,"swrad_");
#else
    insolation_(&lNumProfiles,&myTime,&locallatitude[0],&localswrad[0],&localtau[0]);
#endif                                                 
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localaice,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localaicep,"aice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhice,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localhicep,"hice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhsno,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localhsnop,"hsno_");                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemCyclePeriod,numBiogeochemPeriods,
					                              tdpBiogeochem,&localwindp,"wind_");  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemCyclePeriod,numBiogeochemPeriods,
									              tdpBiogeochem,&localatmospp,"atmosp_");					                              
    if (useEmP) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                    tdpBiogeochem,&localEmPp,"EmP_");                                                      
    }
  }

/* compute global means */
  if (useEmP) {
	for (itr=0; itr<numTracers; itr++) {    
      ierr = VecDot(surfVolFrac,v[itr],&TRglobavg[itr]);CHKERRQ(ierr); /* volume weighted mean surface TR */									              
    }
    EmPglobavg = 0.0; /* set this to zero for now */
  }

  relyr = myTime/secondsPerYear; /* number of years (and fractional years) of model */
  day = myTime/86400.0 - floor(relyr)*daysPerYear; /* relative day number referenced to the beginning of the current year */

  for (ip=0; ip<lNumProfiles; ip++) {
    nzloc=lProfileLength[ip];
    kl=lStartIndices[ip];
	for (itr=0; itr<numTracers; itr++) {    	
	  mobi_copy_data_(&nzloc,&itr,&localTR[itr][kl],&toMobi);
/* 
	  if (ip==0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"MOBI IN: %d, %g\n",itr+1,localTR[itr][kl]);CHKERRQ(ierr);
      }	  

 */
	}  
    mobi_calc_(&nzloc,&locallatitude[ip],&day,&relyr,
               &localTs[kl],&localSs[kl],&TRglobavg[0],
# if defined O_carbon
               &pCO2atm,&localwind[ip],
#endif      
#  if defined O_npzd_nitrogen
               &localsgbathy[kl],
#  endif
#  if defined O_npzd_fe_limitation
               &localFe[kl],
#  endif
#  if defined O_embm
               &localswrad[ip],
#  endif
#  if defined O_ice
#   if !defined O_ice_cpts
               &localaice[ip], &localhice[ip], &localhsno[ip],
#   endif
#  endif
               &localEmP[ip], &EmPglobavg, &debugFlag);

	for (itr=0; itr<numTracers; itr++) {    
	  mobi_copy_data_(&nzloc,&itr,&localJTR[itr][kl],&fromMobi);
/* 
	  if (ip==0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"MOBI OUT: %d, %g\n",itr+1,localJTR[itr][kl]);CHKERRQ(ierr);
      }

 */
	}  

  } /* end loop over profiles */

  if (useSeparateBiogeochemTimeStepping) {  /* return updated tracer field */
//	ierr = VecSetValues(TR1,lSize,gIndices,localTR1,INSERT_VALUES);CHKERRQ(ierr);
//	ierr = VecAssemblyBegin(TR1);CHKERRQ(ierr);
//	ierr = VecAssemblyEnd(TR1);CHKERRQ(ierr);    
  } else {  /* return tracer tendency */
	for (itr=0; itr<numTracers; itr++) {  
	  ierr = VecSetValues(ut[itr],lSize,gIndices,localJTR[itr],INSERT_VALUES);CHKERRQ(ierr);
	  ierr = VecAssemblyBegin(ut[itr]);CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(ut[itr]);CHKERRQ(ierr);    
    }
#if defined (FORSPINUP) || defined (FORJACOBIAN)
/* add relaxation term: ut = ut - lambda*(v-vr) = ut -lambda*v + lambda*vr */
    if (relaxTracer) {
      for (itr=0; itr<numTracers; itr++) {
        ierr = VecAXPY(ut[itr],-relaxLambda[itr],v[itr]);CHKERRQ(ierr); /* ut = ut - lambda*v */
        ierr = VecShift(ut[itr],relaxLambda[itr]*relaxValue[itr]);CHKERRQ(ierr); /* ut = ut + lambda*vr */
      }
    }
#endif

/* 
  if (Iter==0) {
  PetscViewer fd;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"ut.petsc",FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
  for (itr=0; itr<numTracers; itr++) {
    ierr = VecView(ut[itr],fd);CHKERRQ(ierr);
  }
  ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);      
  }

 */
  
/*  Convert to discrete tendency */
	for (itr=0; itr<numTracers; itr++) {
	  ierr = VecScale(ut[itr],DeltaT);CHKERRQ(ierr);
	}
  }
    
  return 0;
}




/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/




#undef __FUNCT__
#define __FUNCT__ "writeExternalForcing"
PetscErrorCode writeExternalForcing(PetscScalar tc, PetscInt iLoop, PetscInt numTracers, Vec *v, Vec *ut)
{

  PetscErrorCode ierr;
  PetscInt ip;
  PetscScalar zero = 0.0, one = 1.0;  

/* Note: tc and iLoop are the time and step at the end of the current time step. */

#ifdef CARBON
  if (useAtmModel) {
/* write instantaneous atmos model state */
    if ((iLoop % atmWriteSteps)==0) {  /*  time to write out */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing atmospheric model output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
      ierr = PetscFPrintf(PETSC_COMM_WORLD,atmfptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);           
      ierr = writeBinaryScalarData("pCO2atm_output.bin",&pCO2atm,1,PETSC_TRUE);
      ierr = writeBinaryScalarData("Foceanint_output.bin",&Foceanint,1,PETSC_TRUE);
      Foceanint = 0.0;
    }
  }
#endif

  return 0;
}




/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/



#undef __FUNCT__
#define __FUNCT__ "finalizeExternalForcing"
PetscErrorCode finalizeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt numTracers)
{

  PetscErrorCode ierr;

/* write final pickup */
#ifdef CARBON
  if (useAtmModel) {
/* write instantaneous atmos model state */
    ierr = writeBinaryScalarData("pickup_pCO2atm.bin",&pCO2atm,1,PETSC_FALSE);
  }
#endif

  ierr = VecDestroy(Ts);CHKERRQ(ierr);
  ierr = VecDestroy(Ss);CHKERRQ(ierr);
  ierr = VecDestroy(Fe);CHKERRQ(ierr);
  
  ierr = PetscFree(gIndices);CHKERRQ(ierr);  

  if (periodicBiogeochemForcing) {    
    ierr = destroyPeriodicVec(&Tsp);CHKERRQ(ierr);
    ierr = destroyPeriodicVec(&Ssp);CHKERRQ(ierr);
    ierr = destroyPeriodicVec(&Fep);CHKERRQ(ierr);    
    ierr = destroyPeriodicArray(&localaicep);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localhicep);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localhsnop);CHKERRQ(ierr);    
    ierr = destroyPeriodicArray(&localwindp);CHKERRQ(ierr);    
    ierr = destroyPeriodicArray(&localatmospp);CHKERRQ(ierr);    
#ifdef READ_SWRAD    
    ierr = destroyPeriodicArray(&localswradp);CHKERRQ(ierr);
#endif    
  }    

#ifdef CARBON
  if (useAtmModel) {
    ierr = PetscFClose(PETSC_COMM_WORLD,atmfptime);CHKERRQ(ierr);
  }
#endif

  return 0;
}





/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/



#undef __FUNCT__
#define __FUNCT__ "reInitializeExternalForcing"
PetscErrorCode reInitializeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt numTracers, Vec *v, Vec *ut)
{
  PetscErrorCode ierr;
  PetscInt ip, kl, nzloc;
  PetscScalar myTime;
  PetscViewer fd;
  PetscInt fp;

  myTime = DeltaT*Iter; /* Iter should start at 0 */

  if (periodicBiogeochemForcing) {   
	ierr = interpPeriodicVector(tc,&Ts,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Ssp,"Ss_");	
	ierr = interpPeriodicVector(tc,&Fe,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Fep,"Fe_");	
#ifdef READ_SWRAD
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localswradp,"swrad_");
#else
    insolation_(&lNumProfiles,&myTime,&locallatitude[0],&localswrad[0],&localtau[0]);
#endif                                                 
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localaice,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localaicep,"aice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhice,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localhicep,"hice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhsno,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localhsnop,"hsno_");                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemCyclePeriod,numBiogeochemPeriods,
					                              tdpBiogeochem,&localwindp,"wind_");   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemCyclePeriod,numBiogeochemPeriods,
									              tdpBiogeochem,&localatmospp,"atmosp_");	
    if (useEmP) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                    tdpBiogeochem,&localEmPp,"EmP_");                                                      
    }
  }
    
  return 0;
}
