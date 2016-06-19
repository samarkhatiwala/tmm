#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#undef READ_SWRAD
#include "petscmat.h"
#include "petsc_matvec_utils.h"
#include "tmm_main.h"
#include "tmm_forcing_utils.h"
#include "tmm_profile_utils.h"
#include "tmm_profile_data.h"
#include "kiel_biogeochem.h"

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

Vec Ts,Ss;
PetscInt *gIndices;
PetscScalar *localTs,*localSs;
PetscScalar *localTR1,*localTR2,*localTR3,*localTR4,*localTR5,*localTR6;
PetscScalar *localJTR1,*localJTR2,*localJTR3,*localJTR4,*localJTR5,*localJTR6;
#ifdef CARBON
PetscScalar *localTR7,*localTR8;
PetscScalar *localJTR7,*localJTR8;
PetscScalar *localph;
PetscBool useEmP = PETSC_FALSE;
PetscScalar *localEmP;
PeriodicArray localEmPp;
Vec surfVolFrac;
#endif
PetscScalar *localwind,*localfice,*localdz,*localatmosp;
PetscScalar *localswrad, *localtau;
#ifndef READ_SWRAD
PetscScalar *locallatitude;
#endif

PetscBool useSeparateBiogeochemTimeStepping = PETSC_FALSE;
PetscInt numBiogeochemStepsPerOceanStep = 1;
PetscInt nzmax,nzeuph;
PetscScalar DeltaT,TheoDeltaT;
PetscScalar *drF;

PeriodicVec Tsp, Ssp;
PeriodicArray localwindp,localficep,localatmospp;
#ifdef READ_SWRAD
PeriodicArray localswradp;
#endif

PetscInt numBiogeochemPeriods;
PetscScalar *tdpBiogeochem; /* arrays for periodic forcing */
PetscBool periodicBiogeochemForcing = PETSC_FALSE;
PetscScalar biogeochemCyclePeriod, biogeochemCycleStep;

PetscBool readBGCParams = PETSC_FALSE;
PetscInt numBGCParams = 0;
char bgcParamsFile[PETSC_MAX_PATH_LEN];
PetscScalar *bgcparams;

PetscScalar *localdA;

#ifdef CARBON
/* atmospheric model variables */
PetscScalar *TpCO2atm_hist, *pCO2atm_hist;
PetscInt numpCO2atm_hist = 0;
PetscBool fixedAtmosCO2 = PETSC_TRUE;
char pCO2atmIniFile[PETSC_MAX_PATH_LEN];  

PetscBool useAtmModel = PETSC_FALSE;
PetscScalar pCO2atm_ini = 280.0; /* default initial value */
PetscScalar pCO2atm = 280.0; /* default initial value */
PetscScalar ppmToPgC=2.1324;
PetscScalar atmModelDeltaT;
PetscScalar secPerYear=86400.0*360.0;
PetscScalar Focean=0.0;
PetscScalar localFocean=0.0;
PetscScalar Foceanint = 0.0;
PetscInt atmModelUpdateTimeSteps=1;

PetscInt atmWriteSteps;
PetscBool atmAppendOutput;
FILE *atmfptime;
PetscViewer atmfd;
PetscInt atmfp;
char atmOutTimeFile[PETSC_MAX_PATH_LEN];  
#endif

#ifdef SEDIMENT
PetscScalar runoff_ini = 0.0;
PetscScalar GRunoff; /* Global runoff, calculated from burial */
PetscScalar *localrunoffvol; /* volume supplied by runoff */
PetscScalar localFburial = 0.0;
PetscScalar Fburial=0.0;
PetscInt burialSumSteps;
char runoffOutTimeFile[PETSC_MAX_PATH_LEN];  
char runoffIniFile[PETSC_MAX_PATH_LEN];  
FILE *runofffptime;
#endif

PetscBool calcDiagnostics = PETSC_FALSE;
PetscInt diagNumTimeSteps, diagStartTimeStep, diagCount;
PetscBool appendDiagnostics = PETSC_FALSE;
/* Add model specific diagnostic variables below */
Vec fbgc1, fbgc2, fbgc3, fbgc4, fbgc5, fbgc1avg, fbgc2avg, fbgc3avg, fbgc4avg, fbgc5avg;
PetscViewer fdfbgc1avg, fdfbgc2avg, fdfbgc3avg, fdfbgc4avg, fdfbgc5avg;
PetscScalar *localfbgc1, *localfbgc2, *localfbgc3, *localfbgc4, *localfbgc5;

#ifdef CARBON
PetscScalar *localco2airseafluxdiag, *localco2airseafluxdiagavg;
#endif

PetscBool TRUEFLAG = PETSC_TRUE, FALSEFLAG = PETSC_FALSE;

#if defined (FORSPINUP) || defined (FORJACOBIAN)
PetscScalar relaxTau[50], relaxLambda[50], relaxValue[50];
PetscBool relaxTracer = PETSC_FALSE;
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
  PetscBool flg;
  PetscInt it;
  PetscScalar myTime;
  PetscScalar zero = 0.0;
  
  PetscScalar DaysPerYear = 360.0;

#if defined (FORSPINUP) || defined (FORJACOBIAN)
  ierr = PetscOptionsHasName(PETSC_NULL,"-relax_tracer",&relaxTracer);CHKERRQ(ierr);
  if (relaxTracer) {  
    PetscInt maxValsToRead, itr;

    maxValsToRead = numTracers;
    ierr = PetscOptionsGetRealArray(PETSC_NULL,"-relax_tau",relaxTau,&maxValsToRead,&flg);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate tracer relaxation tau with the -relax_tau option");
    if (maxValsToRead != numTracers) {
      SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of relaxation tau values specified");
    }

    maxValsToRead = numTracers;
    ierr = PetscOptionsGetRealArray(PETSC_NULL,"-relax_value",relaxValue,&maxValsToRead,&flg);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate relaxation values with the -relax_value option");
    if (maxValsToRead != numTracers) {
      SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of relaxation values specified");
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
  ierr = VecGetArray(TR1,&localTR1);CHKERRQ(ierr);
  ierr = VecGetArray(TR2,&localTR2);CHKERRQ(ierr);
  ierr = VecGetArray(TR3,&localTR3);CHKERRQ(ierr);
  ierr = VecGetArray(TR4,&localTR4);CHKERRQ(ierr);
  ierr = VecGetArray(TR5,&localTR5);CHKERRQ(ierr);
  ierr = VecGetArray(TR6,&localTR6);CHKERRQ(ierr);
#ifdef CARBON
  ierr = VecGetArray(TR7,&localTR7);CHKERRQ(ierr);
  ierr = VecGetArray(TR8,&localTR8);CHKERRQ(ierr);
#endif  

  for (itr=0; itr<numTracers; itr++) {
    ierr = VecSet(ut[itr],zero); CHKERRQ(ierr);
  }
  ierr = VecGetArray(JTR1,&localJTR1);CHKERRQ(ierr);
  ierr = VecGetArray(JTR2,&localJTR2);CHKERRQ(ierr);
  ierr = VecGetArray(JTR3,&localJTR3);CHKERRQ(ierr);
  ierr = VecGetArray(JTR4,&localJTR4);CHKERRQ(ierr);
  ierr = VecGetArray(JTR5,&localJTR5);CHKERRQ(ierr);
  ierr = VecGetArray(JTR6,&localJTR6);CHKERRQ(ierr);
#ifdef CARBON
  ierr = VecGetArray(JTR7,&localJTR7);CHKERRQ(ierr);
  ierr = VecGetArray(JTR8,&localJTR8);CHKERRQ(ierr);
#endif
  ierr = PetscOptionsGetBool(PETSC_NULL,"-separate_biogeochem_time_stepping",&useSeparateBiogeochemTimeStepping,0);CHKERRQ(ierr);
#if defined (FORSPINUP) || defined (FORJACOBIAN)
  if (useSeparateBiogeochemTimeStepping) {
    SETERRQ(PETSC_COMM_WORLD,1,"Cannot use the -separate_biogeochem_time_stepping option with SPINUP or JACOBIAN ");  
  
  }
#endif
  if (useSeparateBiogeochemTimeStepping) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Biogeochem model will be time-stepped independently\n");CHKERRQ(ierr);
  }  
  
  ierr = PetscOptionsGetInt(PETSC_NULL,"-num_biogeochem_steps_per_ocean_step",&numBiogeochemStepsPerOceanStep,&flg);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of biogeochem model time steps per ocean time step = %d\n",numBiogeochemStepsPerOceanStep);CHKERRQ(ierr);

  ierr = PetscOptionsGetInt(PETSC_NULL,"-nzeuph",&nzeuph,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate number of euphotic zone layers with the -nzeuph option");  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of euphotic zone layers is %d \n",nzeuph);CHKERRQ(ierr);

  ierr = PetscOptionsGetReal(PETSC_NULL,"-biogeochem_deltat",&DeltaT,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate biogeochemical time step in seconds with the -biogeochem_deltat option");  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Ocean time step for BGC length is  %12.7f seconds\n",DeltaT);CHKERRQ(ierr);

  TheoDeltaT = DaysPerYear*86400.0*deltaTClock;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Check: using a year length of %12.3f days \n",DaysPerYear);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Theoretical ocean time step length for BGC is then  %12.7f seconds\n",TheoDeltaT);CHKERRQ(ierr);

  ierr = PetscOptionsHasName(PETSC_NULL,"-periodic_biogeochem_forcing",&periodicBiogeochemForcing);CHKERRQ(ierr);

  if (periodicBiogeochemForcing) {    
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic biogeochemical forcing specified\n");CHKERRQ(ierr);

/*  read time data */
/*  IMPORTANT: time units must be the same as that used by the toplevel driver */
    ierr = PetscOptionsGetReal(PETSC_NULL,"-periodic_biogeochem_cycle_period",&biogeochemCyclePeriod,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate biogeochemical forcing cycling time with the -periodic_biogeochem_cycle_period option");
    ierr = PetscOptionsGetReal(PETSC_NULL,"-periodic_biogeochem_cycle_step",&biogeochemCycleStep,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate biogeochemical forcing cycling step with the -periodic_biogeochem_cycle_step option");
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
	ierr = VecLoad(Ts,fd);CHKERRQ(ierr);  /* IntoVector */ 
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Ss.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(Ss,fd);CHKERRQ(ierr);    /* IntoVector */ 
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
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

/* in principle, localdA is only used by the atmospheric exchange - but I might also need it for the sediment,
and who knows what comes next. IK. */
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localdA);CHKERRQ(ierr);
    ierr = readProfileSurfaceScalarData("dA.bin",localdA,1);  

#ifdef CARBON
  ierr = PetscOptionsHasName(PETSC_NULL,"-use_atm_model",&useAtmModel);CHKERRQ(ierr);

  if (useAtmModel) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using interactive atmospheric model\n");CHKERRQ(ierr);  

/* overwrite default value */
	ierr = PetscOptionsGetReal(PETSC_NULL,"-pco2atm_ini",&pCO2atm_ini,&flg);CHKERRQ(ierr); /* read from command line */
    if (!flg) {
      ierr = PetscOptionsGetString(PETSC_NULL,"-pco2atm_ini_file",pCO2atmIniFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
      if (flg) { /* read from binary file */
        ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,pCO2atmIniFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
        ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
        ierr = PetscBinaryRead(fp,&pCO2atm_ini,1,PETSC_SCALAR);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      }
    }
    pCO2atm = pCO2atm_ini;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using initial atmospheric pCO2 of %g ppm\n",pCO2atm);CHKERRQ(ierr);
      
    ierr = PetscOptionsGetInt(PETSC_NULL,"-atm_write_steps",&atmWriteSteps,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate atmospheric model output step with the -atm_write_steps option");

    ierr = PetscOptionsHasName(PETSC_NULL,"-atm_append",&atmAppendOutput);CHKERRQ(ierr);
    if (atmAppendOutput) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output will be appended\n");CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output will overwrite existing file(s)\n");CHKERRQ(ierr);
    }    

/* Output times */
    ierr = PetscOptionsGetString(PETSC_NULL,"-atm_time_file",atmOutTimeFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
    if (!flg) {
      strcpy(atmOutTimeFile,"");
      sprintf(atmOutTimeFile,"%s","atm_output_time.txt");
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output times will be written to %s\n",atmOutTimeFile);CHKERRQ(ierr);

    if (!atmAppendOutput) {
      ierr = PetscFOpen(PETSC_COMM_WORLD,atmOutTimeFile,"w",&atmfptime);CHKERRQ(ierr);  
      ierr = PetscFPrintf(PETSC_COMM_WORLD,atmfptime,"%d   %10.5f\n",Iter0,time0);CHKERRQ(ierr);     
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing atmospheric output at time %10.5f, step %d\n", tc,Iter);CHKERRQ(ierr);  
      ierr = writeBinaryScalarData("pCO2atm_output.bin",&pCO2atm,1,PETSC_FALSE);
      ierr = writeBinaryScalarData("Foceanint_output.bin",&Focean,1,PETSC_FALSE);
    } else {
      ierr = PetscFOpen(PETSC_COMM_WORLD,atmOutTimeFile,"a",&atmfptime);CHKERRQ(ierr);  
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output will be appended. Initial condition will NOT be written\n");CHKERRQ(ierr);      
    }

    ierr = PetscOptionsGetInt(PETSC_NULL,"-atm_update_steps",&atmModelUpdateTimeSteps,&flg);CHKERRQ(ierr);
    if ((maxSteps % atmModelUpdateTimeSteps)!=0) {
      SETERRQ(PETSC_COMM_WORLD,1,"maxSteps not divisible by atmModelUpdateTimeSteps!");
    }
    if ((atmWriteSteps % atmModelUpdateTimeSteps)!=0) {
      SETERRQ(PETSC_COMM_WORLD,1,"atmWriteSteps not divisible by atmModelUpdateTimeSteps!");
    }    
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model will be updated every %d time steps\n",atmModelUpdateTimeSteps);CHKERRQ(ierr);
    if (atmModelUpdateTimeSteps>1) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"WARNING: Focean and pCO2atm diagnostics may not be correct!\n");CHKERRQ(ierr);
    }

    atmModelDeltaT = atmModelUpdateTimeSteps*DeltaT/secPerYear; /* time step in years */    

  } else {  /* not using atm model */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using prescribed atmospheric pCO2\n");CHKERRQ(ierr);
  
    ierr = PetscOptionsGetInt(PETSC_NULL,"-pco2_num_hist",&numpCO2atm_hist,&flg);CHKERRQ(ierr);
    if (flg) { /* Read atmospheric pCO2 history */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading time-dependent atmospheric pCO2 history\n");CHKERRQ(ierr);
    
      fixedAtmosCO2 = PETSC_FALSE;
      ierr = PetscMalloc(numpCO2atm_hist*sizeof(PetscScalar),&TpCO2atm_hist);CHKERRQ(ierr); 
      ierr = PetscMalloc(numpCO2atm_hist*sizeof(PetscScalar),&pCO2atm_hist);CHKERRQ(ierr); 
  
      ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"TpCO2.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
      ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
      ierr = PetscBinaryRead(fp,TpCO2atm_hist,numpCO2atm_hist,PETSC_SCALAR);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  
      ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"pCO2atm.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
      ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
      ierr = PetscBinaryRead(fp,pCO2atm_hist,numpCO2atm_hist,PETSC_SCALAR);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      
      pCO2atm = pCO2atm_hist[0];

    }	else {
      ierr = PetscOptionsGetReal(PETSC_NULL,"-pco2_atm",&pCO2atm,&flg);CHKERRQ(ierr); /* overwrite default value */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Using fixed atmospheric pCO2 of %g ppm\n",pCO2atm);CHKERRQ(ierr);
      
    }    
  }
  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localph);CHKERRQ(ierr);  

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

	ierr = VecDuplicate(TR7,&surfVolFrac);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"surface_volume_fraction.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(surfVolFrac,fd);CHKERRQ(ierr);  /* IntoVector */ 
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
  }

#endif

#ifdef SEDIMENT

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using Burial-Runoff model\n");CHKERRQ(ierr);  

/* Define the interval over which to integrate global burial */
  ierr = PetscOptionsGetInt(PETSC_NULL,"-burial_sum_steps",&burialSumSteps,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate burial integration interval with the -burial_sum_steps option");
  if ((maxSteps % burialSumSteps)!=0) {
    SETERRQ(PETSC_COMM_WORLD,1,"maxSteps not divisible by burialSumSteps!");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Runoff be integrated over and written every %d time steps\n",burialSumSteps);CHKERRQ(ierr);

/* set the name of the runoff time file */
  ierr = PetscOptionsGetString(PETSC_NULL,"-runoff_time_file",runoffOutTimeFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) {
  strcpy(runoffOutTimeFile,"");
  sprintf(runoffOutTimeFile,"%s","runoff_output_time.txt");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Runoff output times will be written to %s\n",runoffOutTimeFile);CHKERRQ(ierr);


/* set inititial runoff: overwrite default value with value from command line*/
  ierr = PetscOptionsGetReal(PETSC_NULL,"-runoff_ini",&runoff_ini,&flg);CHKERRQ(ierr);
/* set inititial runoff: overwrite default value with value from file*/
  if (!flg) {
    ierr = PetscOptionsGetString(PETSC_NULL,"-runoff_ini_file",runoffIniFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
    if (flg) { /* read from binary file */
      ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,runoffIniFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
      ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
      ierr = PetscBinaryRead(fp,&runoff_ini,1,PETSC_SCALAR);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
    }
  }
    
  GRunoff = runoff_ini;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using initial runoff of %g Gmol P/d\n",GRunoff);CHKERRQ(ierr);

/* if run is continued, always append runoff and output times */
  if (Iter0>0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Runoff output will be appended\n");CHKERRQ(ierr);
      ierr = PetscFOpen(PETSC_COMM_WORLD,runoffOutTimeFile,"a",&runofffptime);CHKERRQ(ierr);  
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Initial runoff output will not be written\n");CHKERRQ(ierr);
   } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Runoff output will overwrite existing file(s)\n");CHKERRQ(ierr);
      ierr = PetscFOpen(PETSC_COMM_WORLD,runoffOutTimeFile,"w",&runofffptime);CHKERRQ(ierr);  
      ierr = PetscFPrintf(PETSC_COMM_WORLD,runofffptime,"%d   %10.5f\n",Iter0,time0);CHKERRQ(ierr);     
      ierr = writeBinaryScalarData("Grunoff_output.bin",&GRunoff,1,PETSC_FALSE);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing runoff output at time %10.5f, step %d\n", tc,Iter);CHKERRQ(ierr);  
   }

/* fraction of global river runoff in each box, divided by the box volume (a 3D field) */
/* Note: VecLoadVecIntoArray resides in petsc_matvec_utils.c and is not a generic petsc function*/
  ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localrunoffvol);CHKERRQ(ierr);    
  ierr = VecLoadVecIntoArray(TR1,"runoff_volume_annual.petsc",localrunoffvol);CHKERRQ(ierr);

#endif

/* Grid arrays */
  ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localdz);CHKERRQ(ierr);    
  ierr = VecLoadVecIntoArray(TR1,"dz.petsc",localdz);CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"drF.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fp,&nzmax,1,PETSC_INT);CHKERRQ(ierr);  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of vertical layers is %d \n",nzmax);CHKERRQ(ierr);  
  ierr = PetscMalloc(nzmax*sizeof(PetscScalar),&drF);CHKERRQ(ierr); 
  ierr = PetscBinaryRead(fp,drF,nzmax,PETSC_SCALAR);CHKERRQ(ierr);  
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

/* Forcing fields */  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localswrad);CHKERRQ(ierr);  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localtau);CHKERRQ(ierr);  
#ifdef READ_SWRAD
  if (periodicBiogeochemForcing) {    
    localswradp.firstTime = PETSC_TRUE;
    localswradp.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("swrad.bin",localswrad,1);  
  }
#else
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&locallatitude);CHKERRQ(ierr);  
  ierr = readProfileSurfaceScalarData("latitude.bin",locallatitude,1);  
#endif	  
  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localfice);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localficep.firstTime = PETSC_TRUE;
    localficep.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("fice.bin",localfice,1);  
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
#ifdef READ_SWRAD
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localswradp,"swrad_");
#else
   insolation_(&lNumProfiles,&myTime,&locallatitude[0],&localswrad[0],&localtau[0]);
#endif                                                 
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localficep,"fice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemCyclePeriod,numBiogeochemPeriods,
					                              tdpBiogeochem,&localwindp,"wind_");   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemCyclePeriod,numBiogeochemPeriods,
									              tdpBiogeochem,&localatmospp,"atmosp_");					                              
#ifdef CARBON
    if (useEmP) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                    tdpBiogeochem,&localEmPp,"EmP_");                                                  
    }
#endif									              
  } else {
#ifndef READ_SWRAD
    insolation_(&lNumProfiles,&myTime,&locallatitude[0],&localswrad[0],&localtau[0]);
#endif    
  }
  
  for (ip=0; ip<lNumProfiles; ip++) {
    nzloc=lProfileLength[ip];
    kl=lStartIndices[ip];
    kiel_biogeochem_ini_(&nzloc,&DeltaT,
                         &localTR1[kl],&localTR2[kl],&localTR3[kl],
			 &localTR4[kl],&localTR5[kl],&localTR6[kl],
#ifdef CARBON			   
			 &localTR7[kl],&localTR8[kl],&localph[ip],
#endif
			 &localTs[kl],&localSs[kl],&localdz[kl],&drF[0],&nzmax,&nzeuph,
                         &numBiogeochemStepsPerOceanStep,&TRUEFLAG);
  }

/* Read and overwrite default parameter values here */
  ierr = PetscOptionsGetString(PETSC_NULL,"-bgc_params_file",bgcParamsFile,PETSC_MAX_PATH_LEN-1,&readBGCParams);CHKERRQ(ierr);
  if (readBGCParams) {
    ierr = PetscOptionsGetInt(PETSC_NULL,"-num_bgc_params",&numBGCParams,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate number of BGC parameters to read with the -num_bgc_params option");    

    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,bgcParamsFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
    ierr = PetscMalloc(numBGCParams*sizeof(PetscScalar),&bgcparams);CHKERRQ(ierr); 
    ierr = PetscBinaryRead(fp,bgcparams,numBGCParams,PETSC_SCALAR);CHKERRQ(ierr);  
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
    
    kiel_biogeochem_set_params_(&numBGCParams,&bgcparams[0]);

    myTime = DeltaT*Iter; /* Iter should start at 0 */
    for (ip=0; ip<lNumProfiles; ip++) {
      nzloc=lProfileLength[ip];
      kl=lStartIndices[ip];
      kiel_biogeochem_ini_(&nzloc,&DeltaT,
                           &localTR1[kl],&localTR2[kl],&localTR3[kl],
			   &localTR4[kl],&localTR5[kl],&localTR6[kl],
#ifdef CARBON			   
			   &localTR7[kl],&localTR8[kl],&localph[ip],
#endif
			   &localTs[kl],&localSs[kl],&localdz[kl],&drF[0],&nzmax,&nzeuph,
                           &numBiogeochemStepsPerOceanStep,&FALSEFLAG);
    }
    
  }
  
  ierr = PetscOptionsHasName(PETSC_NULL,"-calc_diagnostics",&calcDiagnostics);CHKERRQ(ierr);
  if (calcDiagnostics) {    
/*Data for diagnostics */
	ierr = PetscOptionsGetInt(PETSC_NULL,"-diag_start_time_step",&diagStartTimeStep,&flg);CHKERRQ(ierr);
	if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate (absolute) time step at which to start storing diagnostics with the -diag_start_time_step flag");
	ierr = PetscOptionsGetInt(PETSC_NULL,"-diag_time_steps",&diagNumTimeSteps,&flg);CHKERRQ(ierr);
	if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate number of time averaging diagnostics time steps with the -diag_time_step flag");
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed starting at (and including) time step: %d\n", diagStartTimeStep);CHKERRQ(ierr);	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed over %d time steps\n", diagNumTimeSteps);CHKERRQ(ierr);	

	ierr = VecDuplicate(TR1,&fbgc1);CHKERRQ(ierr);
	ierr = VecSet(fbgc1,zero);CHKERRQ(ierr);
	ierr = VecGetArray(fbgc1,&localfbgc1);CHKERRQ(ierr);
	ierr = VecDuplicate(TR1,&fbgc1avg);CHKERRQ(ierr);
	ierr = VecSet(fbgc1avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"fbgc1.petsc",FILE_MODE_WRITE,&fdfbgc1avg);CHKERRQ(ierr);

	ierr = VecDuplicate(TR1,&fbgc2);CHKERRQ(ierr);
	ierr = VecSet(fbgc2,zero);CHKERRQ(ierr);
	ierr = VecGetArray(fbgc2,&localfbgc2);CHKERRQ(ierr);
	ierr = VecDuplicate(TR1,&fbgc2avg);CHKERRQ(ierr);
	ierr = VecSet(fbgc2avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"fbgc2.petsc",FILE_MODE_WRITE,&fdfbgc2avg);CHKERRQ(ierr);

	ierr = VecDuplicate(TR1,&fbgc3);CHKERRQ(ierr);
	ierr = VecSet(fbgc3,zero);CHKERRQ(ierr);
	ierr = VecGetArray(fbgc3,&localfbgc3);CHKERRQ(ierr);
	ierr = VecDuplicate(TR1,&fbgc3avg);CHKERRQ(ierr);
	ierr = VecSet(fbgc3avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"fbgc3.petsc",FILE_MODE_WRITE,&fdfbgc3avg);CHKERRQ(ierr);

	ierr = VecDuplicate(TR1,&fbgc4);CHKERRQ(ierr);
	ierr = VecSet(fbgc4,zero);CHKERRQ(ierr);
	ierr = VecGetArray(fbgc4,&localfbgc4);CHKERRQ(ierr);
	ierr = VecDuplicate(TR1,&fbgc4avg);CHKERRQ(ierr);
	ierr = VecSet(fbgc4avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"fbgc4.petsc",FILE_MODE_WRITE,&fdfbgc4avg);CHKERRQ(ierr);

	ierr = VecDuplicate(TR1,&fbgc5);CHKERRQ(ierr);
	ierr = VecSet(fbgc5,zero);CHKERRQ(ierr);
	ierr = VecGetArray(fbgc5,&localfbgc5);CHKERRQ(ierr);
	ierr = VecDuplicate(TR1,&fbgc5avg);CHKERRQ(ierr);
	ierr = VecSet(fbgc5avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"fbgc5.petsc",FILE_MODE_WRITE,&fdfbgc5avg);CHKERRQ(ierr);

#ifdef CARBON
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localco2airseafluxdiag);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localco2airseafluxdiagavg);CHKERRQ(ierr);  

    for (ip=0; ip<lNumProfiles; ip++) {
      localco2airseafluxdiag[ip]=0.0;
      localco2airseafluxdiagavg[ip]=0.0;
    }    
#endif

	diagCount=0;
  }
  
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

#ifdef CARBON
  PetscInt itf;
  PetscScalar alpha;
  PetscScalar DICglobavg = 0.0;
  PetscScalar localco2airseaflux = 0.0;
#endif
  
#ifdef SEDIMENT
  PetscScalar localburial = 0.0;
#endif
  
  myTime = DeltaT*Iter; /* Iter should start at 0 */

  if (periodicBiogeochemForcing) {   
	ierr = interpPeriodicVector(tc,&Ts,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Ssp,"Ss_");	
#ifdef READ_SWRAD
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localswradp,"swrad_");
#else
    insolation_(&lNumProfiles,&myTime,&locallatitude[0],&localswrad[0],&localtau[0]);
#endif                                                 
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localficep,"fice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemCyclePeriod,numBiogeochemPeriods,
					                              tdpBiogeochem,&localwindp,"wind_");  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemCyclePeriod,numBiogeochemPeriods,
									              tdpBiogeochem,&localatmospp,"atmosp_");					                              
#ifdef CARBON
    if (useEmP) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                    tdpBiogeochem,&localEmPp,"EmP_");                                                      
    }
#endif
  }

#ifdef CARBON
  if (useAtmModel) {
  } else {  
/* Interpolate atmospheric pCO2   */
    if (!fixedAtmosCO2) { 
      if (tc>=TpCO2atm_hist[0]) {
        ierr = calcInterpFactor(numpCO2atm_hist,tc,TpCO2atm_hist,&itf,&alpha); CHKERRQ(ierr);
        pCO2atm = alpha*pCO2atm_hist[itf] + (1.0-alpha)*pCO2atm_hist[itf+1];	  
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: time < %10.5f. Assuming pCO2atm=%g\n",TpCO2atm_hist[0],pCO2atm);CHKERRQ(ierr);
      }
    }  
  }

  if (useEmP) {
    ierr = VecDot(surfVolFrac,TR7,&DICglobavg);CHKERRQ(ierr); /* volume weighted mean surface DIC */									              
  }

#endif
  
  for (ip=0; ip<lNumProfiles; ip++) {
    nzloc=lProfileLength[ip];
    kl=lStartIndices[ip];
    kiel_biogeochem_model_(&nzloc,&DeltaT,
                           &localTR1[kl],&localTR2[kl],&localTR3[kl],
			   &localTR4[kl],&localTR5[kl],&localTR6[kl],
#ifdef CARBON			   
			   &localTR7[kl],&localTR8[kl],&DICglobavg,&localEmP[ip],&pCO2atm,
#endif
			   &localTs[kl],&localSs[kl],&localfice[ip],&localswrad[ip],&localtau[ip],&localwind[ip],&localatmosp[ip],&localdz[kl],
                           &localJTR1[kl],&localJTR2[kl],&localJTR3[kl],
			   &localJTR4[kl],&localJTR5[kl],&localJTR6[kl],
#ifdef CARBON			   
			   &localJTR7[kl],&localJTR8[kl],&localph[ip],&localco2airseaflux,
#endif
#ifdef SEDIMENT
                           &localburial,&GRunoff,&localrunoffvol[kl],
#endif			    
			   &useSeparateBiogeochemTimeStepping);                            

#ifdef CARBON			   
    if (useAtmModel) {
      localFocean = localFocean + (localco2airseaflux/DeltaT)*localdA[ip]*(12.0/1.e18)*secPerYear; /* PgC/y */
    }
#endif

#ifdef SEDIMENT
/* integrate burial in sediment over area and all profiles on each processor */
      localFburial = localFburial + localburial*localdA[ip]; 
#endif

	if (calcDiagnostics) {  
	  if (Iter0+iLoop>=diagStartTimeStep) { /* start time averaging (note: diagStartTimeStep is ABSOLUTE time step) */	
        kiel_biogeochem_diagnostics_(&nzloc,&localfbgc1[kl],&localfbgc2[kl],&localfbgc3[kl],&localfbgc4[kl],&localfbgc5[kl]);
#ifdef CARBON        
        localco2airseafluxdiag[ip]=localco2airseaflux;
#endif                
      }
	}
  } /* end loop over profiles */

#ifdef CARBON

  if (useAtmModel) {
    if ((iLoop % atmModelUpdateTimeSteps)==0) {  /*  time to update atmosphere */
  
      Focean = 0.0;
       
      MPI_Reduce(&localFocean, &Focean, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&Focean, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
    
/*    time step atmosphere */
      Focean = Focean/atmModelUpdateTimeSteps; /* average flux over accumulation period Pg C/yr*/
      pCO2atm = pCO2atm + atmModelDeltaT*(-Focean)/ppmToPgC;

/*    reset values */
      localFocean = 0.0; 
      Foceanint = Foceanint + atmModelDeltaT*Focean; /* calculate the time integrated flux */
        
/*      Focean = 0.0; */     
     
    }
  }  

#endif  

#ifdef SEDIMENT

/* sum burial in sediment over all processors, and scale by time step etc.*/
/* do this only once every burialSumSteps , and then take this value for next year's runoff */

    if ((iLoop % burialSumSteps)==0) {

      Fburial = 0.0;

      MPI_Reduce(&localFburial, &Fburial, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD);
      MPI_Bcast(&Fburial, 1, MPI_DOUBLE, 0, PETSC_COMM_WORLD);
      
      GRunoff = Fburial/(1.e12*burialSumSteps)*(86400.0/DeltaT); /* This is Gmol P/day. 
      Note: localrunoff is scaled with 1e12. Note: GRunoff will be scaled with bgc_dt.*/

      localFburial = 0.0;
      }
    

#endif


  if (useSeparateBiogeochemTimeStepping) {  /* return updated tracer field */
	ierr = VecSetValues(TR1,lSize,gIndices,localTR1,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(TR1);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(TR1);CHKERRQ(ierr);    
	ierr = VecSetValues(TR2,lSize,gIndices,localTR2,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(TR2);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(TR2);CHKERRQ(ierr);    
	ierr = VecSetValues(TR3,lSize,gIndices,localTR3,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(TR3);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(TR3);CHKERRQ(ierr);    
	ierr = VecSetValues(TR4,lSize,gIndices,localTR4,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(TR4);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(TR4);CHKERRQ(ierr);    
	ierr = VecSetValues(TR5,lSize,gIndices,localTR5,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(TR5);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(TR5);CHKERRQ(ierr);    
	ierr = VecSetValues(TR6,lSize,gIndices,localTR6,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(TR6);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(TR6);CHKERRQ(ierr); 
#ifdef CARBON
	ierr = VecSetValues(TR7,lSize,gIndices,localTR7,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(TR7);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(TR7);CHKERRQ(ierr);    
	ierr = VecSetValues(TR8,lSize,gIndices,localTR8,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(TR8);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(TR8);CHKERRQ(ierr);    
#endif	   
  } else {  /* return tracer tendency */
	ierr = VecSetValues(JTR1,lSize,gIndices,localJTR1,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(JTR1);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(JTR1);CHKERRQ(ierr);    
	ierr = VecSetValues(JTR2,lSize,gIndices,localJTR2,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(JTR2);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(JTR2);CHKERRQ(ierr);    
	ierr = VecSetValues(JTR3,lSize,gIndices,localJTR3,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(JTR3);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(JTR3);CHKERRQ(ierr);    
	ierr = VecSetValues(JTR4,lSize,gIndices,localJTR4,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(JTR4);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(JTR4);CHKERRQ(ierr);    
	ierr = VecSetValues(JTR5,lSize,gIndices,localJTR5,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(JTR5);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(JTR5);CHKERRQ(ierr);    
	ierr = VecSetValues(JTR6,lSize,gIndices,localJTR6,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(JTR6);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(JTR6);CHKERRQ(ierr);    
#ifdef CARBON
	ierr = VecSetValues(JTR7,lSize,gIndices,localJTR7,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(JTR7);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(JTR7);CHKERRQ(ierr);    
	ierr = VecSetValues(JTR8,lSize,gIndices,localJTR8,INSERT_VALUES);CHKERRQ(ierr);
	ierr = VecAssemblyBegin(JTR8);CHKERRQ(ierr);
	ierr = VecAssemblyEnd(JTR8);CHKERRQ(ierr);    
#endif	   

#if defined (FORSPINUP) || defined (FORJACOBIAN)
/* add relaxation term: ut = ut - lambda*(v-vr) = ut -lambda*v + lambda*vr */
    if (relaxTracer) {
      for (itr=0; itr<numTracers; itr++) {
        ierr = VecAXPY(ut[itr],-relaxLambda[itr],v[itr]);CHKERRQ(ierr); /* ut = ut - lambda*v */
        ierr = VecShift(ut[itr],relaxLambda[itr]*relaxValue[itr]);CHKERRQ(ierr); /* ut = ut + lambda*vr */
      }
    }
#endif
  
/*  Convert to discrete tendency */
	for (itr=0; itr<numTracers; itr++) {
	  ierr = VecScale(ut[itr],DeltaT);CHKERRQ(ierr);
	}
  }
  
  if (calcDiagnostics) {  
	if (Iter0+iLoop>=diagStartTimeStep) { /* start time averaging (note: diagStartTimeStep is ABSOLUTE time step) */  
	  ierr = VecSetValues(fbgc1,lSize,gIndices,localfbgc1,INSERT_VALUES);CHKERRQ(ierr);
	  ierr = VecAssemblyBegin(fbgc1);CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(fbgc1);CHKERRQ(ierr);    

	  ierr = VecSetValues(fbgc2,lSize,gIndices,localfbgc2,INSERT_VALUES);CHKERRQ(ierr);
	  ierr = VecAssemblyBegin(fbgc2);CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(fbgc2);CHKERRQ(ierr);    

	  ierr = VecSetValues(fbgc3,lSize,gIndices,localfbgc3,INSERT_VALUES);CHKERRQ(ierr);
	  ierr = VecAssemblyBegin(fbgc3);CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(fbgc3);CHKERRQ(ierr);    
  
	  ierr = VecSetValues(fbgc4,lSize,gIndices,localfbgc4,INSERT_VALUES);CHKERRQ(ierr);
	  ierr = VecAssemblyBegin(fbgc4);CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(fbgc4);CHKERRQ(ierr);    

	  ierr = VecSetValues(fbgc5,lSize,gIndices,localfbgc5,INSERT_VALUES);CHKERRQ(ierr);
	  ierr = VecAssemblyBegin(fbgc5);CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(fbgc5);CHKERRQ(ierr);      
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

#ifdef SEDIMENT

    if ((iLoop % burialSumSteps)==0) {  /*  time to write out */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing runoff output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
      ierr = PetscFPrintf(PETSC_COMM_WORLD,runofffptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);           
      ierr = writeBinaryScalarData("Grunoff_output.bin",&GRunoff,1,PETSC_TRUE);
    }

#endif

  if (calcDiagnostics) {  
	if (Iter0+iLoop>=diagStartTimeStep) { /* start time averaging (note: diagStartTimeStep is ABSOLUTE time step) */  
  
	  if (diagCount<=diagNumTimeSteps) { /* still within same averaging block so accumulate */
		ierr = VecAXPY(fbgc1avg,one,fbgc1);CHKERRQ(ierr);
		ierr = VecAXPY(fbgc2avg,one,fbgc2);CHKERRQ(ierr);
		ierr = VecAXPY(fbgc3avg,one,fbgc3);CHKERRQ(ierr);
		ierr = VecAXPY(fbgc4avg,one,fbgc4);CHKERRQ(ierr);
		ierr = VecAXPY(fbgc5avg,one,fbgc5);CHKERRQ(ierr);

#ifdef CARBON
        for (ip=0; ip<lNumProfiles; ip++) {
          localco2airseafluxdiagavg[ip]=localco2airseafluxdiag[ip]+localco2airseafluxdiagavg[ip];      
        }	  
#endif

		diagCount = diagCount+1;
	  }
	  if (diagCount==diagNumTimeSteps) { /* time to write averages to file */
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing diagnostics time average at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);                      

		ierr = VecScale(fbgc1avg,1.0/diagCount);CHKERRQ(ierr);
		ierr = VecView(fbgc1avg,fdfbgc1avg);CHKERRQ(ierr);
		ierr = VecSet(fbgc1avg,zero); CHKERRQ(ierr);

		ierr = VecScale(fbgc2avg,1.0/diagCount);CHKERRQ(ierr);
		ierr = VecView(fbgc2avg,fdfbgc2avg);CHKERRQ(ierr);
		ierr = VecSet(fbgc2avg,zero); CHKERRQ(ierr);

		ierr = VecScale(fbgc3avg,1.0/diagCount);CHKERRQ(ierr);
		ierr = VecView(fbgc3avg,fdfbgc3avg);CHKERRQ(ierr);
		ierr = VecSet(fbgc3avg,zero); CHKERRQ(ierr);

		ierr = VecScale(fbgc4avg,1.0/diagCount);CHKERRQ(ierr);
		ierr = VecView(fbgc4avg,fdfbgc4avg);CHKERRQ(ierr);
		ierr = VecSet(fbgc4avg,zero); CHKERRQ(ierr);

		ierr = VecScale(fbgc5avg,1.0/diagCount);CHKERRQ(ierr);
		ierr = VecView(fbgc5avg,fdfbgc5avg);CHKERRQ(ierr);
		ierr = VecSet(fbgc5avg,zero); CHKERRQ(ierr);

#ifdef CARBON
        for (ip=0; ip<lNumProfiles; ip++) {
          localco2airseafluxdiagavg[ip]=localco2airseafluxdiagavg[ip]/diagCount;
        }	  

        ierr = writeProfileSurfaceScalarData("co2airseaflux_surf.bin",localco2airseafluxdiagavg,1,appendDiagnostics);  		

/*      reset diagnostic arrays */
        for (ip=0; ip<lNumProfiles; ip++) {
          localco2airseafluxdiagavg[ip]=0.0;
        }	  
#endif

        appendDiagnostics=PETSC_TRUE;
	diagCount = 0;        

	  }
	}  
  }

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

#ifdef SEDIMENT
  ierr = writeBinaryScalarData("pickup_runoff.bin",&GRunoff,1,PETSC_FALSE);
#endif
  
  ierr = VecDestroy(&Ts);CHKERRQ(ierr);
  ierr = VecDestroy(&Ss);CHKERRQ(ierr);
  ierr = PetscFree(gIndices);CHKERRQ(ierr);  

  if (periodicBiogeochemForcing) {    
    ierr = destroyPeriodicVec(&Tsp);CHKERRQ(ierr);
    ierr = destroyPeriodicVec(&Ssp);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localficep);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localwindp);CHKERRQ(ierr);    
    ierr = destroyPeriodicArray(&localatmospp);CHKERRQ(ierr);    
#ifdef READ_SWRAD    
    ierr = destroyPeriodicArray(&localswradp);CHKERRQ(ierr);
#endif    
  }    

  if (calcDiagnostics) {  
	ierr = VecDestroy(&fbgc1);CHKERRQ(ierr);
	ierr = VecDestroy(&fbgc1avg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fdfbgc1avg);CHKERRQ(ierr);	

	ierr = VecDestroy(&fbgc2);CHKERRQ(ierr);
	ierr = VecDestroy(&fbgc2avg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fdfbgc2avg);CHKERRQ(ierr);	

	ierr = VecDestroy(&fbgc3);CHKERRQ(ierr);
	ierr = VecDestroy(&fbgc3avg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fdfbgc3avg);CHKERRQ(ierr);	

	ierr = VecDestroy(&fbgc4);CHKERRQ(ierr);
	ierr = VecDestroy(&fbgc4avg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fdfbgc4avg);CHKERRQ(ierr);	

	ierr = VecDestroy(&fbgc5);CHKERRQ(ierr);
	ierr = VecDestroy(&fbgc5avg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fdfbgc5avg);CHKERRQ(ierr);	

  }

#ifdef CARBON
  if (useAtmModel) {
    ierr = PetscFClose(PETSC_COMM_WORLD,atmfptime);CHKERRQ(ierr);
  }
#endif

#ifdef SEDIMENT
    ierr = PetscFClose(PETSC_COMM_WORLD,runofffptime);CHKERRQ(ierr);
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
#ifdef READ_SWRAD
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localswradp,"swrad_");
#else
    insolation_(&lNumProfiles,&myTime,&locallatitude[0],&localswrad[0],&localtau[0]);
#endif                                                 
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localficep,"fice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemCyclePeriod,numBiogeochemPeriods,
					                              tdpBiogeochem,&localwindp,"wind_");   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemCyclePeriod,numBiogeochemPeriods,
									              tdpBiogeochem,&localatmospp,"atmosp_");	
#ifdef CARBON
    if (useEmP) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                    tdpBiogeochem,&localEmPp,"EmP_");                                                      
    }
#endif									              
  }

  if (readBGCParams) {
    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,bgcParamsFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
    ierr = PetscMalloc(numBGCParams*sizeof(PetscScalar),&bgcparams);CHKERRQ(ierr); 
    ierr = PetscBinaryRead(fp,bgcparams,numBGCParams,PETSC_SCALAR);CHKERRQ(ierr);  
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
    
    kiel_biogeochem_set_params_(&numBGCParams,&bgcparams[0]);

    for (ip=0; ip<lNumProfiles; ip++) {
      nzloc=lProfileLength[ip];
      kl=lStartIndices[ip];
      kiel_biogeochem_ini_(&nzloc,&DeltaT,
                           &localTR1[kl],&localTR2[kl],&localTR3[kl],
			   &localTR4[kl],&localTR5[kl],&localTR6[kl],
#ifdef CARBON
			   &localTR7[kl],&localTR8[kl],&localph[ip],
#endif
			   &localTs[kl],&localSs[kl],&localdz[kl],&drF[0],&nzmax,&nzeuph,
                           &numBiogeochemStepsPerOceanStep,&FALSEFLAG);
    }    
  }
    
  return 0;
}
