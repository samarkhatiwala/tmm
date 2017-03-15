/* $Header: /Users/ikriest/CVS/mops/external_forcing_mops_biogeochem.c,v 1.3 2016/06/06 09:43:41 ikriest Exp $ */
/* $Name: mops-2_0 $*/

#define EXTERNAL_FORCING

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
#include "tmm_timer.h"
#include "mops_biogeochem.h"
#include "tmm_misfit.h"
#include "mops_biogeochem_misfit_data.h"

/* Macros to map tracer names to vectors */
/* Note: for MOPS, we have the following tracer assignement:*/
/* v[0] PO4; v[1] DOP; v[2] O2; v[3] Phy; v[4] Zoo; v[5] Det; v[6] NO3 */
/* Additionally, or option -DCARBON: v[7] DIC; v[8] Alk */
/* Note that this also affects BGC_PARAMS.h and the runsccript(s) */

#define TR v[0]
#define DIC v[7]
#define ALK v[8]
#define localDIC localTR[7]
#define localALK localTR[8]

Vec Ts,Ss;
PetscScalar *localTs,*localSs;
PetscScalar **localTR, **localJTR;
#ifdef CARBON
PetscScalar *localph;
PetscBool useVirtualFlux = PETSC_FALSE;
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

PetscBool periodicBiogeochemForcing = PETSC_FALSE;
PeriodicTimer biogeochemTimer;

PetscInt toModel = 1; 
PetscInt fromModel = 2;

PetscBool readBGCParams = PETSC_FALSE;
PetscInt numBGCParams = 0;
char bgcParamsFile[PETSC_MAX_PATH_LEN];

#ifdef ASCIIPARAMS
PetscScalar bgcparams[20];
FILE *fpparams;
#else
PetscScalar *bgcparams;
#endif

PetscScalar *localdA;

PetscInt maxValsToRead;

PetscScalar daysPerYear, secondsPerYear;

#ifdef CARBON
/* atmospheric model variables */
char *pCO2atmFiles[2];  
PetscInt numpCO2atm_hist = 0;
PetscScalar *TpCO2atm_hist, *pCO2atm_hist;
PetscBool fixedAtmosCO2 = PETSC_TRUE;
char pCO2atmIniFile[PETSC_MAX_PATH_LEN];  

PetscBool useAtmModel = PETSC_FALSE;
PetscScalar pCO2atm_ini = 280.0; /* default initial value */
PetscScalar pCO2atm = 280.0; /* default initial value */
PetscScalar ppmToPgC=2.1324;
PetscScalar atmModelDeltaT;
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


PetscScalar runoff_ini = 0.0;
PetscScalar GRunoff; /* Global runoff, calculated from burial */
PetscScalar *localrunoffvol; /* volume supplied by runoff */
PetscScalar localFburial = 0.0;
PetscScalar Fburial=0.0;
PetscInt burialSumSteps;
char runoffOutTimeFile[PETSC_MAX_PATH_LEN];  
char runoffIniFile[PETSC_MAX_PATH_LEN];  
char runoffIniOutFile[PETSC_MAX_PATH_LEN];  
FILE *runofffptime;
PetscScalar totalA = 0.0;

PetscBool calcDiagnostics = PETSC_FALSE;
StepTimer diagTimer;
PetscBool appendDiagnostics = PETSC_FALSE;
/* Add model specific diagnostic variables below */
Vec fbgc1, fbgc2, fbgc3, fbgc4, fbgc5, fbgc6, fbgc7, fbgc1avg, fbgc2avg, fbgc3avg, fbgc4avg, fbgc5avg, fbgc6avg, fbgc7avg;
PetscViewer fdfbgc1avg, fdfbgc2avg, fdfbgc3avg, fdfbgc4avg, fdfbgc5avg, fdfbgc6avg, fdfbgc7avg;
PetscScalar *localfbgc1, *localfbgc2, *localfbgc3, *localfbgc4, *localfbgc5, *localfbgc6, *localfbgc7;
char *diagOutFile[7];
PetscInt idiag;
PetscInt numDiag=7;

#ifdef CARBON
PetscScalar *localco2airseafluxdiag, *localco2airseafluxdiagavg;
#endif

PetscBool MYTRUE = PETSC_TRUE, MYFALSE = PETSC_FALSE;

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
  PetscInt ip, kl, nzloc;
  PetscInt itr;
  PetscViewer fd;
  PetscInt fp;
  PetscBool flg;
  PetscInt it,ipar;
  PetscScalar myTime;
  PetscScalar zero = 0.0;
  
#if defined (FORSPINUP) || defined (FORJACOBIAN)
  ierr = PetscOptionsHasName(PETSC_NULL,"-relax_tracer",&relaxTracer);CHKERRQ(ierr);
  if (relaxTracer) {  

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

  for (itr=0; itr<numTracers; itr++) {
    ierr = VecSet(ut[itr],zero); CHKERRQ(ierr);
  }

  ierr = VecGetArrays(v,numTracers,&localTR);CHKERRQ(ierr);
  ierr = VecGetArrays(ut,numTracers,&localJTR);CHKERRQ(ierr);
    
  ierr = PetscOptionsGetBool(PETSC_NULL,"-separate_biogeochem_time_stepping",&useSeparateBiogeochemTimeStepping,0);CHKERRQ(ierr);
#if defined (FORSPINUP) || defined (FORJACOBIAN)
  if (useSeparateBiogeochemTimeStepping) {
    SETERRQ(PETSC_COMM_WORLD,1,"Cannot use the -separate_biogeochem_time_stepping option with SPINUP or JACOBIAN ");  
  
  }
#endif
  if (useSeparateBiogeochemTimeStepping) {
    fromModel = 3;  
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

  ierr = PetscOptionsGetReal(PETSC_NULL,"-days_per_year",&daysPerYear,&flg);CHKERRQ(ierr);
  if (!flg) {
    daysPerYear = 360.0;
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of days per year is %12.7f\n",daysPerYear);CHKERRQ(ierr);
  secondsPerYear = 86400.0*daysPerYear;

  TheoDeltaT = daysPerYear*86400.0*deltaTClock;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Check: using a year length of %12.3f days \n",daysPerYear);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Theoretical ocean time step length for BGC is then  %12.7f seconds\n",TheoDeltaT);CHKERRQ(ierr);

  ierr = PetscOptionsHasName(PETSC_NULL,"-periodic_biogeochem_forcing",&periodicBiogeochemForcing);CHKERRQ(ierr);

  if (periodicBiogeochemForcing) {    
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic biogeochemical forcing specified\n");CHKERRQ(ierr);
    ierr = iniPeriodicTimer("periodic_biogeochem_", &biogeochemTimer);CHKERRQ(ierr);
  }
  
/*   Read T and S */
  ierr = VecDuplicate(TR,&Ts);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&Ss);CHKERRQ(ierr);  
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

/* Need this for atmospheric exchange, river runoff, ... */
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localdA);CHKERRQ(ierr);
  ierr = readProfileSurfaceScalarData("dA.bin",localdA,1);
  PetscScalar localA = 0.0;
  for (ip=0; ip<lNumProfiles; ip++) {
	localA = localA+localdA[ip];
  }  
  MPI_Allreduce(&localA, &totalA, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); /* global surface area */

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

    atmModelDeltaT = atmModelUpdateTimeSteps*DeltaT/secondsPerYear; /* time step in years */

  } else {  /* not using atm model */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using prescribed atmospheric pCO2\n");CHKERRQ(ierr);

    /* prescribed atmospheric CO2 */
  	maxValsToRead = 2;
	pCO2atmFiles[0] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* time file */
	pCO2atmFiles[1] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* atmospheric pCO2 history file */
    ierr = PetscOptionsGetStringArray(PETSC_NULL,"-pco2atm_history",pCO2atmFiles,&maxValsToRead,&flg);CHKERRQ(ierr);
    if (flg) { /* Read atmospheric pCO2 history */
      if (maxValsToRead != 2) {
        SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of file names specified for atmospheric pCO2 history");
      }      
      fixedAtmosCO2 = PETSC_FALSE;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading time-dependent atmospheric pCO2 history\n");CHKERRQ(ierr);      
      /* read time data */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,pCO2atmFiles[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	  ierr = PetscBinaryRead(fp,&numpCO2atm_hist,1,PETSC_INT);CHKERRQ(ierr);  
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of points in atmospheric history file is %d \n",numpCO2atm_hist);CHKERRQ(ierr);  
      ierr = PetscMalloc(numpCO2atm_hist*sizeof(PetscScalar),&TpCO2atm_hist);CHKERRQ(ierr); 
      ierr = PetscBinaryRead(fp,TpCO2atm_hist,numpCO2atm_hist,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      /* read atmospheric pCO2 data */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,pCO2atmFiles[1],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
      ierr = PetscMalloc(numpCO2atm_hist*sizeof(PetscScalar),&pCO2atm_hist);CHKERRQ(ierr); 
      ierr = PetscBinaryRead(fp,pCO2atm_hist,numpCO2atm_hist,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      
      pCO2atm = pCO2atm_hist[0];

    }	else {
      ierr = PetscOptionsGetReal(PETSC_NULL,"-pco2atm",&pCO2atm,&flg);CHKERRQ(ierr); /* overwrite default value */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Using fixed atmospheric pCO2 of %g ppm\n",pCO2atm);CHKERRQ(ierr);
      
    }      
  }
  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localph);CHKERRQ(ierr);  

  ierr = PetscOptionsHasName(PETSC_NULL,"-use_virtual_flux",&useVirtualFlux);CHKERRQ(ierr);

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localEmP);CHKERRQ(ierr);
  if (periodicBiogeochemForcing) {    
	localEmPp.firstTime = PETSC_TRUE;
	localEmPp.arrayLength = lNumProfiles;
  } else {  
	ierr = readProfileSurfaceScalarData("EmP.bin",localEmP,1);  
  }

  if (useVirtualFlux) {
	ierr = VecDuplicate(TR,&surfVolFrac);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"surface_volume_fraction.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(surfVolFrac,fd);CHKERRQ(ierr);  
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
  }

#endif

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using Burial-Runoff model\n");CHKERRQ(ierr);  

/* Define the interval over which to integrate global burial */
  ierr = PetscOptionsGetInt(PETSC_NULL,"-burial_sum_steps",&burialSumSteps,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate burial integration interval with the -burial_sum_steps option");
  if ((maxSteps % burialSumSteps)!=0) {
    SETERRQ(PETSC_COMM_WORLD,1,"maxSteps not divisible by burialSumSteps!");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Runoff will be integrated over every %d time steps\n",burialSumSteps);CHKERRQ(ierr);

#ifdef WRITE_RUNOFF
/* set the name of the runoff time file */
  ierr = PetscOptionsGetString(PETSC_NULL,"-runoff_time_file",runoffOutTimeFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) {
  strcpy(runoffOutTimeFile,"");
  sprintf(runoffOutTimeFile,"%s","runoff_output_time.txt");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Runoff output times will be written to %s\n",runoffOutTimeFile);CHKERRQ(ierr);
#endif

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

  ierr = PetscOptionsGetString(PETSC_NULL,"-pickup_runoff_out",runoffIniOutFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) {
    strcpy(runoffIniOutFile,"");
    sprintf(runoffIniOutFile,"%s","pickup_runoff.bin");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Final runoff output will be written to %s\n",runoffIniOutFile);CHKERRQ(ierr);

#ifdef WRITE_RUNOFF
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
#endif

/* fraction of global river runoff in each box, divided by the box volume (a 3D field) */
/* Note: VecLoadVecIntoArray resides in petsc_matvec_utils.c and is not a generic petsc function*/
  ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localrunoffvol);CHKERRQ(ierr);    

#ifdef RUNOFF
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Runoff will be supplied via rivers %g\n");CHKERRQ(ierr); 
      ierr = VecLoadVecIntoArray(TR,"runoff_volume_annual.petsc",localrunoffvol);CHKERRQ(ierr);
#else
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Runoff will be distributed over total ocean area of %g\n",totalA);CHKERRQ(ierr);
      ierr = VecLoadVecIntoArray(TR,"dz.petsc",localrunoffvol);CHKERRQ(ierr);  
      /* IK: loading dz.petsc is just a dummy for now; runoff will be divided by dz(1) in BGC_MODEL.F */ 
#endif

/* Grid arrays */
  ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localdz);CHKERRQ(ierr);    
  ierr = VecLoadVecIntoArray(TR,"dz.petsc",localdz);CHKERRQ(ierr);

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
	ierr = interpPeriodicVector(tc,&Ts,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Ssp,"Ss_");	
#ifdef READ_SWRAD
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localswradp,"swrad_");
#else
   insolation_(&lNumProfiles,&myTime,&locallatitude[0],&daysPerYear,&localswrad[0],&localtau[0]);
#endif                                                 
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localficep,"fice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localwindp,"wind_");   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localatmospp,"atmosp_");					                              
#ifdef CARBON
	ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localEmPp,"EmP_");                                                  
#endif									              
  } else {
#ifndef READ_SWRAD
    insolation_(&lNumProfiles,&myTime,&locallatitude[0],&daysPerYear,&localswrad[0],&localtau[0]);
#endif    
  }
  
  for (ip=0; ip<lNumProfiles; ip++) {
    nzloc=lProfileLength[ip];
    kl=lStartIndices[ip];
	for (itr=0; itr<numTracers; itr++) {    	
	  mops_biogeochem_copy_data_(&nzloc,&itr,&localTR[itr][kl],&localJTR[itr][kl],&DeltaT,&toModel);
	}  
    
    mops_biogeochem_ini_(&nzloc,&DeltaT,
#ifdef CARBON			   
			 &localph[ip],
#endif
			 &localTs[kl],&localSs[kl],&localdz[kl],&drF[0],&nzmax,&nzeuph,
                         &numBiogeochemStepsPerOceanStep,&MYTRUE);
  }

/* Read and overwrite default parameter values here */
  ierr = PetscOptionsGetString(PETSC_NULL,"-bgc_params_file",bgcParamsFile,PETSC_MAX_PATH_LEN-1,&readBGCParams);CHKERRQ(ierr);

  if (readBGCParams) {
    ierr = PetscOptionsGetInt(PETSC_NULL,"-num_bgc_params",&numBGCParams,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate number of BGC parameters to read with the -num_bgc_params option");    

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading %d parameters from file\n",numBGCParams);CHKERRQ(ierr);

#ifdef ASCIIPARAMS
    fpparams=fopen(bgcParamsFile,"r");
    for(ipar=0;ipar<numBGCParams;ipar++) {
      ierr = fscanf(fpparams,"%lf\n",&bgcparams[ipar]);
    }   
    fclose(fpparams);
#else
    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,bgcParamsFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
    ierr = PetscMalloc(numBGCParams*sizeof(PetscScalar),&bgcparams);CHKERRQ(ierr);
    ierr = PetscBinaryRead(fp,bgcparams,numBGCParams,PETSC_SCALAR);CHKERRQ(ierr); 
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
#endif
 
    for (ipar=0; ipar<numBGCParams; ipar++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Parameter no. %d is %f\n",ipar,bgcparams[ipar]);CHKERRQ(ierr);
    }

    mops_biogeochem_set_params_(&numBGCParams,&bgcparams[0]);

    myTime = DeltaT*Iter; /* Iter should start at 0 */
    for (ip=0; ip<lNumProfiles; ip++) {
      nzloc=lProfileLength[ip];
      kl=lStartIndices[ip];

	  for (itr=0; itr<numTracers; itr++) {    	
		mops_biogeochem_copy_data_(&nzloc,&itr,&localTR[itr][kl],&localJTR[itr][kl],&DeltaT,&toModel);
	  }  
      
      mops_biogeochem_ini_(&nzloc,&DeltaT,
#ifdef CARBON			   
			   &localph[ip],
#endif
			   &localTs[kl],&localSs[kl],&localdz[kl],&drF[0],&nzmax,&nzeuph,
                           &numBiogeochemStepsPerOceanStep,&MYFALSE);
    }

  }
  
  ierr = PetscOptionsHasName(PETSC_NULL,"-calc_diagnostics",&calcDiagnostics);CHKERRQ(ierr);
  if (calcDiagnostics) {    
/*Data for diagnostics */
    ierr = iniStepTimer("diag_", Iter0, &diagTimer);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed starting at (and including) time step: %d\n", diagTimer.startTimeStep);CHKERRQ(ierr);	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed over %d time steps\n", diagTimer.numTimeSteps);CHKERRQ(ierr);	

	for (idiag=0; idiag<numDiag; idiag++) {
	  diagOutFile[idiag] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
	}
	maxValsToRead = numDiag;
	ierr = PetscOptionsGetStringArray(PETSC_NULL,"-diag_files",diagOutFile,&maxValsToRead,&flg);CHKERRQ(ierr);
	if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate file name(s) for writing diagnostics with the -diag_files option");
	if (maxValsToRead != numDiag) {
	  SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of time average file names specified");
	}  
	ierr = VecDuplicate(TR,&fbgc1);CHKERRQ(ierr);
	ierr = VecSet(fbgc1,zero);CHKERRQ(ierr);
	ierr = VecGetArray(fbgc1,&localfbgc1);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&fbgc1avg);CHKERRQ(ierr);
	ierr = VecSet(fbgc1avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,diagOutFile[0],FILE_MODE_WRITE,&fdfbgc1avg);CHKERRQ(ierr);

	ierr = VecDuplicate(TR,&fbgc2);CHKERRQ(ierr);
	ierr = VecSet(fbgc2,zero);CHKERRQ(ierr);
	ierr = VecGetArray(fbgc2,&localfbgc2);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&fbgc2avg);CHKERRQ(ierr);
	ierr = VecSet(fbgc2avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,diagOutFile[1],FILE_MODE_WRITE,&fdfbgc2avg);CHKERRQ(ierr);

	ierr = VecDuplicate(TR,&fbgc3);CHKERRQ(ierr);
	ierr = VecSet(fbgc3,zero);CHKERRQ(ierr);
	ierr = VecGetArray(fbgc3,&localfbgc3);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&fbgc3avg);CHKERRQ(ierr);
	ierr = VecSet(fbgc3avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,diagOutFile[2],FILE_MODE_WRITE,&fdfbgc3avg);CHKERRQ(ierr);

	ierr = VecDuplicate(TR,&fbgc4);CHKERRQ(ierr);
	ierr = VecSet(fbgc4,zero);CHKERRQ(ierr);
	ierr = VecGetArray(fbgc4,&localfbgc4);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&fbgc4avg);CHKERRQ(ierr);
	ierr = VecSet(fbgc4avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,diagOutFile[3],FILE_MODE_WRITE,&fdfbgc4avg);CHKERRQ(ierr);

	ierr = VecDuplicate(TR,&fbgc5);CHKERRQ(ierr);
	ierr = VecSet(fbgc5,zero);CHKERRQ(ierr);
	ierr = VecGetArray(fbgc5,&localfbgc5);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&fbgc5avg);CHKERRQ(ierr);
	ierr = VecSet(fbgc5avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,diagOutFile[4],FILE_MODE_WRITE,&fdfbgc5avg);CHKERRQ(ierr);

	ierr = VecDuplicate(TR,&fbgc6);CHKERRQ(ierr);
	ierr = VecSet(fbgc6,zero);CHKERRQ(ierr);
	ierr = VecGetArray(fbgc6,&localfbgc6);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&fbgc6avg);CHKERRQ(ierr);
	ierr = VecSet(fbgc6avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,diagOutFile[5],FILE_MODE_WRITE,&fdfbgc6avg);CHKERRQ(ierr);

	ierr = VecDuplicate(TR,&fbgc7);CHKERRQ(ierr);
	ierr = VecSet(fbgc7,zero);CHKERRQ(ierr);
	ierr = VecGetArray(fbgc7,&localfbgc7);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&fbgc7avg);CHKERRQ(ierr);
	ierr = VecSet(fbgc7avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,diagOutFile[6],FILE_MODE_WRITE,&fdfbgc7avg);CHKERRQ(ierr);

#ifdef CARBON
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localco2airseafluxdiag);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localco2airseafluxdiagavg);CHKERRQ(ierr);  

    for (ip=0; ip<lNumProfiles; ip++) {
      localco2airseafluxdiag[ip]=0.0;
      localco2airseafluxdiagavg[ip]=0.0;
    }    
#endif

// 	diagCount=0;
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
  PetscScalar DICemp = 0.0, ALKemp = 0.0;
  PetscScalar localco2airseaflux = 0.0;
#endif
  
  PetscScalar localburial = 0.0;
  
  myTime = DeltaT*Iter; /* Iter should start at 0 */

  if (periodicBiogeochemForcing) {   
	ierr = interpPeriodicVector(tc,&Ts,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Ssp,"Ss_");	
#ifdef READ_SWRAD
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localswradp,"swrad_");
#else
    insolation_(&lNumProfiles,&myTime,&locallatitude[0],&daysPerYear,&localswrad[0],&localtau[0]);
#endif                                                 
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localficep,"fice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localwindp,"wind_");  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localatmospp,"atmosp_");					                              
#ifdef CARBON
	ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localEmPp,"EmP_");                                                      
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

  if (useVirtualFlux) { /* use the global surface mean value to calculate E-P contribution */
    ierr = VecDot(surfVolFrac,DIC,&DICemp);CHKERRQ(ierr); /* volume weighted mean surface DIC */									              
    ierr = VecDot(surfVolFrac,ALK,&ALKemp);CHKERRQ(ierr); /* volume weighted mean surface ALK */									                  
  }

#endif
  
  for (ip=0; ip<lNumProfiles; ip++) {
    nzloc=lProfileLength[ip];
    kl=lStartIndices[ip];

#ifdef CARBON
	if (!useVirtualFlux) { /* use the local surface value to calculate E-P contribution */
      DICemp=localDIC[kl];
      ALKemp=localALK[kl];      
    }
#endif

	for (itr=0; itr<numTracers; itr++) {    	
	  mops_biogeochem_copy_data_(&nzloc,&itr,&localTR[itr][kl],&localJTR[itr][kl],&DeltaT,&toModel);
	}  
    
    mops_biogeochem_model_(&nzloc,&DeltaT,
#ifdef CARBON			   
			   &DICemp,&ALKemp,&localEmP[ip],&pCO2atm,
#endif
			   &localTs[kl],&localSs[kl],&localfice[ip],&localswrad[ip],&localtau[ip],&localwind[ip],&localatmosp[ip],&localdz[kl],
#ifdef CARBON			   
			   &localph[ip],&localco2airseaflux,
#endif
                           &localburial,&GRunoff,&localrunoffvol[kl],
			   &useSeparateBiogeochemTimeStepping);                            

	for (itr=0; itr<numTracers; itr++) {    
	  mops_biogeochem_copy_data_(&nzloc,&itr,&localTR[itr][kl],&localJTR[itr][kl],&DeltaT,&fromModel);
	}  

#ifdef CARBON			   
    if (useAtmModel) {
      localFocean = localFocean + (localco2airseaflux/DeltaT)*localdA[ip]*(12.0/1.e18)*secondsPerYear; /* PgC/y */
    }
#endif

/* integrate burial in sediment over area and all profiles on each processor */
      localFburial = localFburial + localburial*localdA[ip]; 

	if (calcDiagnostics) {  
	  if (Iter0+iLoop>=diagTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */
        mops_biogeochem_diagnostics_(&nzloc,&localfbgc1[kl],&localfbgc2[kl],&localfbgc3[kl],&localfbgc4[kl],&localfbgc5[kl],&localfbgc6[kl],&localfbgc7[kl]);
#ifdef CARBON        
        localco2airseafluxdiag[ip]=localco2airseaflux;
#endif                
      }
	}

/* IK : added for misfit function : start */
/* IK : added for misfit function : start */
/* IK : added for misfit function : start */
	
      if (doMisfit) {
		if (Iter0+iLoop>=costTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */	
			mops_biogeochem_misfit_(&nzloc,&localmbgc1[kl],&localmbgc2[kl],&localmbgc3[kl]);
		}
      }

/* IK : added for misfit function : end */
/* IK : added for misfit function : end */
/* IK : added for misfit function : end */

  } /* end loop over profiles */

#ifdef CARBON

  if (useAtmModel) {
    if ((iLoop % atmModelUpdateTimeSteps)==0) {  /*  time to update atmosphere */
  
      Focean = 0.0;
       
      MPI_Allreduce(&localFocean, &Focean, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);    
    
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

/* sum burial in sediment over all processors, and scale by time step etc.*/
/* do this only once every burialSumSteps , and then take this value for next year's runoff */

    if ((iLoop % burialSumSteps)==0) {

      Fburial = 0.0;

      MPI_Allreduce(&localFburial, &Fburial, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);    

#ifdef RUNOFF
      GRunoff = Fburial/(1.e12*burialSumSteps)*(86400.0/DeltaT); /* This is Gmol P/day. 
      Note: localrunoff is scaled with 1e12. Note: GRunoff will be scaled with bgc_dt.*/
#else
      GRunoff = Fburial/(totalA*burialSumSteps)*(86400.0/DeltaT); /* This is mmol P/m2/day. 
      Note: this will later be divided by depth of first layer. Note: GRunoff will be scaled with bgc_dt.*/ 
#endif

      localFburial = 0.0;
      }
    


  if (useSeparateBiogeochemTimeStepping) {  /* return updated tracer field */
	for (itr=0; itr<numTracers; itr++) {  
	  ierr = VecSetValues(v[itr],lSize,gIndices,localTR[itr],INSERT_VALUES);CHKERRQ(ierr);
	  ierr = VecAssemblyBegin(v[itr]);CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(v[itr]);CHKERRQ(ierr);    
	}
  } else {
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
  
/*  Convert to discrete tendency */
	for (itr=0; itr<numTracers; itr++) {
	  ierr = VecScale(ut[itr],DeltaT);CHKERRQ(ierr);
	}
  }
  
  if (calcDiagnostics) {  
	if (Iter0+iLoop>=diagTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */
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

	  ierr = VecSetValues(fbgc6,lSize,gIndices,localfbgc6,INSERT_VALUES);CHKERRQ(ierr);
	  ierr = VecAssemblyBegin(fbgc6);CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(fbgc6);CHKERRQ(ierr);      

	  ierr = VecSetValues(fbgc7,lSize,gIndices,localfbgc7,INSERT_VALUES);CHKERRQ(ierr);
	  ierr = VecAssemblyBegin(fbgc7);CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(fbgc7);CHKERRQ(ierr);      
	}  
  }

/* IK : added for misfit function : start */
/* IK : added for misfit function : start */
/* IK : added for misfit function : start */

      if (doMisfit) {
		if (Iter0+iLoop>=costTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */  
		  ierr = VecSetValues(mbgc1,lSize,gIndices,localmbgc1,INSERT_VALUES);CHKERRQ(ierr);
		  ierr = VecAssemblyBegin(mbgc1);CHKERRQ(ierr);
		  ierr = VecAssemblyEnd(mbgc1);CHKERRQ(ierr);    
  
		  ierr = VecSetValues(mbgc2,lSize,gIndices,localmbgc2,INSERT_VALUES);CHKERRQ(ierr);
		  ierr = VecAssemblyBegin(mbgc2);CHKERRQ(ierr);
		  ierr = VecAssemblyEnd(mbgc2);CHKERRQ(ierr);    
  
		  ierr = VecSetValues(mbgc3,lSize,gIndices,localmbgc3,INSERT_VALUES);CHKERRQ(ierr);
		  ierr = VecAssemblyBegin(mbgc3);CHKERRQ(ierr);
		  ierr = VecAssemblyEnd(mbgc3);CHKERRQ(ierr);    
        }
	  }

/* IK : added for misfit function : end */
/* IK : added for misfit function : end */
/* IK : added for misfit function : end */

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


#ifdef WRITE_RUNOFF
    if ((iLoop % burialSumSteps)==0) {  /*  time to write out */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing runoff output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
      ierr = PetscFPrintf(PETSC_COMM_WORLD,runofffptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);           
      ierr = writeBinaryScalarData("Grunoff_output.bin",&GRunoff,1,PETSC_TRUE);
    }
#endif

  if (calcDiagnostics) {  
	if (Iter0+iLoop>=diagTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */  
  
	  if (diagTimer.count<=diagTimer.numTimeSteps) { /* still within same averaging block so accumulate */
		ierr = VecAXPY(fbgc1avg,one,fbgc1);CHKERRQ(ierr);
		ierr = VecAXPY(fbgc2avg,one,fbgc2);CHKERRQ(ierr);
		ierr = VecAXPY(fbgc3avg,one,fbgc3);CHKERRQ(ierr);
		ierr = VecAXPY(fbgc4avg,one,fbgc4);CHKERRQ(ierr);
		ierr = VecAXPY(fbgc5avg,one,fbgc5);CHKERRQ(ierr);
		ierr = VecAXPY(fbgc6avg,one,fbgc6);CHKERRQ(ierr);
		ierr = VecAXPY(fbgc7avg,one,fbgc7);CHKERRQ(ierr);

#ifdef CARBON
        for (ip=0; ip<lNumProfiles; ip++) {
          localco2airseafluxdiagavg[ip]=localco2airseafluxdiag[ip]+localco2airseafluxdiagavg[ip];      
        }	  
#endif

		diagTimer.count++; // = diagCount+1;
	  }
	  if (diagTimer.count==diagTimer.numTimeSteps) { /* time to write averages to file */
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing diagnostics time average at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);                      

		ierr = VecScale(fbgc1avg,1.0/diagTimer.count);CHKERRQ(ierr);
		ierr = VecView(fbgc1avg,fdfbgc1avg);CHKERRQ(ierr);
		ierr = VecSet(fbgc1avg,zero); CHKERRQ(ierr);

		ierr = VecScale(fbgc2avg,1.0/diagTimer.count);CHKERRQ(ierr);
		ierr = VecView(fbgc2avg,fdfbgc2avg);CHKERRQ(ierr);
		ierr = VecSet(fbgc2avg,zero); CHKERRQ(ierr);

		ierr = VecScale(fbgc3avg,1.0/diagTimer.count);CHKERRQ(ierr);
		ierr = VecView(fbgc3avg,fdfbgc3avg);CHKERRQ(ierr);
		ierr = VecSet(fbgc3avg,zero); CHKERRQ(ierr);

		ierr = VecScale(fbgc4avg,1.0/diagTimer.count);CHKERRQ(ierr);
		ierr = VecView(fbgc4avg,fdfbgc4avg);CHKERRQ(ierr);
		ierr = VecSet(fbgc4avg,zero); CHKERRQ(ierr);

		ierr = VecScale(fbgc5avg,1.0/diagTimer.count);CHKERRQ(ierr);
		ierr = VecView(fbgc5avg,fdfbgc5avg);CHKERRQ(ierr);
		ierr = VecSet(fbgc5avg,zero); CHKERRQ(ierr);

		ierr = VecScale(fbgc6avg,1.0/diagTimer.count);CHKERRQ(ierr);
		ierr = VecView(fbgc6avg,fdfbgc6avg);CHKERRQ(ierr);
		ierr = VecSet(fbgc6avg,zero); CHKERRQ(ierr);

		ierr = VecScale(fbgc7avg,1.0/diagTimer.count);CHKERRQ(ierr);
		ierr = VecView(fbgc7avg,fdfbgc7avg);CHKERRQ(ierr);
		ierr = VecSet(fbgc7avg,zero); CHKERRQ(ierr);

#ifdef CARBON
        for (ip=0; ip<lNumProfiles; ip++) {
          localco2airseafluxdiagavg[ip]=localco2airseafluxdiagavg[ip]/diagTimer.count;
        }	  

        ierr = writeProfileSurfaceScalarData("co2airseaflux_surf.bin",localco2airseafluxdiagavg,1,appendDiagnostics);  		

/*      reset diagnostic arrays */
        for (ip=0; ip<lNumProfiles; ip++) {
          localco2airseafluxdiagavg[ip]=0.0;
        }	  
#endif

        appendDiagnostics=PETSC_TRUE;
        ierr = updateStepTimer("diag_", Iter0+iLoop, &diagTimer);CHKERRQ(ierr);
// 		diagCount = 0;        

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

  ierr = writeBinaryScalarData(runoffIniOutFile,&GRunoff,1,PETSC_FALSE);
  
  ierr = VecDestroy(&Ts);CHKERRQ(ierr);
  ierr = VecDestroy(&Ss);CHKERRQ(ierr);

  if (periodicBiogeochemForcing) {    
    ierr = destroyPeriodicVec(&Tsp);CHKERRQ(ierr);
    ierr = destroyPeriodicVec(&Ssp);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localficep);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localwindp);CHKERRQ(ierr);    
    ierr = destroyPeriodicArray(&localatmospp);CHKERRQ(ierr);    
#ifdef READ_SWRAD    
    ierr = destroyPeriodicArray(&localswradp);CHKERRQ(ierr);
#endif
#ifdef CARBON
	ierr = destroyPeriodicArray(&localEmPp);CHKERRQ(ierr);
#endif	
  }    

#ifdef CARBON
  if (useVirtualFlux) {
    ierr = VecDestroy(&surfVolFrac);CHKERRQ(ierr);  
  }
#endif

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

	ierr = VecDestroy(&fbgc6);CHKERRQ(ierr);
	ierr = VecDestroy(&fbgc6avg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fdfbgc6avg);CHKERRQ(ierr);	

	ierr = VecDestroy(&fbgc7);CHKERRQ(ierr);
	ierr = VecDestroy(&fbgc7avg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fdfbgc7avg);CHKERRQ(ierr);	

  }

#ifdef CARBON
  if (useAtmModel) {
    ierr = PetscFClose(PETSC_COMM_WORLD,atmfptime);CHKERRQ(ierr);
  }
#endif

#ifdef WRITE_RUNOFF
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
  PetscInt ipar;

  myTime = DeltaT*Iter; /* Iter should start at 0 */

  if (periodicBiogeochemForcing) {   
	ierr = interpPeriodicVector(tc,&Ts,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Ssp,"Ss_");	
#ifdef READ_SWRAD
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localswradp,"swrad_");
#else
    insolation_(&lNumProfiles,&myTime,&locallatitude[0],&daysPerYear,&localswrad[0],&localtau[0]);
#endif                                                 
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localficep,"fice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localwindp,"wind_");   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localatmospp,"atmosp_");	
#ifdef CARBON
	ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localEmPp,"EmP_");                                                      
#endif									              
  }

  if (readBGCParams) {

#ifdef ASCIIPARAMS
    fpparams=fopen(bgcParamsFile,"r");
    for(ipar=0;ipar<numBGCParams;ipar++) {
      ierr = fscanf(fpparams,"%lf\n",&bgcparams[ipar]);
    }   
    fclose(fpparams);
#else
    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,bgcParamsFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
    ierr = PetscMalloc(numBGCParams*sizeof(PetscScalar),&bgcparams);CHKERRQ(ierr); 
    ierr = PetscBinaryRead(fp,bgcparams,numBGCParams,PETSC_SCALAR);CHKERRQ(ierr);  
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
#endif
 
    for (ipar=0; ipar<numBGCParams; ipar++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Parameter no. %d is %f\n",ipar,bgcparams[ipar]);CHKERRQ(ierr);
    }

    mops_biogeochem_set_params_(&numBGCParams,&bgcparams[0]);

    for (ip=0; ip<lNumProfiles; ip++) {
      nzloc=lProfileLength[ip];
      kl=lStartIndices[ip];
      mops_biogeochem_ini_(&nzloc,&DeltaT,
#ifdef CARBON
			   &localph[ip],
#endif
			   &localTs[kl],&localSs[kl],&localdz[kl],&drF[0],&nzmax,&nzeuph,
                           &numBiogeochemStepsPerOceanStep,&MYFALSE);
    }    
  }
    
  return 0;
}
