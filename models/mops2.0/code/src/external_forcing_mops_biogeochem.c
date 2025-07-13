/* $Header: /Users/ikriest/CVS/mops/external_forcing_mops_biogeochem.c,v 1.3 2016/06/06 09:43:41 ikriest Exp $ */
/* $Name: mops-2_0 $*/

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#undef READ_SWRAD
#include "petscmat.h"
#include "tmm_petsc_matvec_utils.h"
#include "tmm_timer.h"
#include "tmm_forcing_utils.h"
#include "tmm_profile_utils.h"
#include "tmm.h"
#include "tmm_share.h"
#include "mops_biogeochem_tmm.h"
#include "mops_forcing.h"

/* Macros to map tracer names to vectors */
/* Note: for MOPS, we have the following tracer assignement:*/
/* c[0] PO4; c[1] DOP; c[2] O2; c[3] Phy; c[4] Zoo; c[5] Det; c[6] NO3 */
/* Additionally, or option -DCARBON: c[7] DIC; c[8] Alk */
/* Note that this also affects BGC_PARAMS.h and the runsccript(s) */

#define TR state->c[0]
#define DIC state->c[7]
#define ALK state->c[8]
#define localDIC ef->localTR[7]
#define localALK ef->localTR[8]

static PetscClassId EXTERNALFORCING_CLASSID;

// Common variables
PetscBool debug = PETSC_FALSE;
PetscInt ipdebug = 2278;
PetscScalar ppmToPgC=2.1324;
PetscScalar totalA = 0.0;

PetscBool useSeparateBiogeochemTimeStepping = PETSC_FALSE;
PetscInt numBiogeochemStepsPerOceanStep = 1;
PetscInt nzmax,nzeuph;
PetscScalar DeltaT,TheoDeltaT;
PetscScalar *drF;
PetscScalar *localdz;

PetscScalar *localrunoffvol; /* volume supplied by runoff */

Vec Ts,Ss;
PetscScalar *localTs,*localSs;
#ifdef CARBON
  PetscBool useVirtualFlux;
  PetscScalar *localEmP;
  PeriodicArray localEmPp;
  Vec surfVolFrac;
#endif

PetscScalar *localwind,*localfice,*localatmosp;
PetscScalar *localswrad, *localtau;
#ifndef READ_SWRAD
PetscScalar *locallatitude;
#endif

PeriodicVec Tsp, Ssp;
PeriodicArray localwindp,localficep,localatmospp;
#ifdef READ_SWRAD
PeriodicArray localswradp;
#endif

PetscBool periodicBiogeochemForcing;
PeriodicTimer biogeochemTimer;

PetscInt toModel = 1; 
PetscInt fromModel = 2;

PetscBool readBGCParams = PETSC_FALSE;
PetscInt numBGCParams = 0, lenBGCParam;
char bgcParamsFile[PETSC_MAX_PATH_LEN], numBGCParamsFile[PETSC_MAX_PATH_LEN], bgcNamesFile[PETSC_MAX_PATH_LEN];
FILE *fpparams, *fpnumbgc, *fpnames;
char bgcnames[20][PETSC_MAX_PATH_LEN];
char allowedBGCnames[25][PETSC_MAX_PATH_LEN] = {"rcp", "rnp", "ro2ut", "subox", "subdin", "TempB", "ACmuphy", "ACik", "ACkpo4",
												"AClambda", "AComni", "plambda", "ACmuzoo", "AClambdaz", "AComniz", "ACeff", "zlambda", "graztodop",
												"dlambda", "detlambda", "detmartin", "burdige_fac", "burdige_exp", "ACkbaco2", "ACkbacdin"};
PetscInt maxAllowed = 25, iall;
PetscInt allowedBGC = 0, set_AConmiz = 0;

#ifdef ASCIIPARAMS
PetscScalar bgcparams[20];
#else
PetscScalar *bgcparams;
#endif

PetscScalar *localdA;

PetscInt maxValsToRead;

PetscScalar daysPerYear, secondsPerYear;

PetscInt idiag;

PetscBool MYTRUE = PETSC_TRUE, MYFALSE = PETSC_FALSE;

/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/




#undef __FUNCT__
#define __FUNCT__ "iniExternalForcing"
PetscErrorCode iniExternalForcing(PetscScalar tc, PetscInt Iter, TMMState state, void *userctx)
{
  PetscInt numTracers;
  const char *prefix;

  PetscErrorCode ierr;
  PetscInt ip, kl, nzloc;
  PetscInt itr;
  PetscViewer fd;
  PetscInt fp;
  PetscBool flg;
  PetscInt it,ipar;
  PetscScalar myTime;
  PetscScalar zero = 0.0;
#ifdef CARBON  
  char *pCO2atmFiles[2];  
  char pCO2atmIniFile[PETSC_MAX_PATH_LEN];  
  PetscScalar pCO2atm_ini = 280.0; /* default initial value */
  char atmOutTimeFile[PETSC_MAX_PATH_LEN];  
#endif
  PetscScalar runoff_ini = 0.0;
  char runoffIniFile[PETSC_MAX_PATH_LEN];  
  char *diagOutFile[12];

  static PetscBool registered = PETSC_FALSE;
  static PetscInt efctxId = 0;

  ExternalForcingContext ef;

  MPI_Comm comm = PETSC_COMM_WORLD;
  
  if (!registered) {
    PetscClassIdRegister("ExternalForcing context", &EXTERNALFORCING_CLASSID);
    registered = PETSC_TRUE;
  }
  PetscHeaderCreate(ef, EXTERNALFORCING_CLASSID, "ExternalForcing", "ExternalForcing context", "ExternalForcing", comm, 0, 0); 

  efctxId++;
  ef->efctxId=efctxId;
  ef->stateId=state->stateId;

  PetscContainer ctxcontainer;
  PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
  PetscCall(PetscContainerSetPointer(ctxcontainer, (void*)ef));
  PetscCall(PetscObjectCompose((PetscObject)state, "external forcing ctx", (PetscObject)ctxcontainer));
  state->extforcctxcontainer = ctxcontainer;
  PetscCall(PetscContainerDestroy(&ctxcontainer));

//   if (ef->stateId != state->id) {
//     SETERRQ(PETSC_COMM_WORLD,1,"ERROR: Mismatch between user context and state!");
//   }
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

// Now set problem specific data
// Common data (only initialize/read once)
  if (ef->efctxId==1) {
	ierr = PetscOptionsHasName(NULL,NULL,"-separate_biogeochem_time_stepping",&useSeparateBiogeochemTimeStepping);CHKERRQ(ierr);

	if (useSeparateBiogeochemTimeStepping) {
	  fromModel = 3;  
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Biogeochem model will be time-stepped independently\n");CHKERRQ(ierr);
	}  

	ierr = PetscOptionsGetInt(NULL,NULL,"-num_biogeochem_steps_per_ocean_step",&numBiogeochemStepsPerOceanStep,&flg);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of biogeochem model time steps per ocean time step = %d\n",numBiogeochemStepsPerOceanStep);CHKERRQ(ierr);

	ierr = PetscOptionsGetInt(NULL,NULL,"-nzeuph",&nzeuph,&flg);CHKERRQ(ierr);
	if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate number of euphotic zone layers with the -nzeuph option");  
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of euphotic zone layers is %d \n",nzeuph);CHKERRQ(ierr);

	ierr = PetscOptionsGetReal(NULL,NULL,"-biogeochem_deltat",&DeltaT,&flg);CHKERRQ(ierr);
	if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate biogeochemical time step in seconds with the -biogeochem_deltat option");  
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Ocean time step for BGC length is  %12.7f seconds\n",DeltaT);CHKERRQ(ierr);

	ierr = PetscOptionsGetReal(NULL,NULL,"-days_per_year",&daysPerYear,&flg);CHKERRQ(ierr);
	if (!flg) {
	  daysPerYear = 360.0;
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of days per year is %12.7f\n",daysPerYear);CHKERRQ(ierr);
	secondsPerYear = 86400.0*daysPerYear;

	TheoDeltaT = daysPerYear*86400.0*deltaTClock;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Check: using a year length of %12.3f days \n",daysPerYear);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Theoretical ocean time step length for BGC is then  %12.7f seconds\n",TheoDeltaT);CHKERRQ(ierr);

/* Need this for atmospheric exchange, river runoff, ... */
	ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localdA);CHKERRQ(ierr);
	ierr = readProfileSurfaceScalarData("dA.bin",localdA,1);
	PetscScalar localA = 0.0;
	for (ip=0; ip<lNumProfiles; ip++) {
	  localA = localA+localdA[ip];
	}  
	MPI_Allreduce(&localA, &totalA, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD); /* global surface area */

/* Grid arrays */
	ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localdz);CHKERRQ(ierr);    
	ierr = VecLoadVecIntoArray(TR,"dz.petsc",localdz);CHKERRQ(ierr);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"drF.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscBinaryRead(fp,&nzmax,1,NULL,PETSC_INT);CHKERRQ(ierr);  
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of vertical layers is %d \n",nzmax);CHKERRQ(ierr);  
	ierr = PetscMalloc(nzmax*sizeof(PetscScalar),&drF);CHKERRQ(ierr); 
	ierr = PetscBinaryRead(fp,drF,nzmax,NULL,PETSC_SCALAR);CHKERRQ(ierr);  
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  
	periodicBiogeochemForcing = PETSC_FALSE;

	ierr = PetscOptionsHasName(NULL,NULL,"-periodic_biogeochem_forcing",&periodicBiogeochemForcing);CHKERRQ(ierr);

	if (periodicBiogeochemForcing) {    
	  ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic biogeochemical forcing specified\n");CHKERRQ(ierr);
	  ierr = PeriodicTimerCreate(&biogeochemTimer);CHKERRQ(ierr);
	  ierr = PeriodicTimerIni("periodic_biogeochem_", NULL, NULL, biogeochemTimer);CHKERRQ(ierr);
	}
  
/*   Read T and S */
	ierr = VecDuplicate(TR,&Ts);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&Ss);CHKERRQ(ierr);  
	if (periodicBiogeochemForcing) {    
	  PeriodicVecCreate(&Tsp);
	  PeriodicVecCreate(&Ssp);
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

#ifdef CARBON
// Defaults
	useVirtualFlux = PETSC_FALSE;

	ierr = PetscOptionsHasName(NULL,NULL,"-use_virtual_flux",&useVirtualFlux);CHKERRQ(ierr);

	ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localEmP);CHKERRQ(ierr);
	if (periodicBiogeochemForcing) {
	  ierr = PeriodicArrayCreate(&localEmPp, lNumProfiles);
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

/* fraction of global river runoff in each box, divided by the box volume (a 3D field) */
/* Note: VecLoadVecIntoArray resides in petsc_matvec_utils.c and is not a generic petsc function*/
	ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localrunoffvol);CHKERRQ(ierr);    

#ifdef RUNOFF
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Runoff will be supplied via rivers\n");CHKERRQ(ierr); 
	ierr = VecLoadVecIntoArray(TR,"runoff_volume_annual.petsc",localrunoffvol);CHKERRQ(ierr);
#else
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Runoff will be distributed over total ocean area of %g\n",totalA);CHKERRQ(ierr);
	ierr = VecLoadVecIntoArray(TR,"dz.petsc",localrunoffvol);CHKERRQ(ierr);  
	/* IK: loading dz.petsc is just a dummy for now; runoff will be divided by dz(1) in BGC_MODEL.F */ 
#endif

/* Forcing fields */  
	ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localswrad);CHKERRQ(ierr);  
	ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localtau);CHKERRQ(ierr);  
#ifdef READ_SWRAD
	if (periodicBiogeochemForcing) {
	  ierr = PeriodicArrayCreate(&localswradp, lNumProfiles);CHKERRQ(ierr);
	} else {  
	  ierr = readProfileSurfaceScalarData("swrad.bin",localswrad,1);  
	}
#else
	ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&locallatitude);CHKERRQ(ierr);  
	ierr = readProfileSurfaceScalarData("latitude.bin",locallatitude,1);  
#endif	  
  
	ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localfice);CHKERRQ(ierr);  
	if (periodicBiogeochemForcing) {    
	  ierr = PeriodicArrayCreate(&localficep, lNumProfiles);CHKERRQ(ierr);
	} else {  
	  ierr = readProfileSurfaceScalarData("fice.bin",localfice,1);  
	}

	ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localwind);CHKERRQ(ierr);  
	if (periodicBiogeochemForcing) {    
	  ierr = PeriodicArrayCreate(&localwindp, lNumProfiles);CHKERRQ(ierr);
	} else {  
	  ierr = readProfileSurfaceScalarData("wind.bin",localwind,1);  
	}

	ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localatmosp);CHKERRQ(ierr);  
	if (periodicBiogeochemForcing) {    
	  ierr = PeriodicArrayCreate(&localatmospp, lNumProfiles);CHKERRQ(ierr);
	} else {  
	  ierr = readProfileSurfaceScalarData("atmosp.bin",localatmosp,1);  
	}

  } /* end common data */

  for (itr=0; itr<numTracers; itr++) {
    ierr = VecSet(state->qef[itr],zero); CHKERRQ(ierr);
  }

  ierr = VecGetArrays(state->c,numTracers,&ef->localTR);CHKERRQ(ierr);
  ierr = VecGetArrays(state->qef,numTracers,&ef->localJTR);CHKERRQ(ierr);

// Some defaults
  ef->readBGCParams = PETSC_FALSE;
  ef->calcDiagnostics = PETSC_FALSE;
  ef->appendDiagnostics = PETSC_FALSE;
  ef->numDiag=7;
 
#ifdef CARBON
// Defaults
  ef->numpCO2atm_hist = 0;
  ef->fixedAtmosCO2 = PETSC_TRUE;
  ef->useAtmModel = PETSC_FALSE;
  ef->pCO2atm = 280.0; /* default initial value */
  ef->Focean=0.0;
  ef->Foceanint = 0.0;

  ierr = PetscOptionsHasName(NULL,prefix,"-use_atm_model",&ef->useAtmModel);CHKERRQ(ierr);

  if (ef->useAtmModel) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using interactive atmospheric model\n");CHKERRQ(ierr);  

/* overwrite default value */
	ierr = PetscOptionsGetReal(NULL,prefix,"-pco2atm_ini",&pCO2atm_ini,&flg);CHKERRQ(ierr); /* read from command line */
    if (!flg) {
      ierr = PetscOptionsGetString(NULL,prefix,"-pco2atm_ini_file",pCO2atmIniFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
      if (flg) { /* read from binary file */
        ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,pCO2atmIniFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
        ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
        ierr = PetscBinaryRead(fp,&pCO2atm_ini,1,NULL,PETSC_SCALAR);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      }
    }
    ef->pCO2atm = pCO2atm_ini;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using initial atmospheric pCO2 of %g ppm\n",ef->pCO2atm);CHKERRQ(ierr);
    
    ierr = StepTimerCreate(&ef->atmWriteTimer);CHKERRQ(ierr);
	ierr = StepTimerIni("atm_write_", prefix, Iter0, ef->atmWriteTimer);CHKERRQ(ierr);

    ierr = PetscOptionsHasName(NULL,prefix,"-atm_append",&ef->atmAppendOutput);CHKERRQ(ierr);
    if (ef->atmAppendOutput) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output will be appended\n");CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output will overwrite existing file(s)\n");CHKERRQ(ierr);
    }    

/* Output times */
    ierr = PetscOptionsGetString(NULL,prefix,"-atm_time_file",atmOutTimeFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
    if (!flg) {
      strcpy(atmOutTimeFile,"");
      sprintf(atmOutTimeFile,"%s","atm_output_time.txt");
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output times will be written to %s\n",atmOutTimeFile);CHKERRQ(ierr);

    if (!ef->atmAppendOutput) {
      ierr = PetscFOpen(PETSC_COMM_WORLD,atmOutTimeFile,"w",&ef->atmfptime);CHKERRQ(ierr);  
      ierr = PetscFPrintf(PETSC_COMM_WORLD,ef->atmfptime,"%d   %10.5f\n",Iter0,time0);CHKERRQ(ierr);     
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing atmospheric output at time %10.5f, step %d\n", tc,Iter);CHKERRQ(ierr);  
      ierr = writeBinaryScalarData("pCO2atm_output.bin",&ef->pCO2atm,1,PETSC_FALSE);
      ierr = writeBinaryScalarData("Foceanint_output.bin",&ef->Focean,1,PETSC_FALSE);
    } else {
      ierr = PetscFOpen(PETSC_COMM_WORLD,atmOutTimeFile,"a",&ef->atmfptime);CHKERRQ(ierr);  
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output will be appended. Initial condition will NOT be written\n");CHKERRQ(ierr);      
    }

    ef->atmModelDeltaT = DeltaT/secondsPerYear; /* time step in years */

  } else {  /* not using atm model */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using prescribed atmospheric pCO2\n");CHKERRQ(ierr);

    /* prescribed atmospheric CO2 */
  	maxValsToRead = 2;
	pCO2atmFiles[0] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* time file */
	pCO2atmFiles[1] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* atmospheric pCO2 history file */
    ierr = PetscOptionsGetStringArray(NULL,prefix,"-pco2atm_history",pCO2atmFiles,&maxValsToRead,&flg);CHKERRQ(ierr);
    if (flg) { /* Read atmospheric pCO2 history */
      if (maxValsToRead != 2) {
        SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of file names specified for atmospheric pCO2 history");
      }      
      ef->fixedAtmosCO2 = PETSC_FALSE;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading time-dependent atmospheric pCO2 history\n");CHKERRQ(ierr);      
      /* read time data */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,pCO2atmFiles[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	  ierr = PetscBinaryRead(fp,&ef->numpCO2atm_hist,1,NULL,PETSC_INT);CHKERRQ(ierr);  
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of points in atmospheric history file is %d \n",ef->numpCO2atm_hist);CHKERRQ(ierr);  
      ierr = PetscMalloc(ef->numpCO2atm_hist*sizeof(PetscScalar),&ef->TpCO2atm_hist);CHKERRQ(ierr); 
      ierr = PetscBinaryRead(fp,ef->TpCO2atm_hist,ef->numpCO2atm_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      /* read atmospheric pCO2 data */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,pCO2atmFiles[1],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
      ierr = PetscMalloc(ef->numpCO2atm_hist*sizeof(PetscScalar),&ef->pCO2atm_hist);CHKERRQ(ierr); 
      ierr = PetscBinaryRead(fp,ef->pCO2atm_hist,ef->numpCO2atm_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      
      ef->pCO2atm = ef->pCO2atm_hist[0];

    }	else {
      ierr = PetscOptionsGetReal(NULL,prefix,"-pco2atm",&ef->pCO2atm,&flg);CHKERRQ(ierr); /* overwrite default value */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Using fixed atmospheric pCO2 of %g ppm\n",ef->pCO2atm);CHKERRQ(ierr);
      
    }      
  }
  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&ef->localph);CHKERRQ(ierr);  
#endif

// Defaults
  ef->localFburial = 0.0;
  ef->Fburial=0.0;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using Burial-Runoff model\n");CHKERRQ(ierr);  

/* Define the interval over which to integrate global burial */
  ierr = PetscOptionsGetInt(NULL,prefix,"-burial_sum_steps",&ef->burialSumSteps,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate burial integration interval with the -burial_sum_steps option");
  if ((maxSteps % ef->burialSumSteps)!=0) {
    SETERRQ(PETSC_COMM_WORLD,1,"maxSteps not divisible by burialSumSteps!");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Runoff will be integrated over every %d time steps\n",ef->burialSumSteps);CHKERRQ(ierr);

#ifdef WRITE_RUNOFF
/* set the name of the runoff time file */
  ierr = PetscOptionsGetString(NULL,prefix,"-runoff_time_file",runoffOutTimeFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) {
  strcpy(runoffOutTimeFile,"");
  sprintf(runoffOutTimeFile,"%s","runoff_output_time.txt");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Runoff output times will be written to %s\n",runoffOutTimeFile);CHKERRQ(ierr);
#endif

/* set inititial runoff: overwrite default value with value from command line*/
  ierr = PetscOptionsGetReal(NULL,prefix,"-runoff_ini",&runoff_ini,&flg);CHKERRQ(ierr);
/* set inititial runoff: overwrite default value with value from file*/
  if (!flg) {
    ierr = PetscOptionsGetString(NULL,prefix,"-runoff_ini_file",runoffIniFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
    if (flg) { /* read from binary file */
      ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,runoffIniFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
      ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
      ierr = PetscBinaryRead(fp,&runoff_ini,1,NULL,PETSC_SCALAR);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
    }
  }
    
  ef->GRunoff = runoff_ini;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using initial runoff of %g Gmol P/d\n",ef->GRunoff);CHKERRQ(ierr);

  ierr = PetscOptionsGetString(NULL,prefix,"-pickup_runoff_out",ef->runoffIniOutFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) {
    strcpy(ef->runoffIniOutFile,"");
    sprintf(ef->runoffIniOutFile,"%s","pickup_runoff.bin");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Final runoff output will be written to %s\n",ef->runoffIniOutFile);CHKERRQ(ierr);

#ifdef WRITE_RUNOFF
/* if run is continued, always append runoff and output times */
  if (Iter0>0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Runoff output will be appended\n");CHKERRQ(ierr);
      ierr = PetscFOpen(PETSC_COMM_WORLD,runoffOutTimeFile,"a",&ef->runofffptime);CHKERRQ(ierr);  
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Initial runoff output will not be written\n");CHKERRQ(ierr);
   } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Runoff output will overwrite existing file(s)\n");CHKERRQ(ierr);
      ierr = PetscFOpen(PETSC_COMM_WORLD,runoffOutTimeFile,"w",&ef->runofffptime);CHKERRQ(ierr);  
      ierr = PetscFPrintf(PETSC_COMM_WORLD,runofffptime,"%d   %10.5f\n",Iter0,time0);CHKERRQ(ierr);     
      ierr = writeBinaryScalarData("Grunoff_output.bin",&ef->GRunoff,1,PETSC_FALSE);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing runoff output at time %10.5f, step %d\n", tc,Iter);CHKERRQ(ierr);  
   }
#endif


/* Initialize biogeochem model */
  myTime = DeltaT*Iter; /* Iter should start at 0 */

  if (ef->efctxId==1) {
	if (periodicBiogeochemForcing) {   
	  ierr = PeriodicVecInterp(tc,&Ts,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,Tsp,"Ts_");
	  ierr = PeriodicVecInterp(tc,&Ss,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,Ssp,"Ss_");	
#ifdef READ_SWRAD
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,localswradp,"swrad_");
#else
	 insolation_(&lNumProfiles,&myTime,&locallatitude[0],&daysPerYear,&localswrad[0],&localtau[0]);
#endif                                                 
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,localficep,"fice_");
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,localwindp,"wind_");   
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,localatmospp,"atmosp_");					                              
#ifdef CARBON
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,localEmPp,"EmP_");                                                  
#endif									              
	} else {
#ifndef READ_SWRAD
	  insolation_(&lNumProfiles,&myTime,&locallatitude[0],&daysPerYear,&localswrad[0],&localtau[0]);
#endif    
	}
  }
  
  for (ip=0; ip<lNumProfiles; ip++) {
    nzloc=lProfileLength[ip];
    kl=lStartIndices[ip];
	for (itr=0; itr<numTracers; itr++) {    	
	  mops_biogeochem_copy_data_(&nzloc,&itr,&ef->localTR[itr][kl],&ef->localJTR[itr][kl],&DeltaT,&toModel);
	}  
    
    mops_biogeochem_ini_(&nzloc,&DeltaT,
#ifdef CARBON			   
			 &ef->localph[ip],
#endif
			 &localTs[kl],&localSs[kl],&localdz[kl],&drF[0],&nzmax,&nzeuph,
                         &numBiogeochemStepsPerOceanStep,&MYTRUE);
  }

/* Read and overwrite default parameter values here */
  ierr = PetscOptionsGetString(NULL,NULL,"-bgc_params_file",bgcParamsFile,PETSC_MAX_PATH_LEN-1,&readBGCParams);CHKERRQ(ierr);

  if (readBGCParams) {
    ierr = PetscOptionsGetString(NULL,NULL,"-num_bgc_params_file",numBGCParamsFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must provide a .txt file containing the number of BGC parameters to read with the -num_bgc_params_file option");    
	fpnumbgc=fopen(numBGCParamsFile,"r");
	ierr = fscanf(fpnumbgc,"%d\n",&numBGCParams);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading %d parameters from files\n",numBGCParams);CHKERRQ(ierr);
    
    ierr = PetscOptionsGetString(NULL,NULL,"-bgc_paramnames_file",bgcNamesFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must provide a .txt file containing the names of the BGC parameters to read with the -bgc_paramnames_file option");   
    fpnames=fopen(bgcNamesFile,"r");
    
    for(ipar=0;ipar<numBGCParams;ipar++) {
      allowedBGC = 0;
      ierr = fscanf(fpnames,"%s\n",&bgcnames[ipar]);
      for(iall=0;iall<maxAllowed;iall++) {
      	if (strcmp(bgcnames[ipar], allowedBGCnames[iall]) == 0) {
      		allowedBGC = 1;
      		break;
      	}
      }
      ierr = PetscPrintf(PETSC_COMM_WORLD,"This param is %s\n", bgcnames[ipar]);CHKERRQ(ierr);  
      ierr = PetscPrintf(PETSC_COMM_WORLD,"allowedBGC is %d\n", allowedBGC);CHKERRQ(ierr);
      if (allowedBGC==0) SETERRQ(PETSC_COMM_WORLD,1,"The specified parameters are not one of the allowed parameter options. Check they match those in mops_biogeochem_set_params.F"); 
      
      /* Check rnp is set before AComniz, as when updating parameters it gets multiplied by rnp, so need to update rnp first */
      if (strcmp(bgcnames[ipar], "AConmiz") == 0) {
      	set_AConmiz = 1;
      } else if (strcmp(bgcnames[ipar], "rnp") == 0) {
      	if (set_AConmiz == 1) SETERRQ(PETSC_COMM_WORLD,1,"In the specified parameters files rnp must be ordered before AComniz"); 
      }
    }   
    fclose(fpnames);

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
    ierr = PetscBinaryRead(fp,bgcparams,numBGCParams,NULL,PETSC_SCALAR);CHKERRQ(ierr); 
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
#endif
 
    for (ipar=0; ipar<numBGCParams; ipar++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Parameter no. %d is named %s with the value %f\n",ipar,bgcnames[ipar],bgcparams[ipar]);CHKERRQ(ierr);
    }

	/* Call mops_biogeochem_set_params.F for every parameter */
	for (ipar=0; ipar<numBGCParams; ipar++) {
		lenBGCParam = strlen(bgcnames[ipar]);
    	mops_biogeochem_set_params_(&lenBGCParam,&bgcparams[ipar],&bgcnames[ipar]);
    }

    myTime = DeltaT*Iter; /* Iter should start at 0 */
    for (ip=0; ip<lNumProfiles; ip++) {
      nzloc=lProfileLength[ip];
      kl=lStartIndices[ip];

	  for (itr=0; itr<numTracers; itr++) {    	
		mops_biogeochem_copy_data_(&nzloc,&itr,&ef->localTR[itr][kl],&ef->localJTR[itr][kl],&DeltaT,&toModel);
	  }  
      
      mops_biogeochem_ini_(&nzloc,&DeltaT,
#ifdef CARBON			   
			   &ef->localph[ip],
#endif
			   &localTs[kl],&localSs[kl],&localdz[kl],&drF[0],&nzmax,&nzeuph,
                           &numBiogeochemStepsPerOceanStep,&MYFALSE);
    }

  }
  
  ierr = PetscOptionsHasName(NULL,prefix,"-calc_diagnostics",&ef->calcDiagnostics);CHKERRQ(ierr);
  if (ef->calcDiagnostics) {    
/*Data for diagnostics */
    ierr = StepTimerCreate(&ef->diagTimer);CHKERRQ(ierr);
    ierr = StepTimerIni("diag_", prefix, Iter0, ef->diagTimer);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed starting at (and including) time step: %d\n", ef->diagTimer->startTimeStep);CHKERRQ(ierr);	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed over %d time steps\n", ef->diagTimer->numTimeSteps);CHKERRQ(ierr);	

	for (idiag=0; idiag<ef->numDiag; idiag++) {
	  diagOutFile[idiag] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
	}
	maxValsToRead = ef->numDiag;
	ierr = PetscOptionsGetStringArray(NULL,prefix,"-diag_files",diagOutFile,&maxValsToRead,&flg);CHKERRQ(ierr);
	if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate file name(s) for writing diagnostics with the -diag_files option");
	if (maxValsToRead != ef->numDiag) {
	  SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of time average file names specified");
	}  
	ierr = VecDuplicate(TR,&ef->fbgc1);CHKERRQ(ierr);
	ierr = VecSet(ef->fbgc1,zero);CHKERRQ(ierr);
	ierr = VecGetArray(ef->fbgc1,&ef->localfbgc1);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&ef->fbgc1avg);CHKERRQ(ierr);
	ierr = VecSet(ef->fbgc1avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,diagOutFile[0],FILE_MODE_WRITE,&ef->fdfbgc1avg);CHKERRQ(ierr);

	ierr = VecDuplicate(TR,&ef->fbgc2);CHKERRQ(ierr);
	ierr = VecSet(ef->fbgc2,zero);CHKERRQ(ierr);
	ierr = VecGetArray(ef->fbgc2,&ef->localfbgc2);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&ef->fbgc2avg);CHKERRQ(ierr);
	ierr = VecSet(ef->fbgc2avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,diagOutFile[1],FILE_MODE_WRITE,&ef->fdfbgc2avg);CHKERRQ(ierr);

	ierr = VecDuplicate(TR,&ef->fbgc3);CHKERRQ(ierr);
	ierr = VecSet(ef->fbgc3,zero);CHKERRQ(ierr);
	ierr = VecGetArray(ef->fbgc3,&ef->localfbgc3);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&ef->fbgc3avg);CHKERRQ(ierr);
	ierr = VecSet(ef->fbgc3avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,diagOutFile[2],FILE_MODE_WRITE,&ef->fdfbgc3avg);CHKERRQ(ierr);

	ierr = VecDuplicate(TR,&ef->fbgc4);CHKERRQ(ierr);
	ierr = VecSet(ef->fbgc4,zero);CHKERRQ(ierr);
	ierr = VecGetArray(ef->fbgc4,&ef->localfbgc4);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&ef->fbgc4avg);CHKERRQ(ierr);
	ierr = VecSet(ef->fbgc4avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,diagOutFile[3],FILE_MODE_WRITE,&ef->fdfbgc4avg);CHKERRQ(ierr);

	ierr = VecDuplicate(TR,&ef->fbgc5);CHKERRQ(ierr);
	ierr = VecSet(ef->fbgc5,zero);CHKERRQ(ierr);
	ierr = VecGetArray(ef->fbgc5,&ef->localfbgc5);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&ef->fbgc5avg);CHKERRQ(ierr);
	ierr = VecSet(ef->fbgc5avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,diagOutFile[4],FILE_MODE_WRITE,&ef->fdfbgc5avg);CHKERRQ(ierr);

	ierr = VecDuplicate(TR,&ef->fbgc6);CHKERRQ(ierr);
	ierr = VecSet(ef->fbgc6,zero);CHKERRQ(ierr);
	ierr = VecGetArray(ef->fbgc6,&ef->localfbgc6);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&ef->fbgc6avg);CHKERRQ(ierr);
	ierr = VecSet(ef->fbgc6avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,diagOutFile[5],FILE_MODE_WRITE,&ef->fdfbgc6avg);CHKERRQ(ierr);

	ierr = VecDuplicate(TR,&ef->fbgc7);CHKERRQ(ierr);
	ierr = VecSet(ef->fbgc7,zero);CHKERRQ(ierr);
	ierr = VecGetArray(ef->fbgc7,&ef->localfbgc7);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&ef->fbgc7avg);CHKERRQ(ierr);
	ierr = VecSet(ef->fbgc7avg,zero);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,diagOutFile[6],FILE_MODE_WRITE,&ef->fdfbgc7avg);CHKERRQ(ierr);

#ifdef CARBON
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&ef->localco2airseafluxdiag);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&ef->localco2airseafluxdiagavg);CHKERRQ(ierr);  

    ierr = PetscOptionsGetString(NULL,prefix,"-co2airseaflux_file",ef->co2airseafluxFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
    if (!flg) {
      SETERRQ(PETSC_COMM_WORLD,1,"Must indicate file name for writing co2airsea flux with the -co2airseaflux_file option");
    }

    for (ip=0; ip<lNumProfiles; ip++) {
      ef->localco2airseafluxdiag[ip]=0.0;
      ef->localco2airseafluxdiagavg[ip]=0.0;
    }    
#endif

  }
  
  return 0;
}


/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/




#undef __FUNCT__
#define __FUNCT__ "calcExternalForcing"
PetscErrorCode calcExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *userctx)
{
  PetscInt numTracers;
  const char *prefix;
  void *ctx;

  PetscErrorCode ierr;
  PetscInt itr, ip, nzloc, kl;
  PetscScalar myTime;
  PetscBool debugLoc = PETSC_FALSE;

#ifdef CARBON
  PetscInt itf;
  PetscScalar alpha;
  PetscScalar DICemp = 0.0, ALKemp = 0.0;
  PetscScalar localco2airseaflux = 0.0;
#endif

  PetscScalar localFocean=0.0;
  PetscScalar localburial = 0.0;

  PetscCall(PetscContainerGetPointer(state->extforcctxcontainer, &ctx));
  ExternalForcingContext ef = (ExternalForcingContext)ctx;
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;
  
  myTime = DeltaT*Iter; /* Iter should start at 0 */

  if (ef->efctxId==1) {
	if (periodicBiogeochemForcing) {   
	  ierr = PeriodicVecInterp(tc,&Ts,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,Tsp,"Ts_");
	  ierr = PeriodicVecInterp(tc,&Ss,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,Ssp,"Ss_");	
#ifdef READ_SWRAD
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,localswradp,"swrad_");
#else
	 insolation_(&lNumProfiles,&myTime,&locallatitude[0],&daysPerYear,&localswrad[0],&localtau[0]);
#endif                                                 
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,localficep,"fice_");
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,localwindp,"wind_");   
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,localatmospp,"atmosp_");					                              
#ifdef CARBON
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,localEmPp,"EmP_");                                                  
#endif									              
	} else {
#ifndef READ_SWRAD
	  insolation_(&lNumProfiles,&myTime,&locallatitude[0],&daysPerYear,&localswrad[0],&localtau[0]);
#endif    
	}
  }

#ifdef CARBON
  if (ef->useAtmModel) {
  } else {  
/* Interpolate atmospheric pCO2   */
    if (!ef->fixedAtmosCO2) { 
      if (tc>=ef->TpCO2atm_hist[0]) {
        ierr = calcInterpFactor(ef->numpCO2atm_hist,tc,ef->TpCO2atm_hist,&itf,&alpha); CHKERRQ(ierr);
        ef->pCO2atm = alpha*ef->pCO2atm_hist[itf] + (1.0-alpha)*ef->pCO2atm_hist[itf+1];	  
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: time < %10.5f. Assuming pCO2atm=%g\n",ef->TpCO2atm_hist[0],ef->pCO2atm);CHKERRQ(ierr);
      }
    }  
  }

  if (useVirtualFlux) { /* use the global surface mean value to calculate E-P contribution */
    ierr = VecDot(surfVolFrac,DIC,&DICemp);CHKERRQ(ierr); /* volume weighted mean surface DIC */									              
    ierr = VecDot(surfVolFrac,ALK,&ALKemp);CHKERRQ(ierr); /* volume weighted mean surface ALK */									                  
  }

#endif
  
  for (ip=0; ip<lNumProfiles; ip++) {
    if ((debug) & (ip==ipdebug)) {
      debugLoc = PETSC_TRUE;
    } else {
      debugLoc = PETSC_FALSE;
    }
    nzloc=lProfileLength[ip];
    kl=lStartIndices[ip];

#ifdef CARBON
	if (!useVirtualFlux) { /* use the local surface value to calculate E-P contribution */
      DICemp=localDIC[kl];
      ALKemp=localALK[kl];      
    }
#endif

	for (itr=0; itr<numTracers; itr++) {    	
	  mops_biogeochem_copy_data_(&nzloc,&itr,&ef->localTR[itr][kl],&ef->localJTR[itr][kl],&DeltaT,&toModel);
	}  
    
    mops_biogeochem_model_(&nzloc,&DeltaT,
#ifdef CARBON			   
			   &DICemp,&ALKemp,&localEmP[ip],&ef->pCO2atm,
#endif
			   &localTs[kl],&localSs[kl],&localfice[ip],&localswrad[ip],&localtau[ip],&localwind[ip],&localatmosp[ip],&localdz[kl],
#ifdef CARBON			   
			   &ef->localph[ip],&localco2airseaflux,
#endif
                           &localburial,&ef->GRunoff,&localrunoffvol[kl],
			   &useSeparateBiogeochemTimeStepping);

	for (itr=0; itr<numTracers; itr++) {    
	  mops_biogeochem_copy_data_(&nzloc,&itr,&ef->localTR[itr][kl],&ef->localJTR[itr][kl],&DeltaT,&fromModel);
	}  


#ifdef CARBON			   
    if (ef->useAtmModel) {
      localFocean = localFocean + (localco2airseaflux/DeltaT)*localdA[ip]*(12.0/1.e18)*secondsPerYear; /* PgC/y */
    }
#endif

/* integrate burial in sediment over area and all profiles on each processor */
      ef->localFburial = ef->localFburial + localburial*localdA[ip]; 

	if (ef->calcDiagnostics) {  
	  if (Iter0+iLoop>=ef->diagTimer->startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */
        mops_biogeochem_diagnostics_(&nzloc,&ef->localfbgc1[kl],&ef->localfbgc2[kl],&ef->localfbgc3[kl],&ef->localfbgc4[kl],&ef->localfbgc5[kl],&ef->localfbgc6[kl],&ef->localfbgc7[kl]);
#ifdef CARBON        
        ef->localco2airseafluxdiag[ip]=localco2airseaflux;
#endif                
      }
	}

  } /* end loop over profiles */

#ifdef CARBON

  if (ef->useAtmModel) {
	ef->Focean = 0.0;
	 
	MPI_Allreduce(&localFocean, &ef->Focean, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);    
    
/*  time step atmosphere */
	ef->pCO2atm = ef->pCO2atm + ef->atmModelDeltaT*(-ef->Focean)/ppmToPgC;

/*  reset values */
	ef->Foceanint = ef->Foceanint + ef->atmModelDeltaT*ef->Focean; /* calculate the time integrated flux */
        
  }  

#endif  

/* sum burial in sediment over all processors, and scale by time step etc.*/
/* do this only once every burialSumSteps , and then take this value for next year's runoff */

    if ((iLoop % ef->burialSumSteps)==0) {

      ef->Fburial = 0.0;

      MPI_Allreduce(&ef->localFburial, &ef->Fburial, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);    

#ifdef RUNOFF
      ef->GRunoff = ef->Fburial/(1.e12*ef->burialSumSteps)*(86400.0/DeltaT); /* This is Gmol P/day. 
      Note: localrunoff is scaled with 1e12. Note: GRunoff will be scaled with bgc_dt.*/
#else
      ef->GRunoff = ef->Fburial/(totalA*ef->burialSumSteps)*(86400.0/DeltaT); /* This is mmol P/m2/day. 
      Note: this will later be divided by depth of first layer. Note: GRunoff will be scaled with bgc_dt.*/ 
#endif

      ef->localFburial = 0.0;
      }
    


  if (!useSeparateBiogeochemTimeStepping) { 
/*  Convert to discrete tendency */
	for (itr=0; itr<numTracers; itr++) {
	  ierr = VecScale(state->qef[itr],DeltaT);CHKERRQ(ierr);
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
PetscErrorCode writeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *userctx)
{
  PetscInt numTracers;
  const char *prefix;
  void *ctx;

  PetscErrorCode ierr;
  PetscInt ip;
  PetscScalar zero = 0.0, one = 1.0;  

  PetscCall(PetscContainerGetPointer(state->extforcctxcontainer, &ctx));
  ExternalForcingContext ef = (ExternalForcingContext)ctx;
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

/* Note: tc and iLoop are the time and step at the end of the current time step. */

#ifdef CARBON
  if (ef->useAtmModel) {
/* write instantaneous atmos model state */
	if (Iter0+iLoop>=(ef->atmWriteTimer->startTimeStep)) { /* note: startTimeStep is ABSOLUTE time step */
	  if ((ef->atmWriteTimer->count)<=(ef->atmWriteTimer->numTimeSteps)) {
		ef->atmWriteTimer->count++;
	  }
	  if ((ef->atmWriteTimer->count)==(ef->atmWriteTimer->numTimeSteps)) { /* time to write out */
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing atmospheric model output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
		ierr = PetscFPrintf(PETSC_COMM_WORLD,ef->atmfptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);           
		ierr = writeBinaryScalarData("pCO2atm_output.bin",&ef->pCO2atm,1,PETSC_TRUE);
		ierr = writeBinaryScalarData("Foceanint_output.bin",&ef->Foceanint,1,PETSC_TRUE);
		ef->Foceanint = 0.0;

		ierr = StepTimerUpdate(Iter0+iLoop, ef->atmWriteTimer);CHKERRQ(ierr);

      }
    }
  }
#endif


#ifdef WRITE_RUNOFF
    if ((iLoop % ef->burialSumSteps)==0) {  /*  time to write out */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing runoff output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
      ierr = PetscFPrintf(PETSC_COMM_WORLD,ef->runofffptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);           
      ierr = writeBinaryScalarData("Grunoff_output.bin",&ef->GRunoff,1,PETSC_TRUE);
    }
#endif

  if (ef->calcDiagnostics) {  
	if (Iter0+iLoop>=(ef->diagTimer->startTimeStep)) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */  
  
	  if (ef->diagTimer->count<=ef->diagTimer->numTimeSteps) { /* still within same averaging block so accumulate */
		ierr = VecAXPY(ef->fbgc1avg,one,ef->fbgc1);CHKERRQ(ierr);
		ierr = VecAXPY(ef->fbgc2avg,one,ef->fbgc2);CHKERRQ(ierr);
		ierr = VecAXPY(ef->fbgc3avg,one,ef->fbgc3);CHKERRQ(ierr);
		ierr = VecAXPY(ef->fbgc4avg,one,ef->fbgc4);CHKERRQ(ierr);
		ierr = VecAXPY(ef->fbgc5avg,one,ef->fbgc5);CHKERRQ(ierr);
		ierr = VecAXPY(ef->fbgc6avg,one,ef->fbgc6);CHKERRQ(ierr);
		ierr = VecAXPY(ef->fbgc7avg,one,ef->fbgc7);CHKERRQ(ierr);

#ifdef CARBON
        for (ip=0; ip<lNumProfiles; ip++) {
          ef->localco2airseafluxdiagavg[ip]=ef->localco2airseafluxdiag[ip]+ef->localco2airseafluxdiagavg[ip];      
        }	  
#endif

		ef->diagTimer->count++;
	  }
	  if ((ef->diagTimer->count)==(ef->diagTimer->numTimeSteps)) { /* time to write averages to file */
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing diagnostics time average at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);                      

		ierr = VecScale(ef->fbgc1avg,1.0/ef->diagTimer->count);CHKERRQ(ierr);
		ierr = VecView(ef->fbgc1avg,ef->fdfbgc1avg);CHKERRQ(ierr);
		ierr = VecSet(ef->fbgc1avg,zero); CHKERRQ(ierr);

		ierr = VecScale(ef->fbgc2avg,1.0/ef->diagTimer->count);CHKERRQ(ierr);
		ierr = VecView(ef->fbgc2avg,ef->fdfbgc2avg);CHKERRQ(ierr);
		ierr = VecSet(ef->fbgc2avg,zero); CHKERRQ(ierr);

		ierr = VecScale(ef->fbgc3avg,1.0/ef->diagTimer->count);CHKERRQ(ierr);
		ierr = VecView(ef->fbgc3avg,ef->fdfbgc3avg);CHKERRQ(ierr);
		ierr = VecSet(ef->fbgc3avg,zero); CHKERRQ(ierr);

		ierr = VecScale(ef->fbgc4avg,1.0/ef->diagTimer->count);CHKERRQ(ierr);
		ierr = VecView(ef->fbgc4avg,ef->fdfbgc4avg);CHKERRQ(ierr);
		ierr = VecSet(ef->fbgc4avg,zero); CHKERRQ(ierr);

		ierr = VecScale(ef->fbgc5avg,1.0/ef->diagTimer->count);CHKERRQ(ierr);
		ierr = VecView(ef->fbgc5avg,ef->fdfbgc5avg);CHKERRQ(ierr);
		ierr = VecSet(ef->fbgc5avg,zero); CHKERRQ(ierr);

		ierr = VecScale(ef->fbgc6avg,1.0/ef->diagTimer->count);CHKERRQ(ierr);
		ierr = VecView(ef->fbgc6avg,ef->fdfbgc6avg);CHKERRQ(ierr);
		ierr = VecSet(ef->fbgc6avg,zero); CHKERRQ(ierr);

		ierr = VecScale(ef->fbgc7avg,1.0/ef->diagTimer->count);CHKERRQ(ierr);
		ierr = VecView(ef->fbgc7avg,ef->fdfbgc7avg);CHKERRQ(ierr);
		ierr = VecSet(ef->fbgc7avg,zero); CHKERRQ(ierr);

#ifdef CARBON
        for (ip=0; ip<lNumProfiles; ip++) {
          ef->localco2airseafluxdiagavg[ip]=ef->localco2airseafluxdiagavg[ip]/ef->diagTimer->count;
        }	  

        ierr = writeProfileSurfaceScalarData(ef->co2airseafluxFile,ef->localco2airseafluxdiagavg,1,ef->appendDiagnostics);  		

/*      reset diagnostic arrays */
        for (ip=0; ip<lNumProfiles; ip++) {
          ef->localco2airseafluxdiagavg[ip]=0.0;
        }	  
#endif

        ef->appendDiagnostics=PETSC_TRUE;
        ierr = StepTimerUpdate(Iter0+iLoop, ef->diagTimer);CHKERRQ(ierr);

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
PetscErrorCode finalizeExternalForcing(PetscScalar tc, PetscInt Iter, TMMState state, void *userctx)
{
  PetscInt numTracers;
  const char *prefix;  
  void *ctx;

  PetscErrorCode ierr;
  
  PetscCall(PetscContainerGetPointer(state->extforcctxcontainer, &ctx));
  ExternalForcingContext ef = (ExternalForcingContext)ctx;
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

/* write final pickup */
#ifdef CARBON
  if (ef->useAtmModel) {
/* write instantaneous atmos model state */
    ierr = writeBinaryScalarData("pickup_pCO2atm.bin",&ef->pCO2atm,1,PETSC_FALSE);
  }
#endif

  ierr = writeBinaryScalarData(ef->runoffIniOutFile,&ef->GRunoff,1,PETSC_FALSE);

  if (ef->efctxId==1) {  
	ierr = VecRestoreArray(Ts,&localTs);CHKERRQ(ierr);
	ierr = VecRestoreArray(Ss,&localSs);CHKERRQ(ierr);
  
	ierr = VecDestroy(&Ts);CHKERRQ(ierr);
	ierr = VecDestroy(&Ss);CHKERRQ(ierr);

	if (periodicBiogeochemForcing) {    
	  ierr = PeriodicVecDestroy(&Tsp);CHKERRQ(ierr);
	  ierr = PeriodicVecDestroy(&Ssp);CHKERRQ(ierr);
	  ierr = PeriodicArrayDestroy(&localficep);CHKERRQ(ierr);
	  ierr = PeriodicArrayDestroy(&localwindp);CHKERRQ(ierr);    
	  ierr = PeriodicArrayDestroy(&localatmospp);CHKERRQ(ierr);    
#ifdef READ_SWRAD    
	  ierr = PeriodicArrayDestroy(&localswradp);CHKERRQ(ierr);
#endif
#ifdef CARBON
	  ierr = PeriodicArrayDestroy(&localEmPp);CHKERRQ(ierr);
#endif	
	}    

#ifdef CARBON
	if (useVirtualFlux) {
	  ierr = VecDestroy(&surfVolFrac);CHKERRQ(ierr);  
	}
#endif
  }
  
  if (ef->calcDiagnostics) {  
	ierr = VecDestroy(&ef->fbgc1);CHKERRQ(ierr);
	ierr = VecDestroy(&ef->fbgc1avg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&ef->fdfbgc1avg);CHKERRQ(ierr);	

	ierr = VecDestroy(&ef->fbgc2);CHKERRQ(ierr);
	ierr = VecDestroy(&ef->fbgc2avg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&ef->fdfbgc2avg);CHKERRQ(ierr);	

	ierr = VecDestroy(&ef->fbgc3);CHKERRQ(ierr);
	ierr = VecDestroy(&ef->fbgc3avg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&ef->fdfbgc3avg);CHKERRQ(ierr);	

	ierr = VecDestroy(&ef->fbgc4);CHKERRQ(ierr);
	ierr = VecDestroy(&ef->fbgc4avg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&ef->fdfbgc4avg);CHKERRQ(ierr);	

	ierr = VecDestroy(&ef->fbgc5);CHKERRQ(ierr);
	ierr = VecDestroy(&ef->fbgc5avg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&ef->fdfbgc5avg);CHKERRQ(ierr);	

	ierr = VecDestroy(&ef->fbgc6);CHKERRQ(ierr);
	ierr = VecDestroy(&ef->fbgc6avg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&ef->fdfbgc6avg);CHKERRQ(ierr);	

	ierr = VecDestroy(&ef->fbgc7);CHKERRQ(ierr);
	ierr = VecDestroy(&ef->fbgc7avg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&ef->fdfbgc7avg);CHKERRQ(ierr);	

  }

#ifdef CARBON
  if (ef->useAtmModel) {
    ierr = PetscFClose(PETSC_COMM_WORLD,ef->atmfptime);CHKERRQ(ierr);
  }
#endif

#ifdef WRITE_RUNOFF
    ierr = PetscFClose(PETSC_COMM_WORLD,ef->runofffptime);CHKERRQ(ierr);
#endif

  return 0;
}





/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/



#undef __FUNCT__
#define __FUNCT__ "reInitializeExternalForcing"
PetscErrorCode reInitializeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *userctx)
{
  PetscInt numTracers;
  const char *prefix;  
  void *ctx;

  PetscErrorCode ierr;
  PetscInt ip, kl, nzloc;
  PetscInt itr;
  PetscScalar myTime;
  PetscViewer fd;
  PetscInt fp;
  PetscInt ipar;
  
  PetscCall(PetscContainerGetPointer(state->extforcctxcontainer, &ctx));
  ExternalForcingContext ef = (ExternalForcingContext)ctx;
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

  myTime = DeltaT*Iter; /* Iter should start at 0 */

  if (ef->efctxId==1) {
	if (periodicBiogeochemForcing) {   
	  ierr = PeriodicVecInterp(tc,&Ts,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,Tsp,"Ts_");
	  ierr = PeriodicVecInterp(tc,&Ss,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,Ssp,"Ss_");	
#ifdef READ_SWRAD
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,localswradp,"swrad_");
#else
	 insolation_(&lNumProfiles,&myTime,&locallatitude[0],&daysPerYear,&localswrad[0],&localtau[0]);
#endif                                                 
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,localficep,"fice_");
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,localwindp,"wind_");   
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,localatmospp,"atmosp_");					                              
#ifdef CARBON
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,localEmPp,"EmP_");                                                  
#endif									              
	} else {
#ifndef READ_SWRAD
	  insolation_(&lNumProfiles,&myTime,&locallatitude[0],&daysPerYear,&localswrad[0],&localtau[0]);
#endif    
	}
  }

//   if (readBGCParams) {
// 
//   	fpnames=fopen(bgcNamesFile,"r");
//     for(ipar=0;ipar<numBGCParams;ipar++) {
//       ierr = fscanf(fpnames,"%s\n",&bgcnames[ipar]);
//     }   
//     fclose(fpnames);
// 
// #ifdef ASCIIPARAMS
//     fpparams=fopen(bgcParamsFile,"r");
//     for(ipar=0;ipar<numBGCParams;ipar++) {
//       ierr = fscanf(fpparams,"%lf\n",&bgcparams[ipar]);
//     }   
//     fclose(fpparams);
// #else
//     ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,bgcParamsFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
//     ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
//     ierr = PetscMalloc(numBGCParams*sizeof(PetscScalar),&bgcparams);CHKERRQ(ierr); 
//     ierr = PetscBinaryRead(fp,bgcparams,numBGCParams,NULL,PETSC_SCALAR);CHKERRQ(ierr);  
//     ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
// #endif
//  
//     for (ipar=0; ipar<numBGCParams; ipar++) {
//     ierr = PetscPrintf(PETSC_COMM_WORLD,"Parameter no. %d is %f\n",ipar,bgcparams[ipar]);CHKERRQ(ierr);
//     }
// 
//     /* Call mops_biogeochem_set_params.F for every parameter */
// 	for (ipar=0; ipar<numBGCParams; ipar++) {
// 		lenBGCParam = strlen(bgcnames[ipar]);
//     	mops_biogeochem_set_params_(&lenBGCParam,&bgcparams[ipar],&bgcnames[ipar]);
//     }

    for (ip=0; ip<lNumProfiles; ip++) {
      nzloc=lProfileLength[ip];
      kl=lStartIndices[ip];
      for (itr=0; itr<numTracers; itr++) {    	
        mops_biogeochem_copy_data_(&nzloc,&itr,&ef->localTR[itr][kl],&ef->localJTR[itr][kl],&DeltaT,&toModel);
      }  
      mops_biogeochem_ini_(&nzloc,&DeltaT,
#ifdef CARBON
			   &ef->localph[ip],
#endif
			   &localTs[kl],&localSs[kl],&localdz[kl],&drF[0],&nzmax,&nzeuph,
                           &numBiogeochemStepsPerOceanStep,&MYFALSE);
    }    
//   }
    
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "getExternalForcingContext"
PetscErrorCode getExternalForcingContext(TMMState state, void **ctxret)
{
  void *ctx;
  PetscErrorCode ierr;
//   ExternalForcingContext ef;
  
//   ierr = PetscPrintf(PETSC_COMM_WORLD,"Calling getForcingContext");CHKERRQ(ierr);
  
  PetscCall(PetscContainerGetPointer(state->extforcctxcontainer, &ctx));
//   ef = (ExternalForcingContext)ctx;
//   ierr = PetscPrintf(PETSC_COMM_WORLD,"Runoff is %f\n",ef->GRunoff);CHKERRQ(ierr);
  *ctxret=ctx;
  
  return 0;
}
