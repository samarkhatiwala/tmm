static char help[] = "\n";
/* 
  Include "petscmat.h" so that we can use matrices.
  automatically includes:
     petsc.h       - base PETSc routines   petscvec.h    - vectors
     petscsys.h    - system routines       petscmat.h    - matrices
     petscis.h     - index sets            petscviewer.h - viewers               
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "petsctime.h"

#include "petsc_matvec_utils.h"
#include "tmm_forcing_utils.h"
#include "tmm_external_forcing.h"
#include "tmm_external_bc.h"
#include "tmm_timer.h"
#include "tmm_profile_utils.h"
#include "tmm_profile_data.h"
#include "tmm_monitor.h"
#include "tmm_misfit.h"
#include "tmm_main.h"

PetscScalar deltaTClock, time0;
PetscInt maxSteps, Iter0;
StepTimer writeTimer;
PetscBool appendOutput = PETSC_FALSE;
PetscInt *gIndices, gLow, gHigh;
PetscInt *gBCIndices, lBCSize, gBCSize, gbcLow, gbcHigh;
PetscBool doMisfit = PETSC_FALSE;
PetscBool rescaleForcing = PETSC_FALSE;

extern PetscErrorCode forwardStep(PetscScalar tc, PetscInt iLoop, PetscScalar dt, PetscInt numTracers, 
                                  PetscBool useForcingFromFile, PetscBool useExternalForcing, PetscBool usePrescribedBC, 
                                  Vec *v, Mat Ae, Mat Ai, Mat Be, Mat Bi, Vec *uf, Vec *uef, Vec *bcc, Vec *bcf, Vec *vtmp);

extern PetscErrorCode iniTMMWrite(PetscScalar tc, PetscInt it, PetscInt ntr, Vec *v);
extern PetscErrorCode doTMMWrite(PetscScalar tc, PetscInt it, PetscInt ntr, Vec *v);
extern PetscErrorCode finalizeTMMWrite(PetscScalar tc, PetscInt it, PetscInt ntr);

extern PetscErrorCode waitForSignal(PetscInt waitTime);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  PetscInt numTracers, n;
  Vec templateVec;
  Vec *v, *vtmp;
/* TM's */
  Mat Ae, Ai;
  PeriodicMat Aep, Aip;
  TimeDependentMat Aetd, Aitd;
  char mateFile[PETSC_MAX_PATH_LEN], matiFile[PETSC_MAX_PATH_LEN], rfsFile[PETSC_MAX_PATH_LEN];
  PetscBool periodicMatrix = PETSC_FALSE;
  PetscBool timeDependentMatrix = PETSC_FALSE;
  PeriodicTimer matrixPeriodicTimer;
  TimeDependentTimer matrixTimeDependentTimer;
  
/* Forcing */
  Vec *uf, *uef;
  PeriodicVec up[MAXNUMTRACERS];
  char *forcingFile[MAXNUMTRACERS];
  PeriodicTimer forcingTimer;
  TimeDependentVec utdf[MAXNUMTRACERS];
  TimeDependentTimer forcingTimeDependentTimer;
  PetscInt forcingFromFileCutOffStep = -1;
  PetscInt externalForcingCutOffStep = -1;

/* Rescale forcing */
  Vec Rfs;
  PeriodicVec Rfsp;
  TimeDependentVec Rfstd;
  
/* BCs */
  Vec *bcc, *bcf;
  PeriodicVec bcp[MAXNUMTRACERS];
  char *bcFile[MAXNUMTRACERS];
  PeriodicTimer bcTimer;
  TimeDependentVec bctd[MAXNUMTRACERS];
  TimeDependentTimer bcTimeDependentTimer;          
  PetscInt bcCutOffStep = -1;
  Mat Be, Bi;
  PeriodicMat Bep, Bip;
  TimeDependentMat Betd, Bitd;  
  char matbeFile[PETSC_MAX_PATH_LEN], matbiFile[PETSC_MAX_PATH_LEN];
  Vec bcTemplateVec;
  
/* I/O   */
  char *iniFile[MAXNUMTRACERS];  
  char *outFile[MAXNUMTRACERS];  
  char *bcoutFile[MAXNUMTRACERS];
  char *ufoutFile[MAXNUMTRACERS], *uefoutFile[MAXNUMTRACERS];
  char pickupFile[PETSC_MAX_PATH_LEN];
  char pickupoutFile[PETSC_MAX_PATH_LEN];  
  PetscBool writePickup = PETSC_FALSE;
  StepTimer pickupTimer;
  char outTimeFile[PETSC_MAX_PATH_LEN];  
  PetscFileMode OUTPUT_FILE_MODE;
  PetscBool doWriteBC = PETSC_FALSE;
  PetscBool doWriteUF = PETSC_FALSE;
  PetscBool doWriteUEF = PETSC_FALSE;
  PetscBool pickupFromFile = PETSC_FALSE;
  PetscBool doTimeAverage = PETSC_FALSE;
  StepTimer avgTimer;
  char *avgOutFile[MAXNUMTRACERS];  
  Vec *vavg;
  PetscViewer fdavgout[MAXNUMTRACERS];
  PetscBool avgAppendOutput = PETSC_FALSE;
  PetscFileMode AVG_FILE_MODE;
  FILE *avgfptime;
  char avgOutTimeFile[PETSC_MAX_PATH_LEN];  
  char *bcavgOutFile[MAXNUMTRACERS];
  Vec *bcavg;
  PetscViewer fdbcavgout[MAXNUMTRACERS];
  char *ufavgOutFile[MAXNUMTRACERS], *uefavgOutFile[MAXNUMTRACERS];
  Vec *ufavg, *uefavg;
  PetscViewer fdufavgout[MAXNUMTRACERS], fduefavgout[MAXNUMTRACERS];
  FILE *fptime;
  PetscViewer fd, fdp, fdout[MAXNUMTRACERS];
  PetscViewer fdbcout[MAXNUMTRACERS], fdufout[MAXNUMTRACERS], fduefout[MAXNUMTRACERS];

#if defined (FORSPINUP) || defined (FORJACOBIAN)
  PetscViewer fdin[MAXNUMTRACERS];
  PetscInt itjac;
  int fp;  
#endif
  
/* run time options */
  PetscBool useExternalForcing = PETSC_FALSE;
  PetscBool useForcingFromFile = PETSC_FALSE;
  PetscBool usePrescribedBC = PETSC_FALSE;
  PetscBool applyExternalForcing = PETSC_FALSE;
  PetscBool applyForcingFromFile = PETSC_FALSE;
  PetscBool applyBC = PETSC_FALSE;
  PetscBool periodicForcing = PETSC_FALSE;
  PetscBool timeDependentForcing = PETSC_FALSE;
  PetscBool constantForcing = PETSC_FALSE;
  PetscBool periodicBC = PETSC_FALSE;
  PetscBool timeDependentBC = PETSC_FALSE;
  PetscBool constantBC = PETSC_FALSE;
  PetscBool doCalcBC = PETSC_FALSE;
  PetscBool useMonitor = PETSC_FALSE;

  PetscMPIInt numProcessors, myId;  
  PetscErrorCode ierr;
  PetscBool flg1,flg2;
  PetscScalar t1, t2, tc, tf;
  PetscInt iLoop, Iterc;
  PetscInt itr, maxValsToRead;
  char tmpFile[PETSC_MAX_PATH_LEN];
  PetscScalar zero = 0.0, one = 1.0;
  PetscInt il;
  
  PetscInitialize(&argc,&args,(char *)0,help);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myId);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcessors);CHKERRQ(ierr);  
  myId=myId+1; /* process ID (starting at 1) */

  PetscPushErrorHandler(PetscAbortErrorHandler,NULL); /* force code to abort on error */

/* Some defaults */
  time0=0.0;
  deltaTClock=1.0;
  Iter0=0;  
 
  numTracers=1;
  
/* Process options and load files */  
/* Number of tracers */
  ierr = PetscOptionsGetInt(NULL,NULL,"-numtracers",&numTracers,&flg1);CHKERRQ(ierr);
  if (numTracers>MAXNUMTRACERS) {
   SETERRQ(PETSC_COMM_WORLD,1,"Number of tracers exceeds maximum allowable. Please increase the variable MAXNUMTRACERS and recompile");
  }  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of tracers to be integrated: %d\n", numTracers);CHKERRQ(ierr); 

/* Time step data */
  ierr = PetscOptionsGetReal(NULL,NULL,"-deltat_clock",&deltaTClock,&flg1);CHKERRQ(ierr);

  ierr = PetscOptionsGetReal(NULL,NULL,"-t0",&time0,&flg1);CHKERRQ(ierr);  
  ierr = PetscOptionsGetInt(NULL,NULL,"-iter0",&Iter0,&flg2);CHKERRQ(ierr);
  if ((flg1) && (!flg2)) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate both or neither time0 and Iter0 with the -t0 and -iter0 flags");
  if ((!flg1) && (flg2)) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate both or neither time0 and Iter0 with the -t0 and -iter0 flags");        
  ierr = PetscOptionsGetInt(NULL,NULL,"-max_steps",&maxSteps,&flg1);CHKERRQ(ierr);
  if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate maximum number of steps with the -max_steps option");

// Catch deprecated option
  PetscInt dum;
  ierr = PetscOptionsGetInt(NULL,NULL,"-write_steps",&dum,&flg1);CHKERRQ(ierr);
  if (flg1) SETERRQ(PETSC_COMM_WORLD,1,"ERROR!: The -write_steps option has been deprecated. Use the StepTimer object with prefix 'write'");

  ierr = iniStepTimer("write_", Iter0, &writeTimer);CHKERRQ(ierr);

/*Data for time averaging */
  ierr = PetscOptionsHasName(NULL,NULL,"-time_avg",&doTimeAverage);CHKERRQ(ierr);
  if (doTimeAverage) {
    ierr = iniStepTimer("avg_", Iter0+1, &avgTimer);CHKERRQ(ierr);
	for (itr=0; itr<numTracers; itr++) {
	  avgOutFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
	}
	maxValsToRead = numTracers;
	ierr = PetscOptionsGetStringArray(NULL,NULL,"-avg_files",avgOutFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
	if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate file name(s) for writing time averages with the -avg_files option");
	if (maxValsToRead != numTracers) {
	  SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of time average file names specified");
	}  
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Time averages will be computed starting at and including (absolute) time step: %d\n", avgTimer.startTimeStep);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Time averages will be computed over %d time steps\n", avgTimer.numTimeSteps);CHKERRQ(ierr);	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Time averages will be written to:\n");CHKERRQ(ierr);
	for (itr=0; itr<numTracers; itr++) {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,avgOutFile[itr]);CHKERRQ(ierr);
	}

	ierr = PetscOptionsHasName(NULL,NULL,"-avg_append",&avgAppendOutput);CHKERRQ(ierr);
	if (avgAppendOutput) {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Time averages will be appended\n");CHKERRQ(ierr);
	  AVG_FILE_MODE=FILE_MODE_APPEND;
	} else {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Time averages will overwrite existing file(s)\n");CHKERRQ(ierr);
	  AVG_FILE_MODE=FILE_MODE_WRITE;
	}

/* Output times */
	ierr = PetscOptionsGetString(NULL,NULL,"-avg_time_file",avgOutTimeFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
	if (!flg1) {
	  strcpy(avgOutTimeFile,"");
	  sprintf(avgOutTimeFile,"%s","time_average_output_time.txt");
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Time average output times will be written to %s\n",avgOutTimeFile);CHKERRQ(ierr);

	if (!avgAppendOutput) {
	  ierr = PetscFOpen(PETSC_COMM_WORLD,avgOutTimeFile,"w",&avgfptime);CHKERRQ(ierr);  
	} else {
	  ierr = PetscFOpen(PETSC_COMM_WORLD,avgOutTimeFile,"a",&avgfptime);CHKERRQ(ierr);  
	}
	
  }
  
/* Initialize profile data and create template vector */
  ierr = iniProfileData(myId);CHKERRQ(ierr);
  if (useProfiles) {
    ierr = VecCreate(PETSC_COMM_WORLD,&templateVec);CHKERRQ(ierr);
    ierr = VecSetSizes(templateVec,lSize,PETSC_DECIDE);CHKERRQ(ierr);
    ierr = VecSetFromOptions(templateVec);CHKERRQ(ierr);
    ierr = VecGetSize(templateVec,&n);CHKERRQ(ierr);
  }

#if defined (FORSPINUP) || defined (FORJACOBIAN)
  ierr = waitForSignal(0);CHKERRQ(ierr); /* initialize */
#endif

/* Matrices */
  ierr = PetscOptionsGetString(NULL,NULL,"-me",mateFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
  if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate binary matrix file with the -me option");
  ierr = PetscOptionsGetString(NULL,NULL,"-mi",matiFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
  if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate binary matrix file with the -mi options");

  ierr = MatCreate(PETSC_COMM_WORLD,&Ae);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&Ai);CHKERRQ(ierr);

/*  Set layout information */  
  if ((useProfiles) && (numProcessors>1)) {
	ierr = MatSetSizes(Ae,lSize,lSize,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
	ierr = MatSetSizes(Ai,lSize,lSize,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
  }   
  
  ierr = MatSetType(Ae,MATMPIAIJ);CHKERRQ(ierr);      
  ierr = MatSetFromOptions(Ae);CHKERRQ(ierr);  
  ierr = MatSetType(Ai,MATMPIAIJ);CHKERRQ(ierr);        
  ierr = MatSetFromOptions(Ai);CHKERRQ(ierr);
  
  ierr = PetscOptionsHasName(NULL,NULL,"-periodic_matrix",&periodicMatrix);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-time_dependent_matrix",&timeDependentMatrix);CHKERRQ(ierr);
  if (periodicMatrix & timeDependentMatrix) {
    SETERRQ(PETSC_COMM_WORLD,1,"Cannot use both -periodic_matrix and -time_dependent_matrix together!");
  }
  
  if (periodicMatrix) {    
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic matrices specified\n");CHKERRQ(ierr);
    ierr=PetscStrcat(mateFile,"_");CHKERRQ(ierr);
    ierr=PetscStrcat(matiFile,"_");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Ae basename is %s\n", mateFile);CHKERRQ(ierr); 
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Ai basename is %s\n", matiFile);CHKERRQ(ierr);     
    
/*  read time data */
    ierr = iniPeriodicTimer("matrix_", &matrixPeriodicTimer);CHKERRQ(ierr);
    
/*  Read here to set sparsity pattern */
	strcpy(tmpFile,"");
	sprintf(tmpFile,"%s%02d",mateFile,0);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatLoad(Ae,fd);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

	strcpy(tmpFile,"");
	sprintf(tmpFile,"%s%02d",matiFile,0);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatLoad(Ai,fd);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    

    Aep.firstTime = PETSC_TRUE;
    Aip.firstTime = PETSC_TRUE;

  } else if (timeDependentMatrix) {
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Time-dependent matrices specified\n");CHKERRQ(ierr);
    ierr=PetscStrcat(mateFile,"_");CHKERRQ(ierr);
    ierr=PetscStrcat(matiFile,"_");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Ae basename is %s\n", mateFile);CHKERRQ(ierr); 
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Ai basename is %s\n", matiFile);CHKERRQ(ierr);     
    
/*  read time data */
    ierr = iniTimeDependentTimer("matrix_", &matrixTimeDependentTimer);CHKERRQ(ierr);
    
/*  Read here to set sparsity pattern */
	strcpy(tmpFile,"");
	sprintf(tmpFile,"%s%02d",mateFile,0);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatLoad(Ae,fd);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

	strcpy(tmpFile,"");
	sprintf(tmpFile,"%s%02d",matiFile,0);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatLoad(Ai,fd);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    

    Aetd.firstTime = PETSC_TRUE;
    Aitd.firstTime = PETSC_TRUE;  
    
  } else { /*  not periodic. read matrices here */

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,mateFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading Ae from file %s\n", mateFile);CHKERRQ(ierr);  
	ierr = MatLoad(Ae,fd);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,matiFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading Ai from file %s\n", matiFile);CHKERRQ(ierr);  
	ierr = MatLoad(Ai,fd);CHKERRQ(ierr);    
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
    
  }

/* create template vector here if not using profiles */
  if (!useProfiles) {
    ierr = MatGetSize(Ae,0,&n);CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD,&templateVec);CHKERRQ(ierr);
    ierr = VecSetSizes(templateVec,PETSC_DECIDE,n);CHKERRQ(ierr);
    ierr = VecSetFromOptions(templateVec);CHKERRQ(ierr);
    ierr = VecGetLocalSize(templateVec,&lSize);CHKERRQ(ierr); /* lSize is otherwise set in iniProfiles */
  }  
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix size is %d x %d\n", n,n);CHKERRQ(ierr);  

/*   Compute global indices for local piece of vectors */
  ierr = VecGetOwnershipRange(templateVec,&gLow,&gHigh);CHKERRQ(ierr);
  gHigh = gHigh - 1; /* Note: gHigh is one more than the last local element */
  ierr = PetscMalloc(lSize*sizeof(PetscInt),&gIndices);CHKERRQ(ierr);  
  for (il=0; il<lSize; il++) {
    gIndices[il] = il + gLow;
  }  

/* Output file */
  for (itr=0; itr<numTracers; itr++) {
    outFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
  }
  maxValsToRead = numTracers;
  ierr = PetscOptionsGetStringArray(NULL,NULL,"-o",outFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
  if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate output file name(s) with the -o option");
  if (maxValsToRead != numTracers) {
    SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of output file names specified");
  }  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Output will be written to:\n");CHKERRQ(ierr);
  for (itr=0; itr<numTracers; itr++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,outFile[itr]);CHKERRQ(ierr);
  }  
  ierr = PetscOptionsHasName(NULL,NULL,"-append",&appendOutput);CHKERRQ(ierr);
  if (appendOutput) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Output will be appended\n");CHKERRQ(ierr);
    OUTPUT_FILE_MODE=FILE_MODE_APPEND;
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Output will overwrite existing file(s)\n");CHKERRQ(ierr);
    OUTPUT_FILE_MODE=FILE_MODE_WRITE;
  }    

/* Output times */
  ierr = PetscOptionsGetString(NULL,NULL,"-time_file",outTimeFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
  if (!flg1) {
	strcpy(outTimeFile,"");
    sprintf(outTimeFile,"%s","output_time.txt");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Output times will be written to %s\n",outTimeFile);CHKERRQ(ierr);

/* File name for final pickup */
  ierr = PetscOptionsGetString(NULL,NULL,"-pickup_out",pickupoutFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
  if (!flg1) {
	strcpy(pickupoutFile,"");
    sprintf(pickupoutFile,"%s","pickup.petsc");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Final pickup will be written to %s\n",pickupoutFile);CHKERRQ(ierr);

  ierr = PetscOptionsHasName(NULL,NULL,"-pickup_time_steps",&writePickup);CHKERRQ(ierr);
  if (writePickup) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Intermediate pickups will be written to pickup_ITERATIONNUMBER.petsc\n");CHKERRQ(ierr);
    ierr = iniStepTimer("pickup_", Iter0, &pickupTimer);CHKERRQ(ierr);
  }

/* tracer vectors */
  ierr = VecDuplicateVecs(templateVec,numTracers,&v);CHKERRQ(ierr);
  ierr = VecDuplicateVecs(templateVec,numTracers,&vtmp);CHKERRQ(ierr); /* temporary work space needed by forwardStep */

/* Initial condition     */
  for (itr=0; itr<numTracers; itr++) {
    iniFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
  }
  ierr = PetscOptionsGetString(NULL,NULL,"-pickup",pickupFile,PETSC_MAX_PATH_LEN-1,&pickupFromFile);CHKERRQ(ierr);
  if (pickupFromFile) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Pickup file has been specified\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  Reading initial conditions from %s\n", pickupFile);CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,pickupFile,FILE_MODE_READ,&fdp);CHKERRQ(ierr);
    for (itr=0; itr<numTracers; itr++) {
      ierr = VecLoad(v[itr],fdp);CHKERRQ(ierr); /* IntoVector */
    }
    ierr = PetscViewerDestroy(&fdp);CHKERRQ(ierr);          
  } else {
    maxValsToRead = numTracers;
    ierr = PetscOptionsGetStringArray(NULL,NULL,"-i",iniFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
    if (flg1) {  /* read from file */
      if (maxValsToRead != numTracers) {
        SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of input file names specified");
      }      
#if defined (FORSPINUP) || defined (FORJACOBIAN)
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Waiting for new initial condition ...\n");CHKERRQ(ierr);
      ierr = waitForSignal(10);CHKERRQ(ierr);
#endif	  
      for (itr=0; itr<numTracers; itr++) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading initial condition from file %s\n", itr,iniFile[itr]);CHKERRQ(ierr);
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,iniFile[itr],FILE_MODE_READ,&fd);CHKERRQ(ierr);
        ierr = VecLoad(v[itr],fd);CHKERRQ(ierr);/* IntoVector */
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);          
      }  
    } else {  /* set to zero */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting initial condition to zero\n");CHKERRQ(ierr);
      for (itr=0; itr<numTracers; itr++) {
        VecSet(v[itr],zero);
      }        
    }
  }

/* Forcing/RHS   */
/* The tracer(s) can be forced in 3 ways (any combination of which can be turned on): */
/* 1) Forcing term read from file (can be periodic, constant, or time-dependent) */
/* 2) External forcing computed in S/R calcExternalForcing */
/* 3) Prescribed boundary condition (can be periodic, constant, or time-dependent) */
  ierr = PetscOptionsHasName(NULL,NULL,"-forcing_from_file",&useForcingFromFile);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-prescribed_bc",&usePrescribedBC);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-external_forcing",&useExternalForcing);CHKERRQ(ierr);

  if (useForcingFromFile) {  
  	ierr=PetscPrintf(PETSC_COMM_WORLD,"Forcing from file(s) specified\n");CHKERRQ(ierr);  
	for (itr=0; itr<numTracers; itr++) {
	  forcingFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
	}
	maxValsToRead = numTracers;
	ierr = PetscOptionsGetStringArray(NULL,NULL,"-forcing_files",forcingFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
	if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"MUST specify forcing files with the -forcing_files option");
    if (maxValsToRead != numTracers) {
      SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of forcing file names specified");
    }    
    ierr = VecDuplicateVecs(templateVec,numTracers,&uf);CHKERRQ(ierr);    
/*  There are 3 possibilities: periodic, constant, and time-dependent forcing */
    ierr = PetscOptionsHasName(NULL,NULL,"-periodic_forcing",&periodicForcing);CHKERRQ(ierr);
	ierr = PetscOptionsHasName(NULL,NULL,"-time_dependent_forcing",&timeDependentForcing);CHKERRQ(ierr);
    if ((periodicForcing) && (timeDependentForcing)) SETERRQ(PETSC_COMM_WORLD,1,"Cannot specify both periodicForcing and timeDependentForcing");
    if (periodicForcing) {
      ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic forcing from file(s) specified\n");CHKERRQ(ierr);
/*    read time data */
      ierr = iniPeriodicTimer("forcing_", &forcingTimer);CHKERRQ(ierr);
/*    Forcing is read in interpPeriodicForcing */
      for (itr=0; itr<numTracers; itr++) {
        ierr=PetscStrcat(forcingFile[itr],"_");CHKERRQ(ierr);        
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d periodic forcing basename is %s\n",itr,forcingFile[itr]);CHKERRQ(ierr); 
        up[itr].firstTime = PETSC_TRUE; /* initialize periodic vector */        	    
      }      
    } else if (timeDependentForcing) {      
	  ierr=PetscPrintf(PETSC_COMM_WORLD,"Time dependent forcing specified\n");CHKERRQ(ierr);
/*    read time data */
	  ierr = iniTimeDependentTimer("forcing_", &forcingTimeDependentTimer);CHKERRQ(ierr);
/*    Forcing is read in interpTimeDependentVector */
	  for (itr=0; itr<numTracers; itr++) {   
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: forcing will be read from file %s\n", itr,forcingFile[itr]);CHKERRQ(ierr);	  
		utdf[itr].firstTime = PETSC_TRUE;
	  }        
	} else { /* constant forcing */
	  ierr=PetscPrintf(PETSC_COMM_WORLD,"Constant forcing specified\n");CHKERRQ(ierr);      
	  constantForcing = PETSC_TRUE;
	  for (itr=0; itr<numTracers; itr++) {   
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading forcing from file %s\n", itr,forcingFile[itr]);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,forcingFile[itr],FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = VecLoad(uf[itr],fd);CHKERRQ(ierr); /* IntoVector */
		ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	  }
    }

/*  Initialize data to write out uf */
    if ((periodicForcing) || (timeDependentForcing)) {
	  for (itr=0; itr<numTracers; itr++) {
		ufoutFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
	  }
	  maxValsToRead = numTracers;
	  ierr = PetscOptionsGetStringArray(NULL,NULL,"-ouf",ufoutFile,&maxValsToRead,&doWriteUF);CHKERRQ(ierr);
	  if (doWriteUF) {
		if (maxValsToRead != numTracers) {
		  SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of forcing-from-file (uf) output file names specified");
		}  
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Forcing-from-file (uf) output will be written to:\n");CHKERRQ(ierr);
		for (itr=0; itr<numTracers; itr++) {
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,ufoutFile[itr]);CHKERRQ(ierr);
		}
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Note: forcing-from-file (uf) output is discrete in time. Divide by the appropriate time step to obtain a tendency\n");CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Note: forcing-from-file (uf) output is shifted by one time step relative to the stated output time step\n");CHKERRQ(ierr);
		if (doTimeAverage) {
		  for (itr=0; itr<numTracers; itr++) {
			ufavgOutFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
		  }
		  maxValsToRead = numTracers;
		  ierr = PetscOptionsGetStringArray(NULL,NULL,"-ufavg_files",ufavgOutFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
		  if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate file name(s) for writing time averages with the -ufavg_files option");
		  if (maxValsToRead != numTracers) {
			SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of forcing-from-file (uf) time average file names specified");
		  }  
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Forcing-from-file (uf) time averages will be written to:\n");CHKERRQ(ierr);
		  for (itr=0; itr<numTracers; itr++) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,ufavgOutFile[itr]);CHKERRQ(ierr);
		  }      
		}
	  }
    }

    applyForcingFromFile = PETSC_TRUE;
    
    ierr = PetscOptionsGetInt(NULL,NULL,"-forcing_from_file_cutoff_step",&forcingFromFileCutOffStep,&flg1);CHKERRQ(ierr);
    if (forcingFromFileCutOffStep>0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Forcing from file will be turned off after time step %d\n",forcingFromFileCutOffStep);CHKERRQ(ierr);    
    }
    
  } else {  /* no forcing */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"No forcing from file(s) specified\n");CHKERRQ(ierr);    
  }

  if (useExternalForcing) {  /* external forcing present */  
    ierr = PetscPrintf(PETSC_COMM_WORLD,"External forcing is being used\n");CHKERRQ(ierr);    
    ierr = VecDuplicateVecs(templateVec,numTracers,&uef);CHKERRQ(ierr);        
    ierr = iniExternalForcing(time0,Iter0,numTracers,v,uef);CHKERRQ(ierr);

/*  Initialize data to write out uef */
	for (itr=0; itr<numTracers; itr++) {
	  uefoutFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
	}
	maxValsToRead = numTracers;
	ierr = PetscOptionsGetStringArray(NULL,NULL,"-ouef",uefoutFile,&maxValsToRead,&doWriteUEF);CHKERRQ(ierr);
	if (doWriteUEF) {
	  if (maxValsToRead != numTracers) {
		SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of external forcing (uef) output file names specified");
	  }  
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"External forcing (uef) output will be written to:\n");CHKERRQ(ierr);
	  for (itr=0; itr<numTracers; itr++) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,uefoutFile[itr]);CHKERRQ(ierr);
	  }
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Note: external forcing (uef) output is discrete in time. Divide by the appropriate time step to obtain a tendency\n");CHKERRQ(ierr);	  
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Note: external forcing (uef) output is shifted by one time step relative to the stated output time step\n");CHKERRQ(ierr);
      if (doTimeAverage) {
		for (itr=0; itr<numTracers; itr++) {
		  uefavgOutFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
		}
		maxValsToRead = numTracers;
		ierr = PetscOptionsGetStringArray(NULL,NULL,"-uefavg_files",uefavgOutFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
		if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate file name(s) for writing time averages with the -uefavg_files option");
		if (maxValsToRead != numTracers) {
		  SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of external forcing (uef) time average file names specified");
		}  
		ierr = PetscPrintf(PETSC_COMM_WORLD,"External forcing (uef) time averages will be written to:\n");CHKERRQ(ierr);
		for (itr=0; itr<numTracers; itr++) {
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,uefavgOutFile[itr]);CHKERRQ(ierr);
		}      
	  }
    }

    applyExternalForcing = PETSC_TRUE;

    ierr = PetscOptionsGetInt(NULL,NULL,"-external_forcing_cutoff_step",&externalForcingCutOffStep,&flg1);CHKERRQ(ierr);
    if (externalForcingCutOffStep>0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"External forcing will be turned off after time step %d\n",externalForcingCutOffStep);CHKERRQ(ierr);    
    }
    
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"No external forcing is being used\n");CHKERRQ(ierr);  
  }  

/* Prescribed BCs   */
  if (usePrescribedBC) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Prescribed BCs specified\n");CHKERRQ(ierr);    

	ierr = VecCreate(PETSC_COMM_WORLD,&bcTemplateVec);CHKERRQ(ierr);
    
    ierr = PetscOptionsHasName(NULL,NULL,"-calc_bc",&doCalcBC);CHKERRQ(ierr);
    if (doCalcBC) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"BCs will be calculated\n");CHKERRQ(ierr);    

      if ((useProfiles) && (numProcessors>1)) {
        lBCSize = lNumProfiles;
        ierr = VecSetSizes(bcTemplateVec,lBCSize,PETSC_DECIDE);CHKERRQ(ierr);      
        ierr = VecSetFromOptions(bcTemplateVec);CHKERRQ(ierr);
        ierr = VecGetSize(bcTemplateVec,&gBCSize);CHKERRQ(ierr);
      } else {
        ierr = PetscOptionsGetInt(NULL,NULL,"-bc_vec_size",&gBCSize,&flg1);CHKERRQ(ierr);
        if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate size of BC vector with the -bc_vec_size option");
        ierr = VecSetSizes(bcTemplateVec,PETSC_DECIDE,gBCSize);CHKERRQ(ierr);      
        ierr = VecSetFromOptions(bcTemplateVec);CHKERRQ(ierr);
        ierr = VecGetLocalSize(bcTemplateVec,&lBCSize);CHKERRQ(ierr);    
      }

	  ierr = VecGetOwnershipRange(bcTemplateVec,&gbcLow,&gbcHigh);CHKERRQ(ierr);
	  gbcHigh = gbcHigh - 1; /* Note: gbcHigh is one more than the last local element */
	  ierr = PetscMalloc(lBCSize*sizeof(PetscInt),&gBCIndices);CHKERRQ(ierr);  
	  for (il=0; il<lBCSize; il++) {
		gBCIndices[il] = il + gbcLow;
	  }
      
      ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&bcc);CHKERRQ(ierr);    
      ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&bcf);CHKERRQ(ierr);    

      ierr = iniCalcBC(time0,Iter0,time0+deltaTClock,Iter0+1,numTracers,v,bcc,bcf);CHKERRQ(ierr);
      
    } else { /* read from file */
      for (itr=0; itr<numTracers; itr++) {
        bcFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
      }
      maxValsToRead = numTracers;
      ierr = PetscOptionsGetStringArray(NULL,NULL,"-bc_files",bcFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
      if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"MUST specify BC files with the -bc_files option");
      if (maxValsToRead != numTracers) {
        SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of BC file names specified");
      }    

/*    There are 3 possibilities: periodic, constant, and time-dependent BCs */
      ierr = PetscOptionsHasName(NULL,NULL,"-periodic_bc",&periodicBC);CHKERRQ(ierr);
	  ierr = PetscOptionsHasName(NULL,NULL,"-time_dependent_bc",&timeDependentBC);CHKERRQ(ierr);
	  if ((periodicBC) && (timeDependentBC)) SETERRQ(PETSC_COMM_WORLD,1,"Cannot specify both periodicBC and timeDependentBC");
      if (periodicBC) {
        ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic BC from file(s) specified\n");CHKERRQ(ierr);
/*      read time data */
        ierr = iniPeriodicTimer("bc_", &bcTimer);CHKERRQ(ierr);
        for (itr=0; itr<numTracers; itr++) {
          ierr=PetscStrcat(bcFile[itr],"_");CHKERRQ(ierr);        
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: periodic BC basename is %s\n",itr,bcFile[itr]);CHKERRQ(ierr); 
          bcp[itr].firstTime = PETSC_TRUE; /* initialize periodic vector */        	    
        }        
/*      Load one vector here as a template; BC is read in interpPeriodicVector */
        strcpy(tmpFile,"");
        sprintf(tmpFile,"%s%02d",bcFile[0],0);
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
        ierr = VecLoad(bcTemplateVec,fd);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
        ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&bcc);CHKERRQ(ierr);    
        ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&bcf);CHKERRQ(ierr);    
      } else if (timeDependentBC) {      
		ierr=PetscPrintf(PETSC_COMM_WORLD,"Time dependent BC specified\n");CHKERRQ(ierr);
/*      read time data */
	    ierr = iniTimeDependentTimer("bc_", &bcTimeDependentTimer);CHKERRQ(ierr);
		for (itr=0; itr<numTracers; itr++) {   
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: BC will be read from %s\n",itr,bcFile[itr]);CHKERRQ(ierr); 
		  bctd[itr].firstTime = PETSC_TRUE;
		}
/*      Load one vector here as a template; BC is read in interpTimeDependentVector */
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,bcFile[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = VecLoad(bcTemplateVec,fd);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
		ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&bcc);CHKERRQ(ierr);    
		ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&bcf);CHKERRQ(ierr);    
	  } else { /* constant BC */
		ierr=PetscPrintf(PETSC_COMM_WORLD,"Constant BC specified\n");CHKERRQ(ierr);      
		constantBC = PETSC_TRUE;
/*      Load one vector here as a template */
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,bcFile[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = VecLoad(bcTemplateVec,fd);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
		ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&bcc);CHKERRQ(ierr);    
		ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&bcf);CHKERRQ(ierr);            
/*      Load BCs */
		for (itr=0; itr<numTracers; itr++) {   
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading BC from file %s\n", itr,bcFile[itr]);CHKERRQ(ierr);
		  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,bcFile[itr],FILE_MODE_READ,&fd);CHKERRQ(ierr);
		  ierr = VecLoad(bcc[itr],fd);CHKERRQ(ierr); /* IntoVector */
		  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
		  ierr = VecCopy(bcc[itr],bcf[itr]);CHKERRQ(ierr);		  
		}
      }
    }

    ierr = VecGetLocalSize(bcTemplateVec,&lBCSize);CHKERRQ(ierr); /* just to be safe */

    if ((useProfiles) && (numProcessors>1)) {
      if (lBCSize != lNumProfiles) {
        SETERRQ(PETSC_COMM_WORLD,1,"Problem with partitioning of BC vectors! lNumProfiles must equal lBCSize");
      }
    }
    
    applyBC = PETSC_TRUE;

    ierr = PetscOptionsGetInt(NULL,NULL,"-bc_cutoff_step",&bcCutOffStep,&flg1);CHKERRQ(ierr);
    if (bcCutOffStep>0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Prescribed BC will be turned off after time step %d\n",bcCutOffStep);CHKERRQ(ierr);    
    }

/*  BC output file */
	for (itr=0; itr<numTracers; itr++) {
	  bcoutFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
	}
	maxValsToRead = numTracers;
	ierr = PetscOptionsGetStringArray(NULL,NULL,"-obc",bcoutFile,&maxValsToRead,&doWriteBC);CHKERRQ(ierr);
	if (doWriteBC) {
	  if (maxValsToRead != numTracers) {
		SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of BC output file names specified");
	  }  
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"BC output will be written to:\n");CHKERRQ(ierr);
	  for (itr=0; itr<numTracers; itr++) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,bcoutFile[itr]);CHKERRQ(ierr);
	  }
      if (doTimeAverage) {
		for (itr=0; itr<numTracers; itr++) {
		  bcavgOutFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
		}
		maxValsToRead = numTracers;
		ierr = PetscOptionsGetStringArray(NULL,NULL,"-bcavg_files",bcavgOutFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
		if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate file name(s) for writing time averages with the -bcavg_files option");
		if (maxValsToRead != numTracers) {
		  SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of BC time average file names specified");
		}  
		ierr = PetscPrintf(PETSC_COMM_WORLD,"BC time averages will be written to:\n");CHKERRQ(ierr);
		for (itr=0; itr<numTracers; itr++) {
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,bcavgOutFile[itr]);CHKERRQ(ierr);
		}      
	  }
    }
    
/*  Matrices */
	ierr = PetscOptionsGetString(NULL,NULL,"-mbe",matbeFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
	if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate binary boundary matrix file name with the -mbe option");
	ierr = PetscOptionsGetString(NULL,NULL,"-mbi",matbiFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
	if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate binary boundary matrix file name with the -mbi options");

	ierr = MatCreate(PETSC_COMM_WORLD,&Be);CHKERRQ(ierr);
	ierr = MatCreate(PETSC_COMM_WORLD,&Bi);CHKERRQ(ierr);	

/*  Set layout information */
	if ((useProfiles) && (numProcessors>1)) {
	  ierr = MatSetSizes(Be,lSize,lBCSize,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
	  ierr = MatSetSizes(Bi,lSize,lBCSize,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
	}		
	
	ierr = MatSetType(Be,MATMPIAIJ);CHKERRQ(ierr);      
	ierr = MatSetFromOptions(Be);CHKERRQ(ierr);  
	ierr = MatSetType(Bi,MATMPIAIJ);CHKERRQ(ierr);        
	ierr = MatSetFromOptions(Bi);CHKERRQ(ierr);
  
	if (periodicMatrix) {    
	  ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic boundary matrices specified\n");CHKERRQ(ierr);
	  ierr=PetscStrcat(matbeFile,"_");CHKERRQ(ierr);
	  ierr=PetscStrcat(matbiFile,"_");CHKERRQ(ierr);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Be basename is %s\n", matbeFile);CHKERRQ(ierr); 
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Bi basename is %s\n", matbiFile);CHKERRQ(ierr);     

/*    Read here to set sparsity pattern */  
	  strcpy(tmpFile,"");
	  sprintf(tmpFile,"%s%02d",matbeFile,0);
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = MatLoad(Be,fd);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

	  strcpy(tmpFile,"");
	  sprintf(tmpFile,"%s%02d",matbiFile,0);
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = MatLoad(Bi,fd);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
    
	  Bep.firstTime = PETSC_TRUE;
	  Bip.firstTime = PETSC_TRUE;

    } else if (timeDependentMatrix) {
	  ierr=PetscPrintf(PETSC_COMM_WORLD,"Time-dependent matrices specified\n");CHKERRQ(ierr);
	  ierr=PetscStrcat(matbeFile,"_");CHKERRQ(ierr);
	  ierr=PetscStrcat(matbiFile,"_");CHKERRQ(ierr);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Be basename is %s\n", matbeFile);CHKERRQ(ierr); 
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Bi basename is %s\n", matbiFile);CHKERRQ(ierr);     

/*    Read here to set sparsity pattern */  
	  strcpy(tmpFile,"");
	  sprintf(tmpFile,"%s%02d",matbeFile,0);
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = MatLoad(Be,fd);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
        
	  strcpy(tmpFile,"");
	  sprintf(tmpFile,"%s%02d",matbiFile,0);
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = MatLoad(Bi,fd);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    

	  Betd.firstTime = PETSC_TRUE;
	  Bitd.firstTime = PETSC_TRUE;  
        
    } else { /*  not periodic. read matrices here */

	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,matbeFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading Be from file %s\n", matbeFile);CHKERRQ(ierr);  
	  ierr = MatLoad(Be,fd);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,matbiFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading Bi from file %s\n", matbiFile);CHKERRQ(ierr);  
	  ierr = MatLoad(Bi,fd);CHKERRQ(ierr);    
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
	}
  } else {  /* no BC */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"No prescribed BCs specified\n");CHKERRQ(ierr);  
  }  

  ierr = PetscOptionsGetString(NULL,NULL,"-rescale_forcing_file",rfsFile,PETSC_MAX_PATH_LEN-1,&rescaleForcing);CHKERRQ(ierr);
  if (rescaleForcing) {  
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Forcing will be rescaled\n");CHKERRQ(ierr); 	    
    ierr = VecDuplicate(templateVec,&Rfs);CHKERRQ(ierr);    
	if (periodicMatrix) {    
	  ierr=PetscStrcat(rfsFile,"_");CHKERRQ(ierr);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Rescale forcing factor file basename is %s\n", rfsFile);CHKERRQ(ierr); 	
	  Rfsp.firstTime = PETSC_TRUE;
	  if (constantForcing) {
        SETERRQ(PETSC_COMM_WORLD,1,"Periodic rescaling not supported with constant forcing");
	  }
	} else if (timeDependentMatrix) {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Rescale forcing factor file is %s\n", rfsFile);CHKERRQ(ierr);	
	  Rfstd.firstTime = PETSC_TRUE;
	  if (constantForcing) {
        SETERRQ(PETSC_COMM_WORLD,1,"Time dependent rescaling not supported with constant forcing");
	  }	  
	} else {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading rescale forcing factor from file %s\n", rfsFile);CHKERRQ(ierr);  	
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,rfsFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = VecLoad(Rfs,fd);CHKERRQ(ierr);  
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	  if (constantForcing) {
		for (itr=0; itr<numTracers; itr++) {
		  ierr = VecPointwiseMult(uf[itr],Rfs,uf[itr]);CHKERRQ(ierr);
        }
	  }    
	}  
  }

  ierr = iniTMMWrite(time0,Iter0,numTracers,v);CHKERRQ(ierr);

/* initialize monitor */
  ierr = PetscOptionsHasName(NULL,NULL,"-monitor",&useMonitor);CHKERRQ(ierr);
  if (useMonitor) {  
    ierr = iniMonitor(time0,Iter0,numTracers,v);CHKERRQ(ierr);
  }

/* initialize misfit */
  ierr = PetscOptionsHasName(NULL,NULL,"-calc_misfit",&doMisfit);CHKERRQ(ierr);
  if (doMisfit) {  
    ierr = iniMisfit(time0,Iter0,numTracers,v);CHKERRQ(ierr);
  }

/* Open files for output and optionally write initial conditions */
#if !defined (FORSPINUP) && !defined (FORJACOBIAN)
  for (itr=0; itr<numTracers; itr++) {       
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,outFile[itr],OUTPUT_FILE_MODE,&fdout[itr]);CHKERRQ(ierr);
  }

  if (doWriteUF) {
/*  We only open files here since uf may not yet have been updated */
	for (itr=0; itr<numTracers; itr++) {
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,ufoutFile[itr],OUTPUT_FILE_MODE,&fdufout[itr]);CHKERRQ(ierr);
	}
  }

  if (doWriteUEF) {
/*  We only open files here since uef may not yet have been updated */
	for (itr=0; itr<numTracers; itr++) {
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,uefoutFile[itr],OUTPUT_FILE_MODE,&fduefout[itr]);CHKERRQ(ierr);
	}
  }
    
  if (doWriteBC) {
/*  We only open files here since BCs may not yet have been computed */
	for (itr=0; itr<numTracers; itr++) {
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,bcoutFile[itr],OUTPUT_FILE_MODE,&fdbcout[itr]);CHKERRQ(ierr);
	}
  }

  if (!appendOutput) {
    ierr = PetscFOpen(PETSC_COMM_WORLD,outTimeFile,"w",&fptime);CHKERRQ(ierr);  
	if (Iter0==writeTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
	  ierr = PetscFPrintf(PETSC_COMM_WORLD,fptime,"%d   %10.5f\n",Iter0,time0);CHKERRQ(ierr);     
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing output at time %10.5f, step %d\n", time0,Iter0);CHKERRQ(ierr);  
	  for (itr=0; itr<numTracers; itr++) {       
		ierr = VecView(v[itr],fdout[itr]);CHKERRQ(ierr);
	  }
	}  
  } else {
	ierr = PetscFOpen(PETSC_COMM_WORLD,outTimeFile,"a",&fptime);CHKERRQ(ierr);
	if (Iter0==writeTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Opening file(s) for output. Initial condition will NOT be written\n");CHKERRQ(ierr);  
	}  
  }

  if (writePickup) {
	if (Iter0==pickupTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing pickup at time %10.5f, step %d\n", time0, Iter0);CHKERRQ(ierr);
	  strcpy(tmpFile,"");
	  sprintf(tmpFile,"pickup_%d.petsc",Iter0);
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_WRITE,&fdp);CHKERRQ(ierr);
	  for (itr=0; itr<numTracers; itr++) {
		ierr = VecView(v[itr],fdp);CHKERRQ(ierr);
	  }
	  ierr = PetscViewerDestroy(&fdp);CHKERRQ(ierr);      
	}
  }
  
#endif

  if (doTimeAverage) {  
    ierr = VecDuplicateVecs(templateVec,numTracers,&vavg);CHKERRQ(ierr);
	for (itr=0; itr<numTracers; itr++) {
	  ierr = VecSet(vavg[itr],zero); CHKERRQ(ierr);
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,avgOutFile[itr],AVG_FILE_MODE,&fdavgout[itr]);CHKERRQ(ierr);
    }
	if (doWriteUF) {
      ierr = VecDuplicateVecs(templateVec,numTracers,&ufavg);CHKERRQ(ierr);  	
	  for (itr=0; itr<numTracers; itr++) {
		ierr = VecSet(ufavg[itr],zero); CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,ufavgOutFile[itr],AVG_FILE_MODE,&fdufavgout[itr]);CHKERRQ(ierr);
	  }
	}
	if (doWriteUEF) {
      ierr = VecDuplicateVecs(templateVec,numTracers,&uefavg);CHKERRQ(ierr);  	
	  for (itr=0; itr<numTracers; itr++) {
		ierr = VecSet(uefavg[itr],zero); CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,uefavgOutFile[itr],AVG_FILE_MODE,&fduefavgout[itr]);CHKERRQ(ierr);
	  }
	}
	if (doWriteBC) {
      ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&bcavg);CHKERRQ(ierr);  	
	  for (itr=0; itr<numTracers; itr++) {
		ierr = VecSet(bcavg[itr],zero); CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,bcavgOutFile[itr],AVG_FILE_MODE,&fdbcavgout[itr]);CHKERRQ(ierr);
	  }
	}    
  }

/* reinitialize forcing if required */
#ifdef FORSPINUP
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"itjac.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fp,&itjac,1,NULL,PETSC_INT);CHKERRQ(ierr);  
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting iteration number to: %d\n", itjac);CHKERRQ(ierr);

  tc=time0 + deltaTClock*(itjac-1);  /*  current time (time at beginning of step) */
  Iterc=Iter0+itjac-1;
  if (useExternalForcing) ierr = reInitializeExternalForcing(tc,Iterc,numTracers,v,uef);CHKERRQ(ierr);      
  if (doCalcBC) ierr = reInitializeCalcBC(tc,Iterc,tc+deltaTClock,Iterc+1,numTracers,v,bcc,bcf);CHKERRQ(ierr);        
#endif

#ifdef FORJACOBIAN
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"itjac.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fp,&itjac,1,NULL,PETSC_INT);CHKERRQ(ierr);  
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting iteration number to: %d\n", itjac);CHKERRQ(ierr);

  if (doTimeAverage) {  
    tc=time0 + deltaTClock*(itjac-1);  /*  current time (time at beginning of step) */
    Iterc=Iter0+itjac-1;
    if (useExternalForcing) ierr = reInitializeExternalForcing(tc,Iterc,numTracers,v,uef);CHKERRQ(ierr);           
  } else {
    for (itr=0; itr<numTracers; itr++) {       
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,outFile[itr],FILE_MODE_WRITE,&fdout[itr]);CHKERRQ(ierr);
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,iniFile[itr],FILE_MODE_READ,&fdin[itr]);CHKERRQ(ierr);
    }
  }
#endif

/* Start time stepping loop */
  ierr = PetscTime(&t1);CHKERRQ(ierr); /* start counting wall clock time */  
  for (iLoop = 1; iLoop <= maxSteps; iLoop++) {
/*  iLoop -> iLoop+1 (convention) */  
#if !defined (FORSPINUP) && !defined (FORJACOBIAN)
    tc=time0 + deltaTClock*(iLoop-1);  /*  current time (time at beginning of step) */
    tf=time0 + deltaTClock*iLoop;  /*  future time (time at end of step) */
    Iterc=Iter0+iLoop-1;
#else
    tc=time0 + deltaTClock*(itjac-1);  /*  current time (time at beginning of step) */
    tf=time0 + deltaTClock*itjac;  /*  future time (time at end of step) */
    Iterc=Iter0+itjac-1;
    itjac++;    
#endif

#ifdef FORJACOBIAN
    if (!doTimeAverage) {
      for (itr=0; itr<numTracers; itr++) {
/* reread initial conditions (these were read once already above, but do it again since its easier to code) */    
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading initial condition from file %s\n", itr,iniFile[itr]);CHKERRQ(ierr);
        ierr = VecLoad(v[itr],fdin[itr]);CHKERRQ(ierr); /* IntoVector */
      } 
      if (useExternalForcing) ierr = reInitializeExternalForcing(tc,Iterc,numTracers,v,uef);CHKERRQ(ierr);          
    }
#endif

#ifndef FORJACOBIAN
/*  interpolate Ae,Ai,uf,uef,bcc to current time (tc) and bcf to future time (tf) */
    if (periodicMatrix) {
      ierr = interpPeriodicMatrix(tc,&Ae,matrixPeriodicTimer.cyclePeriod,matrixPeriodicTimer.numPerPeriod,
                                  matrixPeriodicTimer.tdp,&Aep,mateFile);
      ierr = interpPeriodicMatrix(tc,&Ai,matrixPeriodicTimer.cyclePeriod,matrixPeriodicTimer.numPerPeriod,
                                  matrixPeriodicTimer.tdp,&Aip,matiFile);
      if (rescaleForcing) {
		ierr = interpPeriodicVector(tc,&Rfs,matrixPeriodicTimer.cyclePeriod,matrixPeriodicTimer.numPerPeriod,
                                  matrixPeriodicTimer.tdp,&Rfsp,rfsFile);
      }
    } else if (timeDependentMatrix) {
      ierr = interpTimeDependentMatrix(tc,&Ae,matrixTimeDependentTimer.numTimes,matrixTimeDependentTimer.tdt,&Aetd,mateFile);
      ierr = interpTimeDependentMatrix(tc,&Ai,matrixTimeDependentTimer.numTimes,matrixTimeDependentTimer.tdt,&Aitd,matiFile);
      if (rescaleForcing) {
		ierr = interpTimeDependentVector(tc,&Rfs,matrixTimeDependentTimer.numTimes,matrixTimeDependentTimer.tdt,&Rfstd,rfsFile);
      }
    }

/*  Forcing     */
    if (applyForcingFromFile) {
      if ((forcingFromFileCutOffStep>0) && (iLoop==(forcingFromFileCutOffStep+1))) {
        applyForcingFromFile = PETSC_FALSE;
		for (itr=0; itr<numTracers; itr++) {
		  VecSet(uf[itr],zero);
		}
      } else {
		if (periodicForcing) {
		  for (itr=0; itr<numTracers; itr++) {    
			ierr = interpPeriodicVector(tc,&uf[itr],forcingTimer.cyclePeriod,forcingTimer.numPerPeriod,forcingTimer.tdp,&up[itr],forcingFile[itr]);
			if (rescaleForcing) {
				ierr = VecPointwiseMult(uf[itr],Rfs,uf[itr]);CHKERRQ(ierr);
			}
		  }    
		} else if (timeDependentForcing) {
		  for (itr=0; itr<numTracers; itr++) {    
		    ierr = interpTimeDependentVector(tc,&uf[itr],forcingTimeDependentTimer.numTimes,forcingTimeDependentTimer.tdt,&utdf[itr],forcingFile[itr]);
			if (rescaleForcing) {
				ierr = VecPointwiseMult(uf[itr],Rfs,uf[itr]);CHKERRQ(ierr);
			}
		  }
		}
	  }	
	}
#endif

    if (applyExternalForcing) {
      if ((externalForcingCutOffStep>0) && (iLoop==(externalForcingCutOffStep+1))) {
        applyExternalForcing = PETSC_FALSE;
		for (itr=0; itr<numTracers; itr++) {
		  VecSet(uef[itr],zero);
		}
      } else {    
		ierr = calcExternalForcing(tc,Iterc,iLoop,numTracers,v,uef);CHKERRQ(ierr); /* Compute external forcing in uef */
		if (rescaleForcing) {
		  for (itr=0; itr<numTracers; itr++) {    
			ierr = VecPointwiseMult(uef[itr],Rfs,uef[itr]);CHKERRQ(ierr);
		  }  
		}
	  }	
    } 

#ifndef FORJACOBIAN
    if (applyBC) {
      if ((bcCutOffStep>0) && (iLoop==(bcCutOffStep+1))) {
        applyBC = PETSC_FALSE;
        doCalcBC = PETSC_FALSE;
		for (itr=0; itr<numTracers; itr++) {
		  VecSet(bcc[itr],zero);
		  VecSet(bcf[itr],zero);
		}        
      } else {    
        if (periodicMatrix) {
          ierr = interpPeriodicMatrix(tc,&Be,matrixPeriodicTimer.cyclePeriod,matrixPeriodicTimer.numPerPeriod,
                                      matrixPeriodicTimer.tdp,&Bep,matbeFile);
          ierr = interpPeriodicMatrix(tc,&Bi,matrixPeriodicTimer.cyclePeriod,matrixPeriodicTimer.numPerPeriod,
                                      matrixPeriodicTimer.tdp,&Bip,matbiFile);
		} else if (timeDependentMatrix) {
		  ierr = interpTimeDependentMatrix(tc,&Be,matrixTimeDependentTimer.numTimes,matrixTimeDependentTimer.tdt,&Betd,matbeFile);
		  ierr = interpTimeDependentMatrix(tc,&Bi,matrixTimeDependentTimer.numTimes,matrixTimeDependentTimer.tdt,&Bitd,matbiFile);
		}
        if (periodicBC) {
          for (itr=0; itr<numTracers; itr++) {    
            ierr = interpPeriodicVector(tc,&bcc[itr],bcTimer.cyclePeriod,bcTimer.numPerPeriod,bcTimer.tdp,&bcp[itr],bcFile[itr]);
            ierr = interpPeriodicVector(tf,&bcf[itr],bcTimer.cyclePeriod,bcTimer.numPerPeriod,bcTimer.tdp,&bcp[itr],bcFile[itr]);
          }    
        } else if (timeDependentBC) {
		  for (itr=0; itr<numTracers; itr++) {    
			ierr = interpTimeDependentVector(tc,&bcc[itr],bcTimeDependentTimer.numTimes,bcTimeDependentTimer.tdt,&bctd[itr],bcFile[itr]);
			ierr = interpTimeDependentVector(tf,&bcf[itr],bcTimeDependentTimer.numTimes,bcTimeDependentTimer.tdt,&bctd[itr],bcFile[itr]);
		  }
        } else if (doCalcBC) {
          ierr = calcBC(tc,Iterc,tc+deltaTClock,Iterc+1,iLoop,numTracers,v,bcc,bcf);CHKERRQ(ierr); /* Compute BC in bcc and bcf */	  
        }
      }
	  
    }

/*  Write out BC at first time step */
	if (doWriteBC) {
	  if ((iLoop == 1) && (!appendOutput)) {	  
		for (itr=0; itr<numTracers; itr++) {
		  ierr = VecView(bcc[itr],fdbcout[itr]);CHKERRQ(ierr);
		}
	  }	
    }
        
    ierr = forwardStep(tc,Iterc,deltaTClock,numTracers,applyForcingFromFile,applyExternalForcing,applyBC,
                       v,Ae,Ai,Be,Bi,uf,uef,bcc,bcf,vtmp);CHKERRQ(ierr);
#endif    
    tc=time0 + deltaTClock*iLoop;  /*  time at end of step */    

    if (useMonitor) {
      ierr = updateMonitor(tc,iLoop,numTracers,v);CHKERRQ(ierr);
      ierr = writeMonitor(tc,iLoop,numTracers,v);CHKERRQ(ierr);
    }

    if (doMisfit) {
      ierr = calcMisfit(tc,iLoop,numTracers,v);CHKERRQ(ierr);
      ierr = writeMisfit(tc,iLoop,numTracers,v);CHKERRQ(ierr);
    }

/* write output */
#if !defined (FORSPINUP) && !defined (FORJACOBIAN)
	if (Iter0+iLoop>writeTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
/*    We do this here in case writeExternalForcing etc want to use writeTimer */	
	  if (writeTimer.count<writeTimer.numTimeSteps) {
	    writeTimer.count++;
	  }
    }
    if (applyExternalForcing) {
      ierr = writeExternalForcing(tc,iLoop,numTracers,v,uef);CHKERRQ(ierr);
    }
    if (doCalcBC) {
      ierr = writeBC(tc,iLoop,numTracers,v,bcc,bcf);CHKERRQ(ierr); 
    }

	if (Iter0+iLoop>=writeTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
	  if ((writeTimer.count==writeTimer.numTimeSteps) || (Iter0+iLoop==writeTimer.startTimeStep)) { /* time to write out */
	
		 ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
		 ierr = PetscFPrintf(PETSC_COMM_WORLD,fptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);           
		 for (itr=0; itr<numTracers; itr++) {
		   ierr = VecView(v[itr],fdout[itr]);CHKERRQ(ierr);
		 }

		 if (doWriteUF) {
		   for (itr=0; itr<numTracers; itr++) {
			 ierr = VecView(uf[itr],fdufout[itr]);CHKERRQ(ierr);
		   }
		 }

		 if (doWriteUEF) {
		   for (itr=0; itr<numTracers; itr++) {
			 ierr = VecView(uef[itr],fduefout[itr]);CHKERRQ(ierr);
		   }
		 }
	  
		if (doWriteBC) {
		  for (itr=0; itr<numTracers; itr++) {
			ierr = VecView(bcf[itr],fdbcout[itr]);CHKERRQ(ierr);
		  }
		}
		
      }
      
      if (writeTimer.count==writeTimer.numTimeSteps) {
		ierr = updateStepTimer("write_", Iter0+iLoop, &writeTimer);CHKERRQ(ierr);
      }
      
    }
    
    ierr = doTMMWrite(tc,iLoop,numTracers,v);CHKERRQ(ierr);

    if (writePickup) {
	  if (Iter0+iLoop>pickupTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
		if (pickupTimer.count<pickupTimer.numTimeSteps) {
		  pickupTimer.count++;
		}
      }		

	  if (Iter0+iLoop>=pickupTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
        if ((pickupTimer.count==pickupTimer.numTimeSteps) || (Iter0+iLoop==pickupTimer.startTimeStep)) { /* time to write out */
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing pickup at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
		  strcpy(tmpFile,"");
		  sprintf(tmpFile,"pickup_%d.petsc",Iter0+iLoop);
		  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_WRITE,&fdp);CHKERRQ(ierr);
		  for (itr=0; itr<numTracers; itr++) {
			ierr = VecView(v[itr],fdp);CHKERRQ(ierr);
		  }
		  ierr = PetscViewerDestroy(&fdp);CHKERRQ(ierr);      
        }

        if (pickupTimer.count==pickupTimer.numTimeSteps) {
		  ierr = updateStepTimer("pickup_", Iter0+iLoop, &pickupTimer);CHKERRQ(ierr);
		}
	  }
    }

#else
#ifdef FORSPINUP
	if (Iter0+iLoop>=writeTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */
	  if (writeTimer.count<=writeTimer.numTimeSteps) {
	    writeTimer.count++;	  
	  }
	  if (writeTimer.count==writeTimer.numTimeSteps) { /* time to write out */
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
		for (itr=0; itr<numTracers; itr++) {       
		  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,outFile[itr],FILE_MODE_WRITE,&fdout[itr]);CHKERRQ(ierr);
		  ierr = VecView(v[itr],fdout[itr]);CHKERRQ(ierr);
		  ierr = PetscViewerDestroy(&fdout[itr]);CHKERRQ(ierr);		
		}

		ierr = PetscPrintf(PETSC_COMM_WORLD,"Waiting for new initial condition ...\n");CHKERRQ(ierr);
		ierr = waitForSignal(10);CHKERRQ(ierr);

		for (itr=0; itr<numTracers; itr++) {
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading new initial condition from file %s\n", itr,iniFile[itr]);CHKERRQ(ierr);
		  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,iniFile[itr],FILE_MODE_READ,&fd);CHKERRQ(ierr);
		  ierr = VecLoad(v[itr],fd);CHKERRQ(ierr); /* IntoVector */
		  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);          
		}        

		ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"itjac.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
		ierr = PetscBinaryRead(fp,&itjac,1,NULL,PETSC_INT);CHKERRQ(ierr);  
		ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

		ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting iteration number to: %d\n", itjac);CHKERRQ(ierr);

		tc=time0 + deltaTClock*(itjac-1);  /*  current time (time at beginning of step) */
		tf=time0 + deltaTClock*itjac;  /*  future time (time at end of step) */
		Iterc=Iter0+itjac-1;

		if (applyExternalForcing) ierr = reInitializeExternalForcing(tc,Iterc,numTracers,v,uef);CHKERRQ(ierr);
		if (doCalcBC) ierr = reInitializeCalcBC(tc,Iterc,tc+deltaTClock,Iterc+1,numTracers,v,bcc,bcf);CHKERRQ(ierr);
	  }	
    }
#endif
#ifdef FORJACOBIAN
    if (doTimeAverage) {
	  if (Iter0+iLoop>=avgTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */
		if (avgTimer.count<avgTimer.numTimeSteps) { /* still within same averaging block so accumulate */
		  for (itr=0; itr<numTracers; itr++) {
			ierr = VecAXPY(vavg[itr],one,uef[itr]);CHKERRQ(ierr);
		  }          
		  avgTimer.count++;
/*           ierr = PetscPrintf(PETSC_COMM_WORLD,"Accumulating: %d\n", iLoop);CHKERRQ(ierr);                       */
		}
		if (avgTimer.count==avgTimer.numTimeSteps) { /* time to write averages to file */
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing time average q at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);                      
		  for (itr=0; itr<numTracers; itr++) {
			ierr = VecScale(vavg[itr],1.0/avgTimer.count);CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,outFile[itr],FILE_MODE_WRITE,&fdout[itr]);CHKERRQ(ierr);            
			ierr = VecView(vavg[itr],fdout[itr]);CHKERRQ(ierr);              
			ierr = PetscViewerDestroy(&fdout[itr]);CHKERRQ(ierr);		              
			ierr = VecSet(vavg[itr],zero); CHKERRQ(ierr);
		  }        
		  ierr = updateStepTimer("avg_", Iter0+iLoop, &avgTimer);CHKERRQ(ierr);              

		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Waiting for new initial condition ...\n");CHKERRQ(ierr);
		  ierr = waitForSignal(10);CHKERRQ(ierr);

		  for (itr=0; itr<numTracers; itr++) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading new initial condition from file %s\n", itr,iniFile[itr]);CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,iniFile[itr],FILE_MODE_READ,&fd);CHKERRQ(ierr);
			ierr = VecLoad(v[itr],fd);CHKERRQ(ierr); /* IntoVector */
			ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);          
		  }        
  
		  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"itjac.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
		  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
		  ierr = PetscBinaryRead(fp,&itjac,1,NULL,PETSC_INT);CHKERRQ(ierr);  
		  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting iteration number to: %d\n", itjac);CHKERRQ(ierr);
	
		  tc=time0 + deltaTClock*(itjac-1);  /*  current time (time at beginning of step) */
		  tf=time0 + deltaTClock*itjac;  /*  future time (time at end of step) */
		  Iterc=Iter0+itjac-1;
	
		  if (applyExternalForcing) ierr = reInitializeExternalForcing(tc,Iterc,numTracers,v,uef);CHKERRQ(ierr);      
		
		} else {
		
		  for (itr=0; itr<numTracers; itr++) {
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading new initial condition from file %s\n", itr,iniFile[itr]);CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,iniFile[itr],FILE_MODE_READ,&fd);CHKERRQ(ierr);
			ierr = VecLoad(v[itr],fd);CHKERRQ(ierr); /* IntoVector */
			ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);          
		  }        

		  tc=time0 + deltaTClock*(itjac-1);  /* this is now the time at end of time step */
		  Iterc=Iter0+itjac-1;    
		  if (applyExternalForcing) ierr = reInitializeExternalForcing(tc,Iterc,numTracers,v,uef);CHKERRQ(ierr);                  
		}
	  }
    } else { /* no time averaging */
	  if (Iter0+iLoop>=writeTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */
		if (writeTimer.count<=writeTimer.numTimeSteps) {
		  writeTimer.count++;	  
		}
		if (writeTimer.count==writeTimer.numTimeSteps) { /* time to write out */
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing q at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
		  for (itr=0; itr<numTracers; itr++) {       
			ierr = VecView(uef[itr],fdout[itr]);CHKERRQ(ierr);
		  }
		  for (itr=0; itr<numTracers; itr++) {       
			ierr = PetscViewerDestroy(&fdout[itr]);CHKERRQ(ierr);		
			ierr = PetscViewerDestroy(&fdin[itr]);CHKERRQ(ierr);          
		  }
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Waiting for new initial condition ...\n");CHKERRQ(ierr);
		  ierr = waitForSignal(10);CHKERRQ(ierr);
  /*      open files here for I/O */
		  for (itr=0; itr<numTracers; itr++) {       
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,outFile[itr],FILE_MODE_WRITE,&fdout[itr]);CHKERRQ(ierr);
			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,iniFile[itr],FILE_MODE_READ,&fdin[itr]);CHKERRQ(ierr);
		  }

		  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"itjac.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
		  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
		  ierr = PetscBinaryRead(fp,&itjac,1,NULL,PETSC_INT);CHKERRQ(ierr);  
		  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);          
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting iteration number to: %d\n", itjac);CHKERRQ(ierr);
		}		
	  }	
    }	  
#endif
#endif
#ifndef FORJACOBIAN    
    if (doTimeAverage) {
      if (Iter0+iLoop>=avgTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */
        if (avgTimer.count<avgTimer.numTimeSteps) { /* still within same averaging block so accumulate */
/*           ierr = PetscPrintf(PETSC_COMM_WORLD,"Accumulating for time average\n");CHKERRQ(ierr);               */
		  for (itr=0; itr<numTracers; itr++) {
			ierr = VecAXPY(vavg[itr],one,v[itr]);CHKERRQ(ierr);
		  }   
		  if (doWriteUF) {
			for (itr=0; itr<numTracers; itr++) {
				ierr = VecAXPY(ufavg[itr],one,uf[itr]);CHKERRQ(ierr);
			}
		  }
		  if (doWriteUEF) {
			for (itr=0; itr<numTracers; itr++) {
				ierr = VecAXPY(uefavg[itr],one,uef[itr]);CHKERRQ(ierr);
			}
		  }
		  if (doWriteBC) {
			for (itr=0; itr<numTracers; itr++) {
				ierr = VecAXPY(bcavg[itr],one,bcf[itr]);CHKERRQ(ierr);
			}
		  }	
		  avgTimer.count++;
        }
        if (avgTimer.count==avgTimer.numTimeSteps) { /* time to write averages to file */        
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing time average at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
		  ierr = PetscFPrintf(PETSC_COMM_WORLD,avgfptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);
		  for (itr=0; itr<numTracers; itr++) {
			ierr = VecScale(vavg[itr],1.0/avgTimer.count);CHKERRQ(ierr);
            ierr = VecView(vavg[itr],fdavgout[itr]);CHKERRQ(ierr);
            ierr = VecSet(vavg[itr],zero); CHKERRQ(ierr);
		  }
		  if (doWriteUF) {
			for (itr=0; itr<numTracers; itr++) {		  
			  ierr = VecScale(ufavg[itr],1.0/avgTimer.count);CHKERRQ(ierr);
			  ierr = VecView(ufavg[itr],fdufavgout[itr]);CHKERRQ(ierr);
			  ierr = VecSet(ufavg[itr],zero); CHKERRQ(ierr);
			}
		  }
		  if (doWriteUEF) {
			for (itr=0; itr<numTracers; itr++) {		  
			  ierr = VecScale(uefavg[itr],1.0/avgTimer.count);CHKERRQ(ierr);
			  ierr = VecView(uefavg[itr],fduefavgout[itr]);CHKERRQ(ierr);
			  ierr = VecSet(uefavg[itr],zero); CHKERRQ(ierr);
			}
		  }
		  if (doWriteBC) {
			for (itr=0; itr<numTracers; itr++) {		  
			  ierr = VecScale(bcavg[itr],1.0/avgTimer.count);CHKERRQ(ierr);
			  ierr = VecView(bcavg[itr],fdbcavgout[itr]);CHKERRQ(ierr);
			  ierr = VecSet(bcavg[itr],zero); CHKERRQ(ierr);
			}
		  }
		  ierr = updateStepTimer("avg_", Iter0+iLoop, &avgTimer);CHKERRQ(ierr);
        }
      }
    }
#endif    
  }  /* end of time-stepping loop */

#if !defined (FORSPINUP) && !defined (FORJACOBIAN)
  for (itr=0; itr<numTracers; itr++) {
    ierr = PetscViewerDestroy(&fdout[itr]);CHKERRQ(ierr);
  }
  ierr = PetscFClose(PETSC_COMM_WORLD,fptime);CHKERRQ(ierr);

  ierr = finalizeTMMWrite(tc,maxSteps,numTracers);CHKERRQ(ierr);

  if (doWriteUF) {
	for (itr=0; itr<numTracers; itr++) {
	  ierr = PetscViewerDestroy(&fdufout[itr]);CHKERRQ(ierr);
	}  
  }

  if (doWriteUEF) {
	for (itr=0; itr<numTracers; itr++) {
	  ierr = PetscViewerDestroy(&fduefout[itr]);CHKERRQ(ierr);
	}  
  }
  
  if (doWriteBC) {
	for (itr=0; itr<numTracers; itr++) {
	  ierr = PetscViewerDestroy(&fdbcout[itr]);CHKERRQ(ierr);
	}  
  }
#endif

  if (doTimeAverage) {
	ierr = PetscFClose(PETSC_COMM_WORLD,avgfptime);CHKERRQ(ierr);
	for (itr=0; itr<numTracers; itr++) {
	  ierr = PetscViewerDestroy(&fdavgout[itr]);CHKERRQ(ierr);
	}
	ierr = VecDestroyVecs(numTracers,&vavg);CHKERRQ(ierr);
	if (doWriteUF) {
	  for (itr=0; itr<numTracers; itr++) {
		ierr = PetscViewerDestroy(&fdufavgout[itr]);CHKERRQ(ierr);
	  }
	  ierr = VecDestroyVecs(numTracers,&ufavg);CHKERRQ(ierr);
	}
	if (doWriteUEF) {
	  for (itr=0; itr<numTracers; itr++) {
		ierr = PetscViewerDestroy(&fduefavgout[itr]);CHKERRQ(ierr);
	  }
	  ierr = VecDestroyVecs(numTracers,&uefavg);CHKERRQ(ierr);
	}
	if (doWriteBC) {
	  for (itr=0; itr<numTracers; itr++) {
		ierr = PetscViewerDestroy(&fdbcavgout[itr]);CHKERRQ(ierr);
	  }
	  ierr = VecDestroyVecs(numTracers,&bcavg);CHKERRQ(ierr);
	}
  }
  
/* write final pickup */  
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,pickupoutFile,FILE_MODE_WRITE,&fdp);CHKERRQ(ierr);
  for (itr=0; itr<numTracers; itr++) {
    ierr = VecView(v[itr],fdp);CHKERRQ(ierr);
  }
  ierr = PetscViewerDestroy(&fdp);CHKERRQ(ierr);      
  
  ierr=PetscTime(&t2); CHKERRQ(ierr); /* stop counting wall clock time */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Wall clock time: %10.5f\n", t2-t1);CHKERRQ(ierr); 
  
  /* Free data structures */
  ierr = VecDestroyVecs(numTracers,&v);CHKERRQ(ierr);  
  ierr = VecDestroyVecs(numTracers,&vtmp);CHKERRQ(ierr);
  ierr = MatDestroy(&Ae);CHKERRQ(ierr);
  ierr = MatDestroy(&Ai);CHKERRQ(ierr);

  if (periodicMatrix) {
    ierr = destroyPeriodicMat(&Aep);CHKERRQ(ierr);
    ierr = destroyPeriodicMat(&Aip);CHKERRQ(ierr);    
  } else if (timeDependentMatrix) {
	ierr = destroyTimeDependentMat(&Aetd);
	ierr = destroyTimeDependentMat(&Aitd);
  }

  if (rescaleForcing) {
	ierr = VecDestroy(&Rfs);CHKERRQ(ierr);
	if (periodicMatrix) {
      ierr = destroyPeriodicVec(&Rfsp);CHKERRQ(ierr);
    } else if (timeDependentMatrix) {
      ierr = destroyTimeDependentVec(&Rfstd);CHKERRQ(ierr);
    }
  }
  
  if (useMonitor) {
	ierr = finalizeMonitor(tc,maxSteps,numTracers);CHKERRQ(ierr);
  }

  if (doMisfit) {
	ierr = finalizeMisfit(tc,maxSteps,numTracers);CHKERRQ(ierr);
  }
  
  if (useExternalForcing) {
    ierr = VecDestroyVecs(numTracers,&uef);CHKERRQ(ierr);  
    ierr = finalizeExternalForcing(tc,maxSteps,numTracers);CHKERRQ(ierr);
  }  
  
  if (useForcingFromFile) {
    ierr = VecDestroyVecs(numTracers,&uf);CHKERRQ(ierr);  
	if (periodicForcing) {
	  for (itr=0; itr<numTracers; itr++) {  
		ierr = destroyPeriodicVec(&up[itr]);CHKERRQ(ierr);
	  }
	} else if (timeDependentForcing) {
	  for (itr=0; itr<numTracers; itr++) {
	    ierr = destroyTimeDependentVec(&utdf[itr]);
	  }  
	}
  }
  
  if (usePrescribedBC) {
    ierr = MatDestroy(&Be);CHKERRQ(ierr);
    ierr = MatDestroy(&Bi);CHKERRQ(ierr);
	if (periodicMatrix) {
	  ierr = destroyPeriodicMat(&Bep);CHKERRQ(ierr);
	  ierr = destroyPeriodicMat(&Bip);CHKERRQ(ierr);    
	} else if (timeDependentMatrix) {
	  ierr = destroyTimeDependentMat(&Betd);
	  ierr = destroyTimeDependentMat(&Bitd);
	}
    ierr = VecDestroyVecs(numTracers,&bcc);CHKERRQ(ierr);
    ierr = VecDestroyVecs(numTracers,&bcf);CHKERRQ(ierr);  
	if (periodicBC) {
	  for (itr=0; itr<numTracers; itr++) {  
		ierr = destroyPeriodicVec(&bcp[itr]);CHKERRQ(ierr);
	  }
	} else if (timeDependentBC) {
	  for (itr=0; itr<numTracers; itr++) {
	    ierr = destroyTimeDependentVec(&bctd[itr]);
	  }  
	} else if (doCalcBC) {
      ierr = finalizeCalcBC(tc,maxSteps,numTracers);CHKERRQ(ierr);
      ierr = PetscFree(gBCIndices);CHKERRQ(ierr);
	}
  }

  ierr = PetscFree(gIndices);CHKERRQ(ierr);

  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
