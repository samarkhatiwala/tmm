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
#include "tmm_profile_utils.h"
#include "tmm_profile_data.h"
#include "tmm_monitor.h"
#include "tmm_main.h"

PetscScalar deltaTClock, time0;
PetscInt maxSteps, Iter0, writeSteps;
PetscInt *gIndices, gLow, gHigh;

extern PetscErrorCode forwardStep(PetscScalar tc, PetscInt iLoop, PetscScalar dt, PetscInt numTracers, 
                                  PetscBool useForcingFromFile, PetscBool useExternalForcing, PetscBool usePrescribedBC, 
                                  Vec *v, Mat Ae, Mat Ai, Mat Be, Mat Bi, Vec *uf, Vec *uef, Vec *bcc, Vec *bcf, Vec *vtmp);

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
  char mateFile[PETSC_MAX_PATH_LEN], matiFile[PETSC_MAX_PATH_LEN], rfsFile[PETSC_MAX_PATH_LEN];
  PetscScalar *tdpMatrix; /* array for periodic matrix times */
  PetscInt numMatrixPeriods;
  PetscBool periodicMatrix = PETSC_FALSE;
  PetscScalar matrixCyclePeriod, matrixCycleStep;

/* Forcing */
  Vec *uf, *uef;
  PeriodicVec up[MAXNUMTRACERS];
  char *forcingFile[MAXNUMTRACERS];
  PetscInt numForcingPeriods, numForcing;
  PetscScalar *tdpForcing; /* array for periodic forcing */
  Vec **utdf;
  PetscScalar *tdfT; /* array for time dependent (nonperiodic) forcing */
  PetscScalar forcingCyclePeriod, forcingCycleStep;
  PetscScalar tf0, tf1;

/* Rescale forcing */
  Vec Rfs;
  PeriodicVec Rfsp;
  
/* BC's */
  Vec *bcc, *bcf;
  PeriodicVec bcp[MAXNUMTRACERS];
  char *bcFile[MAXNUMTRACERS];
  PetscInt numBCPeriods, numBC;
  PetscScalar *tdpBC; /* array for periodic forcing */  
  Vec **bctd;
  PetscScalar *tdbcT; /* array for time dependent (nonperiodic) forcing */
  PetscScalar bcCyclePeriod, bcCycleStep;  
  PetscScalar tbc0, tbc1;
  PetscInt bcCutOffStep = -1;
  Mat Be, Bi;
  PeriodicMat Bep, Bip;
  char matbeFile[PETSC_MAX_PATH_LEN], matbiFile[PETSC_MAX_PATH_LEN];
  Vec bcTemplateVec;
  PetscInt lBCSize, gBCSize;
  
/* I/O   */
  char *iniFile[MAXNUMTRACERS];  
  char *outFile[MAXNUMTRACERS];  
  char pickupFile[PETSC_MAX_PATH_LEN];
  char pickupoutFile[PETSC_MAX_PATH_LEN];  
  char outTimeFile[PETSC_MAX_PATH_LEN];  
  PetscBool appendOutput = PETSC_FALSE;    
  PetscBool pickupFromFile = PETSC_FALSE;
  PetscBool doTimeAverage = PETSC_FALSE;
  PetscInt avgNumTimeSteps, avgStartTimeStep, avgCount;
  char *avgOutFile[MAXNUMTRACERS];  
  Vec *vavg;
  PetscViewer fdavgout[MAXNUMTRACERS];
  FILE *fptime;
  PetscViewer fd, fdp, fdout[MAXNUMTRACERS], fdin[MAXNUMTRACERS];

/* run time options */
  PetscBool useExternalForcing = PETSC_FALSE;
  PetscBool useForcingFromFile = PETSC_FALSE;
  PetscBool usePrescribedBC = PETSC_FALSE;
  PetscBool applyBC = PETSC_FALSE;  
  PetscBool periodicForcing = PETSC_FALSE;
  PetscBool timeDependentForcing = PETSC_FALSE;
  PetscBool constantForcing = PETSC_FALSE;
  PetscBool periodicBC = PETSC_FALSE;
  PetscBool timeDependentBC = PETSC_FALSE;
  PetscBool constantBC = PETSC_FALSE;
  PetscBool doCalcBC = PETSC_FALSE;
  PetscBool rescaleForcing = PETSC_FALSE;
  PetscBool useMonitor = PETSC_FALSE;

  PetscMPIInt numProcessors, myId;  
  PetscErrorCode ierr;
  PetscBool flg1,flg2;
  PetscScalar t1, t2, tc, tf;
  PetscInt iLoop, Iterc;
  PetscInt it;
  PetscInt itr, maxValsToRead;
  PetscInt itjac;
  PetscInt fp;
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
  ierr = PetscOptionsGetInt(PETSC_NULL,"-numtracers",&numTracers,&flg1);CHKERRQ(ierr);
  if (numTracers>MAXNUMTRACERS) {
   SETERRQ(PETSC_COMM_WORLD,1,"Number of tracers exceeds maximum allowable. Please increase the variable MAXNUMTRACERS and recompile");
  }  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of tracers to be integrated: %d\n", numTracers);CHKERRQ(ierr); 

/* Time step data */
  ierr = PetscOptionsGetReal(PETSC_NULL,"-deltat_clock",&deltaTClock,&flg1);CHKERRQ(ierr);
  ierr = PetscOptionsGetReal(PETSC_NULL,"-t0",&time0,&flg1);CHKERRQ(ierr);  
  ierr = PetscOptionsGetInt(PETSC_NULL,"-iter0",&Iter0,&flg2);CHKERRQ(ierr);
  if ((flg1) && (!flg2)) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate both or neither time0 and Iter0 with the -t0 and -iter0 flags");
  if ((!flg1) && (flg2)) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate both or neither time0 and Iter0 with the -t0 and -iter0 flags");        
  ierr = PetscOptionsGetInt(PETSC_NULL,"-max_steps",&maxSteps,&flg1);CHKERRQ(ierr);
  if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate maximum number of steps with the -max_steps option");
  ierr = PetscOptionsGetInt(PETSC_NULL,"-write_steps",&writeSteps,&flg1);CHKERRQ(ierr);
  if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate output step with the -write_steps option");

/*Data for time averaging */
  ierr = PetscOptionsHasName(PETSC_NULL,"-time_avg",&doTimeAverage);CHKERRQ(ierr);
  if (doTimeAverage) {  
    ierr = PetscOptionsGetInt(PETSC_NULL,"-avg_start_time_step",&avgStartTimeStep,&flg1);CHKERRQ(ierr);
    if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate (absolute) time step at which to start time averaging with the -avg_start_time_step flag");
    ierr = PetscOptionsGetInt(PETSC_NULL,"-avg_time_steps",&avgNumTimeSteps,&flg1);CHKERRQ(ierr);
    if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate number of time averaging time steps with the -avg_time_step flag");
	for (itr=0; itr<numTracers; itr++) {
	  avgOutFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
	}
	maxValsToRead = numTracers;
	ierr = PetscOptionsGetStringArray(PETSC_NULL,"-avg_files",avgOutFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
	if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate file name(s) for writing time averages with the -avg_files option");
	if (maxValsToRead != numTracers) {
	  SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of time average file names specified");
	}  
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Time averages will be computed starting at (and including) time step: %d\n", avgStartTimeStep);CHKERRQ(ierr);	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Time averages will be computed over %d time steps\n", avgNumTimeSteps);CHKERRQ(ierr);	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Time averages will be written to:\n");CHKERRQ(ierr);
	for (itr=0; itr<numTracers; itr++) {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,avgOutFile[itr]);CHKERRQ(ierr);
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
  ierr = PetscOptionsGetString(PETSC_NULL,"-me",mateFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
  if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate binary matrix file with the -me option");
  ierr = PetscOptionsGetString(PETSC_NULL,"-mi",matiFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
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
  
  ierr = PetscOptionsHasName(PETSC_NULL,"-periodic_matrix",&periodicMatrix);CHKERRQ(ierr);  
  if (periodicMatrix) {    
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic matrices specified\n");CHKERRQ(ierr);
    ierr=PetscStrcat(mateFile,"_");CHKERRQ(ierr);
    ierr=PetscStrcat(matiFile,"_");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Ae basename is %s\n", mateFile);CHKERRQ(ierr); 
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Ai basename is %s\n", matiFile);CHKERRQ(ierr);     
    
/*  read time data */
    ierr = PetscOptionsGetReal(PETSC_NULL,"-matrix_cycle_period",&matrixCyclePeriod,&flg2);CHKERRQ(ierr);
    if (!flg2) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate matrix cycling time with the -matrix_cycle_period option");
    ierr = PetscOptionsGetReal(PETSC_NULL,"-matrix_cycle_step",&matrixCycleStep,&flg2);CHKERRQ(ierr);
    if (!flg2) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate cycling step with the -matrix_cycle_step option");
    numMatrixPeriods=matrixCyclePeriod/matrixCycleStep;
/*  array for holding extended time array */
    PetscMalloc((numMatrixPeriods+2)*sizeof(PetscScalar), &tdpMatrix); 
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Periodic matrix specified at times:\n");CHKERRQ(ierr);            
    for (it=0; it<=numMatrixPeriods+1; it++) {
      tdpMatrix[it]=(-matrixCycleStep/2.0) + it*matrixCycleStep;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"tdpMatrix=%10.5f\n", tdpMatrix[it]);CHKERRQ(ierr);        
    }    
    
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
  ierr = PetscOptionsGetStringArray(PETSC_NULL,"-o",outFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
  if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate output file name(s) with the -o option");
  if (maxValsToRead != numTracers) {
    SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of outfile file names specified");
  }  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Output will be written to:\n");CHKERRQ(ierr);
  for (itr=0; itr<numTracers; itr++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,outFile[itr]);CHKERRQ(ierr);
  }  
  ierr = PetscOptionsHasName(PETSC_NULL,"-append",&appendOutput);CHKERRQ(ierr);
  if (appendOutput) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Output will be appended\n");CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Output will overwrite existing file(s)\n");CHKERRQ(ierr);
  }    

/* Output times */
  ierr = PetscOptionsGetString(PETSC_NULL,"-time_file",outTimeFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
  if (!flg1) {
	strcpy(outTimeFile,"");
    sprintf(outTimeFile,"%s","output_time.txt");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Output times will be written to %s\n",outTimeFile);CHKERRQ(ierr);

/* File name for final pickup */
  ierr = PetscOptionsGetString(PETSC_NULL,"-pickup_out",pickupoutFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
  if (!flg1) {
	strcpy(pickupoutFile,"");
    sprintf(pickupoutFile,"%s","pickup.petsc");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Final pickup will be written to %s\n",pickupoutFile);CHKERRQ(ierr);

/* tracer vectors */
  ierr = VecDuplicateVecs(templateVec,numTracers,&v);CHKERRQ(ierr);
  ierr = VecDuplicateVecs(templateVec,numTracers,&vtmp);CHKERRQ(ierr); /* temporary work space needed by forwardStep */

/* Initial condition     */
  for (itr=0; itr<numTracers; itr++) {
    iniFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
  }
  ierr = PetscOptionsGetString(PETSC_NULL,"-pickup",pickupFile,PETSC_MAX_PATH_LEN-1,&pickupFromFile);CHKERRQ(ierr);
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
    ierr = PetscOptionsGetStringArray(PETSC_NULL,"-i",iniFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
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
  ierr = PetscOptionsHasName(PETSC_NULL,"-forcing_from_file",&useForcingFromFile);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL,"-prescribed_bc",&usePrescribedBC);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL,"-external_forcing",&useExternalForcing);CHKERRQ(ierr);

/*   if (useExternalForcing) { */
/*     if ((useForcingFromFile) | (usePrescribedBC)) { */
/*       SETERRQ(PETSC_COMM_WORLD,1,"Cannot use the -external_forcing option with the -forcing_from_file or -prescribed_bc options");   */
/*     }   */
/*   } */

  if (useForcingFromFile) {  
  	ierr=PetscPrintf(PETSC_COMM_WORLD,"Forcing from file(s) specified\n");CHKERRQ(ierr);  
	for (itr=0; itr<numTracers; itr++) {
	  forcingFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
	}
	maxValsToRead = numTracers;
	ierr = PetscOptionsGetStringArray(PETSC_NULL,"-forcing_files",forcingFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
	if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"MUST specify forcing files with the -forcing_files option");
    if (maxValsToRead != numTracers) {
      SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of forcing file names specified");
    }    
    ierr = VecDuplicateVecs(templateVec,numTracers,&uf);CHKERRQ(ierr);    
/*  There are 3 possibilities: periodic, constant, and time-dependent forcing */
    ierr = PetscOptionsHasName(PETSC_NULL,"-periodic_forcing",&periodicForcing);CHKERRQ(ierr);
    if (periodicForcing) {
/*    Read some info. Forcing is read in interpPeriodicForcing. */
      numForcing=-1;
      ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic forcing from file(s) specified\n");CHKERRQ(ierr);
      for (itr=0; itr<numTracers; itr++) {
        ierr=PetscStrcat(forcingFile[itr],"_");CHKERRQ(ierr);        
	    ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d periodic forcing basename is %s\n",itr,forcingFile[itr]);CHKERRQ(ierr); 
        up[itr].firstTime = PETSC_TRUE; /* initialize periodic vector */        	    
      }      

/*    read time data */
      ierr = PetscOptionsGetReal(PETSC_NULL,"-forcing_cycle_period",&forcingCyclePeriod,&flg2);CHKERRQ(ierr);
      if (!flg2) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate forcing cycling time with the -forcing_cycle_period option");
      ierr = PetscOptionsGetReal(PETSC_NULL,"-forcing_cycle_step",&forcingCycleStep,&flg2);CHKERRQ(ierr);
      if (!flg2) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate cycling step with the -forcing_cycle_step option");
      numForcingPeriods=forcingCyclePeriod/forcingCycleStep;
/*    array for holding extended time array */
      PetscMalloc((numForcingPeriods+2)*sizeof(PetscScalar), &tdpForcing); 
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Periodic forcing specified at times:\n");CHKERRQ(ierr);            
      for (it=0; it<=numForcingPeriods+1; it++) {
        tdpForcing[it]=(-forcingCycleStep/2.0) + it*forcingCycleStep;
        ierr = PetscPrintf(PETSC_COMM_WORLD,"tdpForcing=%10.5f\n", tdpForcing[it]);CHKERRQ(ierr);        
      }
    } else { /* constant or (nonperiodic) time dependent forcing */
/*    Read info AND forcing here */    
      ierr = PetscOptionsHasName(PETSC_NULL,"-time_dependent_forcing",&timeDependentForcing);CHKERRQ(ierr);
      if (timeDependentForcing) {      
        ierr=PetscPrintf(PETSC_COMM_WORLD,"Time dependent forcing specified\n");CHKERRQ(ierr);
        ierr = PetscOptionsGetInt(PETSC_NULL,"-number_forcing_vecs",&numForcing,&flg2);CHKERRQ(ierr);
        if (!flg2) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate number of forcing vectors in file with the -number_forcing_vecs option");
/*      Make sure we have all the time info         */
        ierr = PetscOptionsGetReal(PETSC_NULL,"-t0",&time0,&flg2);CHKERRQ(ierr);
        if (!flg2) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate initial time with the -t0 option");
        ierr = PetscOptionsGetReal(PETSC_NULL,"-deltat_clock",&deltaTClock,&flg2);CHKERRQ(ierr);
        if (!flg2) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate time step with the -deltat_clock option");      
        ierr = PetscOptionsGetReal(PETSC_NULL,"-tfini",&tf0,&flg2);CHKERRQ(ierr);
        if (!flg2) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate first forcing time with the -tfini option");
        ierr = PetscOptionsGetReal(PETSC_NULL,"-tfend",&tf1,&flg2);CHKERRQ(ierr);
        if (!flg2) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate last forcing time with the -tfend option");        
        utdf = malloc(numTracers*sizeof(Vec *));
        for (itr=0; itr<numTracers; itr++) {   
          ierr = VecDuplicateVecs(templateVec,numForcing,&utdf[itr]);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading forcing from file %s\n", itr,forcingFile[itr]);CHKERRQ(ierr);
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,forcingFile[itr],FILE_MODE_READ,&fd);CHKERRQ(ierr);
          for (it=0; it<numForcing; it++) {
            ierr = VecLoad(utdf[itr][it],fd);CHKERRQ(ierr); /* IntoVector */
          }
          ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"   Read %d forcing vectors\n", numForcing);CHKERRQ(ierr);            
        }
        PetscMalloc(numForcing*sizeof(PetscScalar), &tdfT); 
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Time dependent forcing specified at times:\n");CHKERRQ(ierr);        
        for (it=0; it<numForcing; it++) {
          tdfT[it] = tf0 + ((tf1-tf0)/(numForcing-1.0))*it;
          ierr = PetscPrintf(PETSC_COMM_WORLD,"tdfT=%10.5f\n", tdfT[it]);CHKERRQ(ierr);        
        }
      } else { /* constant forcing */
        ierr=PetscPrintf(PETSC_COMM_WORLD,"Constant forcing specified\n");CHKERRQ(ierr);      
        constantForcing = PETSC_TRUE;
        numForcing=1;
		for (itr=0; itr<numTracers; itr++) {   
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading forcing from file %s\n", itr,forcingFile[itr]);CHKERRQ(ierr);
		  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,forcingFile[itr],FILE_MODE_READ,&fd);CHKERRQ(ierr);
		  ierr = VecLoad(uf[itr],fd);CHKERRQ(ierr); /* IntoVector */
		  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"   Read %d forcing vectors\n", numForcing);CHKERRQ(ierr);            
		}
      }  /* TD/constant forcing */
    }  /* periodic/nonperiodic forcing */
  } else {  /* no forcing */
    numForcing=0;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"No forcing from file(s) specified\n");CHKERRQ(ierr);    
  }

  if (useExternalForcing) {  /* external forcing present */  
    ierr = PetscPrintf(PETSC_COMM_WORLD,"External forcing is being used\n");CHKERRQ(ierr);    
    ierr = VecDuplicateVecs(templateVec,numTracers,&uef);CHKERRQ(ierr);        
    ierr = iniExternalForcing(time0,Iter0,numTracers,v,uef);CHKERRQ(ierr);
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"No external forcing is being used\n");CHKERRQ(ierr);  
  }  

/* Prescribed BCs   */
  if (usePrescribedBC) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Prescribed BC's specified\n");CHKERRQ(ierr);    

	ierr = VecCreate(PETSC_COMM_WORLD,&bcTemplateVec);CHKERRQ(ierr);
    
    ierr = PetscOptionsHasName(PETSC_NULL,"-calc_bc",&doCalcBC);CHKERRQ(ierr);
    if (doCalcBC) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"BCs will be calculated\n");CHKERRQ(ierr);    

      if ((useProfiles) && (numProcessors>1)) {
        lBCSize = lNumProfiles;
        ierr = VecSetSizes(bcTemplateVec,lBCSize,PETSC_DECIDE);CHKERRQ(ierr);      
        ierr = VecSetFromOptions(bcTemplateVec);CHKERRQ(ierr);        
      } else {
        ierr = PetscOptionsGetInt(PETSC_NULL,"-bc_vec_size",&gBCSize,&flg1);CHKERRQ(ierr);
        if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate size of BC vector with the -bc_vec_size option");
        ierr = VecSetSizes(bcTemplateVec,PETSC_DECIDE,gBCSize);CHKERRQ(ierr);      
        ierr = VecSetFromOptions(bcTemplateVec);CHKERRQ(ierr);
      }
      
      ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&bcc);CHKERRQ(ierr);    
      ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&bcf);CHKERRQ(ierr);    
      
      ierr = iniCalcBC(time0,Iter0,time0+deltaTClock,Iter0+1,numTracers,v,bcc,bcf);CHKERRQ(ierr);
      
    } else { /* read from file */
      for (itr=0; itr<numTracers; itr++) {
        bcFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
      }
      maxValsToRead = numTracers;
      ierr = PetscOptionsGetStringArray(PETSC_NULL,"-bc_files",bcFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
      if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"MUST specify BC files with the -bc_files option");
      if (maxValsToRead != numTracers) {
        SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of BC file names specified");
      }    
/*     ierr = VecDuplicateVecs(templateVec,numTracers,&bcc);CHKERRQ(ierr);     */
/*     ierr = VecDuplicateVecs(templateVec,numTracers,&bcf);CHKERRQ(ierr);     */

/*    There are 3 possibilities: periodic, constant, and time-dependent BCs */
      ierr = PetscOptionsHasName(PETSC_NULL,"-periodic_bc",&periodicBC);CHKERRQ(ierr);
      if (periodicBC) {
/*      Read some info. BC is read in interpPeriodicVector */
        ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic BC from file(s) specified\n");CHKERRQ(ierr);
        for (itr=0; itr<numTracers; itr++) {
          ierr=PetscStrcat(bcFile[itr],"_");CHKERRQ(ierr);        
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d periodic BC basename is %s\n",itr,bcFile[itr]);CHKERRQ(ierr); 
          bcp[itr].firstTime = PETSC_TRUE; /* initialize periodic vector */        	    
        }      

/*      Load one vector here as a template */
        strcpy(tmpFile,"");
        sprintf(tmpFile,"%s%02d",bcFile[0],0);
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
        ierr = VecLoad(bcTemplateVec,fd);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
        ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&bcc);CHKERRQ(ierr);    
        ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&bcf);CHKERRQ(ierr);    

/*      read time data */
        ierr = PetscOptionsGetReal(PETSC_NULL,"-bc_cycle_period",&bcCyclePeriod,&flg2);CHKERRQ(ierr);
        if (!flg2) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate BC cycling time with the -bc_cycle_period option");
        ierr = PetscOptionsGetReal(PETSC_NULL,"-bc_cycle_step",&bcCycleStep,&flg2);CHKERRQ(ierr);
        if (!flg2) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate cycling step with the -bc_cycle_step option");
        numBCPeriods=bcCyclePeriod/bcCycleStep;
/*      array for holding extended time array */
        PetscMalloc((numBCPeriods+2)*sizeof(PetscScalar), &tdpBC); 
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Periodic BC specified at times:\n");CHKERRQ(ierr);            
        for (it=0; it<=numBCPeriods+1; it++) {
          tdpBC[it]=(-bcCycleStep/2.0) + it*bcCycleStep;
          ierr = PetscPrintf(PETSC_COMM_WORLD,"tdpBC=%10.5f\n", tdpBC[it]);CHKERRQ(ierr);        
        }
      } else { /* constant or (nonperiodic) time dependent BC */
/*      Read info AND BC here */    
        ierr = PetscOptionsHasName(PETSC_NULL,"-time_dependent_bc",&timeDependentBC);CHKERRQ(ierr);
        if (timeDependentBC) {      
          ierr=PetscPrintf(PETSC_COMM_WORLD,"Time dependent BC specified\n");CHKERRQ(ierr);
          ierr = PetscOptionsGetInt(PETSC_NULL,"-number_bc_vecs",&numBC,&flg2);CHKERRQ(ierr);
          if (!flg2) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate number of BC vectors in file with the -number_bc_vecs option");
/*        Make sure we have all the time info         */
          ierr = PetscOptionsGetReal(PETSC_NULL,"-t0",&time0,&flg2);CHKERRQ(ierr);
          if (!flg2) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate initial time with the -t0 option");
          ierr = PetscOptionsGetReal(PETSC_NULL,"-deltat_clock",&deltaTClock,&flg2);CHKERRQ(ierr);
          if (!flg2) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate time step with the -deltat_clock option");      
          ierr = PetscOptionsGetReal(PETSC_NULL,"-tbcini",&tbc0,&flg2);CHKERRQ(ierr);
          if (!flg2) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate first BC time with the -tbcini option");
          ierr = PetscOptionsGetReal(PETSC_NULL,"-tbcend",&tbc1,&flg2);CHKERRQ(ierr);
          if (!flg2) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate last BC time with the -tbcend option");        

/*        Load one vector here as a template */
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,bcFile[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
          ierr = VecLoad(bcTemplateVec,fd);CHKERRQ(ierr);
          ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
          ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&bcc);CHKERRQ(ierr);    
          ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&bcf);CHKERRQ(ierr);    
        
          bctd = malloc(numTracers*sizeof(Vec *));
          for (itr=0; itr<numTracers; itr++) {   
            ierr = VecDuplicateVecs(bcTemplateVec,numBC,&bctd[itr]);CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading BC from file %s\n", itr,bcFile[itr]);CHKERRQ(ierr);
            ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,bcFile[itr],FILE_MODE_READ,&fd);CHKERRQ(ierr);
            for (it=0; it<numBC; it++) {
              ierr = VecLoad(bctd[itr][it],fd);CHKERRQ(ierr); /* IntoVector */
            }
            ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"   Read %d BC vectors\n", numBC);CHKERRQ(ierr);            
          }
          PetscMalloc(numBC*sizeof(PetscScalar), &tdbcT); 
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Time dependent BC specified at times:\n");CHKERRQ(ierr);        
          for (it=0; it<numBC; it++) {
            tdbcT[it] = tbc0 + ((tbc1-tbc0)/(numBC-1.0))*it;
            ierr = PetscPrintf(PETSC_COMM_WORLD,"tdbcT=%10.5f\n", tdbcT[it]);CHKERRQ(ierr);        
          }
        } else { /* constant BC */
          ierr=PetscPrintf(PETSC_COMM_WORLD,"Constant BC specified\n");CHKERRQ(ierr);      
          constantBC = PETSC_TRUE;
          numBC=1;
/*        Load one vector here as a template */
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,bcFile[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
          ierr = VecLoad(bcTemplateVec,fd);CHKERRQ(ierr);
          ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
          ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&bcc);CHKERRQ(ierr);    
          ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&bcf);CHKERRQ(ierr);            
          for (itr=0; itr<numTracers; itr++) {   
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading BC from file %s\n", itr,bcFile[itr]);CHKERRQ(ierr);
            ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,bcFile[itr],FILE_MODE_READ,&fd);CHKERRQ(ierr);
            ierr = VecLoad(bcc[itr],fd);CHKERRQ(ierr); /* IntoVector */
            ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
            ierr = PetscPrintf(PETSC_COMM_WORLD,"   Read %d BC vectors\n", numBC);CHKERRQ(ierr);            
            ierr = VecCopy(bcc[itr],bcf[itr]);CHKERRQ(ierr);		  
          }
        }  /* TD/constant BC */
      }  /* periodic/nonperiodic BC */
    }

    ierr = VecGetLocalSize(bcTemplateVec,&lBCSize);CHKERRQ(ierr);

    if ((useProfiles) && (numProcessors>1)) {
      if (lBCSize != lNumProfiles) {
        SETERRQ(PETSC_COMM_WORLD,1,"Problem with partitioning of BC vectors! lNumProfiles must equal lBCSize");
      }
    }
    
    applyBC = PETSC_TRUE;

    ierr = PetscOptionsGetInt(PETSC_NULL,"-bc_cutoff_step",&bcCutOffStep,&flg1);CHKERRQ(ierr);
    if (bcCutOffStep>0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Prescribed BC will be turned off after time step %d\n",bcCutOffStep);CHKERRQ(ierr);    
    }
    
/*  Matrices */
	ierr = PetscOptionsGetString(PETSC_NULL,"-mbe",matbeFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
	if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate binary boundary matrix file name with the -mbe option");
	ierr = PetscOptionsGetString(PETSC_NULL,"-mbi",matbiFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
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
    numBC=0;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"No prescribed BC's specified\n");CHKERRQ(ierr);  
  }  

  ierr = PetscOptionsGetString(PETSC_NULL,"-rescale_forcing_file",rfsFile,PETSC_MAX_PATH_LEN-1,&rescaleForcing);CHKERRQ(ierr);
  if (rescaleForcing) {  
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Forcing will be rescaled\n");CHKERRQ(ierr); 	    
    ierr = VecDuplicate(templateVec,&Rfs);CHKERRQ(ierr);    
	if (periodicMatrix) {    
	  ierr=PetscStrcat(rfsFile,"_");CHKERRQ(ierr);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Rescale forcing factor file basename is %s\n", rfsFile);CHKERRQ(ierr); 	
	  Rfsp.firstTime = PETSC_TRUE;
	} else {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading rescale forcing factor from file %s\n", rfsFile);CHKERRQ(ierr);  	
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,rfsFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = VecLoad(Rfs,fd);CHKERRQ(ierr);  
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
	}  
  }

/* initialize monitor */
  ierr = PetscOptionsHasName(PETSC_NULL,"-monitor",&useMonitor);CHKERRQ(ierr);
  if (useMonitor) {  
    ierr = iniMonitor(time0,Iter0,numTracers,v);CHKERRQ(ierr);
  }

/* Open files for output and optionally write initial conditions */
#if !defined (FORSPINUP) && !defined (FORJACOBIAN)
  if (!appendOutput) {
    ierr = PetscFOpen(PETSC_COMM_WORLD,outTimeFile,"w",&fptime);CHKERRQ(ierr);  
    ierr = PetscFPrintf(PETSC_COMM_WORLD,fptime,"%d   %10.5f\n",Iter0,time0);CHKERRQ(ierr);     
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing output at time %10.5f, step %d\n", time0,Iter0);CHKERRQ(ierr);  
    for (itr=0; itr<numTracers; itr++) {       
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,outFile[itr],FILE_MODE_WRITE,&fdout[itr]);CHKERRQ(ierr);
      ierr = VecView(v[itr],fdout[itr]);CHKERRQ(ierr);
    }
  } else {
	  ierr = PetscFOpen(PETSC_COMM_WORLD,outTimeFile,"a",&fptime);CHKERRQ(ierr);  
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Opening file(s) for output. Initial condition will NOT be written\n");CHKERRQ(ierr);  
      for (itr=0; itr<numTracers; itr++) {       
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,outFile[itr],FILE_MODE_APPEND,&fdout[itr]);CHKERRQ(ierr);
      }
  }
#endif

  if (doTimeAverage) {  
    ierr = VecDuplicateVecs(templateVec,numTracers,&vavg);CHKERRQ(ierr);  
	for (itr=0; itr<numTracers; itr++) {
	  ierr = VecSet(vavg[itr],zero); CHKERRQ(ierr);
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,avgOutFile[itr],FILE_MODE_WRITE,&fdavgout[itr]);CHKERRQ(ierr);
    }    
    avgCount=0;
  }

/* reinitialize forcing if required */
#ifdef FORSPINUP
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"itjac.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fp,&itjac,1,PETSC_INT);CHKERRQ(ierr);  
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
  ierr = PetscBinaryRead(fp,&itjac,1,PETSC_INT);CHKERRQ(ierr);  
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
      ierr = interpPeriodicMatrix(tc,&Ae,matrixCyclePeriod,numMatrixPeriods,
                                  tdpMatrix,&Aep,mateFile);
      ierr = interpPeriodicMatrix(tc,&Ai,matrixCyclePeriod,numMatrixPeriods,
                                  tdpMatrix,&Aip,matiFile);
      if (rescaleForcing) {
		ierr = interpPeriodicVector(tc,&Rfs,matrixCyclePeriod,numMatrixPeriods,tdpMatrix,&Rfsp,rfsFile);
      }
    }

/*  Forcing     */
    if (useForcingFromFile) {
	  if (periodicForcing) {
		for (itr=0; itr<numTracers; itr++) {    
		  ierr = interpPeriodicVector(tc,&uf[itr],forcingCyclePeriod,numForcingPeriods,tdpForcing,&up[itr],forcingFile[itr]);
		}    
	  } else if (timeDependentForcing) {
		ierr = interpTimeDependentVector(tc,uf,numTracers,numForcing,tdfT,utdf);CHKERRQ(ierr);
	  }
	  if (rescaleForcing) {
		for (itr=0; itr<numTracers; itr++) {    
          ierr = VecPointwiseMult(uf[itr],Rfs,uf[itr]);CHKERRQ(ierr);
        }  
	  }
	}
#endif

    if (useExternalForcing) {
      ierr = calcExternalForcing(tc,Iterc,iLoop,numTracers,v,uef);CHKERRQ(ierr); /* Compute external forcing in uef */
	  if (rescaleForcing) {
		for (itr=0; itr<numTracers; itr++) {    
          ierr = VecPointwiseMult(uef[itr],Rfs,uef[itr]);CHKERRQ(ierr);
        }  
	  }      
    }    

#ifndef FORJACOBIAN
    if (applyBC) {
      if ((bcCutOffStep>0) && (iLoop==(bcCutOffStep+1))) {
        applyBC = PETSC_FALSE;
        doCalcBC = PETSC_FALSE;
      } else {    
        if (periodicMatrix) {
          ierr = interpPeriodicMatrix(tc,&Be,matrixCyclePeriod,numMatrixPeriods,
                                      tdpMatrix,&Bep,matbeFile);
          ierr = interpPeriodicMatrix(tc,&Bi,matrixCyclePeriod,numMatrixPeriods,
                                      tdpMatrix,&Bip,matbiFile);
        }    
        if (periodicBC) {
          for (itr=0; itr<numTracers; itr++) {    
            ierr = interpPeriodicVector(tc,&bcc[itr],bcCyclePeriod,numBCPeriods,tdpBC,&bcp[itr],bcFile[itr]);
            ierr = interpPeriodicVector(tf,&bcf[itr],bcCyclePeriod,numBCPeriods,tdpBC,&bcp[itr],bcFile[itr]);
          }    
        } else if (timeDependentBC) {
          ierr = interpTimeDependentVector(tc,bcc,numTracers,numBC,tdbcT,bctd);CHKERRQ(ierr);
          ierr = interpTimeDependentVector(tf,bcf,numTracers,numBC,tdbcT,bctd);CHKERRQ(ierr);
        } else if (doCalcBC) {
          ierr = calcBC(tc,Iterc,tc+deltaTClock,Iterc+1,iLoop,numTracers,v,bcc,bcf);CHKERRQ(ierr); /* Compute BC in bcc and bcf */	  
        }
      }
	  
    }
        
    ierr = forwardStep(tc,Iterc,deltaTClock,numTracers,useForcingFromFile,useExternalForcing,applyBC,
                       v,Ae,Ai,Be,Bi,uf,uef,bcc,bcf,vtmp);CHKERRQ(ierr);
#endif    
    tc=time0 + deltaTClock*iLoop;  /*  time at end of step */    

    if (useMonitor) {
      ierr = updateMonitor(tc,iLoop,numTracers,v);CHKERRQ(ierr);
      ierr = writeMonitor(tc,iLoop,numTracers,v);CHKERRQ(ierr);
    }

/* write output */
#if !defined (FORSPINUP) && !defined (FORJACOBIAN)
    if (useExternalForcing) {
      ierr = writeExternalForcing(tc,iLoop,numTracers,v,uef);CHKERRQ(ierr);
    }
    if (doCalcBC) {
      ierr = writeBC(tc,iLoop,numTracers,v,bcc,bcf);CHKERRQ(ierr);    
    }
    if ((iLoop % writeSteps)==0) {  /*  time to write out */  /* ??? should this be Iter0+iLoop */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
      ierr = PetscFPrintf(PETSC_COMM_WORLD,fptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);           
/*       ierr = PetscViewerASCIIPrintf(fdtime,"\n%d   %10.5f",Iter0+iLoop,tc);CHKERRQ(ierr);        */
      for (itr=0; itr<numTracers; itr++) {
        ierr = VecView(v[itr],fdout[itr]);CHKERRQ(ierr);
      }      
    }
#else
#ifdef FORSPINUP
    if ((iLoop % writeSteps)==0) {  /*  time to write out */  /* ??? should this be Iter0+iLoop */
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
      ierr = PetscBinaryRead(fp,&itjac,1,PETSC_INT);CHKERRQ(ierr);  
      ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

      ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting iteration number to: %d\n", itjac);CHKERRQ(ierr);

      tc=time0 + deltaTClock*(itjac-1);  /*  current time (time at beginning of step) */
      tf=time0 + deltaTClock*itjac;  /*  future time (time at end of step) */
      Iterc=Iter0+itjac-1;

      if (useExternalForcing) ierr = reInitializeExternalForcing(tc,Iterc,numTracers,v,uef);CHKERRQ(ierr);      
      if (doCalcBC) ierr = reInitializeCalcBC(tc,Iterc,tc+deltaTClock,Iterc+1,numTracers,v,bcc,bcf);CHKERRQ(ierr);              
    }
#endif
#ifdef FORJACOBIAN
    if (doTimeAverage) {
      if ((iLoop % writeSteps)==0) {  /*  time to write out */  /* ??? should this be Iter0+iLoop */      
        if (Iter0+iLoop>=avgStartTimeStep) { /* start time averaging (note: avgStartTimeStep is ABSOLUTE time step) */
          if (avgCount<=avgNumTimeSteps) { /* still within same averaging block so accumulate */
            for (itr=0; itr<numTracers; itr++) {
              ierr = VecAXPY(vavg[itr],one,uef[itr]);CHKERRQ(ierr);
            }          
            avgCount = avgCount+1;
/*           ierr = PetscPrintf(PETSC_COMM_WORLD,"Accumulating: %d\n", iLoop);CHKERRQ(ierr);                       */
          }
          if (avgCount==avgNumTimeSteps) { /* time to write averages to file */
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing time average q at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);                      
            for (itr=0; itr<numTracers; itr++) {
              ierr = VecScale(vavg[itr],1.0/avgCount);CHKERRQ(ierr);
              ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,outFile[itr],FILE_MODE_WRITE,&fdout[itr]);CHKERRQ(ierr);            
              ierr = VecView(vavg[itr],fdout[itr]);CHKERRQ(ierr);              
              ierr = PetscViewerDestroy(&fdout[itr]);CHKERRQ(ierr);		              
              ierr = VecSet(vavg[itr],zero); CHKERRQ(ierr);
            }          
            avgCount = 0;     

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
            ierr = PetscBinaryRead(fp,&itjac,1,PETSC_INT);CHKERRQ(ierr);  
            ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      
            ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting iteration number to: %d\n", itjac);CHKERRQ(ierr);
      
            tc=time0 + deltaTClock*(itjac-1);  /*  current time (time at beginning of step) */
            tf=time0 + deltaTClock*itjac;  /*  future time (time at end of step) */
            Iterc=Iter0+itjac-1;
      
            if (useExternalForcing) ierr = reInitializeExternalForcing(tc,Iterc,numTracers,v,uef);CHKERRQ(ierr);      
          
          } else {
          
            for (itr=0; itr<numTracers; itr++) {
              ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading new initial condition from file %s\n", itr,iniFile[itr]);CHKERRQ(ierr);
              ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,iniFile[itr],FILE_MODE_READ,&fd);CHKERRQ(ierr);
              ierr = VecLoad(v[itr],fd);CHKERRQ(ierr); /* IntoVector */
              ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);          
            }        

            tc=time0 + deltaTClock*(itjac-1);  /* this is now the time at end of time step */
            Iterc=Iter0+itjac-1;    
            if (useExternalForcing) ierr = reInitializeExternalForcing(tc,Iterc,numTracers,v,uef);CHKERRQ(ierr);                  
          }
        }
      }
    } else { /* no time averaging */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing q at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
      for (itr=0; itr<numTracers; itr++) {       
        ierr = VecView(uef[itr],fdout[itr]);CHKERRQ(ierr);
      }
      if ((iLoop % writeSteps)==0) {  /*  time to read new set of initial condition(s) */      
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
      }
      ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"itjac.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
      ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
      ierr = PetscBinaryRead(fp,&itjac,1,PETSC_INT);CHKERRQ(ierr);  
      ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);          

      ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting iteration number to: %d\n", itjac);CHKERRQ(ierr);              
    }	  
#endif
#endif
#ifndef FORJACOBIAN    
    if (doTimeAverage) {
      if (Iter0+iLoop>=avgStartTimeStep) { /* start time averaging (note: avgStartTimeStep is ABSOLUTE time step) */
        if (avgCount<=avgNumTimeSteps) { /* still within same averaging block so accumulate */
/*           ierr = PetscPrintf(PETSC_COMM_WORLD,"Accumulating for time average\n");CHKERRQ(ierr);               */
		  for (itr=0; itr<numTracers; itr++) {
			ierr = VecAXPY(vavg[itr],one,v[itr]);CHKERRQ(ierr);
		  }          
		  avgCount = avgCount+1;
        }
        if (avgCount==avgNumTimeSteps) { /* time to write averages to file */
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing time average at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);                      
		  for (itr=0; itr<numTracers; itr++) {
			ierr = VecScale(vavg[itr],1.0/avgCount);CHKERRQ(ierr);
            ierr = VecView(vavg[itr],fdavgout[itr]);CHKERRQ(ierr);
            ierr = VecSet(vavg[itr],zero); CHKERRQ(ierr);
		  }          
          avgCount = 0;        
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
#endif

  if (doTimeAverage) {
	for (itr=0; itr<numTracers; itr++) {
	  ierr = PetscViewerDestroy(&fdavgout[itr]);CHKERRQ(ierr);
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
  }

  if (rescaleForcing) {
	ierr = VecDestroy(&Rfs);CHKERRQ(ierr);
	if (periodicMatrix) {
      ierr = destroyPeriodicVec(&Rfsp);CHKERRQ(ierr);
    }	
  }
  
  if (useMonitor) {
	ierr = finalizeMonitor(tc,maxSteps,numTracers);CHKERRQ(ierr);
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
		ierr = VecDestroyVecs(numForcing,&utdf[itr]);CHKERRQ(ierr);
	  }  
	}
  }
  
  if (usePrescribedBC) {
    ierr = MatDestroy(&Be);CHKERRQ(ierr);
    ierr = MatDestroy(&Bi);CHKERRQ(ierr);
	if (periodicMatrix) {
	  ierr = destroyPeriodicMat(&Bep);CHKERRQ(ierr);
	  ierr = destroyPeriodicMat(&Bip);CHKERRQ(ierr);    
	}
    ierr = VecDestroyVecs(numTracers,&bcc);CHKERRQ(ierr);
    ierr = VecDestroyVecs(numTracers,&bcf);CHKERRQ(ierr);  
	if (periodicBC) {
	  for (itr=0; itr<numTracers; itr++) {  
		ierr = destroyPeriodicVec(&bcp[itr]);CHKERRQ(ierr);
	  }
	} else if (timeDependentBC) {
	  for (itr=0; itr<numTracers; itr++) {
		ierr = VecDestroyVecs(numBC,&bctd[itr]);CHKERRQ(ierr);
	  }  
	} else if (doCalcBC) {
      ierr = finalizeCalcBC(tc,maxSteps,numTracers);CHKERRQ(ierr);	
	}
  }
  
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}
