#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "tmm_petsc_matvec_utils.h"
#include "tmm_timer.h"
#include "tmm_forcing_utils.h"
#include "tmm_profile_utils.h"
#include "tmm.h"
#include "tmm_share.h"
#include "tmm_variables.h"

extern PetscErrorCode TMMWriteExtraInitialize(PetscScalar tc, PetscInt it, TMMState state);
extern PetscErrorCode TMMWriteExtra(PetscScalar tc, PetscInt it, TMMState state);
extern PetscErrorCode TMMWriteExtraFinalize(PetscScalar tc, PetscInt it, TMMState state);

PetscErrorCode TMMOutputInitialize(TMMState state)
{

  PetscErrorCode ierr;
  PetscBool flg1;
  PetscInt itr, maxValsToRead;
  char tmpFile[PETSC_MAX_PATH_LEN];
  PetscScalar zero = 0.0;
  PetscViewer fdp;
  PetscInt numTracers;
  const char *prefix;
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;
  
// Catch deprecated option
  PetscInt dum;
  ierr = PetscOptionsGetInt(NULL,prefix,"-write_steps",&dum,&flg1);CHKERRQ(ierr);
  if (flg1) SETERRQ(PETSC_COMM_WORLD,1,"ERROR!: The -write_steps option has been deprecated. Use the StepTimer object with prefix 'write'");

  ierr = StepTimerCreate(&state->writeTimer);CHKERRQ(ierr);
  ierr = StepTimerIni("write_", prefix, Iter0, state->writeTimer);CHKERRQ(ierr);

/*Data for time averaging */
  ierr = PetscOptionsHasName(NULL,prefix,"-time_avg",&state->doTimeAverage);CHKERRQ(ierr);
  if (state->doTimeAverage) {
    ierr = StepTimerCreate(&state->avgTimer);CHKERRQ(ierr);
    ierr = StepTimerIni("avg_", prefix, Iter0+1, state->avgTimer);CHKERRQ(ierr);
    for (itr=0; itr<numTracers; itr++) {
      state->avgOutFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
    }
    maxValsToRead = numTracers;
    ierr = PetscOptionsGetStringArray(NULL,prefix,"-avg_files",state->avgOutFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
    if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate file name(s) for writing time averages with the -avg_files option");
    if (maxValsToRead != numTracers) {
      SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of time average file names specified");
    }  
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Time averages will be computed starting at and including (absolute) time step: %d\n", state->avgTimer->startTimeStep);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Time averages will be computed over %d time steps\n", state->avgTimer->numTimeSteps);CHKERRQ(ierr);	
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Time averages will be written to:\n");CHKERRQ(ierr);
    for (itr=0; itr<numTracers; itr++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,state->avgOutFile[itr]);CHKERRQ(ierr);
    }

    ierr = PetscOptionsHasName(NULL,prefix,"-avg_append",&state->avgAppendOutput);CHKERRQ(ierr);
    // These are all on the same avg timer but to be safe we use different file mode flags for each
    if (state->avgAppendOutput) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Time averages will be appended\n");CHKERRQ(ierr);
      state->AVG_FILE_MODE=FILE_MODE_APPEND;
      state->QFAVG_FILE_MODE=FILE_MODE_APPEND;
      state->QEFAVG_FILE_MODE=FILE_MODE_APPEND;
      state->BCAVG_FILE_MODE=FILE_MODE_APPEND;
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Time averages will overwrite existing file(s)\n");CHKERRQ(ierr);
      state->AVG_FILE_MODE=FILE_MODE_WRITE;
      state->QFAVG_FILE_MODE=FILE_MODE_WRITE;
      state->QEFAVG_FILE_MODE=FILE_MODE_WRITE;
      state->BCAVG_FILE_MODE=FILE_MODE_WRITE;
    }

/* Output times */
    ierr = PetscOptionsGetString(NULL,prefix,"-avg_time_file",state->avgOutTimeFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
    if (!flg1) {
      strcpy(state->avgOutTimeFile,"");
      sprintf(state->avgOutTimeFile,"%s","time_average_output_time.txt");
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Time average output times will be written to %s\n",state->avgOutTimeFile);CHKERRQ(ierr);

    if (!(state->avgAppendOutput)) {
      ierr = PetscFOpen(PETSC_COMM_WORLD,state->avgOutTimeFile,"w",&state->avgfptime);CHKERRQ(ierr);  
    } else {
      ierr = PetscFOpen(PETSC_COMM_WORLD,state->avgOutTimeFile,"a",&state->avgfptime);CHKERRQ(ierr);  
    }
	
  }


/* Output file */
  for (itr=0; itr<numTracers; itr++) {
    state->outFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
  }
  maxValsToRead = numTracers;
  ierr = PetscOptionsGetStringArray(NULL,prefix,"-o",state->outFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
  if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate output file name(s) with the -o option");
  if (maxValsToRead != numTracers) {
    SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of output file names specified");
  }  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Output will be written to:\n");CHKERRQ(ierr);
  for (itr=0; itr<numTracers; itr++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,state->outFile[itr]);CHKERRQ(ierr);
  }  
  ierr = PetscOptionsHasName(NULL,prefix,"-append",&state->appendOutput);CHKERRQ(ierr);
  // These are all on the same write timer but to be safe we use different file mode flags for each
  if (state->appendOutput) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Output will be appended\n");CHKERRQ(ierr);
    state->OUTPUT_FILE_MODE=FILE_MODE_APPEND;
    state->QF_FILE_MODE=FILE_MODE_APPEND;
    state->QEF_FILE_MODE=FILE_MODE_APPEND;
    state->BC_FILE_MODE=FILE_MODE_APPEND;
    
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Output will overwrite existing file(s)\n");CHKERRQ(ierr);
    state->OUTPUT_FILE_MODE=FILE_MODE_WRITE;
    state->QF_FILE_MODE=FILE_MODE_WRITE;
    state->QEF_FILE_MODE=FILE_MODE_WRITE;
    state->BC_FILE_MODE=FILE_MODE_WRITE;
  }

/* Output times */
  ierr = PetscOptionsGetString(NULL,prefix,"-time_file",state->outTimeFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
  if (!flg1) {
    strcpy(state->outTimeFile,"");
    sprintf(state->outTimeFile,"%s","output_time.txt");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Output times will be written to %s\n",state->outTimeFile);CHKERRQ(ierr);

/* File name for final pickup */
  ierr = PetscOptionsGetString(NULL,prefix,"-pickup_out",state->pickupoutFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
  if (!flg1) {
	   strcpy(state->pickupoutFile,"");
    sprintf(state->pickupoutFile,"%s","pickup.petsc");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Final pickup will be written to %s\n",state->pickupoutFile);CHKERRQ(ierr);

  ierr = PetscOptionsHasName(NULL,prefix,"-pickup_time_steps",&state->writePickup);CHKERRQ(ierr);
  if (state->writePickup) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Intermediate pickups will be written to pickup_ITERATIONNUMBER.petsc\n");CHKERRQ(ierr);
    ierr = StepTimerCreate(&state->pickupTimer);CHKERRQ(ierr);
    ierr = StepTimerIni("pickup_", prefix, Iter0, state->pickupTimer);CHKERRQ(ierr);
  }

// #if defined (FORSPINUP) || defined (FORJACOBIAN)
//   ierr = waitForSignal(0);CHKERRQ(ierr); /* initialize */
// #endif

/* Optionally write initial conditions */
  if (!state->appendOutput) {
    ierr = PetscFOpen(PETSC_COMM_WORLD,state->outTimeFile,"w",&state->fptime);CHKERRQ(ierr);  
    if (Iter0==state->writeTimer->startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
      ierr = PetscFPrintf(PETSC_COMM_WORLD,state->fptime,"%d   %10.5f\n",Iter0,time0);CHKERRQ(ierr);     
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing output at time %10.5f, step %d\n", time0,Iter0);CHKERRQ(ierr);  
      for (itr=0; itr<numTracers; itr++) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing tracer %d to %s\n", itr,state->outFile[itr]);CHKERRQ(ierr);  
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->outFile[itr],state->OUTPUT_FILE_MODE,&state->fdout[itr]);CHKERRQ(ierr);
        ierr = VecView(state->c[itr],state->fdout[itr]);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&state->fdout[itr]);CHKERRQ(ierr);
      }
      state->OUTPUT_FILE_MODE=FILE_MODE_APPEND;
    }
  } else {
    ierr = PetscFOpen(PETSC_COMM_WORLD,state->outTimeFile,"a",&state->fptime);CHKERRQ(ierr);
    if (Iter0==state->writeTimer->startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Initial condition will NOT be written\n");CHKERRQ(ierr);  
    }  
  }

  if (state->writePickup) {
    if (Iter0==state->pickupTimer->startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing pickup at time %10.5f, step %d\n", time0, Iter0);CHKERRQ(ierr);
      strcpy(tmpFile,"");
      sprintf(tmpFile,"pickup_%d.petsc",Iter0);
      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_WRITE,&fdp);CHKERRQ(ierr);
      for (itr=0; itr<numTracers; itr++) {
        ierr = VecView(state->c[itr],fdp);CHKERRQ(ierr);
      }
      ierr = PetscViewerDestroy(&fdp);CHKERRQ(ierr);      
	   }
  }

  ierr = TMMWriteExtraInitialize(time0,Iter0,state);CHKERRQ(ierr);

  if (state->useForcingFromFile) {  
/*  Initialize data to write out qf */
    if ((state->periodicForcing) || (state->timeDependentForcing)) {
   	  for (itr=0; itr<numTracers; itr++) {
      		state->qfoutFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
	     }
      maxValsToRead = numTracers;
      ierr = PetscOptionsGetStringArray(NULL,prefix,"-oqf",state->qfoutFile,&maxValsToRead,&state->doWriteQF);CHKERRQ(ierr);
   	  if (state->doWriteQF) {
        if (maxValsToRead != numTracers) {
          SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of forcing-from-file (qf) output file names specified");
        }  
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Forcing-from-file (qf) output will be written to:\n");CHKERRQ(ierr);
        for (itr=0; itr<numTracers; itr++) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,state->qfoutFile[itr]);CHKERRQ(ierr);
        }
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Note: forcing-from-file (qf) output is discrete in time. Divide by the appropriate time step to obtain a tendency\n");CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Note: forcing-from-file (qf) output is shifted by one time step relative to the stated output time step\n");CHKERRQ(ierr);
        if (state->doTimeAverage) {
          for (itr=0; itr<numTracers; itr++) {
            state->qfavgOutFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
    		    }
          maxValsToRead = numTracers;
          ierr = PetscOptionsGetStringArray(NULL,prefix,"-qfavg_files",state->qfavgOutFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
          if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate file name(s) for writing time averages with the -qfavg_files option");
          if (maxValsToRead != numTracers) {
            SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of forcing-from-file (qf) time average file names specified");
          }  
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Forcing-from-file (qf) time averages will be written to:\n");CHKERRQ(ierr);
      		  for (itr=0; itr<numTracers; itr++) {
         			ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,state->qfavgOutFile[itr]);CHKERRQ(ierr);
		        }      
		      }
	     }
    }
  }

  if (state->useExternalForcing) {  /* external forcing present */  
/*  Initialize data to write out qef */
    for (itr=0; itr<numTracers; itr++) {
      state->qefoutFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
    }
    maxValsToRead = numTracers;
    ierr = PetscOptionsGetStringArray(NULL,prefix,"-oqef",state->qefoutFile,&maxValsToRead,&state->doWriteQEF);CHKERRQ(ierr);
    if (state->doWriteQEF) {
      if (maxValsToRead != numTracers) {
        SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of external forcing (qef) output file names specified");
      }  
      ierr = PetscPrintf(PETSC_COMM_WORLD,"External forcing (qef) output will be written to:\n");CHKERRQ(ierr);
      for (itr=0; itr<numTracers; itr++) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,state->qefoutFile[itr]);CHKERRQ(ierr);
      }
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Note: external forcing (qef) output is discrete in time. Divide by the appropriate time step to obtain a tendency\n");CHKERRQ(ierr);	  
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Note: external forcing (qef) output is shifted by one time step relative to the stated output time step\n");CHKERRQ(ierr);
      if (state->doTimeAverage) {
        for (itr=0; itr<numTracers; itr++) {
          state->qefavgOutFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
        }
        maxValsToRead = numTracers;
        ierr = PetscOptionsGetStringArray(NULL,prefix,"-qefavg_files",state->qefavgOutFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
        if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate file name(s) for writing time averages with the -qefavg_files option");
        if (maxValsToRead != numTracers) {
          SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of external forcing (qef) time average file names specified");
        }  
        ierr = PetscPrintf(PETSC_COMM_WORLD,"External forcing (qef) time averages will be written to:\n");CHKERRQ(ierr);
        for (itr=0; itr<numTracers; itr++) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,state->qefavgOutFile[itr]);CHKERRQ(ierr);
        }      
      }
    }
  }  

/* Prescribed BCs   */
  if (state->usePrescribedBC) {
/*  BC output file */
    for (itr=0; itr<numTracers; itr++) {
      state->bcoutFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
    }
    maxValsToRead = numTracers;
    ierr = PetscOptionsGetStringArray(NULL,prefix,"-obc",state->bcoutFile,&maxValsToRead,&state->doWriteBC);CHKERRQ(ierr);
    if (state->doWriteBC) {
      if (maxValsToRead != numTracers) {
        SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of BC output file names specified");
      }  
      ierr = PetscPrintf(PETSC_COMM_WORLD,"BC output will be written to:\n");CHKERRQ(ierr);
      for (itr=0; itr<numTracers; itr++) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,state->bcoutFile[itr]);CHKERRQ(ierr);
      }
      if (state->doTimeAverage) {
        for (itr=0; itr<numTracers; itr++) {
          state->bcavgOutFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
        }
        maxValsToRead = numTracers;
        ierr = PetscOptionsGetStringArray(NULL,prefix,"-bcavg_files",state->bcavgOutFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
        if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate file name(s) for writing time averages with the -bcavg_files option");
        if (maxValsToRead != numTracers) {
          SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of BC time average file names specified");
        }  
        ierr = PetscPrintf(PETSC_COMM_WORLD,"BC time averages will be written to:\n");CHKERRQ(ierr);
        for (itr=0; itr<numTracers; itr++) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itr,state->bcavgOutFile[itr]);CHKERRQ(ierr);
        }
	     }
    }
  }
    
  if (state->doTimeAverage) {  
    ierr = VecDuplicateVecs(templateVec,numTracers,&state->cavg);CHKERRQ(ierr);
    for (itr=0; itr<numTracers; itr++) {
      ierr = VecSet(state->cavg[itr],zero); CHKERRQ(ierr);
    }
    if (state->doWriteQF) {
      ierr = VecDuplicateVecs(templateVec,numTracers,&state->qfavg);CHKERRQ(ierr);  	
      for (itr=0; itr<numTracers; itr++) {
        ierr = VecSet(state->qfavg[itr],zero); CHKERRQ(ierr);
      }
	   }
    if (state->doWriteQEF) {
      ierr = VecDuplicateVecs(templateVec,numTracers,&state->qefavg);CHKERRQ(ierr);  	
      for (itr=0; itr<numTracers; itr++) {
        ierr = VecSet(state->qefavg[itr],zero); CHKERRQ(ierr);
      }
	   }
    if (state->doWriteBC) {
      ierr = VecDuplicateVecs(bcTemplateVec,numTracers,&state->cbavg);CHKERRQ(ierr);  	
      for (itr=0; itr<numTracers; itr++) {
        ierr = VecSet(state->cbavg[itr],zero); CHKERRQ(ierr);
      }
  	 }    
  }

  return 0;
}

PetscErrorCode TMMOutput(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state)
{

// Note: tc is the time at the end of the time step (when this function is usually called)

  PetscInt itr;
  PetscScalar zero = 0.0, one = 1.0;        
  PetscErrorCode ierr;
  char tmpFile[PETSC_MAX_PATH_LEN];
  PetscViewer fdp;
  PetscInt numTracers;
  const char *prefix;
  
  if (!state->doOutput) {
    return 0;
  }
    
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;
  
/*  Write out BC at first time step */
	if (state->doWriteBC) {
	  if ((iLoop == 1) && (!(state->appendOutput))) {
     ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing BC at first time step\n");CHKERRQ(ierr);
     for (itr=0; itr<numTracers; itr++) {
       ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->bcoutFile[itr],state->BC_FILE_MODE,&state->fdbcout[itr]);CHKERRQ(ierr);
       if (state->applyBC) {
         ierr = VecView(state->cbc[itr],state->fdbcout[itr]);CHKERRQ(ierr);
       } else {
         ierr = VecView(bcTemplateVec,state->fdbcout[itr]);CHKERRQ(ierr);
       }
       ierr = PetscViewerDestroy(&state->fdbcout[itr]);CHKERRQ(ierr);	
     }
     state->BC_FILE_MODE=FILE_MODE_APPEND;
	  }	
 }

/* write output */
// Note: templateVec and bcTemplateVec should have been set to zero elsewhere
	if (Iter0+iLoop>state->writeTimer->startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
/*    We do this here in case writeExternalForcing etc want to use writeTimer */	
	  if (state->writeTimer->count<state->writeTimer->numTimeSteps) {
	    state->writeTimer->count++;
	  }
 }
 if (state->useExternalForcing) { // xxxx should this be applyExternalForcing?
   ierr = TMMComputeExtForcFunction(tc,Iter,iLoop,state,TMM_WRI_FUNC);
 }
 if (state->doCalcBC) { // xxxx should this be applyBC?
   ierr = TMMComputeCalcBCFunction(tc,Iter,-1.,-1,iLoop,state,TMM_WRI_FUNC);
 }

	if (Iter0+iLoop>=state->writeTimer->startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
	  if ((state->writeTimer->count==state->writeTimer->numTimeSteps) || (Iter0+iLoop==state->writeTimer->startTimeStep)) { /* time to write out */
     ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
     ierr = PetscFPrintf(PETSC_COMM_WORLD,state->fptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);           
     for (itr=0; itr<numTracers; itr++) {
       ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->outFile[itr],state->OUTPUT_FILE_MODE,&state->fdout[itr]);CHKERRQ(ierr);
       ierr = VecView(state->c[itr],state->fdout[itr]);CHKERRQ(ierr);
       ierr = PetscViewerDestroy(&state->fdout[itr]);CHKERRQ(ierr);
     }
     state->OUTPUT_FILE_MODE=FILE_MODE_APPEND;
     
     if (state->doWriteQF) {
       for (itr=0; itr<numTracers; itr++) {
         ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->qfoutFile[itr],state->QF_FILE_MODE,&state->fdqfout[itr]);CHKERRQ(ierr);
         if (state->applyForcingFromFile) {
           ierr = VecView(state->qf[itr],state->fdqfout[itr]);CHKERRQ(ierr);
         } else {
           ierr = VecView(templateVec,state->fdqfout[itr]);CHKERRQ(ierr);			 
         }
         ierr = PetscViewerDestroy(&state->fdqfout[itr]);CHKERRQ(ierr);
       }
       state->QF_FILE_MODE=FILE_MODE_APPEND;
     }

     if (state->doWriteQEF) {
       for (itr=0; itr<numTracers; itr++) {
         ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->qefoutFile[itr],state->QEF_FILE_MODE,&state->fdqefout[itr]);CHKERRQ(ierr);
         if (state->applyExternalForcing) {
           ierr = VecView(state->qef[itr],state->fdqefout[itr]);CHKERRQ(ierr);
         } else {
           ierr = VecView(templateVec,state->fdqefout[itr]);CHKERRQ(ierr);			 
         }
         ierr = PetscViewerDestroy(&state->fdqefout[itr]);CHKERRQ(ierr);
       }
       state->QEF_FILE_MODE=FILE_MODE_APPEND;
     }
   
     if (state->doWriteBC) {
       for (itr=0; itr<numTracers; itr++) {
         ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->bcoutFile[itr],state->BC_FILE_MODE,&state->fdbcout[itr]);CHKERRQ(ierr);
         if (state->applyBC) {
           ierr = VecView(state->cbf[itr],state->fdbcout[itr]);CHKERRQ(ierr);
         } else {
           ierr = VecView(bcTemplateVec,state->fdbcout[itr]);CHKERRQ(ierr);
         }
         ierr = PetscViewerDestroy(&state->fdbcout[itr]);CHKERRQ(ierr);
       }
       state->BC_FILE_MODE=FILE_MODE_APPEND;
     }
   }
      
   if (state->writeTimer->count==state->writeTimer->numTimeSteps) {
     ierr = StepTimerUpdate(Iter0+iLoop, state->writeTimer);CHKERRQ(ierr);
   }
 }
    
 ierr = TMMWriteExtra(tc,iLoop,state);CHKERRQ(ierr);

 if (state->writePickup) {
   if (Iter0+iLoop>state->pickupTimer->startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
     if (state->pickupTimer->count<state->pickupTimer->numTimeSteps) {
      state->pickupTimer->count++;
     }
   }		

	  if (Iter0+iLoop>=state->pickupTimer->startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
     if ((state->pickupTimer->count==state->pickupTimer->numTimeSteps) || (Iter0+iLoop==state->pickupTimer->startTimeStep)) { /* time to write out */
       ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing pickup at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
       strcpy(tmpFile,"");
       sprintf(tmpFile,"pickup_%d.petsc",Iter0+iLoop);
       ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_WRITE,&fdp);CHKERRQ(ierr);
       for (itr=0; itr<numTracers; itr++) {
         ierr = VecView(state->c[itr],fdp);CHKERRQ(ierr);
		     }
		     ierr = PetscViewerDestroy(&fdp);CHKERRQ(ierr);      
     }

     if (state->pickupTimer->count==state->pickupTimer->numTimeSteps) {
		     ierr = StepTimerUpdate(Iter0+iLoop, state->pickupTimer);CHKERRQ(ierr);
   		}
	  }
 }

 if (state->doTimeAverage) {
   if (Iter0+iLoop>=state->avgTimer->startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */
     if (state->avgTimer->count<state->avgTimer->numTimeSteps) { /* still within same averaging block so accumulate */
       for (itr=0; itr<numTracers; itr++) {
         ierr = VecAXPY(state->cavg[itr],one,state->c[itr]);CHKERRQ(ierr);
       }   
		     if (state->doWriteQF) {
      			for (itr=0; itr<numTracers; itr++) {
		         if (state->applyForcingFromFile) {
				         ierr = VecAXPY(state->qfavg[itr],one,state->qf[itr]);CHKERRQ(ierr);
			        } else {
         				ierr = VecAXPY(state->qfavg[itr],one,templateVec);CHKERRQ(ierr);			  
      			  }
      			}
		     }
		     if (state->doWriteQEF) {
      			for (itr=0; itr<numTracers; itr++) {
			        if (state->applyExternalForcing) {
         				ierr = VecAXPY(state->qefavg[itr],one,state->qef[itr]);CHKERRQ(ierr);
      			  } else {
         				ierr = VecAXPY(state->qefavg[itr],one,templateVec);CHKERRQ(ierr);			  
      			  }
			      }
		     }
   		  if (state->doWriteBC) {
			      for (itr=0; itr<numTracers; itr++) {
      			  if (state->applyBC) {
         				ierr = VecAXPY(state->cbavg[itr],one,state->cbf[itr]);CHKERRQ(ierr);
      			  } else {
         				ierr = VecAXPY(state->cbavg[itr],one,bcTemplateVec);CHKERRQ(ierr);
      			  }
    			  }
   		  }	
		     state->avgTimer->count++;
     }
     if (state->avgTimer->count==state->avgTimer->numTimeSteps) { /* time to write averages to file */        
       ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing time average at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
   		  ierr = PetscFPrintf(PETSC_COMM_WORLD,state->avgfptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);
		     for (itr=0; itr<numTracers; itr++) {
      			ierr = VecScale(state->cavg[itr],1.0/state->avgTimer->count);CHKERRQ(ierr);
      			ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->avgOutFile[itr],state->AVG_FILE_MODE,&state->fdavgout[itr]);CHKERRQ(ierr);
         ierr = VecView(state->cavg[itr],state->fdavgout[itr]);CHKERRQ(ierr);
         ierr = PetscViewerDestroy(&state->fdavgout[itr]);CHKERRQ(ierr);
         ierr = VecSet(state->cavg[itr],zero); CHKERRQ(ierr);
		     }
		     state->AVG_FILE_MODE=FILE_MODE_APPEND;
		     
		     if (state->doWriteQF) {
      			for (itr=0; itr<numTracers; itr++) {		  
           ierr = VecScale(state->qfavg[itr],1.0/state->avgTimer->count);CHKERRQ(ierr);
           ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->qfavgOutFile[itr],state->QFAVG_FILE_MODE,&state->fdqfavgout[itr]);CHKERRQ(ierr);
           ierr = VecView(state->qfavg[itr],state->fdqfavgout[itr]);CHKERRQ(ierr);
           ierr = PetscViewerDestroy(&state->fdqfavgout[itr]);CHKERRQ(ierr);
           ierr = VecSet(state->qfavg[itr],zero); CHKERRQ(ierr);
         }
         state->QFAVG_FILE_MODE=FILE_MODE_APPEND;
		     }
		     
       if (state->doWriteQEF) {
         for (itr=0; itr<numTracers; itr++) {		  
           ierr = VecScale(state->qefavg[itr],1.0/state->avgTimer->count);CHKERRQ(ierr);
           ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->qefavgOutFile[itr],state->QEFAVG_FILE_MODE,&state->fdqefavgout[itr]);CHKERRQ(ierr);
           ierr = VecView(state->qefavg[itr],state->fdqefavgout[itr]);CHKERRQ(ierr);
           ierr = PetscViewerDestroy(&state->fdqefavgout[itr]);CHKERRQ(ierr);
           ierr = VecSet(state->qefavg[itr],zero); CHKERRQ(ierr);
         }
         state->QEFAVG_FILE_MODE=FILE_MODE_APPEND;
       }
       
       if (state->doWriteBC) {
         for (itr=0; itr<numTracers; itr++) {		  
           ierr = VecScale(state->cbavg[itr],1.0/state->avgTimer->count);CHKERRQ(ierr);
           ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->bcavgOutFile[itr],state->BCAVG_FILE_MODE,&state->fdbcavgout[itr]);CHKERRQ(ierr);
           ierr = VecView(state->cbavg[itr],state->fdbcavgout[itr]);CHKERRQ(ierr);
           ierr = PetscViewerDestroy(&state->fdbcavgout[itr]);CHKERRQ(ierr);
           ierr = VecSet(state->cbavg[itr],zero); CHKERRQ(ierr);
         }
         state->BCAVG_FILE_MODE=FILE_MODE_APPEND;
       }

		     ierr = StepTimerUpdate(Iter0+iLoop, state->avgTimer);CHKERRQ(ierr);
     }
   }
 }

	return 0;

}

#undef __FUNCT__
#define __FUNCT__ "TMMOutputFinalize"
PetscErrorCode TMMOutputFinalize(PetscScalar tc, TMMState state)
// Vec *c, Vec *qf, Vec *qef, Vec *cbc, Vec *cbf)
{

  PetscInt itr;
  PetscScalar zero = 0.0, one = 1.0;        
  PetscErrorCode ierr;
  PetscViewer fdp;
  PetscInt numTracers;

  numTracers=state->numTracers;
  
//   for (itr=0; itr<numTracers; itr++) {
//     ierr = PetscViewerDestroy(&state->fdout[itr]);CHKERRQ(ierr);
//   }
  ierr = PetscFClose(PETSC_COMM_WORLD,state->fptime);CHKERRQ(ierr);

  ierr = TMMWriteExtraFinalize(tc,maxSteps,state);CHKERRQ(ierr);

  if (state->doTimeAverage) {
    ierr = PetscFClose(PETSC_COMM_WORLD,state->avgfptime);CHKERRQ(ierr);
    ierr = VecDestroyVecs(numTracers,&state->cavg);CHKERRQ(ierr);
    if (state->doWriteQF) {
      ierr = VecDestroyVecs(numTracers,&state->qfavg);CHKERRQ(ierr);
    }
   	if (state->doWriteQEF) {
      ierr = VecDestroyVecs(numTracers,&state->qefavg);CHKERRQ(ierr);
    }
    if (state->doWriteBC) {
      ierr = VecDestroyVecs(numTracers,&state->cbavg);CHKERRQ(ierr);
    }
  }
  
/* write final pickup */  
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->pickupoutFile,FILE_MODE_WRITE,&fdp);CHKERRQ(ierr);
  for (itr=0; itr<numTracers; itr++) {
    ierr = VecView(state->c[itr],fdp);CHKERRQ(ierr);
  }
  ierr = PetscViewerDestroy(&fdp);CHKERRQ(ierr);      

  return 0;
}
