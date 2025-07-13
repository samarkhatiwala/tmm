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
// #include "tmm_profile_data.h"
#include "tmm.h"
#include "tmm_share.h"

#undef __FUNCT__
#define __FUNCT__ "TMMWriteExtraInitialize"
PetscErrorCode TMMWriteExtraInitialize(PetscScalar tc, PetscInt Iter, TMMState state)
{

  PetscInt itr;
  PetscBool flg;
  PetscErrorCode ierr;
  char outTimeFile[PETSC_MAX_PATH_LEN];  
  PetscBool appendExtraOutput = PETSC_FALSE;
  const char *prefix;
  
  prefix = ((PetscObject)state)->prefix;

  ierr = PetscOptionsHasName(NULL,prefix,"-write_extra",&state->doExtraWrite);CHKERRQ(ierr);

  if (state->doExtraWrite) {    
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Extra tracer output has been activated\n");CHKERRQ(ierr);
    ierr = StepTimerCreate(&state->extraWriteTimer);CHKERRQ(ierr);
    ierr = StepTimerIni("write_extra_", prefix, Iter0, state->extraWriteTimer);CHKERRQ(ierr);

    state->numExtraTracers=MAXNUMTRACERS;    
    ierr = PetscOptionsGetIntArray(NULL,prefix,"-write_extra_tracer_indices",state->itrExtra,&state->numExtraTracers,&flg);CHKERRQ(ierr);
	   if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate extra tracer indices with the -write_extra_tracer_indices option!");

   /* Output file */
    for (itr=0; itr<state->numExtraTracers; itr++) {
      state->extraOutFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
    }

    ierr = PetscOptionsGetStringArray(NULL,prefix,"-o_extra",state->extraOutFile,&state->numExtraTracers,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate extra output file name(s) with the -o_extra option");

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Extra output will be written to:\n");CHKERRQ(ierr);
    for (itr=0; itr<state->numExtraTracers; itr++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", state->itrExtra[itr],state->extraOutFile[itr]);CHKERRQ(ierr);
    }

    ierr = PetscOptionsHasName(NULL,prefix,"-append_extra",&appendExtraOutput);CHKERRQ(ierr);
    if (appendExtraOutput) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Extra output will be appended\n");CHKERRQ(ierr);
      state->OUTPUT_EXTRA_FILE_MODE=FILE_MODE_APPEND;
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Extra output will overwrite existing file(s)\n");CHKERRQ(ierr);
      state->OUTPUT_EXTRA_FILE_MODE=FILE_MODE_WRITE;
    }    

   /* Output times */
    ierr = PetscOptionsGetString(NULL,prefix,"-time_file_extra",outTimeFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
    if (!flg) {
      strcpy(outTimeFile,"");
      sprintf(outTimeFile,"%s","output_time_extra.txt");
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Extra output times will be written to %s\n",outTimeFile);CHKERRQ(ierr);

    if (!appendExtraOutput) {
      ierr = PetscFOpen(PETSC_COMM_WORLD,outTimeFile,"w",&state->fptimeextra);CHKERRQ(ierr);  
      if (Iter0==state->extraWriteTimer->startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing extra output at time %10.5f, step %d\n", tc, Iter0);CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,state->fptimeextra,"%d   %10.5f\n",Iter0,tc);CHKERRQ(ierr);
        for (itr=0; itr<state->numExtraTracers; itr++) {
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->extraOutFile[itr],state->OUTPUT_EXTRA_FILE_MODE,&state->fdoutextra[itr]);CHKERRQ(ierr);
          ierr = VecView(state->c[state->itrExtra[itr]],state->fdoutextra[itr]);CHKERRQ(ierr);
          ierr = PetscViewerDestroy(&state->fdoutextra[itr]);CHKERRQ(ierr);
        }
        state->OUTPUT_EXTRA_FILE_MODE=FILE_MODE_APPEND;
      }
    } else {
      ierr = PetscFOpen(PETSC_COMM_WORLD,outTimeFile,"a",&state->fptimeextra);CHKERRQ(ierr);
      if (Iter0==state->extraWriteTimer->startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Opening extra output file(s) for output. Initial condition will NOT be written\n");CHKERRQ(ierr);
      }	
    }
    
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TMMWriteExtra"
PetscErrorCode TMMWriteExtra(PetscScalar tc, PetscInt iLoop, TMMState state)
{

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  PetscInt itr;
  PetscErrorCode ierr;
  const char *prefix;
  
  prefix = ((PetscObject)state)->prefix;

  if (state->doExtraWrite) {
    if (Iter0+iLoop>state->extraWriteTimer->startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
      if (state->extraWriteTimer->count<state->extraWriteTimer->numTimeSteps) {
        state->extraWriteTimer->count++;
      }
    }

    if (Iter0+iLoop>=state->extraWriteTimer->startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
      if ((state->extraWriteTimer->count==state->extraWriteTimer->numTimeSteps) || (Iter0+iLoop==state->extraWriteTimer->startTimeStep)) { /* time to write out */
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing extra output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,state->fptimeextra,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);
        for (itr=0; itr<state->numExtraTracers; itr++) {
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->extraOutFile[itr],state->OUTPUT_EXTRA_FILE_MODE,&state->fdoutextra[itr]);CHKERRQ(ierr);
          ierr = VecView(state->c[state->itrExtra[itr]],state->fdoutextra[itr]);CHKERRQ(ierr);
          ierr = PetscViewerDestroy(&state->fdoutextra[itr]);CHKERRQ(ierr);
        }
        state->OUTPUT_EXTRA_FILE_MODE=FILE_MODE_APPEND;
      }
      
      if (state->extraWriteTimer->count==state->extraWriteTimer->numTimeSteps) {
        ierr = StepTimerUpdate(Iter0+iLoop, state->extraWriteTimer);CHKERRQ(ierr);
      }
    }
  }  

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TMMWriteExtraFinalize"
PetscErrorCode TMMWriteExtraFinalize(PetscScalar tc, PetscInt Iter, TMMState state)
{

  PetscInt itr;
  PetscErrorCode ierr;

  if (state->doExtraWrite) {
    ierr = PetscFClose(PETSC_COMM_WORLD,state->fptimeextra);CHKERRQ(ierr);
  }
  return 0;
}
