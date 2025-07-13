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
#define __FUNCT__ "iniMisfit"
PetscErrorCode iniMisfit(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx)
{

  PetscErrorCode ierr;
  PetscBool flg;
  const char *prefix;
  
  prefix = ((PetscObject)state)->prefix;

  ierr = StepTimerCreate(&state->misfitTimer);CHKERRQ(ierr);
  ierr = StepTimerIni("misfit_", prefix, Iter0, state->misfitTimer);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Misfit will be computed starting at (and including) time step: %d\n", state->misfitTimer->startTimeStep);CHKERRQ(ierr);	
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Misfit will be computed over %d time steps\n", state->misfitTimer->numTimeSteps);CHKERRQ(ierr);

/* Add your code here */

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "calcMisfit"
PetscErrorCode calcMisfit(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx)

{

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  PetscErrorCode ierr;

  if (Iter0+iLoop>=state->misfitTimer->startTimeStep) { /* start computing misfit (note: startTimeStep is ABSOLUTE time step) */	
/* Add your code here */
//     ierr = PetscPrintf(PETSC_COMM_WORLD,"Calculating misfit at time step: %d\n", Iter0+iLoop);CHKERRQ(ierr);	
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "writeMisfit"
PetscErrorCode writeMisfit(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx)
{

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  PetscErrorCode ierr;

  if (Iter0+iLoop>=state->misfitTimer->startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */	
/* Add your code here */
  }

  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "finalizeMisfit"
PetscErrorCode finalizeMisfit(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx)
{

  PetscErrorCode ierr;

/* Add your code here */

  return 0;
}
