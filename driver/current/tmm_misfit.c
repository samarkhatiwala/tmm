#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "petsc_matvec_utils.h"
#include "tmm_main.h"
#include "tmm_forcing_utils.h"
#include "tmm_profile_utils.h"
#include "tmm_profile_data.h"
#include "tmm_timer.h"

StepTimer misfitTimer;

#undef __FUNCT__
#define __FUNCT__ "iniMisfit"
PetscErrorCode iniMisfit(PetscScalar tc, PetscInt Iter, PetscInt numTracers, Vec *v)
{

  PetscErrorCode ierr;
  PetscBool flg;

  ierr = iniStepTimer("misfit_", Iter0, &misfitTimer);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Misfit will be computed starting at (and including) time step: %d\n", misfitTimer.startTimeStep);CHKERRQ(ierr);	
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Misfit will be computed over %d time steps\n", misfitTimer.numTimeSteps);CHKERRQ(ierr);

/* Add your code here */

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "calcMisfit"
PetscErrorCode calcMisfit(PetscScalar tc, PetscInt iLoop, PetscInt numTracers, Vec *v)

{

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  PetscErrorCode ierr;

  if (Iter0+iLoop>=misfitTimer.startTimeStep) { /* start computing misfit (note: startTimeStep is ABSOLUTE time step) */	
/* Add your code here */
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "writeMisfit"
PetscErrorCode writeMisfit(PetscScalar tc, PetscInt iLoop, PetscInt numTracers, Vec *v)
{

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  PetscErrorCode ierr;

  if (Iter0+iLoop>=misfitTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */	
/* Add your code here */
  }

  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "finalizeMisfit"
PetscErrorCode finalizeMisfit(PetscScalar tc, PetscInt Iter, PetscInt numTracers)
{

  PetscErrorCode ierr;

/* Add your code here */

  return 0;
}
