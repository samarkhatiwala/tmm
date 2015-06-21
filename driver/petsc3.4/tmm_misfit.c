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

PetscInt misfitStartTimeStep, misfitWriteSteps;

#undef __FUNCT__
#define __FUNCT__ "iniMisfit"
PetscErrorCode iniMisfit(PetscScalar tc, PetscInt Iter, PetscInt numTracers, Vec *v)
{

  PetscErrorCode ierr;
  PetscBool flg;
  
  ierr = PetscOptionsGetInt(PETSC_NULL,"-misfit_start_time_step",&misfitStartTimeStep,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate (absolute) time step at which to start computing misfit with the -misfit_start_time_step flag");
  ierr = PetscOptionsGetInt(PETSC_NULL,"-misfit_write_steps",&misfitWriteSteps,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate frequency at which misfit data is written out with the -misfit_write_steps flag");

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Misfit will be computed starting at time step: %d\n", misfitStartTimeStep);CHKERRQ(ierr);	
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Misfit will be written out every %d time steps\n", misfitWriteSteps);CHKERRQ(ierr);	

/* Add your code here */

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "calcMisfit"
PetscErrorCode calcMisfit(PetscScalar tc, PetscInt iLoop, PetscInt numTracers, Vec *v)

{

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  PetscErrorCode ierr;

  if (Iter0+iLoop>=misfitStartTimeStep) { /* start computing misfit (note: misfitStartTimeStep is ABSOLUTE time step) */	
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

  if (Iter0+iLoop>=misfitStartTimeStep) { /* note: misfitStartTimeStep is ABSOLUTE time step */	
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
