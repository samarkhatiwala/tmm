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

PetscInt monitorStartTimeStep, monitorSteps, monitorWriteSteps;

#undef __FUNCT__
#define __FUNCT__ "iniMonitor"
PetscErrorCode iniMonitor(PetscScalar tc, PetscInt Iter, PetscInt numTracers, Vec *v)
{

  PetscErrorCode ierr;
  PetscBool flg;
  
  ierr = PetscOptionsGetInt(PETSC_NULL,"-monitor_start_time_step",&monitorStartTimeStep,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate (absolute) time step at which to start monitor with the -monitor_start_time_step flag");
  ierr = PetscOptionsGetInt(PETSC_NULL,"-monitor_steps",&monitorSteps,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate frequency with which solution is monitored with the -monitor_steps flag");
  ierr = PetscOptionsGetInt(PETSC_NULL,"-monitor_write_steps",&monitorWriteSteps,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate frequency at which monitor data is written out with the -monitor_write_steps flag");

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Solution will be monitored starting at time step: %d\n", monitorStartTimeStep);CHKERRQ(ierr);	
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Solution will be monitored every %d time steps\n", monitorSteps);CHKERRQ(ierr);	
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Monitor data will be written out every %d time steps\n", monitorWriteSteps);CHKERRQ(ierr);	

/* Add your code here */

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "updateMonitor"
PetscErrorCode updateMonitor(PetscScalar tc, PetscInt iLoop, PetscInt numTracers, Vec *v)

{

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  PetscErrorCode ierr;

  if (Iter0+iLoop>=monitorStartTimeStep) { /* start monitoring solution (note: monitorStartTimeStep is ABSOLUTE time step) */	
/* Add your code here */
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "writeMonitor"
PetscErrorCode writeMonitor(PetscScalar tc, PetscInt iLoop, PetscInt numTracers, Vec *v)
{

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  PetscErrorCode ierr;

  if (Iter0+iLoop>=monitorStartTimeStep) { /* note: monitorStartTimeStep is ABSOLUTE time step */	
/* Add your code here */
  }

  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "finalizeMonitor"
PetscErrorCode finalizeMonitor(PetscScalar tc, PetscInt Iter, PetscInt numTracers)
{

  PetscErrorCode ierr;

/* Add your code here */

  return 0;
}
