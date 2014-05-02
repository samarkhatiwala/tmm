#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "petsc_matvec_utils.h"
#include "forcing_utils.h"

#undef __FUNCT__
#define __FUNCT__ "iniCalcBC"
PetscErrorCode iniCalcBC(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt numTracers, Vec *v, Vec *bcc, Vec *bcf)
{

  PetscErrorCode ierr;
  PetscInt itr;
  PetscScalar one = 1.0;

/* Add your code here */
  for (itr=0; itr<numTracers; itr++) {
    ierr = VecSet(bcc[itr],one); CHKERRQ(ierr);
    ierr = VecSet(bcf[itr],one); CHKERRQ(ierr);    
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "calcBC"
PetscErrorCode calcBC(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt numTracers, Vec *v, Vec *bcc, Vec *bcf)
{

  PetscErrorCode ierr;
  PetscInt itr;
  PetscScalar one = 1.0;

/* Add your code here */
  for (itr=0; itr<numTracers; itr++) {
    ierr = VecSet(bcc[itr],one); CHKERRQ(ierr);
    ierr = VecSet(bcf[itr],one); CHKERRQ(ierr);    
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "finalizeCalcBC"
PetscErrorCode finalizeCalcBC(PetscScalar tc, PetscInt Iter, PetscInt numTracers)
{

  PetscErrorCode ierr;

/* Add your code here */

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "reInitializeCalcBC"
PetscErrorCode reInitializeCalcBC(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt numTracers, Vec *v, Vec *bcc, Vec *bcf)
{

  PetscErrorCode ierr;
  PetscInt itr;
  PetscScalar one = 1.0;

/* Add your code here */
  for (itr=0; itr<numTracers; itr++) {
    ierr = VecSet(bcc[itr],one); CHKERRQ(ierr);
    ierr = VecSet(bcf[itr],one); CHKERRQ(ierr);    
  }

  return 0;
}

