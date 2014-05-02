#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "petsc_matvec_utils.h"
#include "tmm_forcing_utils.h"

#undef __FUNCT__
#define __FUNCT__ "iniCalcBC"
PetscErrorCode iniCalcBC(PetscScalar tc, PetscInt Iter, PetscScalar tf, PetscInt Iterf, PetscInt numTracers, Vec *v, Vec *bc, Vec *bf)
{

  PetscErrorCode ierr;

/* Add your code here */

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "calcBC"
PetscErrorCode calcBC(PetscScalar tc, PetscInt Iter, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, PetscInt numTracers, Vec *v, Vec *bc, Vec *bf)
{

  PetscErrorCode ierr;

/* Add your code here */

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "writeBC"
PetscErrorCode writeBC(PetscScalar tc, PetscInt iLoop, PetscInt numTracers, Vec *v, Vec *bc, Vec *bf)

{

  PetscErrorCode ierr;

/* Add your code here */

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
PetscErrorCode reInitializeCalcBC(PetscScalar tc, PetscInt Iter, PetscScalar tf, PetscInt Iterf, PetscInt numTracers, Vec *v, Vec *bc, Vec *bf)
{

  PetscErrorCode ierr;

/* Add your code here */

  return 0;
}

