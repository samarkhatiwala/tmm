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

static PetscClassId CALCBC_CLASSID;

typedef struct _p_CalcBCCtx *CalcBCContext;
struct _p_CalcBCCtx {
  PETSCHEADER(int);
  PetscInt cbcctxId;
  PetscInt stateId;
/* Add problem-specific variables below */  
};


#undef __FUNCT__
#define __FUNCT__ "iniCalcBC"
PetscErrorCode iniCalcBC(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, TMMState state, void *ctx)
{

  PetscErrorCode ierr;
 
  static PetscBool registered = PETSC_FALSE;
  PetscInt cbcctxId = 0;
 
  CalcBCContext bc;

  MPI_Comm comm = PETSC_COMM_WORLD;
  
  if (!registered) {
    PetscClassIdRegister("CalcBC context", &CALCBC_CLASSID);
    registered = PETSC_TRUE;
  }
  PetscHeaderCreate(bc, CALCBC_CLASSID, "CalcBC", "CalcBC context", "CalcBC", comm, 0, 0);

  cbcctxId++;
  bc->cbcctxId=cbcctxId;
  bc->stateId=state->stateId;

  PetscContainer ctxcontainer;
  PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
  PetscCall(PetscContainerSetPointer(ctxcontainer, (void*)bc));
  PetscCall(PetscObjectCompose((PetscObject)state, "calcbc ctx", (PetscObject)ctxcontainer));
  state->calcbcctxcontainer = ctxcontainer;
  PetscCall(PetscContainerDestroy(&ctxcontainer));

/* Add your code here */

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "calcCalcBC"
PetscErrorCode calcCalcBC(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, TMMState state, void *ctx)
{

  PetscErrorCode ierr;

/* Add your code here */

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "writeCalcBC"
PetscErrorCode writeCalcBC(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, TMMState state, void *ctx)
{

  PetscErrorCode ierr;

/* Add your code here */

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "finalizeCalcBC"
PetscErrorCode finalizeCalcBC(PetscScalar tc, PetscInt Iterc, TMMState state, void *ctx)
{

  PetscErrorCode ierr;

/* Add your code here */

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "reInitializeCalcBC"
PetscErrorCode reInitializeCalcBC(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, TMMState state, void *ctx)
{

  PetscErrorCode ierr;

/* Add your code here */

  return 0;
}

