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

static PetscClassId EXTERNALFORCING_CLASSID;

typedef struct _p_ExternalForcingCtx *ExternalForcingContext;
struct _p_ExternalForcingCtx {
  PETSCHEADER(int);
  PetscInt efctxId;
  PetscInt stateId;
/* Add problem-specific variables below */  
};

#undef __FUNCT__
#define __FUNCT__ "iniExternalForcing"
PetscErrorCode iniExternalForcing(PetscScalar tc, PetscInt Iter, TMMState state, void *userctx)
{

  PetscErrorCode ierr;
  PetscInt numTracers;
  const char *prefix;

  static PetscBool registered = PETSC_FALSE;
  static PetscInt efctxId = 0;

  ExternalForcingContext ef;

  MPI_Comm comm = PETSC_COMM_WORLD;
  
  if (!registered) {
    PetscClassIdRegister("ExternalForcing context", &EXTERNALFORCING_CLASSID);
    registered = PETSC_TRUE;
  }
  PetscHeaderCreate(ef, EXTERNALFORCING_CLASSID, "ExternalForcing", "ExternalForcing context", "ExternalForcing", comm, 0, 0); 

  efctxId++;
  ef->efctxId=efctxId;
  ef->stateId=state->stateId;

  PetscContainer ctxcontainer;
  PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
  PetscCall(PetscContainerSetPointer(ctxcontainer, (void*)ef));
  PetscCall(PetscObjectCompose((PetscObject)state, "external forcing ctx", (PetscObject)ctxcontainer));
  state->extforcctxcontainer = ctxcontainer;
  PetscCall(PetscContainerDestroy(&ctxcontainer));

  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

// Add problem specific code ...
   
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "calcExternalForcing"
PetscErrorCode calcExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *userctx)
{

  PetscErrorCode ierr;
  PetscInt numTracers;
  const char *prefix;

  void *ctx;
  PetscCall(PetscContainerGetPointer(state->extforcctxcontainer, &ctx));
  ExternalForcingContext ef = (ExternalForcingContext)ctx;
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

// Add problem specific code ...
  
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "writeExternalForcing"
PetscErrorCode writeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *userctx)
{

  PetscErrorCode ierr;
  PetscInt numTracers;
  const char *prefix;
  
  void *ctx;
  PetscCall(PetscContainerGetPointer(state->extforcctxcontainer, &ctx));
  ExternalForcingContext ef = (ExternalForcingContext)ctx;
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

// Add problem specific code ...

  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "finalizeExternalForcing"
PetscErrorCode finalizeExternalForcing(PetscScalar tc, PetscInt Iter, TMMState state, void *userctx)
{

  PetscErrorCode ierr;
  PetscInt numTracers;
  const char *prefix;
  
  void *ctx;
  PetscCall(PetscContainerGetPointer(state->extforcctxcontainer, &ctx));
  ExternalForcingContext ef = (ExternalForcingContext)ctx;
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

// Add problem specific code ...

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "reInitializeExternalForcing"
PetscErrorCode reInitializeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *userctx)
{

  PetscErrorCode ierr;
  PetscInt numTracers;
  const char *prefix;
  
  void *ctx;
  PetscCall(PetscContainerGetPointer(state->extforcctxcontainer, &ctx));
  ExternalForcingContext ef = (ExternalForcingContext)ctx;
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

  return 0;
}

