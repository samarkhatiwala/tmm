#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "petsctime.h"

#include "tmm_petsc_matvec_utils.h"
#include "tmm_forcing_utils.h"
#include "tmm_timer.h"
#include "tmm_profile_utils.h"
#include "tmm.h"
#include "tmm_share.h"
#include "tmm_variables.h"

extern PetscErrorCode TMMFinalizeState(PetscScalar tc, TMMState state);
extern PetscErrorCode TMMOutputFinalize(PetscScalar tc, TMMState state);
extern PetscErrorCode TMMForcingInitialize(TMMState state);
extern PetscErrorCode TMMOutputInitialize(TMMState state);

PetscInt stateId = 0;

#undef __FUNCT__
#define __FUNCT__ "TMMDestroy"
PetscErrorCode TMMDestroy(PetscScalar tc, TMMState *state)
{
  PetscInt itr;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  if (!*state) PetscFunctionReturn(PETSC_SUCCESS);
  PetscValidHeaderSpecific((*state), TMM_CLASSID, 1);
  if (--((PetscObject)(*state))->refct > 0) {
    *state = NULL;
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  if ((*state)->doOutput) {
    ierr = TMMOutputFinalize(tc,*state);
  }

  ierr = TMMFinalizeState(tc,*state);

  ierr = VecDestroyVecs((*state)->numTracers,&(*state)->c);CHKERRQ(ierr);  
  ierr = VecDestroyVecs((*state)->numTracers,&(*state)->cwork);CHKERRQ(ierr);
    
  if ((*state)->useExternalForcing) {
    ierr = VecDestroyVecs((*state)->numTracers,&(*state)->qef);CHKERRQ(ierr);  
  }  
  
  if ((*state)->useForcingFromFile) {
    ierr = VecDestroyVecs((*state)->numTracers,&(*state)->qf);CHKERRQ(ierr);  
	if ((*state)->periodicForcing) {
	  for (itr=0; itr<(*state)->numTracers; itr++) {  
		ierr = PeriodicVecDestroy(&(*state)->qp[itr]);CHKERRQ(ierr);
	  }
	} else if ((*state)->timeDependentForcing) {
	  for (itr=0; itr<(*state)->numTracers; itr++) {
	    ierr = TimeDependentVecDestroy(&(*state)->qtdf[itr]);
	  }  
	}
  }

  if ((*state)->relaxTracer) {
    ierr = VecDestroyVecs((*state)->numTracers,&(*state)->qrel);CHKERRQ(ierr);  
  }
  
  if ((*state)->usePrescribedBC) {
    ierr = VecDestroyVecs((*state)->numTracers,&(*state)->cbc);CHKERRQ(ierr);
    ierr = VecDestroyVecs((*state)->numTracers,&(*state)->cbf);CHKERRQ(ierr);  
	if ((*state)->periodicBC) {
	  for (itr=0; itr<(*state)->numTracers; itr++) {  
		ierr = PeriodicVecDestroy(&(*state)->cbp[itr]);CHKERRQ(ierr);
	  }
	} else if ((*state)->timeDependentBC) {
	  for (itr=0; itr<(*state)->numTracers; itr++) {
	    ierr = TimeDependentVecDestroy(&(*state)->cbtd[itr]);
	  }  
	} else if ((*state)->doCalcBC) {
//   Nothing to do	
	}
  }
  PetscCall(PetscHeaderDestroy(state));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "TMMFinalizeState"
PetscErrorCode TMMFinalizeState(PetscScalar tc, TMMState state)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  if (state->useMonitor) {
    ierr = TMMComputeMonitorFunction(tc,maxSteps,-1,state,TMM_FIN_FUNC);
// 	ierr = finalizeMonitor(tc,maxSteps,state);CHKERRQ(ierr);
  }

  if (state->doMisfit) {
    ierr = TMMComputeMisfitFunction(tc,maxSteps,-1,state,TMM_FIN_FUNC);
// 	ierr = finalizeMisfit(tc,maxSteps,state);CHKERRQ(ierr);
  }
  
  if (state->useExternalForcing) {
    ierr = TMMComputeExtForcFunction(tc,maxSteps,-1,state,TMM_FIN_FUNC);
//     ierr = finalizeExternalForcing(tc,maxSteps,state,NULL);CHKERRQ(ierr);
  }  
  
  if (state->usePrescribedBC) {
	if (state->doCalcBC) {
	  ierr = TMMComputeCalcBCFunction(tc, maxSteps, -1.0, -1, -1, state, TMM_FIN_FUNC);
//       ierr = finalizeCalcBC(tc,maxSteps,state,NULL);
	}
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "TMMCreate"
PetscErrorCode TMMCreate(TMMState *state)
{
  static PetscBool registered = PETSC_FALSE;
  TMMState s;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  *state=NULL;

  if (!registered) {
    PetscCall(PetscClassIdRegister("TMM state", &TMM_CLASSID));  
    registered = PETSC_TRUE;
  }
  PetscCall(PetscHeaderCreate(s, TMM_CLASSID, "TMMState", "TMM state", "TMMState", PETSC_COMM_WORLD, 0, 0));

  stateId++;
  if (stateId>MAXNUMSTATES) {
    SETERRQ(PETSC_COMM_WORLD,1,"ERROR!: Maximum number of states exceeded! Increase MAXNUMSTATES in tmmimpl.h");
  }
    
  s->stateId=stateId;

  s->c=NULL;
  s->qf=NULL;
  s->qef=NULL;
  s->qrel=NULL;
  s->cbc=NULL;
  s->cbf=NULL;

  s->isInitialized = PETSC_FALSE;
  
// State defaults
  s->numTracers=1;

/* run time options */
  s->useExternalForcing = PETSC_FALSE;
  s->useForcingFromFile = PETSC_FALSE;
  s->usePrescribedBC = PETSC_FALSE;
  s->applyExternalForcing = PETSC_FALSE;
  s->applyForcingFromFile = PETSC_FALSE;
  s->applyBC = PETSC_FALSE;
  s->periodicForcing = PETSC_FALSE;
  s->timeDependentForcing = PETSC_FALSE;
  s->constantForcing = PETSC_FALSE;
  s->periodicBC = PETSC_FALSE;
  s->timeDependentBC = PETSC_FALSE;
  s->constantBC = PETSC_FALSE;
  s->doCalcBC = PETSC_FALSE;
  s->useMonitor = PETSC_FALSE;
  s->doMisfit = PETSC_FALSE;
  s->relaxTracer = PETSC_FALSE;
  s->isInitializedExternalForcing = PETSC_FALSE;
  s->isInitializedCalcBC = PETSC_FALSE;
  s->isInitializedMonitor = PETSC_FALSE;
  s->isInitializedMisfit = PETSC_FALSE;
  
/* Forcing */
  s->forcingFromFileStartStep = -1;
  s->externalForcingStartStep = -1;
  s->forcingFromFileCutOffStep = -1;
  s->externalForcingCutOffStep = -1;

/* BC's */
  s->bcStartStep = -1;
  s->bcCutOffStep = -1;

/* I/O   */
  s->appendOutput = PETSC_FALSE;
  s->writePickup = PETSC_FALSE;
  s->doWriteBC = PETSC_FALSE;
  s->doWriteQF = PETSC_FALSE;
  s->doWriteQEF = PETSC_FALSE;
  s->pickupFromFile = PETSC_FALSE;
  s->doTimeAverage = PETSC_FALSE;
  s->avgAppendOutput = PETSC_FALSE;
  s->doExtraWrite = PETSC_FALSE;

  s->cbavg=NULL;
  s->qfavg=NULL;
  s->qefavg=NULL;
  s->cavg=NULL;  
  s->cwork=NULL;

/* Default external forcing functions */
//   ierr = TMMSetExternalForcingFunctions(s,iniExternalForcing,calcExternalForcing,writeExternalForcing,
//                                         finalizeExternalForcing,reInitializeExternalForcing,NULL);

/* Default calc BC functions */
//   ierr = TMMSetCalcBCFunctions(s,iniCalcBC,calcCalcBC,writeCalcBC,
//                                         finalizeCalcBC,reInitializeCalcBC,NULL);
  
  *state=s;
  PetscFunctionReturn(PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "TMMSetFromOptions"
PetscErrorCode TMMSetFromOptions(TMMState state, const char pre[], PetscBool doOutput)
{

  PetscErrorCode ierr;
  PetscBool flg1;
  PetscInt itr, maxValsToRead;
  PetscScalar zero = 0.0;
  PetscViewer fd;
  const char *prefix;
  
//----------------------------------------------------------------------

  ierr = TMMSetOptionsPrefix(state, pre);

  state->doOutput = doOutput;

  prefix = ((PetscObject)state)->prefix;

/* tracer vectors */
/* Number of tracers */
  ierr = PetscOptionsGetInt(NULL,prefix,"-numtracers",&state->numTracers,&flg1);CHKERRQ(ierr);
  if (state->numTracers>MAXNUMTRACERS) {
   SETERRQ(PETSC_COMM_WORLD,1,"Number of tracers exceeds maximum allowable. Please increase the variable MAXNUMTRACERS and recompile");
  }  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of tracers to be integrated: %d\n", state->numTracers);CHKERRQ(ierr); 

  ierr = VecDuplicateVecs(templateVec,state->numTracers,&state->c);CHKERRQ(ierr);

/* Initial condition     */
  for (itr=0; itr<state->numTracers; itr++) {
    state->iniFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
  }
  ierr = PetscOptionsGetString(NULL,prefix,"-pickup",state->pickupFile,PETSC_MAX_PATH_LEN-1,&state->pickupFromFile);CHKERRQ(ierr);
  if (state->pickupFromFile) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Pickup file has been specified\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  Reading initial conditions from %s\n", state->pickupFile);CHKERRQ(ierr);
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->pickupFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
    for (itr=0; itr<state->numTracers; itr++) {
      ierr = VecLoad(state->c[itr],fd);CHKERRQ(ierr); /* IntoVector */
    }
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);          
  } else {
    maxValsToRead = state->numTracers;
    ierr = PetscOptionsGetStringArray(NULL,prefix,"-i",state->iniFile,&maxValsToRead,&flg1);CHKERRQ(ierr);
    if (flg1) {  /* read from file */
      if (maxValsToRead != state->numTracers) {
        SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of input file names specified");
      }      
// #if defined (FORSPINUP) || defined (FORJACOBIAN)
// 	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Waiting for new initial condition ...\n");CHKERRQ(ierr);
//       ierr = waitForSignal(10);CHKERRQ(ierr);
// #endif	  
      for (itr=0; itr<state->numTracers; itr++) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading initial condition from file %s\n", itr,state->iniFile[itr]);CHKERRQ(ierr);
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->iniFile[itr],FILE_MODE_READ,&fd);CHKERRQ(ierr);
        ierr = VecLoad(state->c[itr],fd);CHKERRQ(ierr);/* IntoVector */
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);          
      }  
    } else {  /* set to zero */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting initial condition to zero\n");CHKERRQ(ierr);
      for (itr=0; itr<state->numTracers; itr++) {
        VecSet(state->c[itr],zero);
      }        
    }
  }

//SPK Not sure where to put these
/* initialize monitor */
  ierr = PetscOptionsHasName(NULL,prefix,"-use_monitor",&state->useMonitor);CHKERRQ(ierr);
// This is now moved to TMMForcingUpdate because otherwise the functions need to be specified before calling TMMSetFromOptions
//   if (state->useMonitor) {  
//     ierr = iniMonitor(time0,Iter0,state);CHKERRQ(ierr);
//   }

/* initialize misfit */
  ierr = PetscOptionsHasName(NULL,prefix,"-calc_misfit",&state->doMisfit);CHKERRQ(ierr);
// This is now moved to TMMForcingUpdate because otherwise the functions need to be specified before calling TMMSetFromOptions
//   if (state->doMisfit) {  
//     ierr = iniMisfit(time0,Iter0,state);CHKERRQ(ierr);
//   }

  ierr = VecDuplicateVecs(templateVec,state->numTracers,&state->cwork);CHKERRQ(ierr); /* temporary work space needed by forwardStep */

  ierr = TMMForcingInitialize(state);
  
  if (state->doOutput) {
    ierr = TMMOutputInitialize(state);
  }
  
  state->isInitialized = PETSC_TRUE;
  
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "TMMSetOptionsPrefix"
PetscErrorCode TMMSetOptionsPrefix(TMMState state, const char prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  PetscCall(PetscObjectSetOptionsPrefix((PetscObject)state, prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#undef __FUNCT__
#define __FUNCT__ "TMMGetOptionsPrefix"
PetscErrorCode TMMGetOptionsPrefix(TMMState state, const char *prefix[])
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  PetscCall(PetscObjectGetOptionsPrefix((PetscObject)state, prefix));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// PetscErrorCode TMMSetExternalForcingFunctions(TMMState state, 
// TMMExtIniFunctionFn *fini,
// TMMExtCalcFunctionFn *fcalc,
// TMMExtWriFunctionFn *fwri,
// TMMExtFinFunctionFn *ffin,
// TMMExtReiFunctionFn *frei
// )
// {
//   PetscFunctionBegin;
//   PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
//   if (fini) {
// 	state->ops->iniexternalforcing = fini;
// // 	if (ctx) {
// // 	  PetscContainer ctxcontainer;
// // 	  PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
// // 	  PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
// // 	  PetscCall(PetscObjectCompose((PetscObject)state, "ini function ctx", (PetscObject)ctxcontainer));
// // 	  state->iniexternalforcingctxcontainer = ctxcontainer;
// // 	  PetscCall(PetscContainerDestroy(&ctxcontainer));
// // 	}
//   }
//   if (fcalc) {
//     state->ops->calcexternalforcing = fcalc;
// // 	if (ctx) {
// // 	  PetscContainer ctxcontainer;
// // 	  PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
// // 	  PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
// // 	  PetscCall(PetscObjectCompose((PetscObject)state, "calc function ctx", (PetscObject)ctxcontainer));
// // 	  state->calcexternalforcingctxcontainer = ctxcontainer;
// // 	  PetscCall(PetscContainerDestroy(&ctxcontainer));
// // 	}
//   }
//   if (fwri) {
// 	state->ops->wriexternalforcing = fwri;
// // 	if (ctx) {
// // 	  PetscContainer ctxcontainer;
// // 	  PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
// // 	  PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
// // 	  PetscCall(PetscObjectCompose((PetscObject)state, "write function ctx", (PetscObject)ctxcontainer));
// // 	  state->wriexternalforcingctxcontainer = ctxcontainer;
// // 	  PetscCall(PetscContainerDestroy(&ctxcontainer));
// // 	}
//   }
//   if (ffin) {
// 	state->ops->finexternalforcing = ffin;
// // 	if (ctx) {
// // 	  PetscContainer ctxcontainer;
// // 	  PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
// // 	  PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
// // 	  PetscCall(PetscObjectCompose((PetscObject)state, "finalize function ctx", (PetscObject)ctxcontainer));
// // 	  state->finexternalforcingctxcontainer = ctxcontainer;
// // 	  PetscCall(PetscContainerDestroy(&ctxcontainer));
// // 	}
//   }
//   if (frei) {
// 	state->ops->reiexternalforcing = frei;
// // 	if (ctx) {
// // 	  PetscContainer ctxcontainer;
// // 	  PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
// // 	  PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
// // 	  PetscCall(PetscObjectCompose((PetscObject)state, "reinitialize function ctx", (PetscObject)ctxcontainer));
// // 	  state->reiexternalforcingctxcontainer = ctxcontainer;
// // 	  PetscCall(PetscContainerDestroy(&ctxcontainer));
// // 	}
//   }
//   PetscFunctionReturn(PETSC_SUCCESS);
// }

PetscErrorCode TMMSetIniExtForcFunction(TMMState state, TMMExtIniFunctionFn *fini, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (fini) {
//     PetscPrintf(PETSC_COMM_WORLD,"Setting ini ext forcing func\n");
    state->ops->iniexternalforcing = fini;
    if (ctx) {
//       PetscPrintf(PETSC_COMM_WORLD,"Setting ini ext forcing ctx\n");
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "ExtForcIniCtx", (PetscObject)ctxcontainer));
      state->iniextforcuserctxcontainer = ctxcontainer;
//       PetscPrintf(PETSC_COMM_WORLD,"Ini ext forcing ctx=%p\n",state->iniextforcuserctxcontainer);
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMSetCalcExtForcFunction(TMMState state, TMMExtCalcFunctionFn *fcalc, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (fcalc) {
    state->ops->calcexternalforcing = fcalc;
    if (ctx) {
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "ExtForCalcCtx", (PetscObject)ctxcontainer));
      state->calcextforcuserctxcontainer = ctxcontainer;
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMSetWriExtForcFunction(TMMState state, TMMExtWriFunctionFn *fwri, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (fwri) {
    state->ops->wriexternalforcing = fwri;
    if (ctx) {
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "ExtForWriCtx", (PetscObject)ctxcontainer));
      state->wriextforcuserctxcontainer = ctxcontainer;
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMSetFinExtForcFunction(TMMState state, TMMExtFinFunctionFn *ffin, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (ffin) {
    state->ops->finexternalforcing = ffin;
    if (ctx) {
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "ExtForFinCtx", (PetscObject)ctxcontainer));
      state->finextforcuserctxcontainer = ctxcontainer;
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMSetReiExtForcFunction(TMMState state, TMMExtReiFunctionFn *frei, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (frei) {
    state->ops->reiexternalforcing = frei;
    if (ctx) {
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "ExtForReiCtx", (PetscObject)ctxcontainer));
      state->reiextforcuserctxcontainer = ctxcontainer;
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMSetIniCalcBCFunction(TMMState state, TMMCalcBCIniFunctionFn *fini, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (fini) {
    state->ops->inicalcbc = fini;
    if (ctx) {
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "CalcBCIniCtx", (PetscObject)ctxcontainer));
      state->inicalcbcuserctxcontainer = ctxcontainer;
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMSetCalcCalcBCFunction(TMMState state, TMMCalcBCCalcFunctionFn *fcalc, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (fcalc) {
    state->ops->calccalcbc = fcalc;
    if (ctx) {
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "CalcBCCalcCtx", (PetscObject)ctxcontainer));
      state->calccalcbcuserctxcontainer = ctxcontainer;
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMSetWriCalcBCFunction(TMMState state, TMMCalcBCWriFunctionFn *fwri, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (fwri) {
    state->ops->wricalcbc = fwri;
    if (ctx) {
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "CalcBCWriCtx", (PetscObject)ctxcontainer));
      state->wricalcbcuserctxcontainer = ctxcontainer;
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMSetFinCalcBCFunction(TMMState state, TMMCalcBCFinFunctionFn *ffin, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (ffin) {
    state->ops->fincalcbc = ffin;
    if (ctx) {
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "CalcBCFinCtx", (PetscObject)ctxcontainer));
      state->fincalcbcuserctxcontainer = ctxcontainer;
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMSetReiCalcBCFunction(TMMState state, TMMCalcBCReiFunctionFn *frei, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (frei) {
    state->ops->reicalcbc = frei;
    if (ctx) {
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "CalcBCReiCtx", (PetscObject)ctxcontainer));
      state->reicalcbcuserctxcontainer = ctxcontainer;
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMSetIniMonitorFunction(TMMState state, TMMMonitorIniFunctionFn *fini, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (fini) {
    state->ops->inimonitor = fini;
    if (ctx) {
//       PetscPrintf(PETSC_COMM_WORLD,"Setting ini monitor ctx\n");
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "MonitorIniCtx", (PetscObject)ctxcontainer));
      state->inimonitoruserctxcontainer = ctxcontainer;
//       PetscPrintf(PETSC_COMM_WORLD,"Ini monitor ctx=%p\n",state->inimonitoruserctxcontainer);
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMSetCalcMonitorFunction(TMMState state, TMMMonitorCalcFunctionFn *fcalc, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (fcalc) {
    state->ops->calcmonitor = fcalc;
    if (ctx) {
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "MonitorCalcCtx", (PetscObject)ctxcontainer));
      state->calcmonitoruserctxcontainer = ctxcontainer;
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMSetWriMonitorFunction(TMMState state, TMMMonitorWriFunctionFn *fwri, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (fwri) {
    state->ops->wrimonitor = fwri;
    if (ctx) {
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "MonitorWriCtx", (PetscObject)ctxcontainer));
      state->wrimonitoruserctxcontainer = ctxcontainer;
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMSetFinMonitorFunction(TMMState state, TMMMonitorFinFunctionFn *ffin, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (ffin) {
    state->ops->finmonitor = ffin;
    if (ctx) {
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "MonitorFinCtx", (PetscObject)ctxcontainer));
      state->finmonitoruserctxcontainer = ctxcontainer;
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMSetIniMisfitFunction(TMMState state, TMMMisfitIniFunctionFn *fini, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (fini) {
    state->ops->inimisfit = fini;
    if (ctx) {
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "MisfitIniCtx", (PetscObject)ctxcontainer));
      state->inimisfituserctxcontainer = ctxcontainer;
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMSetCalcMisfitFunction(TMMState state, TMMMisfitCalcFunctionFn *fcalc, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (fcalc) {
    state->ops->calcmisfit = fcalc;
    if (ctx) {
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "MisfitCalcCtx", (PetscObject)ctxcontainer));
      state->calcmisfituserctxcontainer = ctxcontainer;
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMSetWriMisfitFunction(TMMState state, TMMMisfitWriFunctionFn *fwri, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (fwri) {
    state->ops->wrimisfit = fwri;
    if (ctx) {
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "MisfitWriCtx", (PetscObject)ctxcontainer));
      state->wrimisfituserctxcontainer = ctxcontainer;
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMSetFinMisfitFunction(TMMState state, TMMMisfitFinFunctionFn *ffin, void *ctx)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
  if (ffin) {
    state->ops->finmisfit = ffin;
    if (ctx) {
      PetscContainer ctxcontainer;
      PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
      PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
      PetscCall(PetscObjectCompose((PetscObject)state, "MisfitFinCtx", (PetscObject)ctxcontainer));
      state->finmisfituserctxcontainer = ctxcontainer;
      PetscCall(PetscContainerDestroy(&ctxcontainer));
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMComputeExtForcFunction(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, TMMFUNCTYPE theFunc)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);

  if (theFunc==TMM_INI_FUNC) {
    PetscCheck(state->ops->iniexternalforcing, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetIniExtForcFunction() before TMMComputeExtForcFunction().");
    if (state->ops->iniexternalforcing) {
      {
        void *ctx;
        TMMExtIniFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
        f = state->ops->iniexternalforcing;
        if (state->iniextforcuserctxcontainer) {
//           PetscPrintf(PETSC_COMM_WORLD,"Fetching iniextforcuserctxcontainer\n");
//           PetscPrintf(PETSC_COMM_WORLD,"ext forcing ctx=%p\n",state->iniextforcuserctxcontainer);
          PetscCall(PetscContainerGetPointer(state->iniextforcuserctxcontainer, &ctx));
//           PetscPrintf(PETSC_COMM_WORLD,"Finished fetching iniextforcuserctxcontainer\n");
        } else {
          ctx = NULL;
        }
        PetscCallBack("TMM callback function", (*f)(tc, Iter, state, ctx));
      }
    }
  } else if (theFunc==TMM_CALC_FUNC) {
    PetscCheck(state->ops->calcexternalforcing, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetCalcExtForcFunction() before TMMComputeExtForcFunction().");
    if (state->ops->calcexternalforcing) {
      {
        void *ctx;
        TMMExtCalcFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *ctx);
        f = state->ops->calcexternalforcing;
        if (state->calcextforcuserctxcontainer) {
          PetscCall(PetscContainerGetPointer(state->calcextforcuserctxcontainer, &ctx));
        } else {
          ctx = NULL;
        }
        PetscCallBack("TMM callback function", (*f)(tc, Iter, iLoop, state, ctx));
      }
    }
  } else if (theFunc==TMM_WRI_FUNC) {
// 	PetscCheck(state->ops->wriexternalforcing, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetWriExtForcFunction() before TMMComputeExtForcFunction().");
    if (state->ops->wriexternalforcing) {
      {
        void *ctx;
        TMMExtWriFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *ctx);
        f = state->ops->wriexternalforcing;
        if (state->wriextforcuserctxcontainer) {
          PetscCall(PetscContainerGetPointer(state->wriextforcuserctxcontainer, &ctx));
        } else {
          ctx = NULL;
        }		
        PetscCallBack("TMM callback function", (*f)(tc, Iter, iLoop, state, ctx));
      }
    }
  } else if (theFunc==TMM_FIN_FUNC) {
    PetscCheck(state->ops->finexternalforcing, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetFinExtForcFunction() before TMMComputeExtForcFunction().");
    if (state->ops->finexternalforcing) {
      {
        void *ctx;
        TMMExtFinFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
        f = state->ops->finexternalforcing;
        if (state->finextforcuserctxcontainer) {
          PetscCall(PetscContainerGetPointer(state->finextforcuserctxcontainer, &ctx));
        } else {
          ctx = NULL;
        }
        PetscCallBack("TMM callback function", (*f)(tc, Iter, state, ctx));
      }
    }
  } else if (theFunc==TMM_REI_FUNC) {
// 	PetscCheck(state->ops->reiexternalforcing, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetReiExtForcFunction() before TMMComputeExtForcFunction().");
    if (state->ops->reiexternalforcing) {
      {
        void *ctx;
        TMMExtReiFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *ctx);
        f = state->ops->reiexternalforcing;
        if (state->reiextforcuserctxcontainer) {
          PetscCall(PetscContainerGetPointer(state->reiextforcuserctxcontainer, &ctx));
        } else {
          ctx = NULL;
        }
        PetscCallBack("TMM callback function", (*f)(tc, Iter, iLoop, state, ctx));
      }
    }
  } else {
    SETERRQ(PETSC_COMM_WORLD,1,"ERROR!: Unknown function type!");
  }  

  PetscFunctionReturn(PETSC_SUCCESS);
}

// PetscErrorCode TMMSetFunction(TMMState state, PetscErrorCode (*f)(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *ctx), void *ctx)
// {
//   PetscFunctionBegin;
//   PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
//   if (f) state->ops->calcexternalforcing = f;
//   if (ctx) {
//     PetscContainer ctxcontainer;
//     PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
//     PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
//     PetscCall(PetscObjectCompose((PetscObject)state, "function ctx", (PetscObject)ctxcontainer));
//     state->calcexternalforcingctxcontainer = ctxcontainer;
//     PetscCall(PetscContainerDestroy(&ctxcontainer));
//   }
//   PetscFunctionReturn(PETSC_SUCCESS);
// }

// PetscErrorCode TMMGetFunction(TMMState state, PetscErrorCode (**f)(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *ctx), void **ctx)
// {
//   PetscFunctionBegin;
//   PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
//   if (f) *f = state->ops->calcexternalforcing;
//   if (ctx) {
//     if (state->calcexternalforcingctxcontainer) PetscCall(PetscContainerGetPointer(state->calcexternalforcingctxcontainer, ctx));
//     else *ctx = NULL;
//   }
//   PetscFunctionReturn(PETSC_SUCCESS);
// }

// PetscErrorCode TMMComputeFunction(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state)
// {
//   PetscFunctionBegin;
//   PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
//   PetscCheck(state->ops->calcexternalforcing, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetFunction() before TMMComputeFunction().");
//   if (state->ops->calcexternalforcing) {
//     {
//       void *ctx;
//       PetscErrorCode (*calcexternalforcing)(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *ctx);
//       PetscCall(TMMGetFunction(state, &calcexternalforcing, &ctx));
//       PetscCallBack("TMM callback function", (*calcexternalforcing)(tc, Iter, iLoop, state, ctx));
//     }
//   }
//   PetscFunctionReturn(PETSC_SUCCESS);
// }


// PetscErrorCode TMMSetCalcBCFunctions(TMMState state, 
// TMMCalcBCIniFunctionFn *fini,
// TMMCalcBCCalcFunctionFn *fcalc,
// TMMCalcBCWriFunctionFn *fwri,
// TMMCalcBCFinFunctionFn *ffin,
// TMMCalcBCReiFunctionFn *frei,
// void *ctx)
// {
//   PetscFunctionBegin;
//   PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
//   if (fini) {
// 	state->ops->inicalcbc = fini;
// 	if (ctx) {
// 	  PetscContainer ctxcontainer;
// 	  PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
// 	  PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
// 	  PetscCall(PetscObjectCompose((PetscObject)state, "ini calcbc function ctx", (PetscObject)ctxcontainer));
// 	  state->inicalcbcctxcontainer = ctxcontainer;
// 	  PetscCall(PetscContainerDestroy(&ctxcontainer));
// 	}
//   }
//   if (fcalc) {
//     state->ops->calccalcbc = fcalc;
// 	if (ctx) {
// 	  PetscContainer ctxcontainer;
// 	  PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
// 	  PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
// 	  PetscCall(PetscObjectCompose((PetscObject)state, "calc calcbc function ctx", (PetscObject)ctxcontainer));
// 	  state->calccalcbcctxcontainer = ctxcontainer;
// 	  PetscCall(PetscContainerDestroy(&ctxcontainer));
// 	}
//   }
//   if (fwri) {
// 	state->ops->wricalcbc = fwri;
// 	if (ctx) {
// 	  PetscContainer ctxcontainer;
// 	  PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
// 	  PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
// 	  PetscCall(PetscObjectCompose((PetscObject)state, "write calcbc function ctx", (PetscObject)ctxcontainer));
// 	  state->wricalcbcctxcontainer = ctxcontainer;
// 	  PetscCall(PetscContainerDestroy(&ctxcontainer));
// 	}
//   }
//   if (ffin) {
// 	state->ops->fincalcbc = ffin;
// 	if (ctx) {
// 	  PetscContainer ctxcontainer;
// 	  PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
// 	  PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
// 	  PetscCall(PetscObjectCompose((PetscObject)state, "finalize calcbc function ctx", (PetscObject)ctxcontainer));
// 	  state->fincalcbcctxcontainer = ctxcontainer;
// 	  PetscCall(PetscContainerDestroy(&ctxcontainer));
// 	}
//   }
//   if (frei) {
// 	state->ops->reicalcbc = frei;
// 	if (ctx) {
// 	  PetscContainer ctxcontainer;
// 	  PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
// 	  PetscCall(PetscContainerSetPointer(ctxcontainer, ctx));
// 	  PetscCall(PetscObjectCompose((PetscObject)state, "reinitialize calcbc function ctx", (PetscObject)ctxcontainer));
// 	  state->reicalcbcctxcontainer = ctxcontainer;
// 	  PetscCall(PetscContainerDestroy(&ctxcontainer));
// 	}
//   }
//   PetscFunctionReturn(PETSC_SUCCESS);
// }

PetscErrorCode TMMComputeCalcBCFunction(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, TMMState state, TMMFUNCTYPE theFunc)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);

  if (theFunc==TMM_INI_FUNC) {
    PetscCheck(state->ops->inicalcbc, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetIniCalcBCFunction() before TMMComputeCalcBCFunction().");
    if (state->ops->inicalcbc) {
      {
        void *ctx;
        TMMCalcBCIniFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, TMMState state, void *ctx);
        f = state->ops->inicalcbc;
        if (state->inicalcbcuserctxcontainer) {
          PetscCall(PetscContainerGetPointer(state->inicalcbcuserctxcontainer, &ctx));
        } else {
          ctx = NULL;
        }
        PetscCallBack("TMM callback function", (*f)(tc, Iterc, tf, Iterf, state, ctx));
      }
    }
  } else if (theFunc==TMM_CALC_FUNC) {
    PetscCheck(state->ops->calccalcbc, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetCalcBCFunctions() before TMMComputeCalcBCFunction().");
    if (state->ops->calccalcbc) {
      {
        void *ctx;
        TMMCalcBCCalcFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, TMMState state, void *ctx);
        f = state->ops->calccalcbc;
        if (state->calccalcbcuserctxcontainer) {
          PetscCall(PetscContainerGetPointer(state->calccalcbcuserctxcontainer, &ctx));
        } else {
          ctx = NULL;
        }
        PetscCallBack("TMM callback function", (*f)(tc, Iterc, tf, Iterf, iLoop, state, ctx));
      }
    }
  } else if (theFunc==TMM_WRI_FUNC) {
// 	PetscCheck(state->ops->wricalcbc, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetWriCalcBCFunction() before TMMComputeCalcBCFunction().");
    if (state->ops->wricalcbc) {
      {
        void *ctx;
        TMMCalcBCWriFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, TMMState state, void *ctx);
        f = state->ops->wricalcbc;
        if (state->wricalcbcuserctxcontainer) {
          PetscCall(PetscContainerGetPointer(state->wricalcbcuserctxcontainer, &ctx));
        } else {
          ctx = NULL;
        }		
        PetscCallBack("TMM callback function", (*f)(tc, Iterc, iLoop, state, ctx));
      }
    }
  } else if (theFunc==TMM_FIN_FUNC) {
    PetscCheck(state->ops->fincalcbc, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetFinCalcBCFunction() before TMMComputeCalcBCFunction().");
    if (state->ops->fincalcbc) {
      {
        void *ctx;
        TMMCalcBCFinFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt Iterc, TMMState state, void *ctx);
        f = state->ops->fincalcbc;
        if (state->fincalcbcuserctxcontainer) {
          PetscCall(PetscContainerGetPointer(state->fincalcbcuserctxcontainer, &ctx));
        } else {
          ctx = NULL;
        }
        PetscCallBack("TMM callback function", (*f)(tc, Iterc, state, ctx));
      }
    }
  } else if (theFunc==TMM_REI_FUNC) {
// 	PetscCheck(state->ops->reicalcbc, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetReiCalcBCFunction() before TMMComputeCalcBCFunction().");
    if (state->ops->reicalcbc) {
      {
        void *ctx;
        TMMCalcBCReiFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, TMMState state, void *ctx);
        f = state->ops->reicalcbc;
        if (state->reicalcbcuserctxcontainer) {
          PetscCall(PetscContainerGetPointer(state->reicalcbcuserctxcontainer, &ctx));
        } else {
          ctx = NULL;
        }
        PetscCallBack("TMM callback function", (*f)(tc, Iterc, tf, Iterf, iLoop, state, ctx));
      }
    }
  } else {
    SETERRQ(PETSC_COMM_WORLD,1,"ERROR!: Unknown function type!");
  }  

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMComputeMonitorFunction(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, TMMFUNCTYPE theFunc)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);

  if (theFunc==TMM_INI_FUNC) {
    PetscCheck(state->ops->inimonitor, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetIniMonitorFunction() before TMMComputeMonitorFunction().");
    if (state->ops->inimonitor) {
      {
        void *ctx;
        TMMMonitorIniFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
        f = state->ops->inimonitor;
        if (state->inimonitoruserctxcontainer) {
          PetscCall(PetscContainerGetPointer(state->inimonitoruserctxcontainer, &ctx));
        } else {
          ctx = NULL;
        }
        PetscCallBack("TMM callback function", (*f)(tc, Iter, state, ctx));
      }
    }
  } else if (theFunc==TMM_CALC_FUNC) {
    PetscCheck(state->ops->calcmonitor, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetCalcMonitorFunction() before TMMComputeMonitorFunction().");
    if (state->ops->calcmonitor) {
      {
        void *ctx;
        TMMMonitorCalcFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx);
        f = state->ops->calcmonitor;
        if (state->calcmonitoruserctxcontainer) {
          PetscCall(PetscContainerGetPointer(state->calcmonitoruserctxcontainer, &ctx));
        } else {
          ctx = NULL;
        }
        PetscCallBack("TMM callback function", (*f)(tc, iLoop, state, ctx));
      }
    }
  } else if (theFunc==TMM_WRI_FUNC) {
// 	PetscCheck(state->ops->wrimonitor, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetWriMonitorFunction() before TMMComputeMonitorFunction().");
    if (state->ops->wrimonitor) {
      {
        void *ctx;
        TMMMonitorWriFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx);
        f = state->ops->wrimonitor;
        if (state->wrimonitoruserctxcontainer) {
          PetscCall(PetscContainerGetPointer(state->wrimonitoruserctxcontainer, &ctx));
        } else {
          ctx = NULL;
        }		
        PetscCallBack("TMM callback function", (*f)(tc, iLoop, state, ctx));
      }
    }
  } else if (theFunc==TMM_FIN_FUNC) {
    PetscCheck(state->ops->finmonitor, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetFinMonitorFunction() before TMMComputeMonitorFunction().");
    if (state->ops->finmonitor) {
      {
        void *ctx;
        TMMMonitorFinFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
        f = state->ops->finmonitor;
        if (state->finmonitoruserctxcontainer) {
          PetscCall(PetscContainerGetPointer(state->finmonitoruserctxcontainer, &ctx));
        } else {
          ctx = NULL;
        }
        PetscCallBack("TMM callback function", (*f)(tc, Iter, state, ctx));
      }
    }
  } else {
    SETERRQ(PETSC_COMM_WORLD,1,"ERROR!: Unknown function type!");
  }  

  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TMMComputeMisfitFunction(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, TMMFUNCTYPE theFunc)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(state, TMM_CLASSID, 1);

  if (theFunc==TMM_INI_FUNC) {
    PetscCheck(state->ops->inimisfit, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetIniMisfitFunction() before TMMComputeMisfitFunction().");
    if (state->ops->inimisfit) {
      {
        void *ctx;
        TMMMisfitIniFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
        f = state->ops->inimisfit;
        if (state->inimisfituserctxcontainer) {
          PetscCall(PetscContainerGetPointer(state->inimisfituserctxcontainer, &ctx));
        } else {
          ctx = NULL;
        }
        PetscCallBack("TMM callback function", (*f)(tc, Iter, state, ctx));
      }
    }
  } else if (theFunc==TMM_CALC_FUNC) {
    PetscCheck(state->ops->calcmisfit, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetCalcMisfitFunction() before TMMComputeMisfitFunction().");
    if (state->ops->calcmisfit) {
      {
        void *ctx;
        TMMMisfitCalcFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx);
        f = state->ops->calcmisfit;
        if (state->calcmisfituserctxcontainer) {
          PetscCall(PetscContainerGetPointer(state->calcmisfituserctxcontainer, &ctx));
        } else {
          ctx = NULL;
        }
        PetscCallBack("TMM callback function", (*f)(tc, iLoop, state, ctx));
      }
    }
  } else if (theFunc==TMM_WRI_FUNC) {
// 	PetscCheck(state->ops->wrimisfit, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetWriMisfitFunction() before TMMComputeMisfitFunction().");
    if (state->ops->wrimisfit) {
      {
        void *ctx;
        TMMMisfitWriFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx);
        f = state->ops->wrimisfit;
        if (state->wrimisfituserctxcontainer) {
          PetscCall(PetscContainerGetPointer(state->wrimisfituserctxcontainer, &ctx));
        } else {
          ctx = NULL;
        }		
        PetscCallBack("TMM callback function", (*f)(tc, iLoop, state, ctx));
      }
    }
  } else if (theFunc==TMM_FIN_FUNC) {
    PetscCheck(state->ops->finmisfit, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetFinMisfitFunction() before TMMComputeMisfitFunction().");
    if (state->ops->finmisfit) {
      {
        void *ctx;
        TMMMisfitFinFunctionFn *f; //PetscErrorCode (*f)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
        f = state->ops->finmisfit;
        if (state->finmisfituserctxcontainer) {
          PetscCall(PetscContainerGetPointer(state->finmisfituserctxcontainer, &ctx));
        } else {
          ctx = NULL;
        }
        PetscCallBack("TMM callback function", (*f)(tc, Iter, state, ctx));
      }
    }
  } else {
    SETERRQ(PETSC_COMM_WORLD,1,"ERROR!: Unknown function type!");
  }  

  PetscFunctionReturn(PETSC_SUCCESS);
}

// void* TMMCreateExternalForcingContext(MPI_Comm comm, TMMState state, TMMExtCtxFunctionFn *fctx)
// {
//   PetscFunctionBegin;
//   void* ctx;  
// 
//   ctx=(*fctx)(comm, state);
//   return ctx;
// //   PetscFunctionReturn(PETSC_SUCCESS);
// }

// extern TMMExtIniFunctionFn* TMMGetIniExternalForcingFunction(TMMState state);
// extern TMMExtCalcFunctionFn TMMGetCalcExternalForcingFunction(TMMState state);
// extern TMMExtWriFunctionFn TMMGetWriExternalForcingFunction(TMMState state);
// extern TMMExtFinFunctionFn TMMGetFinExternalForcingFunction(TMMState state);
// extern TMMExtReiFunctionFn TMMGetReiExternalForcingFunction(TMMState state);
// 
// TMMExtIniFunctionFn* TMMGetIniExternalForcingFunction(TMMState state)
// {
// //   PetscFunctionBegin;
// //   PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
//   TMMExtIniFunctionFn *f=NULL;
//   if (state->ops->iniexternalforcing) {
// 	f=state->ops->iniexternalforcing;
//   }
// 
//   return f;
// }

// void* TMMGetFunction(TMMState state, TMMFUNCTYPE theFunc)
// {
// //   PetscFunctionBegin;
// //   PetscValidHeaderSpecific(state, TMM_CLASSID, 1);
//   void *f=NULL;
//   if (theFunc==TMM_INI_FUNC) {
// // 	PetscCheck(state->ops->iniexternalforcing, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetExternalForcingFunctions() before TMMComputeExtForcFunction().");
// 	if (state->ops->iniexternalforcing) {
//       f=state->ops->iniexternalforcing;
//     }
//   } else if (theFunc==TMM_FIN_FUNC) {
// // 	PetscCheck(state->ops->finexternalforcing, PETSC_COMM_SELF, PETSC_ERR_ARG_WRONGSTATE, "Must call TMMSetExternalForcingFunctions() before TMMComputeExtForcFunction().");
// 	if (state->ops->finexternalforcing) {
//       f=state->ops->finexternalforcing;
//     }  
//   } else {
//     SETERRQ(PETSC_COMM_WORLD,1,"ERROR!: Unknown function type!");
//   }  
// 
//   return f;
// }
