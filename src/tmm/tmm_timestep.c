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
#include "tmm_variables.h"

extern PetscErrorCode TMMForwardStep(PetscScalar tc, PetscInt iLoop, PetscScalar dt, PetscInt numTracers, 
                                  PetscBool useForcingFromFile, PetscBool useExternalForcing, PetscBool relaxTracer, PetscBool usePrescribedBC, 
                                  Vec *c, Mat Ae, Mat Ai, Mat Be, Mat Bi, Vec *qf, Vec *qef, Vec *qrel, Vec *cbc, Vec *cbf, Vec *vtmp);

#undef __FUNCT__
#define __FUNCT__ "TMMTimeStep"
PetscErrorCode TMMTimeStep(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, TMMState state)
{

	PetscErrorCode ierr;

    ierr = TMMForwardStep(tc,Iterc,deltaTClock,state->numTracers,state->applyForcingFromFile,state->applyExternalForcing,state->relaxTracer,state->applyBC,
                       state->c,Ae,Ai,Be,Bi,state->qf,state->qef,state->qrel,state->cbc,state->cbf,state->cwork);CHKERRQ(ierr);

        return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TMMTimeStepPost"
PetscErrorCode TMMTimeStepPost(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, TMMState state)
{

// Note: tc is the time at the end of the time step

    PetscErrorCode ierr;

    if (state->useMonitor) {
      ierr = TMMComputeMonitorFunction(tc,Iterc,iLoop,state,TMM_CALC_FUNC);
      ierr = TMMComputeMonitorFunction(tc,Iterc,iLoop,state,TMM_WRI_FUNC);
//       ierr = updateMonitor(tc,iLoop,state);CHKERRQ(ierr);
//       ierr = writeMonitor(tc,iLoop,state);CHKERRQ(ierr);
    }

    if (state->doMisfit) {
      ierr = TMMComputeMisfitFunction(tc,Iterc,iLoop,state,TMM_CALC_FUNC);
      ierr = TMMComputeMisfitFunction(tc,Iterc,iLoop,state,TMM_WRI_FUNC);
//       ierr = calcMisfit(tc,iLoop,state);CHKERRQ(ierr);
//       ierr = writeMisfit(tc,iLoop,state);CHKERRQ(ierr);
    }

	return 0;
}
