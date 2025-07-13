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

#undef __FUNCT__
#define __FUNCT__ "TMMFinalize"
PetscErrorCode TMMFinalize(PetscScalar tc)
{

  PetscErrorCode ierr;
  
  /* Free data structures */
  ierr = VecDestroy(&templateVec);CHKERRQ(ierr);  
  ierr = MatDestroy(&Ae);CHKERRQ(ierr);
  ierr = MatDestroy(&Ai);CHKERRQ(ierr);

  if (periodicMatrix) {
    ierr = PeriodicMatDestroy(&Aep);CHKERRQ(ierr);
    ierr = PeriodicMatDestroy(&Aip);CHKERRQ(ierr);    
  } else if (timeDependentMatrix) {
	ierr = TimeDependentMatDestroy(&Aetd);
	ierr = TimeDependentMatDestroy(&Aitd);
  }

  if (rescaleForcing) {
	ierr = VecDestroy(&Rfs);CHKERRQ(ierr);
	if (periodicMatrix) {
      ierr = PeriodicVecDestroy(&Rfsp);CHKERRQ(ierr);
    } else if (timeDependentMatrix) {
      ierr = TimeDependentVecDestroy(&Rfstd);CHKERRQ(ierr);
    }
  }

  if (prescribedBCInUse) {
    ierr = VecDestroy(&bcTemplateVec);CHKERRQ(ierr);
    ierr = MatDestroy(&Be);CHKERRQ(ierr);
    ierr = MatDestroy(&Bi);CHKERRQ(ierr);
	if (periodicMatrix) {
	  ierr = PeriodicMatDestroy(&Bep);CHKERRQ(ierr);
	  ierr = PeriodicMatDestroy(&Bip);CHKERRQ(ierr);    
	} else if (timeDependentMatrix) {
	  ierr = TimeDependentMatDestroy(&Betd);
	  ierr = TimeDependentMatDestroy(&Bitd);
	}
	if (calcBCInUse) {
      ierr = PetscFree(gBCIndices);CHKERRQ(ierr);
	}
  }

  ierr = PetscFree(gIndices);CHKERRQ(ierr);

  return 0;
}
