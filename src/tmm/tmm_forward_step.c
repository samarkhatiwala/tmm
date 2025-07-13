#include "petscmat.h"
#include "tmm_petsc_matvec_utils.h"
#include "tmm_forcing_utils.h"

#undef __FUNCT__
#define __FUNCT__ "TMMForwardStep"
PetscErrorCode TMMForwardStep(PetscScalar tc, PetscInt iLoop, PetscScalar dt, PetscInt numTracers, 
                           PetscBool applyForcingFromFile, PetscBool applyExternalForcing, PetscBool relaxTracer, PetscBool applyBC, 
                           Vec *c, Mat Ae, Mat Ai, Mat Be, Mat Bi, Vec *qf, Vec *qef, Vec *qrel, Vec *cbc, Vec *cbf, Vec *vtmp)
{

  PetscInt itr;
  PetscErrorCode ierr;
  PetscScalar one = 1.0;

/* iLoop -> iLoop+1 (convention) */
  for (itr=0; itr<numTracers; itr++) {
	ierr = MatMult(Ae,c[itr],vtmp[itr]);CHKERRQ(ierr); /* vtmp <- Ae*c */
  }

  if (applyForcingFromFile) {
	for (itr=0; itr<numTracers; itr++) {
	  ierr = VecAXPY(vtmp[itr],one,qf[itr]);CHKERRQ(ierr); /* vtmp <- vtmp + qf */
	}
  }

  if (applyExternalForcing) {
	for (itr=0; itr<numTracers; itr++) {
	  ierr = VecAXPY(vtmp[itr],one,qef[itr]);CHKERRQ(ierr); /* vtmp <- vtmp + qef */
	}
  }

  if (relaxTracer) {
	for (itr=0; itr<numTracers; itr++) {
	  ierr = VecAXPY(vtmp[itr],one,qrel[itr]);CHKERRQ(ierr); /* vtmp <- vtmp + qrel */
	}
  }

/* vtmp is now:  vtmp = Ae*c + [qf] + [qef] + [qrel] */

  if (applyBC) {
	for (itr=0; itr<numTracers; itr++) {
	  ierr = MatMultAdd(Be,cbc[itr],vtmp[itr],vtmp[itr]);CHKERRQ(ierr); /* vtmp <- Be*cbc + vtmp = Be*cbc + Ae*c + [qf] + [qef] */
	}
    for (itr=0; itr<numTracers; itr++) {
      ierr = MatMult(Bi,cbf[itr],c[itr]);CHKERRQ(ierr); /* c <- Bi*cbf */
    }
	for (itr=0; itr<numTracers; itr++) {
	  ierr = MatMultAdd(Ai,vtmp[itr],c[itr],c[itr]);CHKERRQ(ierr); /* c <- Ai*vtmp + c = Ai*(Be*cbc + Ae*c + [qf] + [qef]) + Bi*cbf */
	}
  } else {
    for (itr=0; itr<numTracers; itr++) {
      ierr = MatMult(Ai,vtmp[itr],c[itr]);CHKERRQ(ierr); /* c <- Ai*vtmp = Ai*(Ae*c + [qf] + [qef]) */
    }
  }

  return 0;
}
  