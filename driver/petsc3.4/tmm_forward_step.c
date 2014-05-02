#include "petscmat.h"
#include "petsc_matvec_utils.h"
#include "tmm_forcing_utils.h"

#undef __FUNCT__
#define __FUNCT__ "forwardStep"
PetscErrorCode forwardStep(PetscScalar tc, PetscInt iLoop, PetscScalar dt, PetscInt numTracers, 
                           PetscBool useForcingFromFile, PetscBool useExternalForcing, PetscBool usePrescribedBC, 
                           Vec *v, Mat Ae, Mat Ai, Mat Be, Mat Bi, Vec *uf, Vec *uef, Vec *bcc, Vec *bcf, Vec *vtmp)
{

  PetscInt itr;
  PetscErrorCode ierr;
  PetscScalar one = 1.0;
/*   PetscInt bcCutOffStep = 720; */

/* iLoop -> iLoop+1 (convention) */
  for (itr=0; itr<numTracers; itr++) {
	ierr = MatMult(Ae,v[itr],vtmp[itr]);CHKERRQ(ierr); /* vtmp <- Ae*v */
  }

  if (useForcingFromFile) {
	for (itr=0; itr<numTracers; itr++) {
	ierr = VecAXPY(vtmp[itr],one,uf[itr]);CHKERRQ(ierr); /* vtmp <- vtmp + uf */
	}
  }

  if (useExternalForcing) {
	for (itr=0; itr<numTracers; itr++) {
	ierr = VecAXPY(vtmp[itr],one,uef[itr]);CHKERRQ(ierr); /* vtmp <- vtmp + uef */
	}
  }

/* vtmp is now:  vtmp = Ae*v + [uf] + [uef] */

  if (usePrescribedBC) {
	for (itr=0; itr<numTracers; itr++) {
	  ierr = MatMultAdd(Be,bcc[itr],vtmp[itr],vtmp[itr]);CHKERRQ(ierr); /* vtmp <- Be*bcc + vtmp = Be*bcc + Ae*v + [uf] + [uef] */
	}
    for (itr=0; itr<numTracers; itr++) {
      ierr = MatMult(Bi,bcf[itr],v[itr]);CHKERRQ(ierr); /* v <- Bi*bcf */
    }
	for (itr=0; itr<numTracers; itr++) {
	  ierr = MatMultAdd(Ai,vtmp[itr],v[itr],v[itr]);CHKERRQ(ierr); /* v <- Ai*vtmp + v = Ai*(Be*bcc + Ae*v + [uf] + [uef]) + Bi*bcf */
	}
  } else {
    for (itr=0; itr<numTracers; itr++) {
      ierr = MatMult(Ai,vtmp[itr],v[itr]);CHKERRQ(ierr); /* v <- Ai*vtmp = Ai*(Ae*v + [uf] + [uef]) */
    }
  }

  return 0;
}
  