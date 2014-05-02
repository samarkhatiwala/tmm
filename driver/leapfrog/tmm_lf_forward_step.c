#include "petscmat.h"
#include "petsc_matvec_utils.h"
#include "tmm_forcing_utils.h"

#undef __FUNCT__
#define __FUNCT__ "forwardStep"
PetscErrorCode forwardStep(PetscScalar tc, PetscInt iLoop, PetscScalar dt, PetscInt numTracers, 
                           PetscTruth useForcingFromFile, PetscTruth useExternalForcing, PetscTruth usePrescribedBC, 
                           Vec *vold, Vec *vcur, Vec *vnew, Mat At, Mat Ad, Mat Ai, Mat Be, Mat Bi, 
                           Vec *uf, Vec *uef, Vec *bcc, Vec *bcf, Vec *vtmp, PetscScalar deltaT)
{

  PetscInt itr;
  PetscErrorCode ierr;
  PetscScalar one = 1.0;
  
/*   PetscInt bcCutOffStep = 720; */

/* iLoop -> iLoop+1 (convention) */
  for (itr=0; itr<numTracers; itr++) {
	ierr = MatMult(At,vcur[itr],vtmp[itr]);CHKERRQ(ierr); /* vtmp <- Ae*v^(n) */
  }

  for (itr=0; itr<numTracers; itr++) {
	ierr = MatMultAdd(Ad,vold[itr],vtmp[itr],vtmp[itr]);CHKERRQ(ierr); /* vtmp <- Ad*vold + vtmp = Ad*vold + At*vcur */
  }
//MatMultAdd(Mat A,Vec x,Vec y,Vec w);
//w=A∗x+y
//y and w can be identical

  for (itr=0; itr<numTracers; itr++) {
	ierr = VecAYPX(vtmp[itr],2.0*deltaT,vold[itr]);CHKERRQ(ierr); /* vtmp <- v^(n-1) + 2*dt*vtmp */
  }

/* Note: forcing is already multiped by deltaT; so here we only need to multiply by 2 for leapfrog */
  if (useForcingFromFile) {
	for (itr=0; itr<numTracers; itr++) {
	  ierr = VecAXPY(vtmp[itr],2.0,uf[itr]);CHKERRQ(ierr); /* vtmp <- vtmp + 2*uf */
	}
  }

  if (useExternalForcing) {
	for (itr=0; itr<numTracers; itr++) {
	  ierr = VecAXPY(vtmp[itr],2.0,uef[itr]);CHKERRQ(ierr); /* vtmp <- vtmp + 2*uef */
	}
  }

/* vtmp is now:  vtmp = v^(n-1) + 2*dt*(Ae*v^(n) + f^(n)) */

  if (usePrescribedBC) {
	for (itr=0; itr<numTracers; itr++) {
	  ierr = MatMultAdd(Be,bcc[itr],vtmp[itr],vtmp[itr]);CHKERRQ(ierr); /* vtmp <- Be*bcc + vtmp = Be*bcc + Ae*v^(n) + [uf] + [uef] */
	}
    for (itr=0; itr<numTracers; itr++) {
      ierr = MatMult(Bi,bcf[itr],vtmp[itr]);CHKERRQ(ierr); /* v <- Bi*bcf */
    }
	for (itr=0; itr<numTracers; itr++) {
	  ierr = MatMultAdd(Ai,vtmp[itr],vcur[itr],vnew[itr]);CHKERRQ(ierr); /* v <- Ai*vtmp + v = Ai*(Be*bcc + Ae*v^(n) + [uf] + [uef]) + Bi*bcf */
	}
  } else {
    for (itr=0; itr<numTracers; itr++) {
      ierr = MatMult(Ai,vtmp[itr],vnew[itr]);CHKERRQ(ierr); /* v <- Ai*vtmp = Ai*(v^(n-1) + 2*dt*(Ae*v^(n) + f^(n))) */
    }
  }

  return 0;
}
  