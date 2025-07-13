#define DEFINE_JACOBIAN_VARIABLES

#include "petscmat.h"
#include "petsctime.h"

#include "tmm_forcing_utils.h"
// #include "tmm_profile_data.h"
#include "tmm_timer.h"
#include "tmm.h"
#include "tmm_share.h"
#include "tmm_jacobian_share.h"

PetscErrorCode JacobianFunction(void*,Vec,Vec,void*);

extern PetscErrorCode XtoVconvert(PetscInt numTracers, Vec *c, Vec X, PetscScalar *XTovScaleFac);
extern PetscErrorCode VtoXconvert(PetscInt numTracers, Vec *c, Vec X, PetscScalar *vToXScaleFac);

ISColoring    iscoloring;
MatFDColoring fdcoloring;
MatColoring   coloring;
Vec X;
// Mat Q;
// PetscInt lSize;

// PetscInt Itercloc, iLooploc, numTracersloc;
// PetscScalar tcloc, tfloc;

TMMJAC jacstate;

PetscErrorCode doJacobianInitialize(Mat Q, TMMState state)
{

  PetscErrorCode ierr;
  PetscViewer fd;
  PetscInt kx;
  PetscInt itr, maxValsToRead;
  PetscBool flg;  
  PetscInt numTracers;
  const char *prefix;
  
  PetscFunctionBeginUser;

  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;
  
  ierr = PetscMalloc(numTracers*sizeof(PetscScalar),&jacstate.XTovScaleFac);CHKERRQ(ierr);   
  ierr = PetscMalloc(numTracers*sizeof(PetscScalar),&jacstate.vToXScaleFac);CHKERRQ(ierr);   

  for (itr=0; itr<numTracers; itr++) {
    jacstate.XTovScaleFac[itr]=1.0;
    jacstate.vToXScaleFac[itr]=1.0;
  }      

  maxValsToRead = numTracers;
  ierr = PetscOptionsGetRealArray(NULL,prefix,"-vscale_fac",jacstate.XTovScaleFac,&maxValsToRead,&flg);
  if (flg) {
    if (maxValsToRead != numTracers) {
      SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of vscale_fac values specified");
    }
    for (itr=0; itr<numTracers; itr++) {
      if (jacstate.XTovScaleFac[itr]>0.0) {
        jacstate.vToXScaleFac[itr]=1.0/jacstate.XTovScaleFac[itr];
      }      
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d scale factor=%15.11f\n",itr,jacstate.XTovScaleFac[itr]);CHKERRQ(ierr);
    }      
  }

jacstate.numTracers=numTracers;
jacstate.tc=0.0;
jacstate.tf=0.0;
jacstate.Iterc=1;
jacstate.iLoop=1;
jacstate.state=state;

lXSize=lSize*numTracers;

// ierr = VecGetLocalSize(state.c[0],&lSize);CHKERRQ(ierr);

ierr = VecCreate(PETSC_COMM_WORLD,&X);CHKERRQ(ierr);
ierr = VecSetSizes(X,lXSize,PETSC_DECIDE);CHKERRQ(ierr);      
ierr = VecSetFromOptions(X);CHKERRQ(ierr);

/*   Compute global indices for local piece of vectors */
// Don't really need these
  ierr = VecGetOwnershipRange(X,&gXLow,&gXHigh);CHKERRQ(ierr);
  gXHigh = gXHigh - 1; /* Note: gXHigh is one more than the last local element */
  ierr = PetscMalloc(lXSize*sizeof(PetscInt),&gXIndices);CHKERRQ(ierr);  
  for (kx=0; kx<lXSize; kx++) {
    gXIndices[kx] = kx + gXLow;
  }  

// Set the sparsity pattern of Q by reading in a sparsity matrix from file
// ierr = MatCreate(PETSC_COMM_WORLD,&Q);CHKERRQ(ierr);
ierr = MatSetSizes(Q,lXSize,lXSize,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
ierr = MatSetType(Q,MATMPIAIJ);CHKERRQ(ierr);      
ierr = MatSetFromOptions(Q);CHKERRQ(ierr);  
ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"S.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
ierr = MatLoad(Q,fd);CHKERRQ(ierr);
ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

// Color Q
MatColoringCreate(Q, &coloring);
MatColoringSetType(coloring,MATCOLORINGSL); 
// MATCOLORINGSL
// MATCOLORINGJP
MatColoringSetFromOptions(coloring);
MatColoringApply(coloring, &iscoloring);
MatColoringDestroy(&coloring);
MatFDColoringCreate(Q,iscoloring, &fdcoloring);
MatFDColoringSetFromOptions(fdcoloring);
MatFDColoringSetUp(Q,iscoloring,fdcoloring);
MatFDColoringSetFunction(fdcoloring,(PetscErrorCode (*)(void))JacobianFunction, &jacstate);
MatFDColoringView(fdcoloring,PETSC_VIEWER_STDOUT_WORLD);

  PetscFunctionReturn(0);
}

PetscErrorCode doJacobianCalc(PetscScalar tc, PetscScalar tf,PetscInt Iterc, PetscInt iLoop, Mat Q, TMMState state)
{

  PetscErrorCode ierr;
  PetscViewer fd;
  
  PetscFunctionBeginUser;

jacstate.tc=tc;
jacstate.tf=tf;
jacstate.Iterc=Iterc;
jacstate.iLoop=iLoop;
jacstate.state=state;

// c to X
  ierr = VtoXconvert(jacstate.numTracers, jacstate.state->c, X, jacstate.vToXScaleFac);CHKERRQ(ierr);

// Compute Q using finite differences
MatFDColoringApply(Q,fdcoloring,X,NULL);

// Write Q to file
PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Q.petsc",FILE_MODE_WRITE,&fd);
MatView(Q,fd);
ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode doJacobianFinalize()
{

  PetscErrorCode ierr;
  
  PetscFunctionBeginUser;

ISColoringDestroy(&iscoloring);
MatFDColoringDestroy(&fdcoloring);
// MatDestroy(&Q);
VecDestroy(&X);

  PetscFunctionReturn(0);
}

PetscErrorCode JacobianFunction(void *dummy,Vec Xx,Vec F,void *ptr)
{

// #include "tmm_external_forcing.h"

  PetscErrorCode ierr;
TMMJAC         *user = (TMMJAC*)ptr;

//   PetscInt ip, ks, ke, k, itr, kx;

//   PetscScalar **localV, **localUEF;
//   PetscScalar *localX, *localF;


  PetscFunctionBeginUser;

// X -> c{}
  ierr = XtoVconvert(user->numTracers, user->state->c, Xx, user->XTovScaleFac);CHKERRQ(ierr);
//   if (numTracers>1) {
// 	ierr = VecGetArrays(user->c,numTracers,&localV);CHKERRQ(ierr);
// 	ierr = VecGetArrays(user->qef,numTracers,&localUEF);CHKERRQ(ierr);  
// 	ierr = VecGetArray(X,&localX);CHKERRQ(ierr);
// 	ierr = VecGetArray(F,&localF);CHKERRQ(ierr);
//   }

//   X to c
//   if (numTracers==1) {
//     VecCopy(X,user->c[0]);
//   } else {
// 	kx=0;
// 	for (ip=0; ip<lNumProfiles; ip++) {
// 	  ks=lStartIndices[ip];  
// 	  ke=lEndIndices[ip];
// 	  for (itr=0; itr<numTracers; itr++) {
// 		for (k=ks; k<=ke; k++) {
// 		  localV[itr][k]=localX[kx];
// 		  kx++;
// 		}
// 	  }
// 	}
// 	for (itr=0; itr<numTracers; itr++) {  
// 	  ierr = VecSetValues(user->c[itr],lSize,gIndices,localV[itr],INSERT_VALUES);CHKERRQ(ierr);
// 	  ierr = VecAssemblyBegin(user->c[itr]);CHKERRQ(ierr);
// 	  ierr = VecAssemblyEnd(user->c[itr]);CHKERRQ(ierr);    
// 	}
//   }

    ierr = TMMComputeExternalForcingFunction(user->tc,user->Iterc,user->iLoop,user->state,TMM_CALC_FUNC);
//   ierr = calcExternalForcing(user->tc,user->Iterc,user->iLoop,user->numTracers,user->state->c,user->state->qef);CHKERRQ(ierr); /* Compute external forcing in qef */

//   qef{} to F
  ierr = VtoXconvert(user->numTracers, user->state->qef, F, user->vToXScaleFac);CHKERRQ(ierr);

//   if (numTracers==1) {
// 	VecCopy(user->qef[0],F);
//   } else {
// 	kx=0;
// 	for (ip=0; ip<lNumProfiles; ip++) {
// 	  ks=lStartIndices[ip];  
// 	  ke=lEndIndices[ip];
// 	  for (itr=0; itr<numTracers; itr++) {
// 		for (k=ks; k<=ke; k++) {
// 		  localF[kx]=localUEF[itr][k];
// 		  kx++;
// 		}
// 	  }
// 	}
// 	ierr = VecSetValues(F,lXSize,gXIndices,localF,INSERT_VALUES);CHKERRQ(ierr);
// 	ierr = VecAssemblyBegin(F);CHKERRQ(ierr);
// 	ierr = VecAssemblyEnd(F);CHKERRQ(ierr);        
//   }
  
  PetscFunctionReturn(0);
}

PetscErrorCode XtoVconvert(PetscInt numTracers, Vec *c, Vec X, PetscScalar *XTovScaleFac)
{

  PetscErrorCode ierr;

  PetscInt ip, ks, ke, k, itr, kx;

  PetscScalar **localV;
  const PetscScalar *localX;


  PetscFunctionBeginUser;

  if (numTracers==1) {
	VecCopy(X,c[0]);
	ierr = VecScale(c[0],XTovScaleFac[0]);CHKERRQ(ierr);
  } else {
	ierr = VecGetArrays(c,numTracers,&localV);CHKERRQ(ierr);
	ierr = VecGetArrayRead(X,&localX);CHKERRQ(ierr);
	if (useProfiles) {
	  kx=0;
	  for (ip=0; ip<lNumProfiles; ip++) {
		ks=lStartIndices[ip];  
		ke=lEndIndices[ip];
		for (itr=0; itr<numTracers; itr++) {
		  for (k=ks; k<=ke; k++) {
			localV[itr][k]=XTovScaleFac[itr]*localX[kx];
			kx++;
		  }
		}
	  }
	} else {
	  kx=0;
	  for (itr=0; itr<numTracers; itr++) {
		for (k=0; k<lSize; k++) {
		  localV[itr][k]=XTovScaleFac[itr]*localX[kx];
		  kx++;
		}
	  }
	}  
    ierr = VecRestoreArrays(c,numTracers,&localV);CHKERRQ(ierr);
	ierr = VecRestoreArrayRead(X,&localX);CHKERRQ(ierr);
// 	for (itr=0; itr<numTracers; itr++) {  
// 	  ierr = VecSetValues(c[itr],lSize,gIndices,localV[itr],INSERT_VALUES);CHKERRQ(ierr);
// 	  ierr = VecAssemblyBegin(c[itr]);CHKERRQ(ierr);
// 	  ierr = VecAssemblyEnd(c[itr]);CHKERRQ(ierr);    
// 	}
  }
  
  PetscFunctionReturn(0);
}

PetscErrorCode VtoXconvert(PetscInt numTracers, Vec *c, Vec X, PetscScalar *vToXScaleFac)
{

  PetscErrorCode ierr;

  PetscInt ip, ks, ke, k, itr, kx;

  PetscScalar **localV;
  PetscScalar *localX;


  PetscFunctionBeginUser;

//   if (numTracers>1) {
//   }

  if (numTracers==1) {
	VecCopy(c[0],X);
	ierr = VecScale(c[0],vToXScaleFac[0]);CHKERRQ(ierr);	
  } else {
	ierr = VecGetArrays(c,numTracers,&localV);CHKERRQ(ierr);
	ierr = VecGetArray(X,&localX);CHKERRQ(ierr);
	if (useProfiles) {
	  kx=0;
	  for (ip=0; ip<lNumProfiles; ip++) {
		ks=lStartIndices[ip];  
		ke=lEndIndices[ip];
		for (itr=0; itr<numTracers; itr++) {
		  for (k=ks; k<=ke; k++) {
			localX[kx]=vToXScaleFac[itr]*localV[itr][k];
			kx++;
		  }
		}
	  }
	} else {
	  kx=0;
	  for (itr=0; itr<numTracers; itr++) {
		for (k=0; k<lSize; k++) {
		  localX[kx]=vToXScaleFac[itr]*localV[itr][k];
		  kx++;
		}
	  }
	}  
    ierr = VecRestoreArrays(c,numTracers,&localV);CHKERRQ(ierr);
	ierr = VecRestoreArray(X,&localX);CHKERRQ(ierr);
	
// 	ierr = VecSetValues(X,lXSize,gXIndices,localX,INSERT_VALUES);CHKERRQ(ierr);
// 	ierr = VecAssemblyBegin(X);CHKERRQ(ierr);
// 	ierr = VecAssemblyEnd(X);CHKERRQ(ierr);        
  }
  
  PetscFunctionReturn(0);
}

