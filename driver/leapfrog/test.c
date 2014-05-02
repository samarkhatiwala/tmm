static char help[] = "\n";
/* 
  Include "petscmat.h" so that we can use matrices.
  automatically includes:
     petsc.h       - base PETSc routines   petscvec.h    - vectors
     petscsys.h    - system routines       petscmat.h    - matrices
     petscis.h     - index sets            petscviewer.h - viewers               
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"

#include "petsc_matvec_utils.h"
#include "tmm_forcing_utils.h"

PetscErrorCode output(Vec *v, PetscInt numTracers);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  PetscErrorCode ierr;
  PetscInt itr, itl;
  PetscInt numTracers, n;
  Vec templateVec;
  Vec **v;
  Vec *vold, *vcur, *vnew;
  PetscInt iold, icur, inew;
  PetscViewer fd;

  PetscInitialize(&argc,&args,(char *)0,help);

  numTracers=5;
  n=100;
  iold = 0;
  icur = 1;
  inew = 2;
  
  ierr = VecCreate(PETSC_COMM_WORLD,&templateVec);CHKERRQ(ierr);
  ierr = VecSetSizes(templateVec,PETSC_DECIDE,n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(templateVec);CHKERRQ(ierr);
  
/* tracer vectors */
//   v = malloc(numTracers*sizeof(Vec *));
//   vold = malloc(numTracers*sizeof(Vec *));
//   vcur = malloc(numTracers*sizeof(Vec *));
//   vnew = malloc(numTracers*sizeof(Vec *));

  v = malloc(3*sizeof(Vec *));
  vold = malloc(numTracers*sizeof(Vec *));
  vcur = malloc(numTracers*sizeof(Vec *));
  vnew = malloc(numTracers*sizeof(Vec *));

//   for (itr=0; itr<numTracers; itr++) {   
// 	ierr = VecDuplicateVecs(templateVec,3,&v[itr]);CHKERRQ(ierr);
//   }	

  for (itl=0; itl<3; itl++) {
	ierr = VecDuplicateVecs(templateVec,numTracers,&v[itl]);CHKERRQ(ierr);
  }	

  for (itl=0; itl<3; itl++) {
	for (itr=0; itr<numTracers; itr++) {	
      VecSet(v[itl][itr],1.0*(itl+1.0));
    }
  }  

//  for (itr=0; itr<numTracers; itr++) {
//	ierr = VecCopy(v[icur][itr],v[iold][itr]);CHKERRQ(ierr);
//  }  

  for (itr=0; itr<numTracers; itr++) {   
    vold[itr]=v[iold][itr];
    vcur[itr]=v[icur][itr];
    vnew[itr]=v[inew][itr];    
  }

  for (itr=0; itr<numTracers; itr++) {   
    vold[itr]=vnew[itr];
  }

//  for (itr=0; itr<numTracers; itr++) {	
//	VecSet(vold[itr],20.0);
//  }
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing old\n");CHKERRQ(ierr); 
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"old.petsc",FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
  for (itr=0; itr<numTracers; itr++) {
    ierr = VecView(v[iold][itr],fd);CHKERRQ(ierr);
  }
  ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);      

//  ierr = output(v[iold],numTracers);

  ierr = output(vold,numTracers);
//   ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"old1.petsc",FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
//   for (itr=0; itr<numTracers; itr++) {
//     ierr = VecView(vold[itr],fd);CHKERRQ(ierr);
//   }
//   ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);      
  
  ierr = PetscFinalize();CHKERRQ(ierr);
  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "output"
PetscErrorCode output(Vec *v, PetscInt numTracers)
{

  PetscInt itr;
  PetscErrorCode ierr;
  PetscViewer fd;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing old1\n");CHKERRQ(ierr); 
  
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"old1.petsc",FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
  for (itr=0; itr<numTracers; itr++) {
    ierr = VecView(v[itr],fd);CHKERRQ(ierr);
  }
  ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);      


  return 0;
}
  