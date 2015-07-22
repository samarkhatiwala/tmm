#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include "petscmat.h"
#include "petsc_matvec_utils.h"

#undef __FUNCT__
#define __FUNCT__ "MatAXPBYmy"
PetscErrorCode MatAXPBYmy(PetscScalar a,PetscScalar b,Mat X,Mat Y,Mat *Z)
{
/* Compute Z = a*X + b*Y, where X,Y,Z are matrices, and a,b are scalars */
/* If Z doesn't already exist, it is created */

   PetscErrorCode ierr;

   if ((*Z)) {
     ierr=MatCopy(Y,*Z,SAME_NONZERO_PATTERN);CHKERRQ(ierr);  /* Z=Y */
   } else {
     ierr=MatDuplicate(Y,MAT_COPY_VALUES,Z);CHKERRQ(ierr);   
   }
   
   ierr=MatScale(*Z,b);CHKERRQ(ierr);  /* Z=b*Z */
   ierr=MatAXPY(*Z,a,X,SAME_NONZERO_PATTERN);CHKERRQ(ierr); /* Z=a*X+Z */

   return 0;    
}

#undef __FUNCT__
#define __FUNCT__ "VecAXPBYmy"
PetscErrorCode VecAXPBYmy(PetscScalar a,PetscScalar b,Vec x,Vec y,Vec *z)
{
/* Compute z = a*x + b*y, where x,y,z are vectors, and a,b are scalars */
/* If z doesn't already exist, it is created */

   PetscErrorCode ierr;

   if ((*z)) {    
     ierr=VecCopy(y,*z);CHKERRQ(ierr);  /* z=y */
   } else {
     ierr=VecDuplicate(y,z);CHKERRQ(ierr);
   }    

   ierr=VecScale(*z,b);CHKERRQ(ierr);  /* z=b*z */ 
   ierr=VecAXPY(*z,a,x); CHKERRQ(ierr); /* z=a*x+z */

   return 0;    
}

#undef __FUNCT__
#define __FUNCT__ "VecLoadIntoVectorRandomAccess"
PetscErrorCode VecLoadIntoVectorRandomAccess(PetscViewer viewer,Vec vec, PetscInt length, PetscInt iRec)
{
/* Random access version of VecLoadIntoVector */
/* This version takes two additional arguments:  */
/*   length: length of the vector */
/*   iRec: the record to read (iRec=1 is the first record; iRec=numRecs is the last record) */
  
  PetscErrorCode ierr;
  PetscInt fd;
  off_t off,offset;
      
  PetscFunctionBegin;
  off = PETSC_BINARY_SCALAR_SIZE*(length+1)*(iRec-1);  
  ierr = PetscViewerBinaryGetDescriptor(viewer,&fd);CHKERRQ(ierr);
  ierr = PetscBinarySeek(fd,off,PETSC_BINARY_SEEK_SET,&offset);CHKERRQ(ierr);
  ierr = VecLoad(vec,viewer);CHKERRQ(ierr); /* IntoVector */
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecLoadVecIntoArray"
PetscErrorCode VecLoadVecIntoArray(Vec v, const char filename[], PetscScalar *arr)
{
/* Loads vector from file and copies over local piece of the data to an array */
/* You MUST preallocate memory for the array before calling this routine */
/* Input vector v is used as a template so that processor partioning is preserved  */
/* as desired */

  Vec tmpVec;  
  PetscScalar *localtmpVec;
  PetscErrorCode ierr;
  PetscViewer fd;
  PetscInt lDim, i;
      
  PetscFunctionBegin;
/* Load data into vector   */
  ierr = VecDuplicate(v,&tmpVec);CHKERRQ(ierr);  
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = VecLoad(tmpVec,fd);CHKERRQ(ierr); /* IntoVector */
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

/* Get pointer to local data */
  ierr = VecGetArray(tmpVec,&localtmpVec);CHKERRQ(ierr);
  ierr = VecGetLocalSize(tmpVec,&lDim);CHKERRQ(ierr);

/* Copy data to array. First allocate memory   */
/*   ierr = PetscMalloc(lDim*sizeof(PetscScalar),&arr);CHKERRQ(ierr); */
  for (i=0; i<lDim; i++) {
    arr[i]=localtmpVec[i];
  }

/* Destroy objects   */
  ierr = VecRestoreArray(tmpVec,&localtmpVec);CHKERRQ(ierr);
  ierr = VecDestroy(&tmpVec);CHKERRQ(ierr);   

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatGetSizeFromFile"
PetscErrorCode MatGetSizeFromFile(const char filename[], PetscInt *M, PetscInt *N, PetscInt *nnz)
{
/* Function to read matrix size from file filename. */
/* USAGE: MatGetSizeFromFile(filename,&M,&N,&nnz); */
/* You can pass PETSC_NULL for any of the size arguments */

  PetscErrorCode ierr;
  PetscInt       header[4];
  PetscInt       fd;
  PetscViewer    viewer;
      
  PetscFunctionBegin;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(viewer,&fd);CHKERRQ(ierr);

/* Read header */
  ierr = PetscBinaryRead(fd,(char *)header,4,PETSC_INT);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  /* error checking on file */
  if (header[0] != MAT_FILE_CLASSID) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FILE_UNEXPECTED,"not matrix object");

  if (M) *M = header[1];
  if (N) *N = header[2];
  if (nnz) *nnz = header[3];

  PetscFunctionReturn(0);
}