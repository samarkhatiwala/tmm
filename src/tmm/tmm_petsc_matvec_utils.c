#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include "petscmat.h"
#include "tmm_petsc_matvec_utils.h"

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
  int fp;
  off_t off,offset;
      
  PetscFunctionBegin;
  off = PETSC_BINARY_SCALAR_SIZE*(length+1)*(iRec-1);  
  ierr = PetscViewerBinaryGetDescriptor(viewer,&fp);CHKERRQ(ierr);
  ierr = PetscBinarySeek(fp,off,PETSC_BINARY_SEEK_SET,&offset);CHKERRQ(ierr);
  ierr = VecLoad(vec,viewer);CHKERRQ(ierr); /* IntoVector */
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecLoadIntoVectorRandomAccessFromFile"
PetscErrorCode VecLoadIntoVectorRandomAccessFromFile(const char fileName[], Vec vec, PetscInt length, PetscInt iRec)
{
/* Random access version of VecLoadIntoVector */
/* This version takes two additional arguments:  */
/*   length: length of the vector */
/*   iRec: the record to read (iRec=1 is the first record; iRec=numRecs is the last record) */
  
  PetscErrorCode ierr;
  PetscViewer fd;    
      
  PetscFunctionBegin;
	 ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fileName,FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = VecLoadIntoVectorRandomAccess(fd, vec, length, iRec);CHKERRQ(ierr);
	 ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecLoadVecIntoArray"
PetscErrorCode VecLoadVecIntoArray(Vec c, const char filename[], PetscScalar *arr)
{
/* Loads vector from file into supplied array with same parallel layout as supplied Vec c */
/* You MUST preallocate memory for the array before calling this routine */
/* This routine creates a temporary Vec from the array */

  Vec tmpVec;  
  PetscErrorCode ierr;
  PetscViewer fd;
  PetscInt lDim;
      
  PetscFunctionBegin;
// Not sure what to set for block size (the 2nd argument) so I'm just setting it to 1 for now
  ierr = VecGetLocalSize(c,&lDim);CHKERRQ(ierr);
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, lDim, PETSC_DETERMINE, arr, &tmpVec);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = VecLoad(tmpVec,fd);CHKERRQ(ierr); /* IntoVector */
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
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
  int       fp;
  PetscViewer    viewer;
      
  PetscFunctionBegin;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,FILE_MODE_READ,&viewer);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(viewer,&fp);CHKERRQ(ierr);

/* Read header */
  ierr = PetscBinaryRead(fp,(char *)header,4,NULL,PETSC_INT);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);

  /* error checking on file */
  if (header[0] != MAT_FILE_CLASSID) SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_FILE_UNEXPECTED,"not matrix object");

  if (M) *M = header[1];
  if (N) *N = header[2];
  if (nnz) *nnz = header[3];

  PetscFunctionReturn(0);
}

// #undef __FUNCT__
// #define __FUNCT__ "VecWrite"
// PetscErrorCode VecWrite(const char *fileName, Vec c, PetscBool appendToFile)
// {
//   PetscErrorCode ierr;
//   PetscViewer viewer;
// 
//   if (appendToFile) {
// 	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fileName,FILE_MODE_APPEND,&viewer);CHKERRQ(ierr);
//   } else {
// 	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,fileName,FILE_MODE_WRITE,&viewer);CHKERRQ(ierr);
//   }  
// 
//   ierr = VecView(c,viewer);CHKERRQ(ierr);
//   ierr = PetscViewerDestroy(&viewer);CHKERRQ(ierr);      
// 
//   return 0;
// }

#undef __FUNCT__
#define __FUNCT__ "VecLoadIntoArray"
PetscErrorCode VecLoadIntoArray(PetscInt lDim, const char filename[], PetscScalar *arr)
{
/* Loads vector from file into supplied array of local dimension lDim */
/* You MUST preallocate memory for the array before calling this routine */
/* This routine creates a temporary Vec from the array */

  Vec tmpVec;  
  PetscErrorCode ierr;
  PetscViewer fd;
      
  PetscFunctionBegin;
// Not sure what to set for block size (the 2nd argument) so I'm just setting it to 1 for now
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, lDim, PETSC_DETERMINE, arr, &tmpVec);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = VecLoad(tmpVec,fd);CHKERRQ(ierr); /* IntoVector */
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  ierr = VecDestroy(&tmpVec);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecCreateFromLocalSize"
PetscErrorCode VecCreateFromLocalSize(PetscInt lDim, Vec *c)
{
/* Create and return a MPI vector of given local size */
/* Primarily for tmm4py support */

  PetscErrorCode ierr;
      
  PetscFunctionBegin;
  ierr = VecCreate(PETSC_COMM_WORLD,c);CHKERRQ(ierr);
  ierr = VecSetType(*c,VECMPI);CHKERRQ(ierr);
  ierr = VecSetSizes(*c,lDim,PETSC_DETERMINE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecCreateWithArrayFromLocalSize"
PetscErrorCode VecCreateWithArrayFromLocalSize(PetscInt lDim, PetscScalar *arr, Vec *c)
{
/* Create and return a MPI vector using a supplied array to hold the data */
/* Primarily for tmm4py support */

  PetscErrorCode ierr;
      
  PetscFunctionBegin;
// Not sure what to set for block size (the 2nd argument) so I'm just setting it to 1 for now
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, lDim, PETSC_DETERMINE, arr, c);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "VecWriteLocalArrayToVec"
PetscErrorCode VecWriteLocalArrayToVec(PetscInt lDim, PetscScalar *arr, const char filename[], PetscFileMode mode)
{
/* Write an array of local dimension lDim to a PETSc Vec file. This routine creates a temporary  */
/* Vec from the array and then writes it to file */

  PetscErrorCode ierr;
  Vec tmpVec;
  PetscViewer fd;
  
  PetscFunctionBegin;
// Not sure what to set for block size (the 2nd argument) so I'm just setting it to 1 for now
  ierr = VecCreateMPIWithArray(PETSC_COMM_WORLD, 1, lDim, PETSC_DETERMINE, arr, &tmpVec);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,filename,mode,&fd);CHKERRQ(ierr);
  ierr = VecView(tmpVec,fd);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  ierr = VecDestroy(&tmpVec);CHKERRQ(ierr);   

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "dotProdArrays"
PetscErrorCode dotProdArrays(PetscScalar *xarr, PetscScalar *yarr, PetscInt n, PetscScalar *z)
{
/* Function to take a dot product of two arrays. */

  PetscInt i;
  PetscScalar localz = 0.0;

  *z = 0.0;
  for (i=0; i<n; i++) {
	   localz = localz + xarr[i]*yarr[i];
  }
  MPI_Allreduce(&localz, z, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "sumArray"
PetscErrorCode sumArray(PetscScalar *xarr, PetscInt n, PetscScalar *z)
{
/* Function to sum an array. */

  PetscInt i;
  PetscScalar localz = 0.0;

  *z = 0.0;
  for (i=0; i<n; i++) {
	   localz = localz + xarr[i];
  }
  MPI_Allreduce(&localz, z, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "sumScalar"
PetscErrorCode sumScalar(PetscScalar x, PetscScalar *tot)
{
/* Function to sum a scalar. */

  *tot = 0.0;
  MPI_Allreduce(&x, tot, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "Barrier"
PetscErrorCode Barrier(void)
{
/* MPI_barrier primarily for tmm4py support */

  MPI_Barrier(PETSC_COMM_WORLD);

  return 0;
}
