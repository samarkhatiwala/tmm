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

   PetscTruth flg;
   PetscErrorCode ierr;

   MatValid(*Z,&flg);
   if (!flg) {
     ierr=MatDuplicate(Y,MAT_COPY_VALUES,Z);CHKERRQ(ierr);
   } else {
     ierr=MatCopy(Y,*Z,SAME_NONZERO_PATTERN);CHKERRQ(ierr);  /* Z=Y */
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

   PetscTruth flg;
   PetscErrorCode ierr;

   VecValid(*z,&flg);
   if (!flg) {    
     ierr=VecDuplicate(y,z);CHKERRQ(ierr);
   } else {
     ierr=VecCopy(y,*z);CHKERRQ(ierr);  /* z=y */
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
  ierr = VecLoadIntoVector(viewer,vec);CHKERRQ(ierr);  
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
  ierr = VecLoadIntoVector(fd,tmpVec);CHKERRQ(ierr);  
  ierr = PetscViewerDestroy(fd);CHKERRQ(ierr);

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
  ierr = VecDestroy(tmpVec);CHKERRQ(ierr);   

  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatLoadIntoMatrix"
PetscErrorCode MatLoadIntoMatrix(PetscViewer viewer, Mat A)
{
  
  PetscErrorCode ierr;
  Mat tmpMat;
  const PetscInt *cols;
  PetscInt I, Istart, Iend, ncols;
  const PetscScalar *vals;
      
  PetscFunctionBegin;
  ierr = MatLoad(viewer,MATMPIAIJ,&tmpMat);CHKERRQ(ierr);
  ierr = MatGetOwnershipRange(tmpMat,&Istart,&Iend); 
/* Note: Iend is one more than the global index of the last local row */  
  for (I=Istart; I<Iend; I++) {
    ierr = MatGetRow(tmpMat,I,&ncols,&cols,&vals);CHKERRQ(ierr);
    ierr = MatSetValues(A,1,&I,ncols,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatRestoreRow(tmpMat,I,&ncols,&cols,&vals);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatDestroy(tmpMat);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatLoadIntoMatrix2"
PetscErrorCode MatLoadIntoMatrix2(const char filename[], Mat A)
{
/* Function to read a matrix from file filename into a matrix A. */
/* The matrix A must already have been precreated with the desired (or  */
/* default) partitioning across processors */
/* IMPORTANT: This routine preallocates memory for the matrix before  */
/* reading from file. For some unknown reason, too many calls to this routine  */
/* with the same matrix argument results in a crash. Hence, use this routine  */
/* ONLY if calling it ONCE. If you want to load a matrix multiple times  */
/* from file, use MatLoadIntoMatrix3. */

  PetscErrorCode ierr;
  PetscScalar    *aa;
  PetscInt       *jj;
  PetscInt       header[4],M,N,nnz;
  PetscInt       I,Istart,Iend,numLocRows,*numColsPerLocRow,maxNumCols,locNNZ;
  PetscInt       ip,numPrevCols,*allIstart,*allnumNNZ;
  PetscInt       fd1,fd2;
  PetscViewer    viewer1,viewer2;
  off_t          off,offset;
  PetscMPIInt    numProcessors,myId;  
      
  PetscFunctionBegin;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcessors);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myId);CHKERRQ(ierr);
  myId=myId+1;
  
  /* open the file twice so that later we can read entries from two different parts of the
     file at the same time. Note that due to file caching this should not impact performance */  
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,FILE_MODE_READ,&viewer1);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(viewer1,&fd1);CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,FILE_MODE_READ,&viewer2);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(viewer2,&fd2);CHKERRQ(ierr);

/* Read header */
  ierr = PetscBinaryRead(fd1,(char *)header,4,PETSC_INT);CHKERRQ(ierr);

  /* error checking on file */
  if (header[0] != MAT_FILE_COOKIE) SETERRQ(PETSC_ERR_FILE_UNEXPECTED,"not matrix object");

  M = header[1]; N = header[2]; nnz = header[3];

/* Find global indices of rows owned by this process */
  ierr = MatGetOwnershipRange(A,&Istart,&Iend); 
/* Note: Iend is one more than the global index of the last local row */  
  Iend = Iend-1;
  numLocRows = Iend-Istart+1; /* number of rows owned by this process */
  ierr = PetscMalloc(numLocRows*sizeof(PetscInt),&numColsPerLocRow);CHKERRQ(ierr);
 
/* Read in number of columns for each row owned by this process */
  off = PETSC_BINARY_INT_SIZE*(4+Istart); /* remember Istart and Iend are base 0 */ 
  ierr = PetscBinarySeek(fd1,off,PETSC_BINARY_SEEK_SET,&offset);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fd1,numColsPerLocRow,numLocRows,PETSC_INT);CHKERRQ(ierr);

/* Find nnz on this processor and max number of columns owned by any local row */
  maxNumCols = 0;
  locNNZ = 0;
  for (I=0; I<numLocRows; I++) {
    locNNZ = locNNZ + numColsPerLocRow[I];
    maxNumCols = MAX(maxNumCols,numColsPerLocRow[I]);
  }

/* Preallocate memory. This probably overestimates space requirements, but is  */
/* absolutely essential! */
/*   The first call is required to work when not using MPI or with only 1 processor */
  ierr = MatSeqAIJSetPreallocation(A,maxNumCols,PETSC_NULL);CHKERRQ(ierr);
  ierr = MatMPIAIJSetPreallocation(A,maxNumCols,PETSC_NULL,maxNumCols,PETSC_NULL);CHKERRQ(ierr);

/* Allocate space for column indices and values */
  ierr = PetscMalloc(maxNumCols*sizeof(PetscInt),&jj);CHKERRQ(ierr);
  ierr = PetscMalloc(maxNumCols*sizeof(PetscScalar),&aa);CHKERRQ(ierr);

/* Move file pointer to beginning of blocks owned by local rows   */
  ierr = PetscMalloc(numProcessors*sizeof(PetscInt),&allIstart);CHKERRQ(ierr);
  ierr = PetscMalloc(numProcessors*sizeof(PetscInt),&allnumNNZ);CHKERRQ(ierr);  
  ierr = MPI_Allgather(&Istart,1,MPIU_INT,allIstart,1,MPIU_INT,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Allgather(&locNNZ,1,MPIU_INT,allnumNNZ,1,MPIU_INT,PETSC_COMM_WORLD);CHKERRQ(ierr);

  numPrevCols = 0;
  for (ip=0; ip<numProcessors; ip++) {
    if (allIstart[ip]<Istart) numPrevCols = numPrevCols + allnumNNZ[ip];
  }

  off = PETSC_BINARY_INT_SIZE*(4+M+numPrevCols); /* remember Istart and Iend are base 0 */
  ierr = PetscBinarySeek(fd1,off,PETSC_BINARY_SEEK_SET,&offset);CHKERRQ(ierr);
  off = PETSC_BINARY_INT_SIZE*(4+M+nnz) + PETSC_BINARY_SCALAR_SIZE*numPrevCols; /* remember Istart and Iend are base 0 */
  ierr = PetscBinarySeek(fd2,off,PETSC_BINARY_SEEK_SET,&offset);CHKERRQ(ierr);

/* Read column indices and values a row at a time and insert into matrix   */
  for (I=Istart; I<=Iend; I++) {
/*     if (myId==2) { */
/*      printf("[%d] I=%d,numColsPerLocRow=%d\n",myId,I,numColsPerLocRow[I-Istart]); */
/*     } */
    if (numColsPerLocRow[I-Istart]>0) {
      ierr = PetscBinaryRead(fd1,jj,numColsPerLocRow[I-Istart],PETSC_INT);CHKERRQ(ierr);
      ierr = PetscBinaryRead(fd2,aa,numColsPerLocRow[I-Istart],PETSC_SCALAR);CHKERRQ(ierr);
      ierr = MatSetValues(A,1,&I,numColsPerLocRow[I-Istart],jj,aa,INSERT_VALUES);CHKERRQ(ierr);
    }
/*     if (myId==2) { */
/*      for (j=1; j<=numColsPerLocRow[I-Istart]; j++) { */
/*      printf("[%d] col=%d,val=%g\n",myId,jj[j-1],aa[j-1]); */
/*      } */
/*     } */    
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = PetscFree(numColsPerLocRow);CHKERRQ(ierr);
  ierr = PetscFree(allIstart);CHKERRQ(ierr);
  ierr = PetscFree(allnumNNZ);CHKERRQ(ierr);
  ierr = PetscFree(jj);CHKERRQ(ierr);
  ierr = PetscFree(aa);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer1);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer2);CHKERRQ(ierr);  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatLoadIntoMatrix3"
PetscErrorCode MatLoadIntoMatrix3(const char filename[], Mat A)
{
/* Function to read a matrix from file filename into a matrix A.  */
/* The matrix A must already have been precreated with the desired (or  */
/* default) partitioning across processors */
/* IMPORTANT: This routine does NOT preallocates memory for the matrix before  */
/* reading from file. You must do this before calling this routine (or performance  */
/* will be severely degraded). Performance may also be degraded if you try to load  */
/* matrix data with a different sparsity pattern. (It should still work, but is  */
/* probably not advisable.) */
/* Use this routine if you want to load a matrix  */
/* multiple times from file. If you only want to load a matrix once (or a  */
/* few times), use MatLoadIntoMatrix2 (which will preallocate memory). */

  PetscErrorCode ierr;
  PetscScalar    *aa;
  PetscInt       *jj;
  PetscInt       header[4],M,N,nnz;
  PetscInt       I,Istart,Iend,numLocRows,*numColsPerLocRow,maxNumCols,locNNZ;
  PetscInt       ip,numPrevCols,*allIstart,*allnumNNZ;
  PetscInt       fd1,fd2;
  PetscViewer    viewer1,viewer2;
  off_t          off,offset;
  PetscMPIInt    numProcessors,myId;  
      
  PetscFunctionBegin;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcessors);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myId);CHKERRQ(ierr);
  myId=myId+1;
  
  /* open the file twice so that later we can read entries from two different parts of the
     file at the same time. Note that due to file caching this should not impact performance */  
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,FILE_MODE_READ,&viewer1);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(viewer1,&fd1);CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,FILE_MODE_READ,&viewer2);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(viewer2,&fd2);CHKERRQ(ierr);

/* Read header */
  ierr = PetscBinaryRead(fd1,(char *)header,4,PETSC_INT);CHKERRQ(ierr);

  /* error checking on file */
  if (header[0] != MAT_FILE_COOKIE) SETERRQ(PETSC_ERR_FILE_UNEXPECTED,"not matrix object");

  M = header[1]; N = header[2]; nnz = header[3];

/* Find global indices of rows owned by this process */
  ierr = MatGetOwnershipRange(A,&Istart,&Iend); 
/* Note: Iend is one more than the global index of the last local row */  
  Iend = Iend-1;
  numLocRows = Iend-Istart+1; /* number of rows owned by this process */
  ierr = PetscMalloc(numLocRows*sizeof(PetscInt),&numColsPerLocRow);CHKERRQ(ierr);
 
/* Read in number of columns for each row owned by this process */
  off = PETSC_BINARY_INT_SIZE*(4+Istart); /* remember Istart and Iend are base 0 */ 
  ierr = PetscBinarySeek(fd1,off,PETSC_BINARY_SEEK_SET,&offset);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fd1,numColsPerLocRow,numLocRows,PETSC_INT);CHKERRQ(ierr);

/* Find nnz on this processor and max number of columns owned by any local row */
  maxNumCols = 0;
  locNNZ = 0;
  for (I=0; I<numLocRows; I++) {
    locNNZ = locNNZ + numColsPerLocRow[I];
    maxNumCols = MAX(maxNumCols,numColsPerLocRow[I]);
  }

/* Preallocate memory. This probably overestimates space requirements, but is  */
/* absolutely essential! */
/*   The first call is required to work when not using MPI or with only 1 processor */
/*   ierr = MatSeqAIJSetPreallocation(A,maxNumCols,PETSC_NULL);CHKERRQ(ierr); */
/*   ierr = MatMPIAIJSetPreallocation(A,maxNumCols,PETSC_NULL,maxNumCols,PETSC_NULL);CHKERRQ(ierr); */

/* Allocate space for column indices and values */
  ierr = PetscMalloc(maxNumCols*sizeof(PetscInt),&jj);CHKERRQ(ierr);
  ierr = PetscMalloc(maxNumCols*sizeof(PetscScalar),&aa);CHKERRQ(ierr);

/* Move file pointer to beginning of blocks owned by local rows   */
  ierr = PetscMalloc(numProcessors*sizeof(PetscInt),&allIstart);CHKERRQ(ierr);
  ierr = PetscMalloc(numProcessors*sizeof(PetscInt),&allnumNNZ);CHKERRQ(ierr);  
  ierr = MPI_Allgather(&Istart,1,MPIU_INT,allIstart,1,MPIU_INT,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Allgather(&locNNZ,1,MPIU_INT,allnumNNZ,1,MPIU_INT,PETSC_COMM_WORLD);CHKERRQ(ierr);

  numPrevCols = 0;
  for (ip=0; ip<numProcessors; ip++) {
    if (allIstart[ip]<Istart) numPrevCols = numPrevCols + allnumNNZ[ip];
  }

  off = PETSC_BINARY_INT_SIZE*(4+M+numPrevCols); /* remember Istart and Iend are base 0 */
  ierr = PetscBinarySeek(fd1,off,PETSC_BINARY_SEEK_SET,&offset);CHKERRQ(ierr);
  off = PETSC_BINARY_INT_SIZE*(4+M+nnz) + PETSC_BINARY_SCALAR_SIZE*numPrevCols; /* remember Istart and Iend are base 0 */
  ierr = PetscBinarySeek(fd2,off,PETSC_BINARY_SEEK_SET,&offset);CHKERRQ(ierr);

/* Read column indices and values a row at a time and insert into matrix   */
  for (I=Istart; I<=Iend; I++) {
/*     if (myId==2) { */
/*      printf("[%d] I=%d,numColsPerLocRow=%d\n",myId,I,numColsPerLocRow[I-Istart]); */
/*     } */
    if (numColsPerLocRow[I-Istart]>0) {
      ierr = PetscBinaryRead(fd1,jj,numColsPerLocRow[I-Istart],PETSC_INT);CHKERRQ(ierr);
      ierr = PetscBinaryRead(fd2,aa,numColsPerLocRow[I-Istart],PETSC_SCALAR);CHKERRQ(ierr);
      ierr = MatSetValues(A,1,&I,numColsPerLocRow[I-Istart],jj,aa,INSERT_VALUES);CHKERRQ(ierr);
    }
/*     if (myId==2) { */
/*      for (j=1; j<=numColsPerLocRow[I-Istart]; j++) { */
/*      printf("[%d] col=%d,val=%g\n",myId,jj[j-1],aa[j-1]); */
/*      } */
/*     } */    
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = PetscFree(numColsPerLocRow);CHKERRQ(ierr);
  ierr = PetscFree(allIstart);CHKERRQ(ierr);
  ierr = PetscFree(allnumNNZ);CHKERRQ(ierr);
  ierr = PetscFree(jj);CHKERRQ(ierr);
  ierr = PetscFree(aa);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer1);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer2);CHKERRQ(ierr);  
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "MatLoadIntoMatrix4"
PetscErrorCode MatLoadIntoMatrix4(const char filename[], Mat A, MatLayout *user)
{
/* Function to read a matrix from file filename into a matrix A.  */
/* The matrix A must already have been precreated with the desired (or  */
/* default) partitioning across processors */
/* This version of MatLoadIntoMatrix* stores matrix layout data in a structure  */
/* 'user' of type MatLayout (defined in petsc_matvec_utils.h). A pointer to this  */
/* structure must be passed as the 3d argument to MatLoadIntoMatrix4. You MUST  */
/* initialize this structure as follows before calling MatLoadIntoMatrix4 for  */
/* the first time (for a given matrix).  */
/*     user.firstTime = PETSC_TRUE; */
/*     user.lastTime = PETSC_FALSE; */
/* The first time you call this S/R, it will figure out the layout of the matrix  */
/* and return the data in the structure 'user'. Subsequent calls will only load  */
/* data from file using this layout information. Consequently, you must ONLY  */
/* call this S/R if you are certain that every file has matrix data with the  */
/* exact same sparsity pattern as that used for the initialization.  */
/* IMPORTANT: This routine does NOT preallocates memory for the matrix before  */
/* reading from file. You must do this before calling this routine (or performance  */
/* will be severely degraded). Use this routine if you want to load a matrix  */
/* multiple times from file, and want something more efficient than  */
/* MatLoadIntoMatrix3. (If you are not sure whether each file has matrix data  */
/* with the same sparsity pattern, use MatLoadIntoMatrix3.) */
/* If you only want to load a matrix once (or a few times), use MatLoadIntoMatrix2 (which  */
/* will preallocate memory). */

  PetscErrorCode ierr;
  PetscInt       header[4],M,N,nnz;
  PetscInt       I,Istart,Iend,numLocRows,maxNumCols,locNNZ;
  PetscInt       ip,numPrevCols;
  PetscMPIInt    numProcessors,myId;
  PetscInt       fd1,fd2;
  PetscViewer    viewer1,viewer2;
  off_t          off,offset;
      
  PetscFunctionBegin;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcessors);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myId);CHKERRQ(ierr);
  myId=myId+1;
  
/* open the file twice so that later we can read entries from two different parts of the
   file at the same time. Note that due to file caching this should not impact performance */  
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,FILE_MODE_READ,&viewer1);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(viewer1,&fd1);CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,filename,FILE_MODE_READ,&viewer2);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(viewer2,&fd2);CHKERRQ(ierr);

/* Read header */
  ierr = PetscBinaryRead(fd1,(char *)header,4,PETSC_INT);CHKERRQ(ierr);

  /* error checking on file */
  if (header[0] != MAT_FILE_COOKIE) SETERRQ(PETSC_ERR_FILE_UNEXPECTED,"not matrix object");

  M = header[1]; N = header[2]; nnz = header[3];

  if (user->firstTime) { /* Initialize matrix layout data */
/*  Find global indices of rows owned by this process */
	ierr = MatGetOwnershipRange(A,&user->Istart,&user->Iend); 
/*  Note: Iend is one more than the global index of the last local row */  
	user->Iend = (user->Iend)-1;
	numLocRows = (user->Iend)-(user->Istart)+1; /* number of rows owned by this process */
	ierr = PetscMalloc(numLocRows*sizeof(PetscInt),&user->numColsPerLocRow);CHKERRQ(ierr);
 
/*  Read in number of columns for each row owned by this process */
	off = PETSC_BINARY_INT_SIZE*(4+(user->Istart)); /* remember Istart and Iend are base 0 */ 
	ierr = PetscBinarySeek(fd1,off,PETSC_BINARY_SEEK_SET,&offset);CHKERRQ(ierr);
	ierr = PetscBinaryRead(fd1,user->numColsPerLocRow,numLocRows,PETSC_INT);CHKERRQ(ierr);

/*  Find nnz on this processor and max number of columns owned by any local row */
	maxNumCols = 0;
	locNNZ = 0;
	for (I=0; I<numLocRows; I++) {
	  locNNZ = locNNZ + (user->numColsPerLocRow[I]);
	  maxNumCols = MAX(maxNumCols,(user->numColsPerLocRow[I]));
	}

/* Preallocate memory. This probably overestimates space requirements, but is  */
/* absolutely essential! */
/*   The first call is required to work when not using MPI or with only 1 processor */
/*   ierr = MatSeqAIJSetPreallocation(A,maxNumCols,PETSC_NULL);CHKERRQ(ierr); */
/*   ierr = MatMPIAIJSetPreallocation(A,maxNumCols,PETSC_NULL,maxNumCols,PETSC_NULL);CHKERRQ(ierr); */

/*  Allocate space for column indices and values */
	ierr = PetscMalloc(maxNumCols*sizeof(PetscInt),&user->jj);CHKERRQ(ierr);
	ierr = PetscMalloc(maxNumCols*sizeof(PetscScalar),&user->aa);CHKERRQ(ierr);

	ierr = PetscMalloc(numProcessors*sizeof(PetscInt),&user->allIstart);CHKERRQ(ierr);
	ierr = PetscMalloc(numProcessors*sizeof(PetscInt),&user->allnumNNZ);CHKERRQ(ierr);  
	ierr = MPI_Allgather(&user->Istart,1,MPIU_INT,user->allIstart,1,MPIU_INT,PETSC_COMM_WORLD);CHKERRQ(ierr);
	ierr = MPI_Allgather(&locNNZ,1,MPIU_INT,user->allnumNNZ,1,MPIU_INT,PETSC_COMM_WORLD);CHKERRQ(ierr);
  
	numPrevCols = 0;
	for (ip=0; ip<numProcessors; ip++) {
	  if (user->allIstart[ip]<user->Istart) numPrevCols = numPrevCols + (user->allnumNNZ[ip]);
	}

	user->off1 = PETSC_BINARY_INT_SIZE*(4+M+numPrevCols); /* remember Istart and Iend are base 0 */
	user->off2 = PETSC_BINARY_INT_SIZE*(4+M+nnz) + PETSC_BINARY_SCALAR_SIZE*numPrevCols; /* remember Istart and Iend are base 0 */
	
	user->firstTime = PETSC_FALSE;
  }  /* end of initialization */

  Istart = user->Istart;
  Iend = user->Iend;

/*  Move file pointer to beginning of blocks owned by local rows   */
  ierr = PetscBinarySeek(fd1,user->off1,PETSC_BINARY_SEEK_SET,&offset);CHKERRQ(ierr);
  ierr = PetscBinarySeek(fd2,user->off2,PETSC_BINARY_SEEK_SET,&offset);CHKERRQ(ierr);

/* Read column indices and values a row at a time and insert into matrix   */
  for (I=Istart; I<=Iend; I++) {
    if (user->numColsPerLocRow[I-Istart]>0) {
      ierr = PetscBinaryRead(fd1,user->jj,user->numColsPerLocRow[I-Istart],PETSC_INT);CHKERRQ(ierr);
      ierr = PetscBinaryRead(fd2,user->aa,user->numColsPerLocRow[I-Istart],PETSC_SCALAR);CHKERRQ(ierr);
      ierr = MatSetValues(A,1,&I,user->numColsPerLocRow[I-Istart],user->jj,user->aa,INSERT_VALUES);CHKERRQ(ierr);
    }
  }
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = PetscViewerDestroy(viewer1);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(viewer2);CHKERRQ(ierr);  

  if (user->lastTime) { /* destroy layout data for this matrix */
	ierr = PetscFree(user->numColsPerLocRow);CHKERRQ(ierr);
	ierr = PetscFree(user->allIstart);CHKERRQ(ierr);
	ierr = PetscFree(user->allnumNNZ);CHKERRQ(ierr);
	ierr = PetscFree(user->jj);CHKERRQ(ierr);
	ierr = PetscFree(user->aa);CHKERRQ(ierr);
  }
  
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
  ierr = PetscViewerDestroy(viewer);CHKERRQ(ierr);

  /* error checking on file */
  if (header[0] != MAT_FILE_COOKIE) SETERRQ(PETSC_ERR_FILE_UNEXPECTED,"not matrix object");

  if (M) *M = header[1];
  if (N) *N = header[2];
  if (nnz) *nnz = header[3];

  PetscFunctionReturn(0);
}