#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "tmm_forcing_utils.h"
PetscInt *gNumProfiles, *gStartIndices, *gEndIndices, *lStartIndices, *lEndIndices;
PetscInt *lProfileLength, lNumProfiles, lSize, numPrevProfiles, totalNumProfiles;
PetscBool useProfiles;
#include "tmm_profile_utils.h"

#undef __FUNCT__
#define __FUNCT__ "iniProfileData"
PetscErrorCode iniProfileData(PetscInt myId)
{
  PetscMPIInt numProcessors;
  PetscErrorCode ierr;
  PetscInt ipro, ip;
  PetscViewer fd;
  PetscInt fp;
  PetscInt dum;
  Vec templateVec;

  useProfiles = PETSC_FALSE;
  ierr = PetscOptionsHasName(PETSC_NULL,"-use_profiles",&useProfiles);CHKERRQ(ierr);

  if (useProfiles) {
  
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcessors);CHKERRQ(ierr);
  
/*   Read in total number of profiles */
  ierr = PetscMalloc(totalNumProfiles*sizeof(PetscInt),&gStartIndices);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"gStartIndices.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fp,&totalNumProfiles,1,PETSC_INT);CHKERRQ(ierr);
  if (totalNumProfiles<=0) SETERRQ(PETSC_COMM_WORLD,1,"Invalid total number of profiles! Must be >0");

/*   Read in starting and ending global indices of profiles. NOTE: these have a base 1 index. */
  ierr = PetscMalloc(totalNumProfiles*sizeof(PetscInt),&gStartIndices);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fp,gStartIndices,totalNumProfiles,PETSC_INT);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Done reading gStartIndices.bin\n");CHKERRQ(ierr);

  ierr = PetscMalloc(totalNumProfiles*sizeof(PetscInt),&gEndIndices);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"gEndIndices.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fp,&dum,1,PETSC_INT);CHKERRQ(ierr);
  if (dum != totalNumProfiles) SETERRQ(PETSC_COMM_WORLD,1,"Total number of profiles don't match!");
  ierr = PetscBinaryRead(fp,gEndIndices,totalNumProfiles,PETSC_INT);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Done reading gEndIndices.bin\n");CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Total number of profiles specified: %d\n",totalNumProfiles);CHKERRQ(ierr);

/* Figure out optimum partitioning of profiles over processors */
  ierr = PetscMalloc(numProcessors*sizeof(PetscInt),&gNumProfiles);CHKERRQ(ierr);
/* Let PETSc do this for us */
  ierr = VecCreate(PETSC_COMM_WORLD,&templateVec);CHKERRQ(ierr);
  ierr = VecSetSizes(templateVec,PETSC_DECIDE,totalNumProfiles);CHKERRQ(ierr);
  ierr = VecSetFromOptions(templateVec);CHKERRQ(ierr);
  ierr = VecGetLocalSize(templateVec,&lNumProfiles);CHKERRQ(ierr);
  ierr = VecDestroy(&templateVec);CHKERRQ(ierr);

/* NOTE: myId starts at 1 */  
  gNumProfiles[myId-1]=lNumProfiles;
  MPI_Allgather(&lNumProfiles, 1, MPI_INT, gNumProfiles, 1, MPI_INT, PETSC_COMM_WORLD);

  for (ipro=1; ipro<=numProcessors; ipro++) {
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Number of profiles on processor %d = %d\n",ipro-1,gNumProfiles[ipro-1]);CHKERRQ(ierr);
  }
  
/* Compute total number of profiles upto (but not including) current processor */
/* NOTE: myId starts at 1 */
  numPrevProfiles=0;
  for (ipro=1; ipro<=myId-1; ipro++) {
    numPrevProfiles = numPrevProfiles + gNumProfiles[ipro-1];
  }
  
/*   Compute starting and ending LOCAL indices of profiles on current processor. NOTE: These have a base 0 index. */
  lNumProfiles = gNumProfiles[myId-1];
/*   ierr=PetscPrintf(PETSC_COMM_WORLD,"lNumProfiles = %d\n",lNumProfiles);CHKERRQ(ierr);   */
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscInt),&lStartIndices);CHKERRQ(ierr);
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscInt),&lEndIndices);CHKERRQ(ierr);
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscInt),&lProfileLength);CHKERRQ(ierr);
  lSize=0; /* local size of vectors */
  for (ip=1; ip<=lNumProfiles; ip++) {
    lStartIndices[ip-1]=gStartIndices[numPrevProfiles+ip-1]-gStartIndices[numPrevProfiles+1-1];
    lEndIndices[ip-1]=gEndIndices[numPrevProfiles+ip-1]-gStartIndices[numPrevProfiles+1-1];
    lProfileLength[ip-1]=lEndIndices[ip-1]-lStartIndices[ip-1]+1;
/*     ierr=PetscPrintf(PETSC_COMM_WORLD,"lStartIndices=%d, lEndIndices=%d, lProfileLength=%d\n",lStartIndices[ip-1],lEndIndices[ip-1],lProfileLength[ip-1]);CHKERRQ(ierr);     */
    lSize=lSize+lProfileLength[ip-1];
  }
  
  } /* useProfiles */
  
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "readProfileSurfaceIntData"
PetscErrorCode readProfileSurfaceIntData(char *fileName, PetscInt *arr, PetscInt numValsPerProfile)
{
  PetscErrorCode ierr;
/*   PetscInt *tmpArr; */
/*   PetscInt ip; */
/*   size_t m1, m2; */
  off_t  off, offset;  
  PetscViewer fd;
  PetscInt fp;
  PetscInt iShift;

/*   m1 = totalNumProfiles*sizeof(PetscInt); */
/*   m2 = lNumProfiles*sizeof(PetscInt); */

/*   ierr = PetscMalloc(m1,&tmpArr);CHKERRQ(ierr); */
/*   ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,fileName,FILE_MODE_READ,&fd);CHKERRQ(ierr); */
/*   ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr); */
/*   ierr = PetscBinaryRead(fp,tmpArr,totalNumProfiles,PETSC_INT);CHKERRQ(ierr); */
/*   ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr); */

/* Shift file pointer to start of data owned by local process */
  iShift = numValsPerProfile*numPrevProfiles;
  off = PETSC_BINARY_INT_SIZE*iShift;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,fileName,FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
  ierr = PetscBinarySeek(fp,off,PETSC_BINARY_SEEK_SET,&offset);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fp,arr,numValsPerProfile*lNumProfiles,PETSC_INT);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

/*   for (ip=0; ip<totalNumProfiles; ip++) { */
/*     ierr = PetscPrintf(PETSC_COMM_WORLD,"ip=%d,kc=%d\n",ip,tmpArr[ip]);CHKERRQ(ierr); */
/*   }  */
/*   ierr = PetscMalloc(m2,&arr);CHKERRQ(ierr); */
/*   for (ip=1; ip<=lNumProfiles; ip++) {   */
/*     arr[ip-1]=tmpArr[numPrevProfiles+ip-1]; */
/*     ierr = PetscPrintf(PETSC_COMM_WORLD,"ip=%d,kc=%d\n",ip,arr[ip-1]);CHKERRQ(ierr);     */
/*   } */

/*   ierr = PetscFree(tmpArr);CHKERRQ(ierr); */
    
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "readProfileSurfaceScalarData"
PetscErrorCode readProfileSurfaceScalarData(char *fileName, PetscScalar *arr, PetscInt numValsPerProfile)
{
  PetscErrorCode ierr;
/*   PetscScalar *tmpArr; */
/*   PetscInt ip; */
/*   size_t m1, m2; */
  off_t  off, offset;  
  PetscViewer fd;
  PetscInt fp;
  PetscInt iShift;
  PetscMPIInt numProcessors, myId;

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myId);CHKERRQ(ierr);  
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcessors);CHKERRQ(ierr);
  
/*   m1 = numValsPerProfile*totalNumProfiles*sizeof(PetscScalar); */
/*   m2 = numValsPerProfile*lNumProfiles*sizeof(PetscScalar); */

/*Read all data into temporary array */
/*   ierr = PetscMalloc(m1,&tmpArr);CHKERRQ(ierr); */
/*   ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,fileName,FILE_MODE_READ,&fd);CHKERRQ(ierr); */
/*   ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr); */
/*   ierr = PetscBinaryRead(fp,tmpArr,numValsPerProfile*totalNumProfiles,PETSC_SCALAR);CHKERRQ(ierr); */
/*   ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr); */

/* Shift file pointer to start of data owned by local process */
  iShift = numValsPerProfile*numPrevProfiles;
/*   printf("ipro=%d,iShift=%d\n",myId,iShift); */
  off = PETSC_BINARY_SCALAR_SIZE*iShift;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,fileName,FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
  ierr = PetscBinarySeek(fp,off,PETSC_BINARY_SEEK_SET,&offset);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fp,arr,numValsPerProfile*lNumProfiles,PETSC_SCALAR);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  
/*   ierr = PetscMalloc(m2,&arr);CHKERRQ(ierr); */
/*   for (ip=1; ip<=lNumProfiles; ip++) {   */
/*     arr[ip-1]=tmpArr[numPrevProfiles+ip-1]; */
/*   } */

/*   ierr = PetscFree(tmpArr);CHKERRQ(ierr); */
    
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "readProfileSurfaceScalarDataRecord"
PetscErrorCode readProfileSurfaceScalarDataRecord(char *fileName, PetscScalar *arr, PetscInt numValsPerProfile, PetscInt iRec)
{
/* Random access version of readProfileSurfaceScalarData */
/* This version takes 1 additional argument:  */
/*   iRec: the record to read (iRec=1 is the first record) */
  PetscErrorCode ierr;
  off_t  off, offset;  
  PetscViewer fd;
  PetscInt fp;
  PetscInt iShift;
  PetscMPIInt numProcessors, myId;

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myId);CHKERRQ(ierr);  
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcessors);CHKERRQ(ierr);
  

/* Shift file pointer to start of data owned by local process */
  iShift = (iRec-1)*numValsPerProfile*totalNumProfiles + numValsPerProfile*numPrevProfiles;
  off = PETSC_BINARY_SCALAR_SIZE*iShift;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,fileName,FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
  ierr = PetscBinarySeek(fp,off,PETSC_BINARY_SEEK_SET,&offset);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fp,arr,numValsPerProfile*lNumProfiles,PETSC_SCALAR);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "readProfileSurfaceData"
/* PetscErrorCode readProfileSurfaceData(char *fileName, void *arr, PetscDataType type) */
/* { */
/*   PetscErrorCode ierr; */
/*   void *tmpArr; */
/*   PetscInt ip; */
/*   size_t m1, m2; */
/*   PetscViewer fd; */
/*   PetscInt fp; */
/*  */
/*   if (type == PETSC_INT) { */
/*     m1 = totalNumProfiles*sizeof(PetscInt); */
/*     m2 = lNumProfiles*sizeof(PetscInt); */
/*   } else if (type == PETSC_SCALAR) { */
/*     m1 = totalNumProfiles*sizeof(PetscScalar); */
/*     m2 = lNumProfiles*sizeof(PetscScalar);     */
/*   } else { */
/*     SETERRQ(PETSC_COMM_WORLD,PETSC_ERR_ARG_OUTOFRANGE,"Unknown type"); */
/*   } */
/*  */
/*   ierr = PetscMalloc(m1,&tmpArr);CHKERRQ(ierr); */
/*   ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,fileName,FILE_MODE_READ,&fd);CHKERRQ(ierr); */
/*   ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr); */
/*   ierr = PetscBinaryRead(fp,tmpArr,totalNumProfiles,type);CHKERRQ(ierr); */
/*   ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr); */
/*    */
/*   ierr = PetscMalloc(m2,&arr);CHKERRQ(ierr); */
/*   for (ip=1; ip<=lNumProfiles; ip++) {   */
/*     arr[ip-1]=tmpArr[numPrevProfiles+ip-1]; */
/*   } */

/*   ierr = PetscFree(tmpArr);CHKERRQ(ierr); */
/*      */
/*   return 0; */
/* } */

#undef __FUNCT__
#define __FUNCT__ "interpPeriodicProfileSurfaceScalarData"
PetscErrorCode interpPeriodicProfileSurfaceScalarData(PetscScalar tc, PetscScalar *uarr, PetscScalar cyclePeriod,
                                    PetscInt numPerPeriod, PetscScalar *tdp, 
                                    PeriodicArray *user, char *fileName)
{
/* Function to interpolate an array that is periodic in time with period cyclePeriod.  */
/* tc is the current time and numPerPeriod is the number of instances per period   */
/* at which data are available (to be read from files). */
/* IMPORTANT: Arrays u0 and u1 MUST have been created and preallocated before  */
/* calling this routine.  */

#include <math.h>
#include "petsc_matvec_utils.h"

  PetscScalar t,t1;
  PetscInt im,it0,it1;
/*   static PetscInt iCurrTimeReadLast=-1; */
  PetscErrorCode ierr;
  PetscScalar alpha[2];  
  char tmpFile[PETSC_MAX_PATH_LEN];
  PetscInt ip;

  if (user->firstTime) {
    if (numPerPeriod>MAX_FORCING_NUM_PER_PERIOD) {
      SETERRQ(PETSC_COMM_WORLD,1,"Number of allowable arrays in PeriodicArray struct exceeded by requested number ! Increase MAX_FORCING_NUM_PER_PERIOD.");
    }      
    user->numPerPeriod = numPerPeriod;  
    for (im=0; im<numPerPeriod; im++) {
      ierr = PetscMalloc((user->arrayLength)*sizeof(PetscScalar),&user->up[im]);CHKERRQ(ierr);    
      strcpy(tmpFile,"");
      sprintf(tmpFile,"%s%02d",fileName,im);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading data from file %s\n", tmpFile);CHKERRQ(ierr);  
	  ierr = readProfileSurfaceScalarData(tmpFile,user->up[im],1);    
    }
    user->firstTime = PETSC_FALSE;  
  }

  t=tc; /* current time */
  if (t<0.) t=cyclePeriod+t;
  t1=t-cyclePeriod*floor(t/cyclePeriod);
  ierr=calcPeriodicInterpFactor(numPerPeriod,t1,tdp,&it0,&it1,&alpha[0],&alpha[1]);  CHKERRQ(ierr);  
/*   ierr = PetscPrintf(PETSC_COMM_WORLD,"tc=%lf,t1=%lf,it0=%d,it1=%d,a1=%17.16lf,a2=%17.16lf\n",tc,t1,it0,it1,alpha[0],alpha[1]);CHKERRQ(ierr);   */
  
/* interpolate to current time   */
  for (ip=0;ip<lNumProfiles;ip++) {
    uarr[ip]=alpha[0]*user->up[it0][ip]+alpha[1]*user->up[it1][ip];
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "writeProfileSurfaceScalarData"
PetscErrorCode writeProfileSurfaceScalarData(char *fileName, PetscScalar *arr, PetscInt numValsPerProfile, PetscBool appendToFile)
{
  PetscErrorCode ierr;
  PetscScalar *tmpArr;
  PetscInt *displs, *rcounts, cumpro;
  PetscInt ipro;
  size_t m1, m2;
/*   off_t  off, offset;   */
  PetscViewer fd;
  PetscInt fp;
/*   PetscInt iShift; */
  PetscMPIInt numProcessors, myId;

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myId);CHKERRQ(ierr);  
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcessors);CHKERRQ(ierr);

  m1 = numValsPerProfile*totalNumProfiles*sizeof(PetscScalar);
  m2 = numProcessors*sizeof(PetscInt);  
/*   Allocate memory for temporary arrays */
  ierr = PetscMalloc(m1,&tmpArr);CHKERRQ(ierr);
  ierr = PetscMalloc(m2,&displs);CHKERRQ(ierr);
  ierr = PetscMalloc(m2,&rcounts);CHKERRQ(ierr);

  cumpro=0;
  for (ipro=1; ipro<=numProcessors; ipro++) {
    displs[ipro-1]=numValsPerProfile*cumpro;
    rcounts[ipro-1]=numValsPerProfile*gNumProfiles[ipro-1];
    cumpro = cumpro + gNumProfiles[ipro-1];
/*     ierr=PetscPrintf(PETSC_COMM_WORLD,"cumpro=%d, displs=%d\n",cumpro,displs[ipro-1],rcounts[ipro-1]);CHKERRQ(ierr);         */
  }
  
  MPI_Gatherv(arr,numValsPerProfile*lNumProfiles,MPI_DOUBLE,tmpArr,rcounts,displs,MPI_DOUBLE,0, PETSC_COMM_WORLD); 

  if (myId==0) { /* this shouldn't really be necessary, but without it, all processors seem to be writing in append mode */
	if (appendToFile) {
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,fileName,FILE_MODE_APPEND,&fd);CHKERRQ(ierr);
	} else {
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,fileName,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	}  
  
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscBinaryWrite(fp,tmpArr,numValsPerProfile*totalNumProfiles,PETSC_SCALAR,PETSC_TRUE);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  }
  
  ierr = PetscFree(tmpArr);CHKERRQ(ierr);
  ierr = PetscFree(displs);CHKERRQ(ierr);
  ierr = PetscFree(rcounts);CHKERRQ(ierr);
    
  return 0;
}
