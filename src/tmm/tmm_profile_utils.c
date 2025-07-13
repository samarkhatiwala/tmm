#define DEFINE_PROFILES

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "tmm_forcing_utils.h"
#include "tmm_share.h"
#include "tmm_profile_utils.h"

#undef __FUNCT__
#define __FUNCT__ "iniProfileData"
PetscErrorCode iniProfileData(PetscInt myId)
{
  PetscMPIInt numProcessors;
  PetscErrorCode ierr;
  PetscInt ipro, ip;
  PetscViewer fd;
  int fp;
  PetscInt dum;
  PetscBool useProfileNumberPartitioning = PETSC_FALSE;

  useProfiles = PETSC_FALSE;
  ierr = PetscOptionsHasName(NULL,NULL,"-use_profiles",&useProfiles);CHKERRQ(ierr);

  if (useProfiles) {
  
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcessors);CHKERRQ(ierr);
  
/*  Read in total number of profiles */
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"gStartIndices.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscBinaryRead(fp,&totalNumProfiles,1,NULL,PETSC_INT);CHKERRQ(ierr);
	if (totalNumProfiles<=0) SETERRQ(PETSC_COMM_WORLD,1,"Invalid total number of profiles! Must be >0");

/*  Read in starting and ending global indices of profiles. NOTE: these have a base 1 index. */
	ierr = PetscMalloc(totalNumProfiles*sizeof(PetscInt),&gStartIndices);CHKERRQ(ierr);
	ierr = PetscBinaryRead(fp,gStartIndices,totalNumProfiles,NULL,PETSC_INT);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Done reading gStartIndices.bin\n");CHKERRQ(ierr);

	ierr = PetscMalloc(totalNumProfiles*sizeof(PetscInt),&gEndIndices);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"gEndIndices.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscBinaryRead(fp,&dum,1,NULL,PETSC_INT);CHKERRQ(ierr);
	if (dum != totalNumProfiles) SETERRQ(PETSC_COMM_WORLD,1,"Total number of profiles don't match!");
	ierr = PetscBinaryRead(fp,gEndIndices,totalNumProfiles,NULL,PETSC_INT);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Done reading gEndIndices.bin\n");CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Total number of profiles specified: %d\n",totalNumProfiles);CHKERRQ(ierr);

/*  Figure out optimum partitioning of profiles over processors */
	ierr = PetscMalloc(numProcessors*sizeof(PetscInt),&gNumProfiles);CHKERRQ(ierr);
	ierr = PetscMalloc(numProcessors*sizeof(PetscInt),&gSizes);CHKERRQ(ierr);
	ierr = PetscMalloc(totalNumProfiles*sizeof(PetscInt),&gProfileLengths);CHKERRQ(ierr);
	for (ip=1; ip<=totalNumProfiles; ip++) {
	  gProfileLengths[ip-1]=gEndIndices[ip-1]-gStartIndices[ip-1]+1;
	}

	for (ipro=0; ipro<numProcessors; ipro++) {
	  gNumProfiles[ipro] = 0; /* initialize all with 0 */
	  gSizes[ipro] = 0; /* initialize all with 0 */
	}

    ierr = PetscOptionsHasName(NULL,NULL,"-partition_by_number_of_profiles",&useProfileNumberPartitioning);CHKERRQ(ierr);
    if (useProfileNumberPartitioning) {
  /* Let PETSc do this for us but this allocates the same number of profiles (rather than boxes) per process */
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Partitioning by number of profiles\n");CHKERRQ(ierr);
	  Vec templateVec;
	  ierr = VecCreate(PETSC_COMM_WORLD,&templateVec);CHKERRQ(ierr);
	  ierr = VecSetSizes(templateVec,PETSC_DECIDE,totalNumProfiles);CHKERRQ(ierr);
	  ierr = VecSetFromOptions(templateVec);CHKERRQ(ierr);
	  ierr = VecGetLocalSize(templateVec,&lNumProfiles);CHKERRQ(ierr);
	  ierr = VecDestroy(&templateVec);CHKERRQ(ierr);

  /* NOTE: myId starts at 1 */  
	  gNumProfiles[myId-1]=lNumProfiles;
	  MPI_Allgather(&lNumProfiles, 1, MPI_INT, gNumProfiles, 1, MPI_INT, PETSC_COMM_WORLD);

    } else {
  /* Alternatively, try to allocate the same number of boxes per process */
  /*  Code from Iris Kriest */
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Partitioning by number of boxes\n");CHKERRQ(ierr);  
	  PetscInt totalNumBoxes;
	  PetscScalar avgNumBoxesPerProcessor;
	  totalNumBoxes = gEndIndices[totalNumProfiles-1];
	  avgNumBoxesPerProcessor = totalNumBoxes/(PetscScalar)numProcessors; /* also the ideal number of boxes per processor, if all profiles were of equal length */

	  for (ip=0; ip<totalNumProfiles; ip++) {
		ipro = (PetscInt)floor((gStartIndices[ip]-1+0.5*gProfileLengths[ip])/avgNumBoxesPerProcessor); /* determine process Id */
		gNumProfiles[ipro] =  gNumProfiles[ipro] + 1; /* add one profile to this processor */
		gSizes[ipro] = gSizes[ipro] + gProfileLengths[ip];
	  }

  /*  Check */
	  PetscInt tot=0;
	  PetscInt totp=0;
	  tot=0;
	  totp=0; 
	  for (ipro=1; ipro<=numProcessors; ipro++) {
		tot=tot+gSizes[ipro-1];
		totp=totp+gNumProfiles[ipro-1];
	  }
 
	  if (tot != totalNumBoxes) {
		ierr=PetscPrintf(PETSC_COMM_WORLD,"Number of boxes calculated during load balancing (%d) not equal to actual number of boxes (%d) \n",tot,totalNumBoxes);CHKERRQ(ierr);   
		SETERRQ(PETSC_COMM_WORLD,1,"Problem distributing profiles across processors!");
	  }

	  if (totp != totalNumProfiles) {
		ierr=PetscPrintf(PETSC_COMM_WORLD,"Number of profiles calculated during load balancing (%d) not equal to actual number of profiles (%d) \n",totp,totalNumProfiles);CHKERRQ(ierr);   
		SETERRQ(PETSC_COMM_WORLD,1,"Problem distributing profiles across processors!");
	  }
 
	  lNumProfiles=gNumProfiles[myId-1];
    } 

/* Compute total number of profiles upto (but not including) current processor */
/* NOTE: myId starts at 1 */
	numPrevProfiles=0;
	for (ipro=1; ipro<=myId-1; ipro++) {
	  numPrevProfiles = numPrevProfiles + gNumProfiles[ipro-1];
	}
  
/*  Compute starting and ending LOCAL indices of profiles on current processor. NOTE: These have a base 0 index. */
/*  ierr=PetscPrintf(PETSC_COMM_WORLD,"lNumProfiles = %d\n",lNumProfiles);CHKERRQ(ierr);   */
	ierr = PetscMalloc(lNumProfiles*sizeof(PetscInt),&lStartIndices);CHKERRQ(ierr);
	ierr = PetscMalloc(lNumProfiles*sizeof(PetscInt),&lEndIndices);CHKERRQ(ierr);
	ierr = PetscMalloc(lNumProfiles*sizeof(PetscInt),&lProfileLength);CHKERRQ(ierr);
	lSize=0; /* local size of vectors */
	for (ip=1; ip<=lNumProfiles; ip++) {
	  lStartIndices[ip-1]=gStartIndices[numPrevProfiles+ip-1]-gStartIndices[numPrevProfiles+1-1];
	  lEndIndices[ip-1]=gEndIndices[numPrevProfiles+ip-1]-gStartIndices[numPrevProfiles+1-1];
	  lProfileLength[ip-1]=lEndIndices[ip-1]-lStartIndices[ip-1]+1;
	  lSize=lSize+lProfileLength[ip-1];
	}

    if (useProfileNumberPartitioning) {
	  MPI_Allgather(&lSize, 1, MPI_INT, gSizes, 1, MPI_INT, PETSC_COMM_WORLD);
    }
    
/*  Check */
/*  NOTE: myId starts at 1 */  
	if (lSize != gSizes[myId-1]) SETERRQ(PETSC_COMM_WORLD,1,"Problem distributing profiles across processors!");

	for (ipro=1; ipro<=numProcessors; ipro++) {
	  ierr=PetscPrintf(PETSC_COMM_WORLD,"Number of profiles (total length) on processor %d = %d (%d)\n",ipro-1,gNumProfiles[ipro-1],gSizes[ipro-1]);CHKERRQ(ierr);
	}  
  
  } /* useProfiles */
  
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "readProfileSurfaceIntData"
PetscErrorCode readProfileSurfaceIntData(const char *fileName, PetscInt *arr, PetscInt numValsPerProfile)
{
  PetscErrorCode ierr;
/*   PetscInt *tmpArr; */
/*   PetscInt ip; */
/*   size_t m1, m2; */
  off_t  off, offset;  
  PetscViewer fd;
  int fp;
  off_t iShift;

/*   m1 = totalNumProfiles*sizeof(PetscInt); */
/*   m2 = lNumProfiles*sizeof(PetscInt); */

/*   ierr = PetscMalloc(m1,&tmpArr);CHKERRQ(ierr); */
/*   ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,fileName,FILE_MODE_READ,&fd);CHKERRQ(ierr); */
/*   ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr); */
/*   ierr = PetscBinaryRead(fp,tmpArr,totalNumProfiles,NULL,PETSC_INT);CHKERRQ(ierr); */
/*   ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr); */

  if (!useProfiles) SETERRQ(PETSC_COMM_WORLD,1,"You must switch on profiles with the -use_profiles option to use the *ProfileSurfaceScalar* routines!");
  
/* Shift file pointer to start of data owned by local process */
  iShift = (off_t)numValsPerProfile*(off_t)numPrevProfiles;
  off = PETSC_BINARY_INT_SIZE*iShift;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,fileName,FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
  ierr = PetscBinarySeek(fp,off,PETSC_BINARY_SEEK_SET,&offset);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fp,arr,numValsPerProfile*lNumProfiles,NULL,PETSC_INT);CHKERRQ(ierr);
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
PetscErrorCode readProfileSurfaceScalarData(const char *fileName, PetscScalar *arr, PetscInt numValsPerProfile)
{
  PetscErrorCode ierr;
/*   PetscScalar *tmpArr; */
/*   PetscInt ip; */
/*   size_t m1, m2; */
  off_t  off, offset;  
  PetscViewer fd;
  int fp;
  off_t iShift;
  PetscMPIInt numProcessors, myId;

  if (!useProfiles) SETERRQ(PETSC_COMM_WORLD,1,"You must switch on profiles with the -use_profiles option to use the *ProfileSurfaceScalar* routines!");

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myId);CHKERRQ(ierr);  
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcessors);CHKERRQ(ierr);
  
/* Shift file pointer to start of data owned by local process */
  iShift = (off_t)numValsPerProfile*(off_t)numPrevProfiles;
/*   printf("ipro=%d,iShift=%d\n",myId,iShift); */
  off = PETSC_BINARY_SCALAR_SIZE*iShift;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,fileName,FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
  ierr = PetscBinarySeek(fp,off,PETSC_BINARY_SEEK_SET,&offset);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fp,arr,numValsPerProfile*lNumProfiles,NULL,PETSC_SCALAR);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "readProfileSurfaceScalarDataRecord"
PetscErrorCode readProfileSurfaceScalarDataRecord(const char *fileName, PetscScalar *arr, PetscInt numValsPerProfile, PetscInt iRec)
{
/* Random access version of readProfileSurfaceScalarData */
/* This version takes 1 additional argument:  */
/*   iRec: the record to read (iRec=1 is the first record) */
  PetscErrorCode ierr;
  off_t  off, offset;  
  PetscViewer fd;
  int fp;
  off_t iShift;
  PetscMPIInt numProcessors, myId;

  if (!useProfiles) SETERRQ(PETSC_COMM_WORLD,1,"You must switch on profiles with the -use_profiles option to use the *ProfileSurfaceScalar* routines!");

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myId);CHKERRQ(ierr);  
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcessors);CHKERRQ(ierr);
  

/* Shift file pointer to start of data owned by local process */
  iShift = ((off_t)iRec-1)*(off_t)numValsPerProfile*(off_t)totalNumProfiles + (off_t)numValsPerProfile*(off_t)numPrevProfiles;
  off = PETSC_BINARY_SCALAR_SIZE*iShift;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,fileName,FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
  ierr = PetscBinarySeek(fp,off,PETSC_BINARY_SEEK_SET,&offset);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fp,arr,numValsPerProfile*lNumProfiles,NULL,PETSC_SCALAR);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "interpPeriodicProfileSurfaceScalarData"
PetscErrorCode interpPeriodicProfileSurfaceScalarData(PetscScalar tc, PetscScalar *uarr, PetscScalar cyclePeriod,
                                    PetscInt numPerPeriod, PetscScalar *tdp, 
                                    PeriodicArray user, const char *fileName)
{
/* Function to interpolate an array that is periodic in time with period cyclePeriod.  */
/* tc is the current time and numPerPeriod is the number of instances per period   */
/* at which data are available (to be read from files). */

  if (!useProfiles) SETERRQ(PETSC_COMM_WORLD,1,"You must switch on profiles with the -use_profiles option to use the *ProfileSurfaceScalar* routines!");

  PetscScalar t,t1;
  PetscInt im,it0,it1;
/*   static PetscInt iCurrTimeReadLast=-1; */
  PetscErrorCode ierr;
  PetscScalar alpha[2];  
  char tmpFile[PETSC_MAX_PATH_LEN];
  PetscInt ip;

  if (user->firstTime) {
//     if (numPerPeriod>MAX_FORCING_NUM_PER_PERIOD) {
//       SETERRQ(PETSC_COMM_WORLD,1,"Number of allowable arrays in PeriodicArray struct exceeded by requested number ! Increase MAX_FORCING_NUM_PER_PERIOD.");
//     }
    if (((user->arrayLength) % lNumProfiles) != 0) {
      SETERRQ(PETSC_COMM_WORLD,1,"arrayLength for PeriodicArray is not divisible by lNumProfiles!");
    }
    user->numValsPerProfile = (user->arrayLength)/lNumProfiles;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Initializing PeriodicArray object %s with %d value(s) per profile\n",fileName,user->numValsPerProfile);CHKERRQ(ierr);    
    user->numPerPeriod = numPerPeriod;  
    user->qp = (PetscScalar **) malloc(numPerPeriod*sizeof(PetscScalar *));
    for (im=0; im<numPerPeriod; im++) {
      ierr = PetscMalloc((user->arrayLength)*sizeof(PetscScalar),&user->qp[im]);CHKERRQ(ierr);    
      strcpy(tmpFile,"");
      sprintf(tmpFile,"%s%02d",fileName,im);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading data from file %s\n", tmpFile);CHKERRQ(ierr);  
	  ierr = readProfileSurfaceScalarData(tmpFile,user->qp[im],user->numValsPerProfile);    
    }
    user->firstTime = PETSC_FALSE;  
  }

  t=tc; /* current time */
  if (t<0.) t=cyclePeriod+t;
  t1=t-cyclePeriod*floor(t/cyclePeriod);
  ierr=calcPeriodicInterpFactor(numPerPeriod,t1,tdp,&it0,&it1,&alpha[0],&alpha[1]);  CHKERRQ(ierr);  
/*   ierr = PetscPrintf(PETSC_COMM_WORLD,"tc=%lf,t1=%lf,it0=%d,it1=%d,a1=%17.16lf,a2=%17.16lf\n",tc,t1,it0,it1,alpha[0],alpha[1]);CHKERRQ(ierr);   */
  
/* interpolate to current time   */
  for (ip=0;ip<(user->arrayLength);ip++) {
    uarr[ip]=alpha[0]*user->qp[it0][ip]+alpha[1]*user->qp[it1][ip];
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "interpTimeDependentProfileSurfaceScalarData"
PetscErrorCode interpTimeDependentProfileSurfaceScalarData(PetscScalar tc, PetscScalar *uarr, PetscInt numTimes, PetscScalar *tdt,
                                    TimeDependentArray user, const char *fileName)
{
/* Function to interpolate a time-dependent array.  */
/* tc is the current time and numTimes are the number of time slices (in tdt) at */
/* at which data are available. */

  PetscErrorCode ierr;
  PetscInt itc;  
  PetscScalar alpha;
  PetscInt ip;

  if (!useProfiles) SETERRQ(PETSC_COMM_WORLD,1,"You must switch on profiles with the -use_profiles option to use the *ProfileSurfaceScalar* routines!");
  
  if (user->firstTime) {
    if (((user->arrayLength) % lNumProfiles) != 0) {
      SETERRQ(PETSC_COMM_WORLD,1,"arrayLength for TimeDependentArray is not divisible by lNumProfiles!");
    }
    user->numValsPerProfile = (user->arrayLength)/lNumProfiles;
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Initializing TimeDependentArray object %s with %d value(s) per profile\n",fileName,user->numValsPerProfile);CHKERRQ(ierr);    
	ierr = PetscMalloc((user->arrayLength)*sizeof(PetscScalar),&user->utd[0]);CHKERRQ(ierr);
	ierr = PetscMalloc((user->arrayLength)*sizeof(PetscScalar),&user->utd[1]);CHKERRQ(ierr);
	user->itcurr=-1;
    user->firstTime = PETSC_FALSE;
  }

  if ((tc<tdt[0]) || (tc>tdt[numTimes-1])) {
    SETERRQ(PETSC_COMM_WORLD,1,"Error in interpTimeDependentProfileSurfaceScalarData: time out of bound");
  }

  ierr = calcInterpFactor(numTimes,tc,tdt,&itc,&alpha); CHKERRQ(ierr);
  if (itc != user->itcurr) { /* time to read new bracketing slices: itc uses 0-based, while readProfileSurfaceScalarDataRecord uses 1-based indexing*/
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading new bracketing slices for array %s at time = %g: %d and %d\n",fileName,tc,itc+1,itc+2);CHKERRQ(ierr);
	ierr = readProfileSurfaceScalarDataRecord(fileName,user->utd[0],user->numValsPerProfile,itc+1);
	ierr = readProfileSurfaceScalarDataRecord(fileName,user->utd[1],user->numValsPerProfile,itc+2);	
	user->itcurr=itc;
  }
/* interpolate to current time   */
  for (ip=0;ip<(user->arrayLength);ip++) {
	uarr[ip] = alpha*user->utd[0][ip] + (1.0-alpha)*user->utd[1][ip];
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "writeProfileSurfaceScalarData"
PetscErrorCode writeProfileSurfaceScalarData(const char *fileName, PetscScalar *arr, PetscInt numValsPerProfile, PetscBool appendToFile)
{
  PetscErrorCode ierr;
  PetscScalar *tmpArr;
  PetscInt *displs, *rcounts, cumpro;
  PetscInt ipro;
  size_t m1, m2;
/*   off_t  off, offset;   */
  PetscViewer fd;
  int fp;
/*   PetscInt iShift; */
  PetscMPIInt numProcessors, myId;

  if (!useProfiles) SETERRQ(PETSC_COMM_WORLD,1,"You must switch on profiles with the -use_profiles option to use the *ProfileSurfaceScalar* routines!");

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
	ierr = PetscBinaryWrite(fp,tmpArr,numValsPerProfile*totalNumProfiles,PETSC_SCALAR);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  }
  
  ierr = PetscFree(tmpArr);CHKERRQ(ierr);
  ierr = PetscFree(displs);CHKERRQ(ierr);
  ierr = PetscFree(rcounts);CHKERRQ(ierr);
    
  return 0;
}

PetscErrorCode dotProdProfileSurfaceScalarData(PetscScalar *xarr, PetscScalar *yarr, PetscScalar *z)
{
/* Function to take a dot product of two profile surface arrays. */
/* I think this only works if numValsPerProfile=1 */

  PetscErrorCode ierr;
  PetscInt ip;
  PetscScalar localz = 0.0;

  if (!useProfiles) SETERRQ(PETSC_COMM_WORLD,1,"You must switch on profiles with the -use_profiles option to use the *ProfileSurfaceScalar* routines!");
  
  *z = 0.0;
  for (ip=0; ip<lNumProfiles; ip++) {
	   localz = localz + xarr[ip]*yarr[ip];
  }
  MPI_Allreduce(&localz, z, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

  return 0;
}
