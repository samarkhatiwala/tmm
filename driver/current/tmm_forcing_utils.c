#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "petsc_matvec_utils.h"
#include "tmm_forcing_utils.h"

#undef __FUNCT__
#define __FUNCT__ "calcInterpFactor"
PetscErrorCode calcInterpFactor( PetscInt n, PetscScalar t, PetscScalar tarr[], PetscInt *itf, PetscScalar *alpha)
{
   PetscInt it1, it2;
   
   it1=findindex(tarr,n,t);
   
   if (it1<0) SETERRQ(PETSC_COMM_WORLD,1,"Error in findindex: time out of bound");

   it1=MIN(it1,n-2);  /*   it1=n-2 if t==tarr[n-1] */
   it2=it1+1;
   *itf=it1;
   *alpha=(tarr[it2]-t)/(tarr[it2]-tarr[it1]);
   return 0;
}

#undef __FUNCT__
#define __FUNCT__ "calcPeriodicInterpFactor"
PetscInt calcPeriodicInterpFactor(PetscInt n,PetscScalar t,PetscScalar tparr[],PetscInt *itf1,PetscInt *itf2,PetscScalar *al1,PetscScalar *al2)
/* Compute interpolation factor and indices for interpolation on a periodic grid */
/* tparr is a periodic array of length n+2 */

{
  PetscInt it1, it2;

  /* 0<=t<=Tc, tparr[n] < Tc < tparr[n+1] */  
  it1=findindex(tparr,n+2,t);   
  if (it1<0) SETERRQ(PETSC_COMM_WORLD,1,"Error in findindex: time out of bound");

  it1=MIN(it1,n); /* it1=n, if t==tparr[n+1]. Not really necessary. */  
  it2=it1+1;  
  /* it1,it2 are referenced to tparr (tparr[0] tparr[1] ... tparr[n] tparr[n+1]) */
  /* 0<=it1<=n, 1<=it2<=n+1 */
  /* Data are given at tparr[1] ... tparr[n] */
  *al1=(tparr[it2]-t)/(tparr[it2]-tparr[it1]);
  *al2=1.0-(*al1);
  

  /* now index to actual data arrays */
  it1=it1-1;
  it2=it2-1;
  
  if (it1==-1) {
    it1=n-1;
    it2=0;
  }
  if (it1==n-1) it2=0;

  *itf1=it1;
  *itf2=it2;
   
  return 0;
}

PetscInt findindex(PetscScalar x[], PetscInt n, PetscScalar xf)
/*      Function to return index i such that x(i)<=xf<x(i+1) */
/*      Inputs: n,x[0:n-1],xf */
/*      Output: findindex */
/*      If xf<x[0], findindex=-1 */
/*      If xf>x[n-1], findindex=-2 */
/*      If x(i)<=xf<x(i+1), findindex=i */
/*      if xf=x[n-1], findindex=n-1 */
{
   PetscInt i;
 
   if (xf<x[0]) {
/*      findindex=-1; */
     return -1;
   }

   if (xf>x[n-1]) {
/*      findindex=-2 */
     return -2;
   }

   for (i=0; i<n-1; i++) {
     if ((x[i]<=xf) && (x[i+1]>xf)) {
/*        findindex=i */
       return i;
     }
   }

/*       findindex=n  !  x(n)==xf */
   return n-1;
}

#undef __FUNCT__
#define __FUNCT__ "interpPeriodicMatrix"
PetscErrorCode interpPeriodicMatrix(PetscScalar tc, Mat *A, PetscScalar cyclePeriod,
                                    PetscInt numPerPeriod, PetscScalar *tdp, 
                                    PeriodicMat *user, const char *fileName)
{
/* Function to interpolate a matrix that is periodic in time with period cyclePeriod.  */
/* tc is the current time and numPerPeriod is the number of instances per period   */
/* at which data are available (to be read from files). */
/* IMPORTANT: Matrix *A MUST have been created and preallocated before  */
/* calling this routine. All other matrices will be created using *A as a template.  */

  PetscScalar t,t1;
  PetscInt im,it0,it1;
/*   static PetscInt iCurrTimeReadLast=-1; */
  PetscErrorCode ierr;
  PetscScalar alpha[2];  
  char tmpFile[PETSC_MAX_PATH_LEN];
  PetscViewer fd;

  if (user->firstTime) {
    if (numPerPeriod>MAX_MATRIX_NUM_PER_PERIOD) {
      SETERRQ(PETSC_COMM_WORLD,1,"Number of allowable matrices in PeriodicMat struct exceeded by requested number ! Increase MAX_MATRIX_NUM_PER_PERIOD.");
    }
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Initializing PeriodicMat object %s\n",fileName);CHKERRQ(ierr);        
    user->numPerPeriod = numPerPeriod;
    for (im=0; im<numPerPeriod; im++) {
      ierr = MatDuplicate(*A,MAT_DO_NOT_COPY_VALUES,&user->Ap[im]);CHKERRQ(ierr);        
	  strcpy(tmpFile,"");
	  sprintf(tmpFile,"%s%02d",fileName,im);	  
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading matrix from file %s\n", tmpFile);CHKERRQ(ierr);        
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = MatLoad(user->Ap[im],fd);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
    }
    user->firstTime = PETSC_FALSE;
  }
  
  t=tc; /* current time */
  if (t<0.) t=cyclePeriod+t;
  t1=t-cyclePeriod*floor(t/cyclePeriod);
  ierr=calcPeriodicInterpFactor(numPerPeriod,t1,tdp,&it0,&it1,&alpha[0],&alpha[1]);  CHKERRQ(ierr);  
/*   ierr = PetscPrintf(PETSC_COMM_WORLD,"tc=%lf,t1=%lf,it0=%d,it1=%d,a1=%17.16lf,a2=%17.16lf\n",tc,t1,it0,it1,alpha[0],alpha[1]);CHKERRQ(ierr);   */
  
/* interpolate to current time   */
  ierr = MatAXPBYmy(alpha[0],alpha[1],user->Ap[it0],user->Ap[it1],A);CHKERRQ(ierr);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "interpPeriodicVector"
PetscErrorCode interpPeriodicVector(PetscScalar tc, Vec *u, PetscScalar cyclePeriod,
                                    PetscInt numPerPeriod, PetscScalar *tdp, 
                                    PeriodicVec *user, const char *fileName)
{
/* Function to interpolate a vector that is periodic in time with period cyclePeriod.  */
/* tc is the current time and numPerPeriod is the number of instances per period   */
/* at which data are available (to be read from files). */
/* IMPORTANT: Vectors u0 and u1 MUST have been created and preallocated before  */
/* calling this routine. Use VecDuplicate to do this.  */

  PetscScalar t,t1;
  PetscInt im,it0,it1;
/*   static PetscInt iCurrTimeReadLast=-1; */
  PetscErrorCode ierr;
  PetscScalar alpha[2];  
  char tmpFile[PETSC_MAX_PATH_LEN];
  PetscViewer fd;

  if (user->firstTime) {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Initializing PeriodicVec object %s\n",fileName);CHKERRQ(ierr);
    user->numPerPeriod = numPerPeriod;  
    ierr = VecDuplicateVecs(*u,numPerPeriod,&user->up);CHKERRQ(ierr);    
    for (im=0; im<numPerPeriod; im++) {
	  strcpy(tmpFile,"");
	  sprintf(tmpFile,"%s%02d",fileName,im);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading vector from file %s\n", tmpFile);CHKERRQ(ierr);  
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
      ierr = VecLoad(user->up[im],fd);CHKERRQ(ierr); /* IntoVector */
      ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);          
    }
    user->firstTime = PETSC_FALSE;
  }

  t=tc; /* current time */
  if (t<0.) t=cyclePeriod+t;
  t1=t-cyclePeriod*floor(t/cyclePeriod);
  ierr=calcPeriodicInterpFactor(numPerPeriod,t1,tdp,&it0,&it1,&alpha[0],&alpha[1]);  CHKERRQ(ierr);  
/*   ierr = PetscPrintf(PETSC_COMM_WORLD,"tc=%lf,t1=%lf,it0=%d,it1=%d,a1=%17.16lf,a2=%17.16lf\n",tc,t1,it0,it1,alpha[0],alpha[1]);CHKERRQ(ierr);   */
  
/* interpolate to current time   */
  ierr = VecAXPBYmy(alpha[0],alpha[1],user->up[it0],user->up[it1],u);CHKERRQ(ierr);  

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "destroyPeriodicVec"
PetscErrorCode destroyPeriodicVec(PeriodicVec *user)
{
/* Function to destroy Vec's in a PeriodicVec struct */

  PetscErrorCode ierr;

  ierr = VecDestroyVecs(user->numPerPeriod,&(user->up));CHKERRQ(ierr);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "destroyPeriodicMat"
/* Function to destroy Mats's in a PeriodicMat struct */

PetscErrorCode destroyPeriodicMat(PeriodicMat *user)
{
  PetscErrorCode ierr;
  PetscInt im;
  
  for (im=0; im<user->numPerPeriod; im++) {
    ierr = MatDestroy(&(user->Ap[im]));CHKERRQ(ierr);
  }
  
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "destroyPeriodicArray"
PetscErrorCode destroyPeriodicArray(PeriodicArray *user)
{
/* Function to destroy arrays in a PeriodicArray struct */

  PetscErrorCode ierr;
  PetscInt im;
  
  for (im=0; im<user->numPerPeriod; im++) {  
    ierr = PetscFree(user->up[im]);CHKERRQ(ierr);    
  }
  
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "interpTimeDependentMatrix"
PetscErrorCode interpTimeDependentMatrix(PetscScalar tc, Mat *A, PetscInt numTimes, PetscScalar *tdt, 
                                    TimeDependentMat *user, const char *fileName)
{
/* Function to interpolate a time-dependent matrix. */
/* tc is the current time and numTimes are the number of time slices (in tdt) at */
/* at which data are available. */
/* IMPORTANT: Matrix *A MUST have been created and preallocated before  */
/* calling this routine. All other matrices will be created using *A as a template.  */

  PetscErrorCode ierr;
  PetscInt itc;  
  PetscScalar alpha;
  char tmpFile[PETSC_MAX_PATH_LEN];
  PetscViewer fd;

  if (user->firstTime) {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Initializing TimeDependentMat object %s\n",fileName);CHKERRQ(ierr);
	ierr = MatDuplicate(*A,MAT_DO_NOT_COPY_VALUES,&user->Atd[0]);CHKERRQ(ierr);
	ierr = MatDuplicate(*A,MAT_DO_NOT_COPY_VALUES,&user->Atd[1]);CHKERRQ(ierr);
	user->itcurr=-1;
    user->firstTime = PETSC_FALSE;
  }

  if ((tc<tdt[0]) || (tc>tdt[numTimes-1])) {
    SETERRQ(PETSC_COMM_WORLD,1,"Error in interpTimeDependentMatrix: time out of bound");
  }

  ierr = calcInterpFactor(numTimes,tc,tdt,&itc,&alpha); CHKERRQ(ierr);
  if (itc != user->itcurr) { /* time to read new bracketing slices */
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading new bracketing slices for matrix %s at time = %g: %d and %d\n",fileName,tc,itc,itc+1);CHKERRQ(ierr);
	strcpy(tmpFile,"");
	sprintf(tmpFile,"%s%02d",fileName,itc);	  
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading matrix from file %s\n", tmpFile);CHKERRQ(ierr);        
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatLoad(user->Atd[0],fd);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

	strcpy(tmpFile,"");
	sprintf(tmpFile,"%s%02d",fileName,itc+1);	  
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading matrix from file %s\n", tmpFile);CHKERRQ(ierr);        
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatLoad(user->Atd[1],fd);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

	user->itcurr=itc;
  }
/* interpolate to current time   */
  ierr = MatAXPBYmy(alpha,1.0-alpha,user->Atd[0],user->Atd[1],A);CHKERRQ(ierr);

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "destroyTimeDependentMat"
/* Function to destroy Mats's in a TimeDependentMat struct */

PetscErrorCode destroyTimeDependentMat(TimeDependentMat *user)
{
  PetscErrorCode ierr;
  
  ierr = MatDestroy(&(user->Atd[0]));CHKERRQ(ierr);
  ierr = MatDestroy(&(user->Atd[1]));CHKERRQ(ierr);
  
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "destroyTimeDependentArray"
/* Function to destroy arrays in a TimeDependentArray struct */

PetscErrorCode destroyTimeDependentArray(TimeDependentArray *user)
{
  PetscErrorCode ierr;

  ierr = PetscFree(user->utd[0]);CHKERRQ(ierr);    
  ierr = PetscFree(user->utd[1]);CHKERRQ(ierr);    
  
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "interpTimeDependentVector"
PetscErrorCode interpTimeDependentVector(PetscScalar tc, Vec *u, PetscInt numTracers, 
                                      PetscInt nt, PetscScalar *t, Vec **ut)
{

  PetscInt itf, itr;
  PetscErrorCode ierr;
  PetscScalar alpha[2];  
  PetscScalar zero=0.0;

  for (itr=0; itr<numTracers; itr++) {
    ierr = VecSet(u[itr],zero); CHKERRQ(ierr);
  }
  if (tc>=t[0]) {
    ierr = calcInterpFactor(nt,tc,t,&itf,&alpha[0]); CHKERRQ(ierr);
    alpha[1]=1.0-alpha[0];
    for (itr=0; itr<numTracers; itr++) {
        VecMAXPY(u[itr],2,alpha,&ut[itr][itf]);
    }
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: time < %10.5f. Assuming u=0\n", t[0]);CHKERRQ(ierr);
  }
    
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "writeBinaryScalarData"
PetscErrorCode writeBinaryScalarData(const char *fileName, PetscScalar *arr, PetscInt N, PetscBool appendToFile)
{
  PetscErrorCode ierr;
  PetscViewer fd;
  int fp;
  PetscMPIInt myId;

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myId);CHKERRQ(ierr);  

  if (myId==0) { /* this shouldn't really be necessary, but without it, all processors seem to be writing in append mode */
	if (appendToFile) {
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,fileName,FILE_MODE_APPEND,&fd);CHKERRQ(ierr);
	} else {
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,fileName,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	}  
  
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscBinaryWrite(fp,arr,N,PETSC_SCALAR);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  }
  
  return 0;
}
