#define DEFINE_VARIABLES

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "petsctime.h"

#include "tmm_petsc_matvec_utils.h"
#include "tmm_forcing_utils.h"
#include "tmm_timer.h"
#include "tmm_profile_utils.h"
#include "tmm.h"
#include "tmm_share.h"
#include "tmm_variables.h"

extern PetscErrorCode TMMInitializeTMs();

#undef __FUNCT__
#define __FUNCT__ "TMMInitialize"
PetscErrorCode TMMInitialize(PetscInt *Iter0ret, PetscInt *maxStepsret, PetscScalar *time0ret, PetscScalar *deltaTClockret)
{

//   TMMState s;
  PetscErrorCode ierr;
  PetscBool flg1,flg2;
  PetscInt il;
  PetscInt n;  
  PetscMPIInt myId;  
  PetscScalar zero = 0.0;
  
//----------------------------------------------------------------------

//   PetscCall(PetscObjectGetComm((PetscObject)state, &comm));
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myId);CHKERRQ(ierr);
  myId=myId+1; /* process ID (starting at 1) */

  PetscPushErrorHandler(PetscAbortErrorHandler,NULL); /* force code to abort on error */

/* Some global defaults */
  time0=0.0;
  deltaTClock=1.0;
  Iter0=0;
    
/* To be safe we also set these */
  nb=-1;
  lSize=-1;
  gBCSize=-1;
  lBCSize=-1;
 
/* TM's */
  periodicMatrix = PETSC_FALSE;
  timeDependentMatrix = PETSC_FALSE;

  rescaleForcing = PETSC_FALSE;
  prescribedBCInUse = PETSC_FALSE;
  calcBCInUse = PETSC_FALSE;

/* Process options and load files */  
/* Time step data */
  ierr = PetscOptionsGetReal(NULL,NULL,"-deltat_clock",&deltaTClock,&flg1);CHKERRQ(ierr);

  ierr = PetscOptionsGetReal(NULL,NULL,"-t0",&time0,&flg1);CHKERRQ(ierr);  
  ierr = PetscOptionsGetInt(NULL,NULL,"-iter0",&Iter0,&flg2);CHKERRQ(ierr);
  if ((flg1) && (!flg2)) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate both or neither time0 and Iter0 with the -t0 and -iter0 flags");
  if ((!flg1) && (flg2)) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate both or neither time0 and Iter0 with the -t0 and -iter0 flags");        
  ierr = PetscOptionsGetInt(NULL,NULL,"-max_steps",&maxSteps,&flg1);CHKERRQ(ierr);
  if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate maximum number of steps with the -max_steps option");

/* Initialize profile data and create template vector */
  ierr = iniProfileData(myId);CHKERRQ(ierr);
  if (useProfiles) {
/*  If useProfiles is switched on then global and local sizes are computed in iniProfileData */
/*  and we use them here. Otherwise, we set these in doTMMInitializeTMs by reading in a matrix. */
    ierr = VecCreate(PETSC_COMM_WORLD,&templateVec);CHKERRQ(ierr);
    ierr = VecSetSizes(templateVec,lSize,PETSC_DECIDE);CHKERRQ(ierr);
    ierr = VecSetFromOptions(templateVec);CHKERRQ(ierr);
    ierr = VecGetSize(templateVec,&n);CHKERRQ(ierr);
  }

  ierr = TMMInitializeTMs();
  
/* create template vector here if not using profiles */
  if (!useProfiles) {
/*  Use layout info computed in TMMInitializeTMs */
    ierr = MatGetSize(Ae,0,&n);CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_WORLD,&templateVec);CHKERRQ(ierr);
    ierr = VecSetSizes(templateVec,PETSC_DECIDE,n);CHKERRQ(ierr);
    ierr = VecSetFromOptions(templateVec);CHKERRQ(ierr);
    ierr = VecGetLocalSize(templateVec,&lSize);CHKERRQ(ierr); /* lSize is otherwise set in iniProfiles */
  }  

// We now know both the global and local sizes of the vectors
  nb = n;
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer vectors are of length %d\n", n);CHKERRQ(ierr);  

/*   Compute global indices for local piece of vectors */
  ierr = VecGetOwnershipRange(templateVec,&gLow,&gHigh);CHKERRQ(ierr);
  gHigh = gHigh - 1; /* Note: gHigh is one more than the last local element */
  ierr = PetscMalloc(lSize*sizeof(PetscInt),&gIndices);CHKERRQ(ierr);  
  for (il=0; il<lSize; il++) {
    gIndices[il] = il + gLow;
  }  

  ierr = VecSet(templateVec,zero);CHKERRQ(ierr);
  
  *Iter0ret=Iter0;
  *maxStepsret=maxSteps;
  *time0ret=time0;
  *deltaTClockret=deltaTClock;

  return 0;
}

