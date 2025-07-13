#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "petsctime.h"

#include "tmm_petsc_matvec_utils.h"
#include "tmm_forcing_utils.h"
#include "tmm_timer.h"
#include "tmm_profile_utils.h"
// #include "tmm_profile_data.h"
#include "tmm.h"
// #include "tmm_external_forcing.h"
// #include "tmm_external_bc.h"
#include "tmm_share.h"
#include "tmm_variables.h"

#undef __FUNCT__
#define __FUNCT__ "TMMInitializeTMs"
PetscErrorCode TMMInitializeTMs()
{

  PetscErrorCode ierr;
  PetscBool flg;
  char tmpFile[PETSC_MAX_PATH_LEN];
  PetscMPIInt numProcessors;  
  PetscViewer fd;
  
//----------------------------------------------------------------------

  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcessors);CHKERRQ(ierr);  

/* Matrices */
  ierr = PetscOptionsGetString(NULL,NULL,"-me",mateFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate binary matrix file with the -me option");
  ierr = PetscOptionsGetString(NULL,NULL,"-mi",matiFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate binary matrix file with the -mi options");

  ierr = MatCreate(PETSC_COMM_WORLD,&Ae);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&Ai);CHKERRQ(ierr);

/*  Set layout information */ 
/*  If useProfiles is switched on then this has been previously computed in iniProfileData and  */
/*  we use it here. Otherwise, read in a matrix later below to let PETSc determine both the */
/*  global size of the matrix (which we don't know yet if useProfiles is off) and local size. */
/*  Note: PETSc uses a simple formula to figure out the partitioning (see PetscSplitOwnership). */ 
  if ((useProfiles) && (numProcessors>1)) {
	ierr = MatSetSizes(Ae,lSize,lSize,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
	ierr = MatSetSizes(Ai,lSize,lSize,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
  }   
  
  ierr = MatSetType(Ae,MATMPIAIJ);CHKERRQ(ierr);      
  ierr = MatSetFromOptions(Ae);CHKERRQ(ierr);  
  ierr = MatSetType(Ai,MATMPIAIJ);CHKERRQ(ierr);        
  ierr = MatSetFromOptions(Ai);CHKERRQ(ierr);
  
  ierr = PetscOptionsHasName(NULL,NULL,"-periodic_matrix",&periodicMatrix);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-time_dependent_matrix",&timeDependentMatrix);CHKERRQ(ierr);
  if (periodicMatrix & timeDependentMatrix) {
    SETERRQ(PETSC_COMM_WORLD,1,"Cannot use both -periodic_matrix and -time_dependent_matrix together!");
  }
  
  if (periodicMatrix) {    
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic matrices specified\n");CHKERRQ(ierr);
    ierr=PetscStrcat(mateFile,"_");CHKERRQ(ierr);
    ierr=PetscStrcat(matiFile,"_");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Ae basename is %s\n", mateFile);CHKERRQ(ierr); 
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Ai basename is %s\n", matiFile);CHKERRQ(ierr);     
    
/*  read time data */
    ierr = PeriodicTimerCreate(&matrixPeriodicTimer);CHKERRQ(ierr);
    ierr = PeriodicTimerIni("matrix_", NULL, NULL, matrixPeriodicTimer);CHKERRQ(ierr);
    
/*  Read here to set sparsity pattern */
	strcpy(tmpFile,"");
	sprintf(tmpFile,"%s%02d",mateFile,0);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatLoad(Ae,fd);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

	strcpy(tmpFile,"");
	sprintf(tmpFile,"%s%02d",matiFile,0);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatLoad(Ai,fd);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    

    Aep.firstTime = PETSC_TRUE;
    Aip.firstTime = PETSC_TRUE;

  } else if (timeDependentMatrix) {
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Time-dependent matrices specified\n");CHKERRQ(ierr);
    ierr=PetscStrcat(mateFile,"_");CHKERRQ(ierr);
    ierr=PetscStrcat(matiFile,"_");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Ae basename is %s\n", mateFile);CHKERRQ(ierr); 
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Ai basename is %s\n", matiFile);CHKERRQ(ierr);     
    
/*  read time data */
    ierr = TimeDependentTimerCreate(&matrixTimeDependentTimer);CHKERRQ(ierr);
    ierr = TimeDependentTimerIni("matrix_", NULL, NULL, matrixTimeDependentTimer);CHKERRQ(ierr);
    
/*  Read here to set sparsity pattern */
	strcpy(tmpFile,"");
	sprintf(tmpFile,"%s%02d",mateFile,0);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatLoad(Ae,fd);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

	strcpy(tmpFile,"");
	sprintf(tmpFile,"%s%02d",matiFile,0);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = MatLoad(Ai,fd);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    

    Aetd.firstTime = PETSC_TRUE;
    Aitd.firstTime = PETSC_TRUE;  
    
  } else { /*  not periodic. read matrices here */

    constantMatrix = PETSC_TRUE;

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,mateFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading Ae from file %s\n", mateFile);CHKERRQ(ierr);  
	ierr = MatLoad(Ae,fd);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,matiFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading Ai from file %s\n", matiFile);CHKERRQ(ierr);  
	ierr = MatLoad(Ai,fd);CHKERRQ(ierr);    
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
    
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TMMUpdateTMs"
PetscErrorCode TMMUpdateTMs(PetscScalar tc)
{

	PetscErrorCode ierr;

/*  interpolate Ae,Ai,Rfs to current time (tc) */
    if (periodicMatrix) {
      ierr = PeriodicMatInterp(tc,&Ae,matrixPeriodicTimer->cyclePeriod,matrixPeriodicTimer->numPerPeriod,
                                  matrixPeriodicTimer->tdp,&Aep,mateFile);
      ierr = PeriodicMatInterp(tc,&Ai,matrixPeriodicTimer->cyclePeriod,matrixPeriodicTimer->numPerPeriod,
                                  matrixPeriodicTimer->tdp,&Aip,matiFile);
      if (rescaleForcing) {
		ierr = PeriodicVecInterp(tc,&Rfs,matrixPeriodicTimer->cyclePeriod,matrixPeriodicTimer->numPerPeriod,
                                  matrixPeriodicTimer->tdp,Rfsp,rfsFile);
      }
    } else if (timeDependentMatrix) {
      ierr = TimeDependentMatInterp(tc,&Ae,matrixTimeDependentTimer->numTimes,matrixTimeDependentTimer->tdt,&Aetd,mateFile);
      ierr = TimeDependentMatInterp(tc,&Ai,matrixTimeDependentTimer->numTimes,matrixTimeDependentTimer->tdt,&Aitd,matiFile);
      if (rescaleForcing) {
		ierr = TimeDependentVecInterp(tc,&Rfs,matrixTimeDependentTimer->numTimes,matrixTimeDependentTimer->tdt,Rfstd,rfsFile);
      }
    }

        return 0;
}

