#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "petsc_matvec_utils.h"
#include "tmm_timer.h"
#include "tmm_forcing_utils.h"
#include "tmm_profile_utils.h"
#include "tmm_profile_data.h"
#include "tmm_main.h"

PetscBool doExtraWrite = PETSC_FALSE;
StepTimer extraWriteTimer;
char *outFile[MAXNUMTRACERS];
PetscFileMode OUTPUT_FILE_MODE;
char outTimeFile[PETSC_MAX_PATH_LEN];  
PetscBool appendOutput = PETSC_FALSE;
FILE *fptime;
PetscViewer fdout[MAXNUMTRACERS];
PetscInt numExtraTracers;
PetscInt itrExtra[MAXNUMTRACERS];

#undef __FUNCT__
#define __FUNCT__ "iniTMMWrite"
PetscErrorCode iniTMMWrite(PetscScalar tc, PetscInt Iter, PetscInt numTracers, Vec *v, PetscBool append)
{

  PetscInt itr;
  PetscBool flg;
  PetscErrorCode ierr;

  ierr = PetscOptionsHasName(PETSC_NULL,"-write_extra",&doExtraWrite);CHKERRQ(ierr);

  if (doExtraWrite) {    
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Extra tracer output has been activated\n");CHKERRQ(ierr);
    ierr = iniStepTimer("write_extra_", Iter0, &extraWriteTimer);CHKERRQ(ierr);

    numExtraTracers=MAXNUMTRACERS;    
    ierr = PetscOptionsGetIntArray(PETSC_NULL,"-write_extra_tracer_indices",itrExtra,&numExtraTracers,&flg);CHKERRQ(ierr);
	if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate extra tracer indices with the -write_extra_tracer_indices option!");

/* Output file */
	for (itr=0; itr<numExtraTracers; itr++) {
	  outFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
	}

	ierr = PetscOptionsGetStringArray(PETSC_NULL,"-o_extra",outFile,&numExtraTracers,&flg);CHKERRQ(ierr);
	if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate extra output file name(s) with the -o_extra option");

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Extra output will be written to:\n");CHKERRQ(ierr);
	for (itr=0; itr<numExtraTracers; itr++) {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"   Tracer %d: %s\n", itrExtra[itr],outFile[itr]);CHKERRQ(ierr);
	}

	ierr = PetscOptionsHasName(PETSC_NULL,"-append_extra",&appendOutput);CHKERRQ(ierr);
	if (appendOutput) {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Extra output will be appended\n");CHKERRQ(ierr);
	  OUTPUT_FILE_MODE=FILE_MODE_APPEND;
	} else {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Extra output will overwrite existing file(s)\n");CHKERRQ(ierr);
	  OUTPUT_FILE_MODE=FILE_MODE_WRITE;
	}    

/* Output times */
	ierr = PetscOptionsGetString(PETSC_NULL,"-time_file_extra",outTimeFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
	if (!flg) {
	  strcpy(outTimeFile,"");
	  sprintf(outTimeFile,"%s","output_time_extra.txt");
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Extra output times will be written to %s\n",outTimeFile);CHKERRQ(ierr);

	for (itr=0; itr<numExtraTracers; itr++) {       
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,outFile[itr],OUTPUT_FILE_MODE,&fdout[itr]);CHKERRQ(ierr);
	}

	if (!appendOutput) {
	  ierr = PetscFOpen(PETSC_COMM_WORLD,outTimeFile,"w",&fptime);CHKERRQ(ierr);  
	  ierr = PetscFPrintf(PETSC_COMM_WORLD,fptime,"%d   %10.5f\n",Iter0,time0);CHKERRQ(ierr);     
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing extra output at time %10.5f, step %d\n", time0,Iter0);CHKERRQ(ierr);  
	  for (itr=0; itr<numExtraTracers; itr++) {       
		ierr = VecView(v[itrExtra[itr]],fdout[itr]);CHKERRQ(ierr);
	  }
	} else {
		ierr = PetscFOpen(PETSC_COMM_WORLD,outTimeFile,"a",&fptime);CHKERRQ(ierr);  
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Opening extra output file(s) for output. Initial condition will NOT be written\n");CHKERRQ(ierr);
	}

  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "doTMMWrite"
PetscErrorCode doTMMWrite(PetscScalar tc, PetscInt iLoop, PetscInt numTracers, Vec *v)
{

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  PetscInt itr;
  PetscErrorCode ierr;

  if (doExtraWrite) {
	if (Iter0+iLoop>=extraWriteTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
	  if (extraWriteTimer.count<=extraWriteTimer.numTimeSteps) {
		extraWriteTimer.count++;
	  }
	  if (extraWriteTimer.count==extraWriteTimer.numTimeSteps) { /* time to write out */

		ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing extra output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
		ierr = PetscFPrintf(PETSC_COMM_WORLD,fptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);
		for (itr=0; itr<numExtraTracers; itr++) {		
		  ierr = VecView(v[itrExtra[itr]],fdout[itr]);CHKERRQ(ierr);
        }
        
		ierr = updateStepTimer("write_extra_", Iter0+iLoop, &extraWriteTimer);CHKERRQ(ierr);

	  }
	}
  }  

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "finalizeTMMWrite"
PetscErrorCode finalizeTMMWrite(PetscScalar tc, PetscInt Iter, PetscInt numTracers)
{

  PetscInt itr;
  PetscErrorCode ierr;

  if (doExtraWrite) {

	for (itr=0; itr<numExtraTracers; itr++) {
	  ierr = PetscViewerDestroy(&fdout[itr]);CHKERRQ(ierr);
	}
	ierr = PetscFClose(PETSC_COMM_WORLD,fptime);CHKERRQ(ierr);

  }
  return 0;
}
