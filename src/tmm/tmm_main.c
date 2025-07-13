static char help[] = "\n";

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "petsctime.h"

#include "tmm_forcing_utils.h"
#include "tmm_timer.h"
#include "tmm.h"
#include "tmm_external_forcing.h"
#include "tmm_external_bc.h"
#include "tmm_monitor.h"
#include "tmm_misfit.h"

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{

  PetscScalar deltaTClock, time0;
  PetscInt maxSteps, Iter0;

  PetscErrorCode ierr;
  PetscScalar t1, t2;
  PetscScalar tc=0.0;
  PetscScalar tf=0.0;
  PetscInt iLoop=1, Iterc=1;
  
  TMMState states[5];
  PetscBool doOutput = PETSC_TRUE;
  
  PetscInt is;
  PetscInt numStates=MAXNUMSTATES; /* This will be overwritten by the actual number of states below */
  char *prefixes[10];
  
  PetscBool flg;
  
  PetscInitialize(&argc,&args,(char *)0,help);

  ierr = PetscTime(&t1);CHKERRQ(ierr); /* start counting wall clock time */  
  
  for (is=0; is<MAXNUMSTATES; is++) {
    prefixes[is] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
    prefixes[is] = NULL;
  }

  doOutput = PETSC_TRUE;

  ierr = TMMInitialize(&Iter0,&maxSteps,&time0,&deltaTClock);

  ierr = PetscOptionsGetStringArray(NULL,NULL,"-prefixes",prefixes,&numStates,&flg);CHKERRQ(ierr);
  if (flg) {
    if (numStates>MAXNUMSTATES) {
      SETERRQ(PETSC_COMM_WORLD,1,"ERROR!: Maximum number of states exceeded! Increase MAXNUMSTATES in tmmimpl.h");
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Prefixes have been specified for %d states\n",numStates); CHKERRQ(ierr);
    for (is=0; is<numStates; is++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  Prefix for state %d is %s\n", is,prefixes[is]); CHKERRQ(ierr);
    }  
  } else {
    numStates=1;
  }

  for (is=0; is<numStates; is++) {
    ierr = TMMCreate(&states[is]);
    ierr = TMMSetFromOptions(states[is],prefixes[is],doOutput);
    if (states[is]->useExternalForcing) {
      ierr = TMMSetIniExtForcFunction(states[is],iniExternalForcing,NULL);
      ierr = TMMSetCalcExtForcFunction(states[is],calcExternalForcing,NULL);
      ierr = TMMSetWriExtForcFunction(states[is],writeExternalForcing,NULL);
      ierr = TMMSetFinExtForcFunction(states[is],finalizeExternalForcing,NULL);
      ierr = TMMSetReiExtForcFunction(states[is],reInitializeExternalForcing,NULL);
    }

    if (states[is]->doCalcBC) {
      ierr = TMMSetIniCalcBCFunction(states[is],iniCalcBC,NULL);
      ierr = TMMSetCalcCalcBCFunction(states[is],calcCalcBC,NULL);
      ierr = TMMSetWriCalcBCFunction(states[is],writeCalcBC,NULL);
      ierr = TMMSetFinCalcBCFunction(states[is],finalizeCalcBC,NULL);
      ierr = TMMSetReiCalcBCFunction(states[is],reInitializeCalcBC,NULL);
    }
    
    if (states[is]->useMonitor) {
      ierr = TMMSetIniMonitorFunction(states[is],iniMonitor,NULL);
      ierr = TMMSetCalcMonitorFunction(states[is],calcMonitor,NULL);
      ierr = TMMSetWriMonitorFunction(states[is],writeMonitor,NULL);
      ierr = TMMSetFinMonitorFunction(states[is],finalizeMonitor,NULL);
    }
 
    if (states[is]->doMisfit) {
      ierr = TMMSetIniMisfitFunction(states[is],iniMisfit,NULL);
      ierr = TMMSetCalcMisfitFunction(states[is],calcMisfit,NULL);
      ierr = TMMSetWriMisfitFunction(states[is],writeMisfit,NULL);
      ierr = TMMSetFinMisfitFunction(states[is],finalizeMisfit,NULL);
    }  
  }
  
/* Start time stepping loop */
  for (iLoop = 1; iLoop <= maxSteps; iLoop++) {
/*  iLoop -> iLoop+1 (convention) */  
    tc=time0 + deltaTClock*(iLoop-1);  /*  current time (time at beginning of step) */
    tf=time0 + deltaTClock*iLoop;  /*  future time (time at end of step) */
    Iterc=Iter0+iLoop-1;

    ierr = TMMUpdateTMs(tc);

    for (is=0; is<numStates; is++) {
      ierr = TMMForcingUpdate(tc,Iterc,iLoop,states[is]);
      ierr = TMMTimeStep(tc,Iterc,iLoop,states[is]);
      ierr = TMMTimeStepPost(tf,Iterc,iLoop,states[is]);
      ierr = TMMOutput(tf,Iterc,iLoop,states[is]);
    }
    	
  } /* end of time-stepping loop */

  tc=time0 + deltaTClock*maxSteps; /*  time at end of integration */

  for (is=0; is<numStates; is++) {
    ierr = TMMDestroy(tc,&states[is]);
  }

  ierr = TMMFinalize(tc);

  ierr=PetscTime(&t2); CHKERRQ(ierr); /* stop counting wall clock time */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Wall clock time: %10.5f\n", t2-t1); CHKERRQ(ierr);

  ierr = PetscFinalize(); CHKERRQ(ierr);
        
  return 0;
}
