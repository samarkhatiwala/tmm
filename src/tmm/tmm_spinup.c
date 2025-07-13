#define DEFINE_SPINUP_VARIABLES

#include <petscsnes.h>
#include "petsctime.h"

#include "tmm_forcing_utils.h"
// #include "tmm_profile_data.h"
#include "tmm_timer.h"
#include "tmm.h"
#include "tmm_spinup_share.h"
#include "tmm_signal_utils.h"

// NGMRES debug
// #include "/Volumes/data2/spk/petsc-3.16.2/include/petsc/private/petscimpl.h"
// #include "/Volumes/data2/spk/petsc-3.16.2/include/petsc/private/snesimpl.h"
// #include "/Volumes/data2/spk/petsc-3.16.2/src/snes/impls/ngmres/snesngmres.h"

PetscErrorCode spinupFunction(SNES,Vec,Vec,void*);
PetscErrorCode integrateOnePeriod(PetscInt maxSteps, PetscInt Iter0, PetscScalar time0, PetscScalar deltaTClock, TMMState state, PetscBool doOutput, PetscBool runExternalModel);
PetscErrorCode spinupMonitor(SNES,PetscInt its,PetscReal fnorm,void *);
PetscErrorCode checkConvergence(SNES snes,PetscInt its,PetscReal xnorm,PetscReal snorm,PetscReal fnorm,SNESConvergedReason *reason,void *);

extern PetscErrorCode XtoVconvert(PetscInt numTracers, Vec *c, Vec X, PetscScalar *XTovScaleFac);
extern PetscErrorCode VtoXconvert(PetscInt numTracers, Vec *c, Vec X, PetscScalar *vToXScaleFac);

PetscBool reinitializeForcing = PETSC_FALSE;

PetscErrorCode spinupCalc(PetscInt maxSteps, PetscInt Iter0, PetscScalar time0, PetscScalar deltaTClock, TMMState state, PetscInt *neval, PetscReal *nm)
{

  PetscErrorCode ierr;
  TMMSPINUP spinupstate;  
  Vec X, r;
  SNES snes;         /* nonlinear solver context */
  PetscInt itr, maxValsToRead;
  PetscBool flg, flg1;
  PetscBool doMonitor = PETSC_FALSE;
  PetscInt numTracers;
  const char *prefix;
  
  PetscFunctionBeginUser;

  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;
  
  ierr = PetscMalloc(numTracers*sizeof(PetscScalar),&spinupstate.XTovScaleFac);CHKERRQ(ierr);   
  ierr = PetscMalloc(numTracers*sizeof(PetscScalar),&spinupstate.vToXScaleFac);CHKERRQ(ierr);   

  for (itr=0; itr<numTracers; itr++) {
    spinupstate.XTovScaleFac[itr]=1.0;
    spinupstate.vToXScaleFac[itr]=1.0;
  }      

  maxValsToRead = numTracers;
  ierr = PetscOptionsGetRealArray(NULL,prefix,"-vscale_fac",spinupstate.XTovScaleFac,&maxValsToRead,&flg);
  if (flg) {
    if (maxValsToRead != numTracers) {
      SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of vscale_fac values specified");
    }
    for (itr=0; itr<numTracers; itr++) {
      if (spinupstate.XTovScaleFac[itr]>0.0) {
        spinupstate.vToXScaleFac[itr]=1.0/spinupstate.XTovScaleFac[itr];
      }      
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d scale factor=%15.11f\n",itr,spinupstate.XTovScaleFac[itr]);CHKERRQ(ierr);
    }      
  }

spinupstate.maxSteps=maxSteps;
spinupstate.Iter0=Iter0;
spinupstate.time0=time0;
spinupstate.deltaTClock=deltaTClock;
spinupstate.itf=0;
spinupstate.checkpointFreq=-1;
spinupstate.state=state;

  ierr = SNESCreate(PETSC_COMM_WORLD,&snes);CHKERRQ(ierr);

//   ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished creating SNES\n"); CHKERRQ(ierr);

  ierr = SNESSetType(snes,SNESNGMRES);CHKERRQ(ierr);
  ierr = SNESSetTolerances(snes,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT,100);

  ierr = VecCreate(PETSC_COMM_WORLD,&X);CHKERRQ(ierr);
  ierr = VecSetSizes(X,lSize*numTracers,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(X);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&r);CHKERRQ(ierr);
  ierr = VecDuplicate(X,&spinupstate.Xini);CHKERRQ(ierr);

//   ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished creating Xini\n"); CHKERRQ(ierr);
  
    ierr = SNESSetFunction(snes,r,spinupFunction,&spinupstate);CHKERRQ(ierr);
    
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,prefix,"-spinup_monitor",&doMonitor);CHKERRQ(ierr);
  if (doMonitor) {
	ierr = SNESMonitorSet(snes,spinupMonitor,&spinupstate,NULL);  CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL,prefix,"-spinup_checkpoint_freq",&spinupstate.checkpointFreq,&flg1);CHKERRQ(ierr);
    if (flg1) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"SNES solution will be written out every %d iterations\n", spinupstate.checkpointFreq); CHKERRQ(ierr);    
    }
	ierr = PetscPrintf(PETSC_COMM_WORLD,"SNES history will be written to %s\n","snes_spinup_log.txt");CHKERRQ(ierr);
	ierr = PetscFOpen(PETSC_COMM_WORLD,"snes_spinup_log.txt","w",&spinupstate.logfp);CHKERRQ(ierr);
  }

  ierr = PetscOptionsHasName(NULL,prefix,"-run_external_model",&spinupstate.runExternalModel);CHKERRQ(ierr);
  if (spinupstate.runExternalModel) {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"External model will be run\n");CHKERRQ(ierr);
  }


  ierr = PetscOptionsHasName(NULL,prefix,"-check_external_convergence",&flg);CHKERRQ(ierr);
  if (flg) {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"External convergence will be checked\n");CHKERRQ(ierr);
	spinupstate.convergedValue = 1; // default value to 
	ierr = PetscOptionsGetInt(NULL,prefix,"-external_converged_value",&spinupstate.convergedValue,&flg1);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"External converged value is set to: %d\n", spinupstate.convergedValue); CHKERRQ(ierr);    
	SNESSetConvergenceTest(snes,checkConvergence,&spinupstate,NULL);
  }	
  
SNESHasNPC(snes, &flg);
  if (flg) {
  ierr = PetscPrintf(PETSC_COMM_WORLD,"NPC is set\n"); CHKERRQ(ierr);

} else {
  ierr = PetscPrintf(PETSC_COMM_WORLD,"NPC is NOT set\n"); CHKERRQ(ierr);
}

// c to X
  ierr = VtoXconvert(spinupstate.state->numTracers, spinupstate.state->c, X, spinupstate.vToXScaleFac);CHKERRQ(ierr);

//   ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished converting c to X\n"); CHKERRQ(ierr);

  ierr = SNESSolve(snes,NULL,X);CHKERRQ(ierr);
//   if (flg) {
//     Vec f;
//     ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//     ierr = SNESGetFunction(snes,&f,0,0);CHKERRQ(ierr);
//     ierr = VecView(r,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
//   }

  ierr = XtoVconvert(spinupstate.state->numTracers, spinupstate.state->c, X, spinupstate.XTovScaleFac);CHKERRQ(ierr);

  *neval=spinupstate.itf;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Checking equilibrium ...\n"); CHKERRQ(ierr);

  ierr = checkEquilibrium(maxSteps,Iter0,time0,deltaTClock,state,PETSC_TRUE,spinupstate.runExternalModel,nm);

//   ierr = VecDuplicateVecs(state->c[0],numTracers,&vini);CHKERRQ(ierr);
//   for (itr=0; itr<numTracers; itr++) {
//     ierr = VecCopy(state->c[itr],vini[itr]);CHKERRQ(ierr);
//   }  
// 
// 
// 
//   for (itr=0; itr<numTracers; itr++) {
//     ierr = VecAXPY(vini[itr],minusone,state->c[itr]);CHKERRQ(ierr);
//     ierr = VecNorm(vini[itr],NORM_2,&nm[itr]);
//   }  
//   

  if (doMonitor) {
	ierr = PetscFClose(PETSC_COMM_WORLD,spinupstate.logfp);CHKERRQ(ierr);
  }
  
  ierr = VecDestroy(&X);CHKERRQ(ierr); 
  ierr = VecDestroy(&r);CHKERRQ(ierr);
//   ierr = MatDestroy(&J);CHKERRQ(ierr); 
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
ierr = VecDestroy(&spinupstate.Xini);CHKERRQ(ierr);
//   ierr = VecDestroyVecs(numTracers,&vini);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

// PetscErrorCode doJacobianFinalize()
// {
// 
//   PetscErrorCode ierr;
//   
//   PetscFunctionBeginUser;
// 
// 
//   PetscFunctionReturn(0);
// }

PetscErrorCode spinupFunction(SNES snes,Vec Xx,Vec F,void *ptr)
{

  PetscErrorCode ierr;
TMMSPINUP         *user = (TMMSPINUP*)ptr;

PetscScalar minusone=-1.0;
//   PetscScalar tc, tf;
//   PetscInt iLoop, Iterc;

  PetscFunctionBeginUser;

// X -> c{}

//   ierr = PetscPrintf(PETSC_COMM_WORLD,"Starting copying X to c\n"); CHKERRQ(ierr);

  ierr = XtoVconvert(user->state->numTracers, user->state->c, Xx, user->XTovScaleFac);CHKERRQ(ierr);
//   ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished copying X to c\n"); CHKERRQ(ierr);
  
VecCopy(Xx,user->Xini);

//   ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished copying Xini\n"); CHKERRQ(ierr);


/* Start time stepping loop */
//   for (iLoop = 1; iLoop <= user->maxSteps; iLoop++) {
// /*  iLoop -> iLoop+1 (convention) */  
//     tc=user->time0 + (user->deltaTClock)*(iLoop-1);  /*  current time (time at beginning of step) */
//     tf=user->time0 + (user->deltaTClock)*iLoop;  /*  future time (time at end of step) */
//     Iterc=(user->Iter0)+iLoop-1;
// 
// 	ierr = TMMUpdateTMs(tc);
// 	ierr = TMMUpdateForcing(tc,tf,Iterc,iLoop,user->state);
// 
// // 	ierr = TMMTimeStep(tc,tf,Iterc,iLoop,c,qf,qef,cbc,cbf);
// 	ierr = TMMTimeStep(tc,tf,Iterc,iLoop,user->state);
// 
// 	tc=user->time0 + (user->deltaTClock)*iLoop; /*  time at end of step */
// 
// 	ierr = TMMTimeStepPost(tc,tf,Iterc,iLoop,user->state);
// // 	ierr = TMMOutput(tc,tf,Iterc,iLoop,user->state);
// 
//   } /* end of time-stepping loop */

  ierr = integrateOnePeriod(user->maxSteps, user->Iter0, user->time0, user->deltaTClock, user->state, PETSC_FALSE, user->runExternalModel);CHKERRQ(ierr);

  ierr = VtoXconvert(user->state->numTracers, user->state->c, F, user->vToXScaleFac);CHKERRQ(ierr);

//   ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished copying c to F\n"); CHKERRQ(ierr);

  ierr = VecAXPY(F,minusone,user->Xini);CHKERRQ(ierr);

  user->itf++; // increment function evaluation counter

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished computing function (%d)\n",user->itf); CHKERRQ(ierr);
  
  PetscFunctionReturn(0);
}

PetscErrorCode integrateOnePeriod(PetscInt maxSteps, PetscInt Iter0, PetscScalar time0, PetscScalar deltaTClock, TMMState state, PetscBool doOutput, PetscBool runExternalModel)
{

  PetscErrorCode ierr;

  PetscScalar tc, tf;
  PetscInt iLoop, Iterc;
  
  FILE *fp;
  PetscViewer fd;
  PetscInt itr;
  PetscBool isExternalCommand = PETSC_FALSE;
  char extcommand[PETSC_MAX_PATH_LEN]; 
  PetscInt numTracers;
  const char *prefix;
  
  PetscFunctionBeginUser;

  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;
  
  if (runExternalModel) {
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"trini.petsc",FILE_MODE_WRITE,&fd);CHKERRQ(ierr);  
	for (itr=0; itr<numTracers; itr++) {
	  ierr = VecView(state->c[itr],fd);CHKERRQ(ierr);
	}
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	ierr = PetscOptionsGetString(NULL,prefix,"-external_model_command",extcommand,sizeof(extcommand),&isExternalCommand);CHKERRQ(ierr);
	if (!isExternalCommand) {
      SETERRQ(PETSC_COMM_WORLD,1,"Must provide external model command with the -external_model_command option!");	
	}

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Running external model: %s\n",extcommand);CHKERRQ(ierr);    
	ierr = PetscPOpen(PETSC_COMM_WORLD,PETSC_NULL,extcommand,"w",&fp);CHKERRQ(ierr);  
// 	ierr = PetscPOpen(PETSC_COMM_WORLD,PETSC_NULL,"/Volumes/data2/spk/UVic_OSU_Matrix/LGM_WindPerturbation_Experiments/no_embm_awind2/picdefault/SpinupC14/SNES/run_time_stepper trini.petsc trend.petsc","w",&fp);CHKERRQ(ierr);  
	ierr = PetscPClose(PETSC_COMM_WORLD,fp);CHKERRQ(ierr);  
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished running external model\n");CHKERRQ(ierr);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"trend.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	for (itr=0; itr<numTracers; itr++) {
	  ierr = VecLoad(state->c[itr],fd);CHKERRQ(ierr);
	}
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);          
  } else {
// 	 if (reinitializeForcing) {
// 	   iLoop = 1;
// 	   tc=time0 + (deltaTClock)*(iLoop-1);  /*  current time (time at beginning of step) */
// 	   tf=time0 + (deltaTClock)*iLoop;  /*  future time (time at end of step) */
// 	   Iterc=Iter0+iLoop-1;
// 	   TMMForcingReinitialize(tc,tf,Iterc,iLoop,state);
// 	 }
 /*   Start time stepping loop */
	 for (iLoop = 1; iLoop <= maxSteps; iLoop++) {
   /*  iLoop -> iLoop+1 (convention) */  
	   tc=time0 + (deltaTClock)*(iLoop-1);  /*  current time (time at beginning of step) */
	   tf=time0 + (deltaTClock)*iLoop;  /*  future time (time at end of step) */
	   Iterc=Iter0+iLoop-1;

	   ierr = TMMUpdateTMs(tc);
	   ierr = TMMForcingUpdate(tc,tf,Iterc,iLoop,state);

	   ierr = TMMTimeStep(tc,tf,Iterc,iLoop,state);

	   tc=time0 + (deltaTClock)*iLoop; /*  time at end of step */

	   ierr = TMMTimeStepPost(tc,Iterc,iLoop,state);
	   if (doOutput) {
		 ierr = TMMOutput(tc,Iterc,iLoop,state);
	   }
     } /* end of time-stepping loop */

// 	 if (reinitializeForcing) {
// 	   iLoop = 1;
// 	   tc=time0 + (deltaTClock)*(iLoop-1);  /*  current time (time at beginning of step) */
// 	   tf=time0 + (deltaTClock)*iLoop;  /*  future time (time at end of step) */
// 	   Iterc=Iter0+iLoop-1;
// 	   doTMMReinitializeForcing(tc,tf,Iterc,iLoop,state);
// 	 }
  }
  
  PetscFunctionReturn(0);
}

PetscErrorCode checkEquilibrium(PetscInt maxSteps, PetscInt Iter0, PetscScalar time0, PetscScalar deltaTClock, TMMState state, PetscBool doOutput, PetscBool runExternalModel, PetscReal *nm)
{

  PetscErrorCode ierr;
  PetscInt itr;
  PetscScalar minusone=-1.0;
  Vec *vini;
  
  PetscFunctionBeginUser;
  
  ierr = VecDuplicateVecs(state->c[0],state->numTracers,&vini);CHKERRQ(ierr);
  for (itr=0; itr<state->numTracers; itr++) {
    ierr = VecCopy(state->c[itr],vini[itr]);CHKERRQ(ierr);
  }  

  ierr = integrateOnePeriod(maxSteps, Iter0, time0, deltaTClock, state, doOutput, runExternalModel);CHKERRQ(ierr);

  for (itr=0; itr<state->numTracers; itr++) {
    ierr = VecAXPY(vini[itr],minusone,state->c[itr]);CHKERRQ(ierr);
    ierr = VecNorm(vini[itr],NORM_2,&nm[itr]);
  }  

// Restore initial state
//   for (itr=0; itr<numTracers; itr++) {
//     ierr = VecCopy(vini[itr],state->c[itr]);CHKERRQ(ierr);
//   }  

  ierr = VecDestroyVecs(state->numTracers,&vini);CHKERRQ(ierr);

  PetscFunctionReturn(0);
}

PetscErrorCode runOnExternalSignal(PetscInt maxSteps, PetscInt Iter0, PetscScalar time0, PetscScalar deltaTClock, TMMState state)
{

  PetscErrorCode ierr;
  PetscInt itr, maxValsToRead;
  char *iniFile[100]; 
  char *outFile[100]; 
  char pickupFile[PETSC_MAX_PATH_LEN];
  char pickupOutFile[PETSC_MAX_PATH_LEN];
  PetscBool pickupFromFile = PETSC_FALSE;
  PetscBool readFromIniFiles = PETSC_FALSE;
  PetscBool writeToPickupFile = PETSC_FALSE;  
  PetscBool flg;
  PetscInt extSignal;
  PetscInt waitTime = 5;
  PetscScalar zero = 0.0;
  PetscViewer fd;
  char command[PETSC_MAX_PATH_LEN], extcommand[PETSC_MAX_PATH_LEN]; 
  PetscMPIInt rank;
  PetscInt numTracers;
  const char *prefix;
  
  PetscFunctionBeginUser;

  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

//   ierr = PetscOptionsGetString(NULL,"run_on_external_signal_","-extcommand",extcommand,sizeof(extcommand),&flg);CHKERRQ(ierr);
//   ierr = PetscPrintf(PETSC_COMM_WORLD,"External command to be executed: %s\n",extcommand);CHKERRQ(ierr);

  ierr = PetscOptionsHasName(NULL,prefix,"-reinitialize_forcing",&reinitializeForcing);CHKERRQ(ierr);
  if (reinitializeForcing) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Forcing will be reinitialized after each integration\n");CHKERRQ(ierr);
  }  

/* Initial condition     */
  ierr = PetscOptionsGetString(NULL,prefix,"-run_on_external_signal_pickup",pickupFile,PETSC_MAX_PATH_LEN-1,&pickupFromFile);CHKERRQ(ierr);
  if (pickupFromFile) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Pickup file has been specified for run on external signal\n");CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"  Initial conditions will be read from %s\n", pickupFile);CHKERRQ(ierr);
  } else {
	for (itr=0; itr<numTracers; itr++) {
	  iniFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
	}
    maxValsToRead = numTracers;
    ierr = PetscOptionsGetStringArray(NULL,prefix,"-run_on_external_signal_i",iniFile,&maxValsToRead,&readFromIniFiles);CHKERRQ(ierr);
    if (readFromIniFiles) {  /* read from file */
      if (maxValsToRead != numTracers) {
        SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of input file names specified");
      }
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Initial condition files have been specified for run on external signal\n");CHKERRQ(ierr);
    } else {  /* set to zero */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Initial condition will be set to zero for run on external signal\n");CHKERRQ(ierr);
    }
  }

  for (itr=0; itr<numTracers; itr++) {
    outFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
  }
  maxValsToRead = numTracers;
  ierr = PetscOptionsGetStringArray(NULL,prefix,"-run_on_external_signal_o",outFile,&maxValsToRead,&flg);CHKERRQ(ierr);
  if (flg) {
	if (maxValsToRead != numTracers) {
	  SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of output file names specified");
    }  
  } else {
    ierr = PetscOptionsGetString(NULL,prefix,"-run_on_external_signal_pickup_out",pickupOutFile,PETSC_MAX_PATH_LEN-1,&writeToPickupFile);CHKERRQ(ierr);
    if (!writeToPickupFile) {
      SETERRQ(PETSC_COMM_WORLD,1,"Must indicate one of output file name(s) or pickup out name!");
    }
  }  

//   strcpy(command,"");
//   sprintf(command,"%s %d",extcommand,1);
//   ierr = PetscPrintf(PETSC_COMM_WORLD,"Command: %s\n",command);CHKERRQ(ierr);    
//   extSignal = getExternalSignal(command);
  extSignal = getExternalSignalFile(1,waitTime);
  while (extSignal != 0) {
//     PetscPrintf(PETSC_COMM_WORLD,"extSignal=%d\n",extSignal);CHKERRQ(ierr);
//     PetscSynchronizedPrintf(PETSC_COMM_WORLD,"rank %d: extSignal=%d\n",rank,extSignal);
// 	PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
	MPI_Barrier(PETSC_COMM_WORLD);
    
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading initial conditions for run on external signal\n");CHKERRQ(ierr);

	if (pickupFromFile) {
// 	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Pickup file has been specified\n");CHKERRQ(ierr);
// 	  ierr = PetscPrintf(PETSC_COMM_WORLD,"  Reading initial conditions from %s\n", pickupFile);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,pickupFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  for (itr=0; itr<numTracers; itr++) {
		ierr = VecLoad(state->c[itr],fd);CHKERRQ(ierr); /* IntoVector */
	  }
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);          
	} else if (readFromIniFiles) {
	  for (itr=0; itr<numTracers; itr++) {
// 		ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading initial condition from file %s\n", itr,iniFile[itr]);CHKERRQ(ierr);
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,iniFile[itr],FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = VecLoad(state->c[itr],fd);CHKERRQ(ierr);/* IntoVector */
		ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);          
	  }  
	} else {  /* set to zero */
// 	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Setting initial condition to zero\n");CHKERRQ(ierr);
	  for (itr=0; itr<numTracers; itr++) {
		VecSet(state->c[itr],zero);
	  }        
	}

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Running model for run on external signal\n");CHKERRQ(ierr);

// ??? Not quite sure where to call this from: here or from integrateOnePeriod?
	if (reinitializeForcing) {
	  TMMForcingReinitialize(time0,time0+deltaTClock,Iter0,1,state);
	}

    ierr = integrateOnePeriod(maxSteps, Iter0, time0, deltaTClock, state, PETSC_FALSE, PETSC_FALSE);CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished running model for run on external signal\n");CHKERRQ(ierr);

    ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing output for run on external signal\n");CHKERRQ(ierr);

    if (writeToPickupFile) {
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,pickupOutFile,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
	  for (itr=0; itr<numTracers; itr++) {
		ierr = VecView(state->c[itr],fd);CHKERRQ(ierr);
	  }
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);          
    } else {
	  for (itr=0; itr<numTracers; itr++) {       
		ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,outFile[itr],FILE_MODE_WRITE,&fd);CHKERRQ(ierr);
		ierr = VecView(state->c[itr],fd);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	  }
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Sending finished signal for run on external signal\n");CHKERRQ(ierr);

//     strcpy(command,"");
// 	sprintf(command,"%s %d",extcommand,2);
//     ierr = PetscPrintf(PETSC_COMM_WORLD,"Command: %s\n",command);CHKERRQ(ierr);	
// 	extSignal = getExternalSignal(command);
	extSignal = getExternalSignalFile(2,waitTime);
//  we ignore the previous signal
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Sending waiting to run signal for run on external signal\n");CHKERRQ(ierr);
// 	strcpy(command,"");
// 	sprintf(command,"%s %d",extcommand,1);
// 	ierr = PetscPrintf(PETSC_COMM_WORLD,"Command: %s\n",command);CHKERRQ(ierr);    
// 	extSignal = getExternalSignal(command);
	extSignal = getExternalSignalFile(1,waitTime);
  }


  PetscFunctionReturn(0);
}

PetscErrorCode spinupMonitor(SNES snes,PetscInt its,PetscReal fnorm,void *ptr)
{

  PetscErrorCode ierr;
TMMSPINUP         *user = (TMMSPINUP*)ptr;
PetscInt itr;
PetscReal nm, totalnm;
  Vec F;
  Vec X;
SNESConvergedReason reason;
  
// NGMRES debug
//   SNES_NGMRES    *ngmres = (SNES_NGMRES*) snes->data;
//   PetscInt       i;
//   Vec            *Fdot      = ngmres->Fdot;
//   Vec            *Xdot      = ngmres->Xdot;
//   PetscScalar    *beta      = ngmres->beta;
//   PetscScalar    *xi        = ngmres->xi;
//   PetscInt       msize      = ngmres->msize;
//   
//   PetscScalar    alph_total = 0.;
  PetscViewer fd;
  char tmpFile[PETSC_MAX_PATH_LEN];

  PetscFunctionBeginUser;

// NGMRES debug
//   for (i=0; i<msize; i++) {
//     ierr = PetscPrintf(PETSC_COMM_WORLD,"Beta[%d] at SNES iteration %d: %g\n", i, its, beta[i]); CHKERRQ(ierr);
//     alph_total += beta[i];
//   }   
//   ierr = PetscPrintf(PETSC_COMM_WORLD,"Total alpha: %g\n", alph_total); CHKERRQ(ierr);
// 
//   strcpy(tmpFile,"");
//   sprintf(tmpFile,"ngmres_solutions_%d.petsc",its);
//   ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);  
//   for (i=0; i<msize; i++) {
//     ierr = XtoVconvert(user->numTracers, user->state->c, Xdot[i], user->XTovScaleFac);CHKERRQ(ierr);
//     for (itr=0; itr<user->numTracers; itr++) {
//   	  ierr = VecView(user->state->c[itr],fd);CHKERRQ(ierr);
//     }
//   }   
//   ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      


SNESGetFunction(snes,&F,0,0);
ierr = XtoVconvert(user->state->numTracers, user->state->c, F, user->XTovScaleFac);CHKERRQ(ierr);
totalnm=0.0;
  ierr = PetscFPrintf(PETSC_COMM_WORLD,user->logfp,"%d  %d  ",its, user->itf);CHKERRQ(ierr);
  for (itr=0; itr<user->state->numTracers; itr++) {
    ierr = VecNorm(user->state->c[itr],NORM_2,&nm);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm for tracer %d at SNES iteration %d, function eval %d: %g\n", itr, its, user->itf, nm); CHKERRQ(ierr);
	totalnm = totalnm + pow(nm,2.0);
	ierr = PetscFPrintf(PETSC_COMM_WORLD,user->logfp,"%g  ",nm);CHKERRQ(ierr);
  }
  totalnm = sqrt(totalnm);
// Note: fnorm is the norm computed by SNES and could be different if tracer scale factors are not 1
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Total norm at SNES iteration %d, function eval %d: %g\n", its, user->itf, totalnm); CHKERRQ(ierr);
  ierr = PetscFPrintf(PETSC_COMM_WORLD,user->logfp,"%g  %g\n",totalnm,fnorm);CHKERRQ(ierr);

  if ((user->checkpointFreq)>0) {
    ierr = SNESGetConvergedReason(snes,&reason);CHKERRQ(ierr);  // I think this doesn't help as it will always be SNES_CONVERGED_ITERATING
//  Write checkpoint if either a multiple of frequency or SNESSolve has stopped (either converged or diverged)
    if (((its % (user->checkpointFreq)) ==  0) || (reason != SNES_CONVERGED_ITERATING)) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing out solution at SNES iteration %d, function evals %d\n", its,user->itf); CHKERRQ(ierr);
      SNESGetSolution(snes,&X);
      ierr = XtoVconvert(user->state->numTracers, user->state->c, X, user->XTovScaleFac);CHKERRQ(ierr);
	  strcpy(tmpFile,"");
	  sprintf(tmpFile,"snes_spinup_checkpoint_%d.petsc",its);
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_WRITE,&fd);CHKERRQ(ierr);  
	  for (itr=0; itr<user->state->numTracers; itr++) {
		ierr = VecView(user->state->c[itr],fd);CHKERRQ(ierr);
	  }
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
 	}  
  }  
  
PetscFunctionReturn(0);
}

PetscErrorCode checkConvergence(SNES snes,PetscInt its,PetscReal xnorm,PetscReal snorm,PetscReal fnorm,SNESConvergedReason *reason,void *ptr)
{

  PetscErrorCode ierr;
  TMMSPINUP         *user = (TMMSPINUP*)ptr;

  PetscViewer fd;
  int fp;
  PetscInt converged;
  
  PetscFunctionBeginUser;

ierr = SNESConvergedDefault(snes,its,xnorm,snorm,fnorm,reason,ptr);CHKERRQ(ierr);

// Note: this check uses information from the last function evaluation (model run). Depending on how SNES 
// is doing thing, this may or may not coincide with the solution being passed to this routine, potentially 
// introducing a small offset. I think it will coincide as SNES uses a line search to calculate the step size. 

	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"convergence.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	  ierr = PetscBinaryRead(fp,&converged,1,NULL,PETSC_INT);CHKERRQ(ierr);  
// 	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of points in emission files is %d \n",numEmission_hist);CHKERRQ(ierr);  
//       ierr = PetscMalloc(numEmission_hist*sizeof(PetscScalar),&Tem_hist);CHKERRQ(ierr); 
//       ierr = PetscBinaryRead(fp,Tem_hist,numEmission_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

// Overwrite default convergence reason
//==(user->convergedValue)

if (converged==0) {
  if (((*reason) > 0) & ((*reason) != SNES_CONVERGED_ITS)) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Overwriting default convergence reason as external convergence criteria have not been met: %d\n",converged);CHKERRQ(ierr);    
    *reason = SNES_CONVERGED_ITERATING;
  }  
} else if (converged>0) {
  *reason = SNES_CONVERGED_FNORM_ABS;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"External convergence criteria have been met: %d\n",converged);CHKERRQ(ierr);    
} else {
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Divergence indicated by external convergence check: %d\n",converged);CHKERRQ(ierr);    
  *reason = SNES_DIVERGED_FNORM_NAN;
}

// if (converged>0) {
//   *reason = SNES_CONVERGED_FNORM_ABS;
//   ierr = PetscPrintf(PETSC_COMM_WORLD,"External convergence criteria have been met: %d\n",converged);CHKERRQ(ierr);    
// } else {
//   if ((converged==0) & ((*reason) > 0) & ((*reason) != SNES_CONVERGED_ITS)) {
//     ierr = PetscPrintf(PETSC_COMM_WORLD,"Overwriting default convergence reason as external convergence criteria have not been met: %d\n",converged);CHKERRQ(ierr);    
//     *reason = SNES_CONVERGED_ITERATING;
//   }  
// }

  PetscFunctionReturn(0);
}

// PetscErrorCode ViewData(SNES snes, void *ptr)
// {
// 
//   PetscErrorCode ierr;
// 
//   SNES_NGMRES    *ngmres = (SNES_NGMRES*) snes->data;
//   PetscInt       i;
//   Vec            *Fdot      = ngmres->Fdot;
//   Vec            *Xdot      = ngmres->Xdot;
//   PetscScalar    *beta      = ngmres->beta;
//   PetscScalar    *xi        = ngmres->xi;
//   PetscInt       msize      = ngmres->msize;
// 
//   TMMSPINUP         *user = (TMMSPINUP*)ptr;
//   
//   PetscFunctionBeginUser;
// 
//   for (i=0; i<msize; i++) {
//     ierr = XtoVconvert(user->numTracers, user->state->c, Xdot[i], user->XTovScaleFac);CHKERRQ(ierr);
//     for (itr=0; itr<user->numTracers; itr++) {
//       VecView(user->state->c[itr]);
//     }
//   }   
// 
// PetscFunctionReturn(0);
// }
