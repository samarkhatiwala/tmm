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

#undef __FUNCT__
#define __FUNCT__ "TMMForcingInitialize"
PetscErrorCode TMMForcingInitialize(TMMState state)
{


  PetscErrorCode ierr;
  PetscBool flg;
  PetscInt itr, maxValsToRead;
  char tmpFile[PETSC_MAX_PATH_LEN];
  PetscInt il;
  PetscViewer fd;
  PetscScalar zero = 0.0;
  PetscMPIInt numProcessors;  
// These static variables are there to prevent recalculation of global variables such as 
// bcTemplateVec, gBCIndices, etc. Clunky but seems to work.
  static PetscBool firstTimeBC = PETSC_TRUE;
  static PetscBool firstTimeForcing = PETSC_TRUE;
  const char *prefix;
  
  prefix = ((PetscObject)state)->prefix;
  
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&numProcessors);CHKERRQ(ierr);  

/* Forcing/RHS   */
/* The tracer(s) can be forced in 3 ways (any combination of which can be turned on): */
/* 1) Forcing term read from file (can be periodic, constant, or time-dependent) */
/* 2) External forcing computed in S/R calcExternalForcing */
/* 3) Prescribed boundary condition (can be periodic, constant, or time-dependent) */
  ierr = PetscOptionsHasName(NULL,prefix,"-forcing_from_file",&state->useForcingFromFile);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,prefix,"-prescribed_bc",&state->usePrescribedBC);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,prefix,"-external_forcing",&state->useExternalForcing);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,prefix,"-relax_tracer",&state->relaxTracer);CHKERRQ(ierr);

  if (state->useForcingFromFile) {  
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Forcing from file(s) specified\n");CHKERRQ(ierr);  
    for (itr=0; itr<state->numTracers; itr++) {
      state->forcingFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
    }
    maxValsToRead = state->numTracers;
    ierr = PetscOptionsGetStringArray(NULL,prefix,"-forcing_files",state->forcingFile,&maxValsToRead,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"MUST specify forcing files with the -forcing_files option");
     if (maxValsToRead != state->numTracers) {
       SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of forcing file names specified");
     }
     ierr = VecDuplicateVecs(templateVec,state->numTracers,&state->qf);CHKERRQ(ierr);
/*   There are 3 possibilities: periodic, constant, and time-dependent forcing */
     ierr = PetscOptionsHasName(NULL,prefix,"-periodic_forcing",&state->periodicForcing);CHKERRQ(ierr);
    	ierr = PetscOptionsHasName(NULL,prefix,"-time_dependent_forcing",&state->timeDependentForcing);CHKERRQ(ierr);
     if ((state->periodicForcing) && (state->timeDependentForcing)) SETERRQ(PETSC_COMM_WORLD,1,"Cannot specify both periodicForcing and timeDependentForcing");
     if (state->periodicForcing) {
       ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic forcing from file(s) specified\n");CHKERRQ(ierr);
/*     read time data */
       ierr = PeriodicTimerCreate(&state->forcingTimer);CHKERRQ(ierr);
       ierr = PeriodicTimerIni("forcing_", prefix, NULL, state->forcingTimer);CHKERRQ(ierr);
/*     Forcing is read in interpPeriodicForcing */
       for (itr=0; itr<state->numTracers; itr++) {
         ierr=PetscStrcat(state->forcingFile[itr],"_");CHKERRQ(ierr);        
         ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d periodic forcing basename is %s\n",itr,state->forcingFile[itr]);CHKERRQ(ierr); 
         ierr = PeriodicVecCreate(&state->qp[itr]);CHKERRQ(ierr);
//         state->qp[itr].firstTime = PETSC_TRUE; /* initialize periodic vector */        	    
       }
     } else if (state->timeDependentForcing) {
       ierr=PetscPrintf(PETSC_COMM_WORLD,"Time dependent forcing specified\n");CHKERRQ(ierr);
/*     read time data */
       ierr = TimeDependentTimerCreate(&state->forcingTimeDependentTimer);CHKERRQ(ierr);
       ierr = TimeDependentTimerIni("forcing_", prefix, NULL, state->forcingTimeDependentTimer);CHKERRQ(ierr);
/*     Forcing is read in TimeDependentVecInterp */
       for (itr=0; itr<state->numTracers; itr++) {   
         ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: forcing will be read from file %s\n", itr,state->forcingFile[itr]);CHKERRQ(ierr);	  
         ierr = TimeDependentVecCreate(&state->qtdf[itr]);CHKERRQ(ierr);
// 		state->qtdf[itr].firstTime = PETSC_TRUE;
       }
	    } else { /* constant forcing */
       ierr=PetscPrintf(PETSC_COMM_WORLD,"Constant forcing specified\n");CHKERRQ(ierr);
       state->constantForcing = PETSC_TRUE;
       for (itr=0; itr<state->numTracers; itr++) {   
         ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading forcing from file %s\n", itr,state->forcingFile[itr]);CHKERRQ(ierr);
         ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->forcingFile[itr],FILE_MODE_READ,&fd);CHKERRQ(ierr);
         ierr = VecLoad(state->qf[itr],fd);CHKERRQ(ierr); /* IntoVector */
         ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
       }
    }

    state->applyForcingFromFile = PETSC_FALSE; // we switch this on in TMMForcingUpdate
    state->forcingFromFileStartStep = Iter0 + 1; // default is to start at the first time step

    ierr = PetscOptionsGetInt(NULL,prefix,"-forcing_from_file_start_step",&state->forcingFromFileStartStep,&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Forcing from file will be switched on at (absolute) time step %d\n",state->forcingFromFileStartStep);CHKERRQ(ierr);
    }
    
    ierr = PetscOptionsGetInt(NULL,prefix,"-forcing_from_file_cutoff_step",&state->forcingFromFileCutOffStep,&flg);CHKERRQ(ierr);
    if (state->forcingFromFileCutOffStep>0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Forcing from file will be turned off after (absolute) time step %d\n",state->forcingFromFileCutOffStep);CHKERRQ(ierr);
    }
    
  } else {  /* no forcing from file */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"No forcing from file(s) specified\n");CHKERRQ(ierr);
  }

  if (state->useExternalForcing) {  /* external forcing present */  
    ierr = PetscPrintf(PETSC_COMM_WORLD,"External forcing is being used\n");CHKERRQ(ierr);
    ierr = VecDuplicateVecs(templateVec,state->numTracers,&state->qef);CHKERRQ(ierr);
// I'm moving this to update because otherwise the functions need to be specified before calling TMMSetFromOptions    
//     ierr = TMMComputeExtForcFunction(time0,Iter0,-1,state,TMM_INI_FUNC);
//     ierr = iniExternalForcing(time0,Iter0,state,NULL);CHKERRQ(ierr);
//     ierr = iniExternalForcing(time0,Iter0,state->numTracers,state->c,state->qef);CHKERRQ(ierr);

    state->applyExternalForcing = PETSC_FALSE; // we switch this on in TMMForcingUpdate
    state->externalForcingStartStep = Iter0 + 1; // default is to start at the first time step

    ierr = PetscOptionsGetInt(NULL,prefix,"-external_forcing_start_step",&state->externalForcingStartStep,&flg);CHKERRQ(ierr);
    if (flg) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"External forcing will be turned on at (absolute) time step %d\n",state->externalForcingStartStep);CHKERRQ(ierr);
    }

    ierr = PetscOptionsGetInt(NULL,prefix,"-external_forcing_cutoff_step",&state->externalForcingCutOffStep,&flg);CHKERRQ(ierr);
    if (state->externalForcingCutOffStep>0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"External forcing will be turned off after (absolute) time step %d\n",state->externalForcingCutOffStep);CHKERRQ(ierr);
    }
    
  } else { /* no external forcing */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"No external forcing is being used\n");CHKERRQ(ierr);  
  }  

  if (state->relaxTracer) {
  	 ierr=PetscPrintf(PETSC_COMM_WORLD,"Tracer relaxation specified\n");CHKERRQ(ierr);  
    
    ierr = PetscMalloc(state->numTracers*sizeof(PetscScalar),&state->relaxTracerValue);CHKERRQ(ierr);   
    ierr = PetscMalloc(state->numTracers*sizeof(PetscScalar),&state->relaxTracerLambda);CHKERRQ(ierr);   

    maxValsToRead = state->numTracers;
    ierr = PetscOptionsGetRealArray(NULL,prefix,"-relax_tracer_lambda",state->relaxTracerLambda,&maxValsToRead,&flg);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate tracer relaxation lambda with the -relax_tracer_lambda option");
    if (maxValsToRead != state->numTracers) {
      SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of relaxation lambda values specified");
    }

    maxValsToRead = state->numTracers;
    ierr = PetscOptionsGetRealArray(NULL,prefix,"-relax_tracer_value",state->relaxTracerValue,&maxValsToRead,&flg);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate relaxation values with the -relax_tracer_value option");
    if (maxValsToRead != state->numTracers) {
      SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of tracer relaxation values specified");
    }

    ierr = VecDuplicateVecs(templateVec,state->numTracers,&state->qrel);CHKERRQ(ierr);
    
    for (itr=0; itr<state->numTracers; itr++) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d relaxation lambda=%15.11f, relaxation value=%10.8f\n",itr,state->relaxTracerLambda[itr],state->relaxTracerValue[itr]);CHKERRQ(ierr);
    }
  }

/* Prescribed BCs   */
  if (state->usePrescribedBC) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Prescribed BCs specified\n");CHKERRQ(ierr);

    if (firstTimeBC) {
      prescribedBCInUse = PETSC_TRUE;
      ierr = VecCreate(PETSC_COMM_WORLD,&bcTemplateVec);CHKERRQ(ierr);
   //    We assign sizes below depending on whether we are using profiles or not	  
    }  

//  There are two cases: calculated BCs and BCs from file    
    ierr = PetscOptionsHasName(NULL,prefix,"-calc_bc",&state->doCalcBC);CHKERRQ(ierr);
    if (state->doCalcBC) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"BCs will be calculated\n");CHKERRQ(ierr);

      if (firstTimeBC) {
        calcBCInUse = PETSC_TRUE;
        if ((useProfiles) && (numProcessors>1)) {
          lBCSize = lNumProfiles;
          ierr = VecSetSizes(bcTemplateVec,lBCSize,PETSC_DECIDE);CHKERRQ(ierr);
          ierr = VecSetFromOptions(bcTemplateVec);CHKERRQ(ierr);
          ierr = VecGetSize(bcTemplateVec,&gBCSize);CHKERRQ(ierr);
        } else {
      //        There seems to be no way of setting this information automatically short of reading Be or Bi		
          ierr = PetscOptionsGetInt(NULL,NULL,"-bc_vec_size",&gBCSize,&flg);CHKERRQ(ierr);
          if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate size of BC vector with the -bc_vec_size option");
          ierr = VecSetSizes(bcTemplateVec,PETSC_DECIDE,gBCSize);CHKERRQ(ierr);
          ierr = VecSetFromOptions(bcTemplateVec);CHKERRQ(ierr);
          ierr = VecGetLocalSize(bcTemplateVec,&lBCSize);CHKERRQ(ierr);    
        }

//      gBCSize and lBCSize are now known when doCalcBC
      }

      ierr = VecDuplicateVecs(bcTemplateVec,state->numTracers,&state->cbc);CHKERRQ(ierr);
      ierr = VecDuplicateVecs(bcTemplateVec,state->numTracers,&state->cbf);CHKERRQ(ierr);

// I'm moving this to update because otherwise the functions need to be specified before calling TMMSetFromOptions    
//       ierr = TMMComputeCalcBCFunction(time0,Iter0,time0+deltaTClock,Iter0+1,-1,state,TMM_INI_FUNC);
//       ierr = iniCalcBC(time0,Iter0,time0+deltaTClock,Iter0+1,state,NULL);CHKERRQ(ierr);
      
    } else { /* read from file */
      for (itr=0; itr<state->numTracers; itr++) {
        state->bcFile[itr] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
      }
      maxValsToRead = state->numTracers;
      ierr = PetscOptionsGetStringArray(NULL,prefix,"-bc_files",state->bcFile,&maxValsToRead,&flg);CHKERRQ(ierr);
      if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"MUST specify BC files with the -bc_files option");
      if (maxValsToRead != state->numTracers) {
        SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of BC file names specified");
      }

/*    There are 3 possibilities: periodic, constant, and time-dependent BCs */
      ierr = PetscOptionsHasName(NULL,prefix,"-periodic_bc",&state->periodicBC);CHKERRQ(ierr);
      ierr = PetscOptionsHasName(NULL,prefix,"-time_dependent_bc",&state->timeDependentBC);CHKERRQ(ierr);
      if ((state->periodicBC) && (state->timeDependentBC)) SETERRQ(PETSC_COMM_WORLD,1,"Cannot specify both periodicBC and timeDependentBC");
      if (state->periodicBC) {
        ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic BC from file(s) specified\n");CHKERRQ(ierr);
 /*      read time data */
        ierr = PeriodicTimerCreate(&state->bcTimer);CHKERRQ(ierr);
        ierr = PeriodicTimerIni("bc_", prefix, NULL, state->bcTimer);CHKERRQ(ierr);
        for (itr=0; itr<state->numTracers; itr++) {
          ierr=PetscStrcat(state->bcFile[itr],"_");CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: periodic BC basename is %s\n",itr,state->bcFile[itr]);CHKERRQ(ierr); 
          ierr = PeriodicVecCreate(&state->cbp[itr]);CHKERRQ(ierr);
 //           state->cbp[itr].firstTime = PETSC_TRUE; /* initialize periodic vector */
        }
        if (firstTimeBC) {
 /*        Load one vector here as a template; BC is read in PeriodicVecInterp */
          strcpy(tmpFile,"");
          sprintf(tmpFile,"%s%02d",state->bcFile[0],0);
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
          ierr = VecLoad(bcTemplateVec,fd);CHKERRQ(ierr);
          ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
        }  
        ierr = VecDuplicateVecs(bcTemplateVec,state->numTracers,&state->cbc);CHKERRQ(ierr);
        ierr = VecDuplicateVecs(bcTemplateVec,state->numTracers,&state->cbf);CHKERRQ(ierr);
      } else if (state->timeDependentBC) {
        ierr=PetscPrintf(PETSC_COMM_WORLD,"Time dependent BC specified\n");CHKERRQ(ierr);
/*      read time data */
        ierr = TimeDependentTimerCreate(&state->bcTimeDependentTimer);CHKERRQ(ierr);
        ierr = TimeDependentTimerIni("bc_", prefix, NULL, state->bcTimeDependentTimer);CHKERRQ(ierr);
        for (itr=0; itr<state->numTracers; itr++) {   
                ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: BC will be read from %s\n",itr,state->bcFile[itr]);CHKERRQ(ierr); 
                ierr = TimeDependentVecCreate(&state->cbtd[itr]);CHKERRQ(ierr);
      // 		  state->cbtd[itr].firstTime = PETSC_TRUE;
        }
        if (firstTimeBC) {
      /*        Load one vector here as a template; BC is read in TimeDependentVectorCreate */
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->bcFile[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
          ierr = VecLoad(bcTemplateVec,fd);CHKERRQ(ierr);
          ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
        }  
        ierr = VecDuplicateVecs(bcTemplateVec,state->numTracers,&state->cbc);CHKERRQ(ierr);    
        ierr = VecDuplicateVecs(bcTemplateVec,state->numTracers,&state->cbf);CHKERRQ(ierr); 
		
	     } else { /* constant BC */
        ierr=PetscPrintf(PETSC_COMM_WORLD,"Constant BC specified\n");CHKERRQ(ierr);
        state->constantBC = PETSC_TRUE;
        if (firstTimeBC) {
      /*        Load one vector here as a template */
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->bcFile[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
          ierr = VecLoad(bcTemplateVec,fd);CHKERRQ(ierr);
          ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
        }  
        ierr = VecDuplicateVecs(bcTemplateVec,state->numTracers,&state->cbc);CHKERRQ(ierr);
        ierr = VecDuplicateVecs(bcTemplateVec,state->numTracers,&state->cbf);CHKERRQ(ierr);
/*      Load BCs */
        for (itr=0; itr<state->numTracers; itr++) {   
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: reading BC from file %s\n", itr,state->bcFile[itr]);CHKERRQ(ierr);
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->bcFile[itr],FILE_MODE_READ,&fd);CHKERRQ(ierr);
          ierr = VecLoad(state->cbc[itr],fd);CHKERRQ(ierr); /* IntoVector */
          ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
          ierr = VecCopy(state->cbc[itr],state->cbf[itr]);CHKERRQ(ierr);		  
        }
      }
    }

    if (firstTimeBC) {
      /* These are not known for the case of BC from file */
      ierr = VecGetSize(bcTemplateVec,&gBCSize);CHKERRQ(ierr);    
      ierr = VecGetLocalSize(bcTemplateVec,&lBCSize);CHKERRQ(ierr);

//    We only need the following information when BCs are calculated but do it anyway
      ierr = VecGetOwnershipRange(bcTemplateVec,&gbcLow,&gbcHigh);CHKERRQ(ierr);
      gbcHigh = gbcHigh - 1; /* Note: gbcHigh is one more than the last local element */
      ierr = PetscMalloc(lBCSize*sizeof(PetscInt),&gBCIndices);CHKERRQ(ierr);  
      for (il=0; il<lBCSize; il++) {
        gBCIndices[il] = il + gbcLow;
      }

      if ((useProfiles) && (numProcessors>1)) {
        if (lBCSize != lNumProfiles) {
          SETERRQ(PETSC_COMM_WORLD,1,"Problem with partitioning of BC vectors! lNumProfiles must equal lBCSize");
        }
      }
    }
        
    state->applyBC = PETSC_FALSE; // we switch this on in TMMForcingUpdate
    state->bcStartStep = Iter0 + 1; // default is to start at the first time step

    ierr = PetscOptionsGetInt(NULL,prefix,"-bc_start_step",&state->bcStartStep,&flg);CHKERRQ(ierr);
  	 ierr = PetscPrintf(PETSC_COMM_WORLD,"Prescribed BC will be turned on at (absolute) time step %d\n",state->bcStartStep);CHKERRQ(ierr);    

    ierr = PetscOptionsGetInt(NULL,prefix,"-bc_cutoff_step",&state->bcCutOffStep,&flg);CHKERRQ(ierr);
    if (state->bcCutOffStep>0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Prescribed BC will be turned off after (absolute) time step %d\n",state->bcCutOffStep);CHKERRQ(ierr);    
    }

    if (firstTimeBC) {
/*    Matrices */
      ierr = PetscOptionsGetString(NULL,NULL,"-mbe",matbeFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
      if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate binary boundary matrix file name with the -mbe option");
      ierr = PetscOptionsGetString(NULL,NULL,"-mbi",matbiFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
      if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate binary boundary matrix file name with the -mbi options");

      ierr = MatCreate(PETSC_COMM_WORLD,&Be);CHKERRQ(ierr);
      ierr = MatCreate(PETSC_COMM_WORLD,&Bi);CHKERRQ(ierr);	

/*    Set layout information */
//    Why is this in an if statement? Both lSize and lBCSize are known at this point. I'm going to remove the if for now.
// 	  if ((useProfiles) && (numProcessors>1)) {
      ierr = MatSetSizes(Be,lSize,lBCSize,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
      ierr = MatSetSizes(Bi,lSize,lBCSize,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
// 	  }		
	
      ierr = MatSetType(Be,MATMPIAIJ);CHKERRQ(ierr);      
      ierr = MatSetFromOptions(Be);CHKERRQ(ierr);  
      ierr = MatSetType(Bi,MATMPIAIJ);CHKERRQ(ierr);        
      ierr = MatSetFromOptions(Bi);CHKERRQ(ierr);
  
      if (periodicMatrix) {    
        ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic boundary matrices specified\n");CHKERRQ(ierr);
        ierr=PetscStrcat(matbeFile,"_");CHKERRQ(ierr);
        ierr=PetscStrcat(matbiFile,"_");CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Be basename is %s\n", matbeFile);CHKERRQ(ierr); 
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Bi basename is %s\n", matbiFile);CHKERRQ(ierr);     

   /*      Read here to set sparsity pattern */  
        strcpy(tmpFile,"");
        sprintf(tmpFile,"%s%02d",matbeFile,0);
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
        ierr = MatLoad(Be,fd);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

        strcpy(tmpFile,"");
        sprintf(tmpFile,"%s%02d",matbiFile,0);
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
        ierr = MatLoad(Bi,fd);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
 
        Bep.firstTime = PETSC_TRUE;
        Bip.firstTime = PETSC_TRUE;

      } else if (timeDependentMatrix) {
        ierr=PetscPrintf(PETSC_COMM_WORLD,"Time-dependent matrices specified\n");CHKERRQ(ierr);
        ierr=PetscStrcat(matbeFile,"_");CHKERRQ(ierr);
        ierr=PetscStrcat(matbiFile,"_");CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Be basename is %s\n", matbeFile);CHKERRQ(ierr); 
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Bi basename is %s\n", matbiFile);CHKERRQ(ierr);

/*      Read here to set sparsity pattern */  
        strcpy(tmpFile,"");
        sprintf(tmpFile,"%s%02d",matbeFile,0);
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
        ierr = MatLoad(Be,fd);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  
        strcpy(tmpFile,"");
        sprintf(tmpFile,"%s%02d",matbiFile,0);
        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmpFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
        ierr = MatLoad(Bi,fd);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

        Betd.firstTime = PETSC_TRUE;
        Bitd.firstTime = PETSC_TRUE;  
		
      } else { /*  not periodic. read matrices here */

        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,matbeFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading Be from file %s\n", matbeFile);CHKERRQ(ierr);  
        ierr = MatLoad(Be,fd);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

        ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,matbiFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading Bi from file %s\n", matbiFile);CHKERRQ(ierr);  
        ierr = MatLoad(Bi,fd);CHKERRQ(ierr);    
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
      }
	     ierr = VecSet(bcTemplateVec,zero);CHKERRQ(ierr);
	   }
	
    if (firstTimeBC) firstTimeBC = PETSC_FALSE;
	
  } else {  /* no BC forcing */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"No prescribed BCs specified\n");CHKERRQ(ierr);  
  }  

//?????  
// Rfs depends on the physical circulation and should be applied to every state with a relevant forcing term
  if ((state->useForcingFromFile) || (state->useExternalForcing)) {
    if (firstTimeForcing) ierr = PetscOptionsGetString(NULL,NULL,"-rescale_forcing_file",rfsFile,PETSC_MAX_PATH_LEN-1,&rescaleForcing);CHKERRQ(ierr);
   	if (rescaleForcing) {  
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Forcing will be rescaled\n");CHKERRQ(ierr);
      if (firstTimeForcing) {
        ierr = VecDuplicate(templateVec,&Rfs);CHKERRQ(ierr);    
        if (periodicMatrix) {    
          if (state->constantForcing) {
            SETERRQ(PETSC_COMM_WORLD,1,"Periodic rescaling not supported with constant forcing");
          }
          ierr=PetscStrcat(rfsFile,"_");CHKERRQ(ierr);
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Rescale forcing factor file basename is %s\n", rfsFile);CHKERRQ(ierr); 	
          ierr = PeriodicVecCreate(&Rfsp);CHKERRQ(ierr);
        } else if (timeDependentMatrix) {
          if (state->constantForcing) {
            SETERRQ(PETSC_COMM_WORLD,1,"Time dependent rescaling not supported with constant forcing");
          }
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Rescale forcing factor file is %s\n", rfsFile);CHKERRQ(ierr);	
          ierr = TimeDependentVecCreate(&Rfstd);CHKERRQ(ierr);
        } else {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading rescale forcing factor from file %s\n", rfsFile);CHKERRQ(ierr);  	
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,rfsFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
          ierr = VecLoad(Rfs,fd);CHKERRQ(ierr);  
          ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
		      }
	     }
      if ((constantMatrix) & (state->constantForcing)) {
        for (itr=0; itr<state->numTracers; itr++) {
          ierr = VecPointwiseMult(state->qf[itr],Rfs,state->qf[itr]);CHKERRQ(ierr);
        }
      }
      if (firstTimeForcing) firstTimeForcing = PETSC_FALSE;
	   }
  }
  
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TMMForcingUpdate"
PetscErrorCode TMMForcingUpdate(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, TMMState state)
//                              Vec *c, Vec *qf, Vec *qef, Vec *cbc, Vec *cbf)
{

 PetscScalar tf;
	PetscInt itr;
	PetscErrorCode ierr;

 tf = tc + deltaTClock;
 
/*  interpolate/update qf,qef,cbc to current time (tc) and cbf to future time (tf) */

/*  Forcing     */
    if ((state->useForcingFromFile) && ((Iter0+iLoop==state->forcingFromFileStartStep))) {
      if (!state->applyForcingFromFile) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Switching on forcing from file at (absolute) time step %d\n",Iter0+iLoop);CHKERRQ(ierr);    
      }
      state->applyForcingFromFile = PETSC_TRUE;
    }
    
    if (state->applyForcingFromFile) {
      if ((state->forcingFromFileCutOffStep>0) && ((Iter0+iLoop)>=(state->forcingFromFileCutOffStep+1))) { /* Note: the >= deals with the situation when Iter0 is > 0 */
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Switching off forcing from file at (absolute) time step %d\n",Iter0+iLoop);CHKERRQ(ierr);
        state->applyForcingFromFile = PETSC_FALSE;
      } else {
        if (state->periodicForcing) {
          for (itr=0; itr<state->numTracers; itr++) {    
            ierr = PeriodicVecInterp(tc,&state->qf[itr],state->forcingTimer->cyclePeriod,state->forcingTimer->numPerPeriod,state->forcingTimer->tdp,state->qp[itr],state->forcingFile[itr]);
            if (rescaleForcing) {
             ierr = VecPointwiseMult(state->qf[itr],Rfs,state->qf[itr]);CHKERRQ(ierr);
            }
          }    
        } else if (state->timeDependentForcing) {
          for (itr=0; itr<state->numTracers; itr++) {    
            ierr = TimeDependentVecInterp(tc,&state->qf[itr],state->forcingTimeDependentTimer->numTimes,state->forcingTimeDependentTimer->tdt,state->qtdf[itr],state->forcingFile[itr]);
            if (rescaleForcing) {
             ierr = VecPointwiseMult(state->qf[itr],Rfs,state->qf[itr]);CHKERRQ(ierr);
            }
          }
        }
      }	
	   }

    if ((state->useExternalForcing) && ((Iter0+iLoop==state->externalForcingStartStep))) {
      if (!state->applyExternalForcing) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Switching on external forcing at (absolute) time step %d\n",Iter0+iLoop);CHKERRQ(ierr);
      }
      if (!state->isInitializedExternalForcing) {
        ierr = TMMComputeExtForcFunction(time0,Iter0,-1,state,TMM_INI_FUNC);
        state->isInitializedExternalForcing = PETSC_TRUE;
      }
	     state->applyExternalForcing = PETSC_TRUE;
    }

    if (state->applyExternalForcing) {
      if ((state->externalForcingCutOffStep>0) && ((Iter0+iLoop)>=(state->externalForcingCutOffStep+1))) { /* Note: the >= deals with the situation when Iter0 is > 0 */
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Switching off external forcing at (absolute) time step %d\n",Iter0+iLoop);CHKERRQ(ierr);
        state->applyExternalForcing = PETSC_FALSE;
      } else {
//         ierr = PetscPrintf(PETSC_COMM_WORLD,"Calling TMMComputeExtForcFunction-calc\n");CHKERRQ(ierr);        
        ierr = TMMComputeExtForcFunction(tc,Iterc,iLoop,state,TMM_CALC_FUNC);
        if (rescaleForcing) {
          for (itr=0; itr<state->numTracers; itr++) {    
            ierr = VecPointwiseMult(state->qef[itr],Rfs,state->qef[itr]);CHKERRQ(ierr);
          }  
        }
	     }
    } 

    if ((state->usePrescribedBC) && ((Iter0+iLoop==state->bcStartStep))) {
      if (!state->applyBC) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Switching on prescribed BC at (absolute) time step %d\n",Iter0+iLoop);CHKERRQ(ierr);
      }
      if ((state->doCalcBC) && (!state->isInitializedCalcBC)) {
        ierr = TMMComputeCalcBCFunction(time0,Iter0,time0+deltaTClock,Iter0+1,-1,state,TMM_INI_FUNC);
        state->isInitializedCalcBC = PETSC_TRUE;
      }
      state->applyBC = PETSC_TRUE;
    }

    if (state->applyBC) {
      if ((state->bcCutOffStep>0) && ((Iter0+iLoop)>=(state->bcCutOffStep+1))) { /* Note: the >= deals with the situation when Iter0 is > 0 */
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Switching off prescribed BC at (absolute) time step %d\n",Iter0+iLoop);CHKERRQ(ierr);
        state->applyBC = PETSC_FALSE;
      } else {
        if (periodicMatrix) {
          ierr = PeriodicMatInterp(tc,&Be,matrixPeriodicTimer->cyclePeriod,matrixPeriodicTimer->numPerPeriod,
                 matrixPeriodicTimer->tdp,&Bep,matbeFile);
          ierr = PeriodicMatInterp(tc,&Bi,matrixPeriodicTimer->cyclePeriod,matrixPeriodicTimer->numPerPeriod,
                 matrixPeriodicTimer->tdp,&Bip,matbiFile);
        } else if (timeDependentMatrix) {
          ierr = TimeDependentMatInterp(tc,&Be,matrixTimeDependentTimer->numTimes,matrixTimeDependentTimer->tdt,&Betd,matbeFile);
          ierr = TimeDependentMatInterp(tc,&Bi,matrixTimeDependentTimer->numTimes,matrixTimeDependentTimer->tdt,&Bitd,matbiFile);
        }
        if (state->periodicBC) {
          for (itr=0; itr<state->numTracers; itr++) {
            ierr = PeriodicVecInterp(tc,&state->cbc[itr],state->bcTimer->cyclePeriod,state->bcTimer->numPerPeriod,state->bcTimer->tdp,state->cbp[itr],state->bcFile[itr]);
            ierr = PeriodicVecInterp(tf,&state->cbf[itr],state->bcTimer->cyclePeriod,state->bcTimer->numPerPeriod,state->bcTimer->tdp,state->cbp[itr],state->bcFile[itr]);
          }
        } else if (state->timeDependentBC) {
          for (itr=0; itr<state->numTracers; itr++) {
            ierr = TimeDependentVecInterp(tc,&state->cbc[itr],state->bcTimeDependentTimer->numTimes,state->bcTimeDependentTimer->tdt,state->cbtd[itr],state->bcFile[itr]);
            ierr = TimeDependentVecInterp(tf,&state->cbf[itr],state->bcTimeDependentTimer->numTimes,state->bcTimeDependentTimer->tdt,state->cbtd[itr],state->bcFile[itr]);
          }
        } else if (state->doCalcBC) {
          ierr = TMMComputeCalcBCFunction(tc,Iterc,tc+deltaTClock,Iterc+1,iLoop,state,TMM_CALC_FUNC);
		      }
	     }	  
	   }

    if (state->relaxTracer) {
      for (itr=0; itr<state->numTracers; itr++) {
        ierr = VecSet(state->qrel[itr],(state->relaxTracerLambda[itr])*(state->relaxTracerValue[itr]));CHKERRQ(ierr); /* qrel = lambda*vrel */
        ierr = VecAXPY(state->qrel[itr],-(state->relaxTracerLambda[itr]),state->c[itr]);CHKERRQ(ierr); /* qrel = qrel - lambda*c */
      }
    }

/* initialize monitor */
    if (state->useMonitor) {
      if (!state->isInitializedMonitor) {
        ierr = TMMComputeMonitorFunction(time0,Iter0,-1,state,TMM_INI_FUNC);
        state->isInitializedMonitor = PETSC_TRUE;
      }
    }

/* initialize misfit */
    if (state->doMisfit) {  
      if (!state->isInitializedMisfit) {
        ierr = TMMComputeMisfitFunction(time0,Iter0,-1,state,TMM_INI_FUNC);
        state->isInitializedMisfit = PETSC_TRUE;
      }
    }

    return 0;
}

#undef __FUNCT__
#define __FUNCT__ "TMMForcingReinitialize"
PetscErrorCode TMMForcingReinitialize(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, TMMState state)
{

   PetscScalar tf;
   PetscInt itr;
   PetscViewer fd;	
   PetscErrorCode ierr;

//    ierr = PetscPrintf(PETSC_COMM_WORLD,"Reinitializing forcing ...\n");CHKERRQ(ierr);

   tf = tc + deltaTClock;
   
   if ((state->useForcingFromFile) & (state->constantForcing)) { /* Reinitialize (constant) forcing from file */
     for (itr=0; itr<state->numTracers; itr++) {   
       ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d: re-reading forcing from file %s\n", itr,state->forcingFile[itr]);CHKERRQ(ierr);
       ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,state->forcingFile[itr],FILE_MODE_READ,&fd);CHKERRQ(ierr);
       ierr = VecLoad(state->qf[itr],fd);CHKERRQ(ierr); /* IntoVector */
       ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
     }
   }
   
   if (state->useExternalForcing) {
     ierr = TMMComputeExtForcFunction(tc,Iterc,iLoop,state,TMM_REI_FUNC);
   }

   if (state->doCalcBC) {
     ierr = TMMComputeCalcBCFunction(tc,Iterc,tf,Iterc+1,iLoop,state,TMM_REI_FUNC);
   }

   return 0;
}
