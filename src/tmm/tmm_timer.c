#define MAX_NUM_INTERVALS 1000

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "tmm_timer.h"

// Note: There are two prefixes in the *TimerIni functions: pre1 and pre2pre. The full option name will be of 
// the form: -[pre1][pre2pre]name, e.g., -state1_write_start_time_steps, if pre1 is "write_" and pre2pre is "state1_". 
// You can pass "" for pre1 and NULL for pre2pre if either of those are absent, in which case the option name would 
// simply be -name, e.g., -start_time_step.

PetscErrorCode StepTimerCreate(StepTimer *timer)
{
  static PetscBool registered = PETSC_FALSE;
  StepTimer s;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  *timer=NULL;

  if (!registered) {
    PetscCall(PetscClassIdRegister("StepTimer", &StepTimer_CLASSID));  
    registered = PETSC_TRUE;
  }
  PetscCall(PetscHeaderCreate(s, StepTimer_CLASSID, "StepTimer", "StepTimer", "StepTimer", PETSC_COMM_WORLD, 0, 0));

  s->isInitialized = PETSC_FALSE;
  
/* run time options */
  s->fixedStep = PETSC_TRUE;
  s->haveResetStartTimeStep = PETSC_FALSE;
  s->count = -1;
  s->startTimeStep = 0;
  s->startTimeStepResetFreq = -1;
  s->numTimeSteps = -1;
  s->maxNumIntervals = -1;
  s->currInterval = -1;
  s->timeIntervals = NULL;
    
  *timer=s;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PeriodicTimerCreate(PeriodicTimer *timer)
{
  static PetscBool registered = PETSC_FALSE;
  PeriodicTimer s;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  *timer=NULL;

  if (!registered) {
    PetscCall(PetscClassIdRegister("PeriodicTimer", &PeriodicTimer_CLASSID));  
    registered = PETSC_TRUE;
  }
  PetscCall(PetscHeaderCreate(s, PeriodicTimer_CLASSID, "PeriodicTimer", "PeriodicTimer", "PeriodicTimer", PETSC_COMM_WORLD, 0, 0));

  s->isInitialized = PETSC_FALSE;
  
/* run time options */
  s->tdp = NULL;
  s->cyclePeriod = -1;
  s->numPerPeriod = -1;
    
  *timer=s;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode TimeDependentTimerCreate(TimeDependentTimer *timer)
{
  static PetscBool registered = PETSC_FALSE;
  TimeDependentTimer s;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;

  *timer=NULL;

  if (!registered) {
    PetscCall(PetscClassIdRegister("TimeDependentTimer", &TimeDependentTimer_CLASSID));  
    registered = PETSC_TRUE;
  }
  PetscCall(PetscHeaderCreate(s, TimeDependentTimer_CLASSID, "TimeDependentTimer", "TimeDependentTimer", "TimeDependentTimer", PETSC_COMM_WORLD, 0, 0));

  s->isInitialized = PETSC_FALSE;
  
/* run time options */
  s->tdt = NULL;
  s->numTimes = -1;
  
  *timer=s;
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode PeriodicTimerIni(const char pre1[], const char pre2pre[], PetscScalar *fromtdp, PeriodicTimer thetimer)
{

  PetscErrorCode ierr;
  PetscBool flg, flg1;
  PetscScalar cycleStep;
  PetscInt it;
  PetscViewer fd;
  int fp;
  char timeFile[PETSC_MAX_PATH_LEN];

// Note: pass NULL if the fromtdp argument is to be ignored

  if (pre2pre && pre2pre[0]) {
    strcpy(thetimer->pre,pre2pre);
    ierr=PetscStrcat(thetimer->pre,pre1);CHKERRQ(ierr);  
  } else {
    strcpy(thetimer->pre,pre1);    
  }

/*  read time data */
  ierr = PetscOptionsGetReal(NULL,thetimer->pre,"-cycle_period",&thetimer->cyclePeriod,&flg);CHKERRQ(ierr);
#if PETSC_VERSION_LT(3,17,0)
  if (!flg) SETERRQ1(PETSC_COMM_WORLD,1,"Must indicate cycling time with the -%scycle_period option",thetimer->pre);
#else
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate cycling time with the -%scycle_period option",thetimer->pre);
#endif
	 ierr = PetscOptionsGetInt(NULL,thetimer->pre,"-num_per_period",&thetimer->numPerPeriod,&flg1);CHKERRQ(ierr);
#if PETSC_VERSION_LT(3,17,0)      
 	if (!flg1) SETERRQ1(PETSC_COMM_WORLD,1,"Must indicate number of fields per period with the -%snum_per_period option",thetimer->pre);
#else
	 if (!flg1) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate number of fields per period with the -%snum_per_period option",thetimer->pre);
#endif
/*    array for holding extended time array */
	 PetscMalloc((thetimer->numPerPeriod+2)*sizeof(PetscScalar), &thetimer->tdp);
  if (fromtdp) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using provided array to set times for periodic timer %s\n",thetimer->pre);CHKERRQ(ierr);
    for (it=1; it<=thetimer->numPerPeriod; it++) {
      thetimer->tdp[it]=fromtdp[it-1];
    }
    thetimer->tdp[0]=thetimer->tdp[thetimer->numPerPeriod]-thetimer->cyclePeriod;
    thetimer->tdp[thetimer->numPerPeriod+1]=thetimer->tdp[1]+thetimer->cyclePeriod;
  } else {	  
    ierr = PetscOptionsGetString(NULL,thetimer->pre,"-periodic_times_file",timeFile,PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
    if (flg1) {
      ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,timeFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
      ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
      ierr = PetscBinaryRead(fp,&thetimer->tdp[1],thetimer->numPerPeriod,NULL,PETSC_SCALAR);CHKERRQ(ierr);  
      ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      thetimer->tdp[0]=thetimer->tdp[thetimer->numPerPeriod]-thetimer->cyclePeriod;
      thetimer->tdp[thetimer->numPerPeriod+1]=thetimer->tdp[1]+thetimer->cyclePeriod;
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Assuming equally-spaced fields for periodic object %s\n",thetimer->pre);CHKERRQ(ierr);            
      cycleStep=thetimer->cyclePeriod/thetimer->numPerPeriod;
      for (it=0; it<=thetimer->numPerPeriod+1; it++) {
        thetimer->tdp[it]=(-cycleStep/2.0) + it*cycleStep;
      } 	
    }
  }	
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Periodic object %s specified at times:\n",thetimer->pre);CHKERRQ(ierr);            
  for (it=0; it<=thetimer->numPerPeriod+1; it++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"tdp=%10.5f\n", thetimer->tdp[it]);CHKERRQ(ierr);        
  }    
    
  return 0;
    
}

PetscErrorCode StepTimerIni(const char pre1[], const char pre2pre[], PetscInt startTimeStep, StepTimer thetimer)
{

  PetscErrorCode ierr;
  PetscBool flg;
  PetscInt it;
  PetscInt tmparr[MAX_NUM_INTERVALS];

  if (pre2pre && pre2pre[0]) {
 	  strcpy(thetimer->pre,pre2pre);
    ierr=PetscStrcat(thetimer->pre,pre1);CHKERRQ(ierr);  
  } else {
    strcpy(thetimer->pre,pre1);    
  }
    
/*  read time step data */
  thetimer->startTimeStep = startTimeStep; /* by default we start at first time step */
  ierr = PetscOptionsGetInt(NULL,thetimer->pre,"-start_time_step",&thetimer->startTimeStep,&flg);CHKERRQ(ierr);
  if (thetimer->startTimeStep < startTimeStep) {
#if PETSC_VERSION_LT(3,17,0)
    SETERRQ2(PETSC_COMM_WORLD,1,"Specified start time step for StepTimer object %s is < %d",thetimer->pre,startTimeStep);
#else
    SETERRQ(PETSC_COMM_WORLD,1,"Specified start time step for StepTimer object %s is < %d",thetimer->pre,startTimeStep);
#endif
  }  
 	ierr = PetscPrintf(PETSC_COMM_WORLD,"Start time step (absolute) for StepTimer object %s is %d\n", thetimer->pre, thetimer->startTimeStep);CHKERRQ(ierr);	  

  thetimer->maxNumIntervals=MAX_NUM_INTERVALS;
  ierr = PetscOptionsGetIntArray(NULL,thetimer->pre,"-time_steps",tmparr,&thetimer->maxNumIntervals,&flg);CHKERRQ(ierr);
#if PETSC_VERSION_LT(3,17,0)
  if (!flg) SETERRQ1(PETSC_COMM_WORLD,1,"Must indicate number of step timer time steps with the -%stime_steps flag",thetimer->pre);
#else
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate number of step timer time steps with the -%stime_steps flag",thetimer->pre);
#endif
  if (thetimer->maxNumIntervals==1) {
    thetimer->fixedStep=PETSC_TRUE;
    thetimer->currInterval=0; /* Not used but we set it anyway to be safe */
    thetimer->numTimeSteps=tmparr[0];
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Fixed interval of %d specified for StepTimer object %s\n", thetimer->numTimeSteps, thetimer->pre);CHKERRQ(ierr);	  
  } else {
    thetimer->fixedStep=PETSC_FALSE;      
    PetscMalloc(thetimer->maxNumIntervals*sizeof(PetscInt), &thetimer->timeIntervals);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Variable number of intervals specified for StepTimer object %s\n", thetimer->pre);CHKERRQ(ierr);	  
    for (it=0; it<thetimer->maxNumIntervals; it++) {
      thetimer->timeIntervals[it] = tmparr[it];
      ierr = PetscPrintf(PETSC_COMM_WORLD,"  Interval #%d=%d\n", it+1,thetimer->timeIntervals[it]);CHKERRQ(ierr);        
    }	  
    thetimer->currInterval=0;
    thetimer->numTimeSteps=thetimer->timeIntervals[thetimer->currInterval];
  }

  thetimer->startTimeStepResetFreq=-1;
  ierr = PetscOptionsGetInt(NULL,thetimer->pre,"-start_time_step_reset_freq",&thetimer->startTimeStepResetFreq,&flg);CHKERRQ(ierr);
  if (flg) {
    PetscInt tmp=0;
    if (!thetimer->fixedStep) {
      for (it=0; it<thetimer->maxNumIntervals; it++) {
        tmp=tmp+(thetimer->timeIntervals[it]);
      } 
    } else {
      tmp=thetimer->numTimeSteps;
    }
    if (tmp > thetimer->startTimeStepResetFreq) {
#if PETSC_VERSION_LT(3,17,0)
      SETERRQ1(PETSC_COMM_WORLD,1,"Start time reset frequency less than total number of timer steps for StepTimer object %s",thetimer->pre);
#else
      SETERRQ(PETSC_COMM_WORLD,1,"Start time reset frequency less than total number of timer steps for StepTimer object %s",thetimer->pre);
#endif
    } 
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Start time will be reset every %d steps for StepTimer object %s\n", thetimer->startTimeStepResetFreq, thetimer->pre);CHKERRQ(ierr);	  
  }    
    
  thetimer->haveResetStartTimeStep=PETSC_FALSE;
	 thetimer->count=0;
	    
  return 0;
    
}

PetscErrorCode StepTimerUpdate(PetscInt Iter, StepTimer thetimer)
{

  PetscBool endOfSequence = PETSC_TRUE;
  PetscErrorCode ierr;

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Updating StepTimer object %s at (absolute) time step %d\n", thetimer->pre,Iter);CHKERRQ(ierr);        

  thetimer->count=0; /* reset counter */
  thetimer->haveResetStartTimeStep=PETSC_FALSE;
    
  if (!thetimer->fixedStep) {
    thetimer->currInterval++;
    if (thetimer->currInterval==thetimer->maxNumIntervals) {
//    We're now at the end of the sequence      
      thetimer->currInterval=0;
    } else {
//    Still within sequence      
      endOfSequence = PETSC_FALSE;
    }       
    thetimer->numTimeSteps=thetimer->timeIntervals[thetimer->currInterval];
    ierr = PetscPrintf(PETSC_COMM_WORLD,"New interval for StepTimer object %s at (absolute) time step %d is %d\n", thetimer->pre,Iter,thetimer->numTimeSteps);CHKERRQ(ierr);
  }

  if ((thetimer->startTimeStepResetFreq > 0) && (endOfSequence)) {
    thetimer->startTimeStep = thetimer->startTimeStep + thetimer->startTimeStepResetFreq;
    thetimer->haveResetStartTimeStep=PETSC_TRUE;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"New start time step (absolute) for StepTimer object %s is %d\n", thetimer->pre, thetimer->startTimeStep);CHKERRQ(ierr);
  }

  return 0;
    
}

PetscErrorCode TimeDependentTimerIni(const char pre1[], const char pre2pre[], PetscScalar *fromtdt, TimeDependentTimer thetimer)
{

  PetscErrorCode ierr;
  PetscBool flg;
  PetscInt it;    
  PetscViewer fd;
  int fp;
  char timeFile[PETSC_MAX_PATH_LEN];

// Note: pass NULL if the fromtdt argument is to be ignored

  if (pre2pre && pre2pre[0]) {
	   strcpy(thetimer->pre,pre2pre);
    ierr=PetscStrcat(thetimer->pre,pre1);CHKERRQ(ierr);  
  } else {
    strcpy(thetimer->pre,pre1);    
  }

/*  read time data */
  ierr = PetscOptionsGetInt(NULL,thetimer->pre,"-num_times",&thetimer->numTimes,&flg);CHKERRQ(ierr);
#if PETSC_VERSION_LT(3,17,0)
  if (!flg) SETERRQ1(PETSC_COMM_WORLD,1,"Must indicate number of time slices with the -%snum_times option",thetimer->pre);
#else
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate number of time slices with the -%snum_times option",thetimer->pre);
#endif
	 ierr = PetscMalloc((thetimer->numTimes)*sizeof(PetscScalar), &thetimer->tdt);CHKERRQ(ierr);
	 if (fromtdt) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using provided array to set times for time dependent timer %s\n",thetimer->pre);CHKERRQ(ierr);
    for (it=0; it<thetimer->numTimes; it++) {
      thetimer->tdt[it]=fromtdt[it];
    }    
  } else {
    ierr = PetscOptionsGetString(NULL,thetimer->pre,"-times_file",timeFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
#if PETSC_VERSION_LT(3,17,0)
    if (!flg) SETERRQ1(PETSC_COMM_WORLD,1,"Must indicate name of time file with the -%stimes_file option",thetimer->pre);
#else
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate name of time file with the -%stimes_file option",thetimer->pre);
#endif	
    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,timeFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
    ierr = PetscBinaryRead(fp,&thetimer->tdt[0],thetimer->numTimes,NULL,PETSC_SCALAR);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Time-dependent object %s specified at times:\n",thetimer->pre);CHKERRQ(ierr);            
  for (it=0; it<thetimer->numTimes; it++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"tdt=%10.5f\n", thetimer->tdt[it]);CHKERRQ(ierr);        
  }    
    
  return 0;
    
}
