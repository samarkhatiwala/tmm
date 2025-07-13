#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "tmm_petsc_matvec_utils.h"
#include "tmm_timer.h"
#include "tmm_forcing_utils.h"
#include "tmm_profile_utils.h"
#include "tmm.h"
#include "tmm_share.h"

static PetscClassId MONITOR_CLASSID;

typedef struct _p_MonitorCtx *MonitorContext;
struct _p_MonitorCtx {
  PETSCHEADER(int);
  PetscInt monctxId;
  PetscInt stateId;
/* Add problem-specific variables below */  
};

#undef __FUNCT__
#define __FUNCT__ "iniMonitor"
PetscErrorCode iniMonitor(PetscScalar tc, PetscInt Iter, TMMState state, void *userctx)
{

  PetscErrorCode ierr;
  PetscInt numTracers;
  const char *prefix;
  PetscBool flg;

  static PetscBool registered = PETSC_FALSE;
  static PetscInt monctxId = 0;

  MonitorContext mon;

  MPI_Comm comm = PETSC_COMM_WORLD;
      
  if (!registered) {
    PetscClassIdRegister("Monitor context", &MONITOR_CLASSID);
    registered = PETSC_TRUE;
  }
  PetscHeaderCreate(mon, MONITOR_CLASSID, "Monitor", "Monitor context", "Monitor", comm, 0, 0); 

  monctxId++;
  mon->monctxId=monctxId;
  mon->stateId=state->stateId;

  PetscContainer ctxcontainer;
  PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
  PetscCall(PetscContainerSetPointer(ctxcontainer, (void*)mon));
  PetscCall(PetscObjectCompose((PetscObject)state, "monitor ctx", (PetscObject)ctxcontainer));
  state->monitorctxcontainer = ctxcontainer;
  PetscCall(PetscContainerDestroy(&ctxcontainer));

  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;
  
// Add problem specific code ...
  
  ierr = PetscOptionsGetInt(NULL,prefix,"-monitor_start_time_step",&state->monitorStartTimeStep,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate (absolute) time step at which to start monitor with the -monitor_start_time_step flag");
  ierr = PetscOptionsGetInt(NULL,prefix,"-monitor_steps",&state->monitorSteps,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate frequency with which solution is monitored with the -monitor_steps flag");
  ierr = PetscOptionsGetInt(NULL,prefix,"-monitor_write_steps",&state->monitorWriteSteps,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate frequency at which monitor data is written out with the -monitor_write_steps flag");

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Solution will be monitored starting at time step: %d\n", state->monitorStartTimeStep);CHKERRQ(ierr);	
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Solution will be monitored every %d time steps\n", state->monitorSteps);CHKERRQ(ierr);	
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Monitor data will be written out every %d time steps\n", state->monitorWriteSteps);CHKERRQ(ierr);	

/* Add your code here */

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "calcMonitor"
PetscErrorCode calcMonitor(PetscScalar tc, PetscInt iLoop, TMMState state, void *userctx)

{

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  PetscErrorCode ierr;
  PetscInt numTracers;
  const char *prefix;

  void *ctx;
  PetscCall(PetscContainerGetPointer(state->monitorctxcontainer, &ctx));
  MonitorContext mon = (MonitorContext)ctx;
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

// Add problem specific code ...

  if (Iter0+iLoop>=state->monitorStartTimeStep) { /* start monitoring solution (note: monitorStartTimeStep is ABSOLUTE time step) */	
/* Add your code here */
//     ierr = PetscPrintf(PETSC_COMM_WORLD,"Monitoring solution at time step: %d\n", Iter0+iLoop);CHKERRQ(ierr);	
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "writeMonitor"
PetscErrorCode writeMonitor(PetscScalar tc, PetscInt iLoop, TMMState state, void *userctx)
{

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  PetscErrorCode ierr;
  PetscInt numTracers;
  const char *prefix;

  void *ctx;
  PetscCall(PetscContainerGetPointer(state->monitorctxcontainer, &ctx));
  MonitorContext mon = (MonitorContext)ctx;
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

// Add problem specific code ...

  if (Iter0+iLoop>=state->monitorStartTimeStep) { /* note: monitorStartTimeStep is ABSOLUTE time step */	
/* Add your code here */
//     ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing monitor at time step: %d\n", Iter0+iLoop);CHKERRQ(ierr);	
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "finalizeMonitor"
PetscErrorCode finalizeMonitor(PetscScalar tc, PetscInt Iter, TMMState state, void *userctx)
{

  PetscErrorCode ierr;
  PetscInt numTracers;
  const char *prefix;

  void *ctx;
  PetscCall(PetscContainerGetPointer(state->monitorctxcontainer, &ctx));
  MonitorContext mon = (MonitorContext)ctx;
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

// Add problem specific code ...
//   ierr = PetscPrintf(PETSC_COMM_WORLD,"Finalizing monitor\n");CHKERRQ(ierr);	

  return 0;
}
