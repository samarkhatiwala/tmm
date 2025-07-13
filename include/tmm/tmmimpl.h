#define MAXNUMTRACERS 100
#define MAXNUMSTATES 10

#include <petsc/private/petscimpl.h>

static PetscClassId TMM_CLASSID;

typedef enum {
  TMM_INI_FUNC,
  TMM_CALC_FUNC,
  TMM_WRI_FUNC,
  TMM_FIN_FUNC,
  TMM_REI_FUNC
} TMMFUNCTYPE;

typedef struct _p_TMMState *TMMState;

typedef struct _TMMOps *TMMOps;
struct _TMMOps {
  PetscErrorCode (*iniexternalforcing)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
  PetscErrorCode (*calcexternalforcing)(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *ctx);
  PetscErrorCode (*wriexternalforcing)(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *ctx);
  PetscErrorCode (*finexternalforcing)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
  PetscErrorCode (*reiexternalforcing)(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *ctx);

  PetscErrorCode (*inicalcbc)(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, TMMState state, void *ctx);
  PetscErrorCode (*calccalcbc)(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, TMMState state, void *ctx);
  PetscErrorCode (*wricalcbc)(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, TMMState state, void *ctx);
  PetscErrorCode (*fincalcbc)(PetscScalar tc, PetscInt Iterc, TMMState state, void *ctx);
  PetscErrorCode (*reicalcbc)(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, TMMState state, void *ctx);

  PetscErrorCode (*inimonitor)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
  PetscErrorCode (*calcmonitor)(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx);
  PetscErrorCode (*wrimonitor)(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx);
  PetscErrorCode (*finmonitor)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);

  PetscErrorCode (*inimisfit)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
  PetscErrorCode (*calcmisfit)(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx);
  PetscErrorCode (*wrimisfit)(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx);
  PetscErrorCode (*finmisfit)(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);

};

struct _p_TMMState {
  PETSCHEADER(struct _TMMOps);
// Containers to hold data (these are created in the inExternalForcing, iniCalcBC etc routines)
  PetscContainer extforcctxcontainer;
  PetscContainer calcbcctxcontainer;
  PetscContainer monitorctxcontainer;
  PetscContainer misfitctxcontainer;

// The following are really only needed for python support. 
// You can use these to store additional user context data when 
// calling the TMMSet* functions from C. The TMMCompute* functions will fetch 
// and pass it to the actual functions as the final user context argument. But 
// application data is now stored in a container within the state struct so 
// there's really no need for it as far as C is concerned. They are however 
// needed by the python interface to store pointers to the actual python 
// callback function, arguments to it, and the python state object (this is done 
// in the setIni* functions in tmm_interface.pyx). The TMMCompute* functions will 
// pass the user context to the *ForcingFn callback functions, which in turn will 
// unpack the pointers and call the actual python callback function with the 
// arguments and python state object.
  PetscContainer iniextforcuserctxcontainer;
  PetscContainer calcextforcuserctxcontainer;
  PetscContainer wriextforcuserctxcontainer;
  PetscContainer finextforcuserctxcontainer;
  PetscContainer reiextforcuserctxcontainer;

  PetscContainer inicalcbcuserctxcontainer;
  PetscContainer calccalcbcuserctxcontainer;
  PetscContainer wricalcbcuserctxcontainer;
  PetscContainer fincalcbcuserctxcontainer;
  PetscContainer reicalcbcuserctxcontainer;

  PetscContainer inimonitoruserctxcontainer;
  PetscContainer calcmonitoruserctxcontainer;
  PetscContainer wrimonitoruserctxcontainer;
  PetscContainer finmonitoruserctxcontainer;

  PetscContainer inimisfituserctxcontainer;
  PetscContainer calcmisfituserctxcontainer;
  PetscContainer wrimisfituserctxcontainer;
  PetscContainer finmisfituserctxcontainer;

  PetscInt stateId;
  
  Vec *c;
  Vec *qf;
  Vec *qef;
  Vec *qrel;
  Vec *cbc;
  Vec *cbf;

  PetscInt numTracers;
  
  PetscBool isInitialized;
  PetscBool useExternalForcing;
  PetscBool useForcingFromFile;
  PetscBool usePrescribedBC;
  PetscBool applyExternalForcing;
  PetscBool applyForcingFromFile;
  PetscBool applyBC;
  PetscBool periodicForcing;
  PetscBool timeDependentForcing;
  PetscBool constantForcing;
  PetscBool periodicBC;
  PetscBool timeDependentBC;
  PetscBool constantBC;
  PetscBool doCalcBC;
  PetscBool useMonitor;
  PetscBool doMisfit;
  PetscBool relaxTracer;
  PetscBool isInitializedExternalForcing;
  PetscBool isInitializedCalcBC;
  PetscBool isInitializedMonitor;
  PetscBool isInitializedMisfit;
  PetscBool doOutput;
  PetscBool appendOutput;
  PetscBool writePickup;
  PetscBool doWriteBC;
  PetscBool doWriteQF;
  PetscBool doWriteQEF;
  PetscBool pickupFromFile;
  PetscBool doTimeAverage;
  PetscBool avgAppendOutput;
  PetscBool doExtraWrite;

/* Forcing */
  PeriodicVec qp[MAXNUMTRACERS];
  char *forcingFile[MAXNUMTRACERS];
  PeriodicTimer forcingTimer;
  TimeDependentVec qtdf[MAXNUMTRACERS];
  TimeDependentTimer forcingTimeDependentTimer;
  PetscInt forcingFromFileStartStep;
  PetscInt externalForcingStartStep;
  PetscInt forcingFromFileCutOffStep;
  PetscInt externalForcingCutOffStep;

  PetscScalar *relaxTracerLambda;
  PetscScalar *relaxTracerValue;

/* BCs */
  PeriodicVec cbp[MAXNUMTRACERS];
  char *bcFile[MAXNUMTRACERS];
  PeriodicTimer bcTimer;
  TimeDependentVec cbtd[MAXNUMTRACERS];
  TimeDependentTimer bcTimeDependentTimer;
  PetscInt bcStartStep;
  PetscInt bcCutOffStep;

/* I/O   */
  char *iniFile[MAXNUMTRACERS]; 
  
  StepTimer writeTimer;
  char *outFile[MAXNUMTRACERS]; 
  char outTimeFile[PETSC_MAX_PATH_LEN]; 
  PetscFileMode OUTPUT_FILE_MODE;
  FILE *fptime;
  PetscViewer fdout[MAXNUMTRACERS];

  char *bcoutFile[MAXNUMTRACERS];
  PetscViewer fdbcout[MAXNUMTRACERS];
  char *bcavgOutFile[MAXNUMTRACERS];
  PetscFileMode BC_FILE_MODE;
  PetscViewer fdbcavgout[MAXNUMTRACERS];
  PetscFileMode BCAVG_FILE_MODE;
  Vec *cbavg;

  char *qfoutFile[MAXNUMTRACERS];
  PetscViewer fdqfout[MAXNUMTRACERS];  
  PetscFileMode QF_FILE_MODE;
  PetscViewer fdqfavgout[MAXNUMTRACERS]; 
  char *qfavgOutFile[MAXNUMTRACERS];
  PetscFileMode QFAVG_FILE_MODE;
  Vec *qfavg;

  char *qefoutFile[MAXNUMTRACERS];
  PetscViewer fdqefout[MAXNUMTRACERS];
  PetscFileMode QEF_FILE_MODE;
  char *qefavgOutFile[MAXNUMTRACERS];
  PetscViewer fdqefavgout[MAXNUMTRACERS];
  PetscFileMode QEFAVG_FILE_MODE;
  Vec *qefavg;

  StepTimer pickupTimer;
  char pickupFile[PETSC_MAX_PATH_LEN];
  char pickupoutFile[PETSC_MAX_PATH_LEN]; 

  StepTimer avgTimer;
  char *avgOutFile[MAXNUMTRACERS];
  PetscViewer fdavgout[MAXNUMTRACERS];
  PetscFileMode AVG_FILE_MODE;
  FILE *avgfptime;
  char avgOutTimeFile[PETSC_MAX_PATH_LEN]; 
  Vec *cavg;

  StepTimer extraWriteTimer;
  char *extraOutFile[MAXNUMTRACERS];
  FILE *fptimeextra;
  PetscViewer fdoutextra[MAXNUMTRACERS];
  PetscFileMode OUTPUT_EXTRA_FILE_MODE;
  PetscInt numExtraTracers;
  PetscInt itrExtra[MAXNUMTRACERS];
  
  StepTimer misfitTimer;

  PetscInt monitorStartTimeStep;
  PetscInt monitorSteps;
  PetscInt monitorWriteSteps;
  
  Vec *cwork;
};
