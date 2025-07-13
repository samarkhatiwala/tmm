#include <petsc/private/petscimpl.h>

static PetscClassId StepTimer_CLASSID;
static PetscClassId PeriodicTimer_CLASSID;
static PetscClassId TimeDependentTimer_CLASSID;

typedef struct _p_StepTimer *StepTimer;

struct _p_StepTimer {
	PetscBool isInitialized;
	char pre[PETSC_MAX_PATH_LEN];
    PetscBool fixedStep;
    PetscBool haveResetStartTimeStep;
    PetscInt count;
    PetscInt startTimeStep;
    PetscInt startTimeStepResetFreq;
    PetscInt numTimeSteps;
    PetscInt maxNumIntervals;
    PetscInt currInterval;
    PetscInt *timeIntervals;
};

typedef struct _p_PeriodicTimer *PeriodicTimer;

struct _p_PeriodicTimer {
	PetscBool isInitialized;
	char pre[PETSC_MAX_PATH_LEN];
    PetscScalar *tdp;    
    PetscScalar cyclePeriod;
    PetscInt numPerPeriod;
};

typedef struct _p_TimeDependentTimer *TimeDependentTimer;

struct _p_TimeDependentTimer {
	PetscBool isInitialized;
	char pre[PETSC_MAX_PATH_LEN];
    PetscScalar *tdt;
    PetscInt numTimes;
};

extern PetscErrorCode StepTimerCreate(StepTimer *timer);
extern PetscErrorCode PeriodicTimerCreate(PeriodicTimer *timer);
extern PetscErrorCode TimeDependentTimerCreate(TimeDependentTimer *timer);
// extern PetscErrorCode TMMDestroy(PetscScalar tc, TMMState *state);

extern PetscErrorCode StepTimerIni(const char pre1[], const char pre2pre[], PetscInt Iter0, StepTimer thetimer);
extern PetscErrorCode StepTimerUpdate(PetscInt Iter, StepTimer thetimer);
extern PetscErrorCode PeriodicTimerIni(const char pre1[], const char pre2pre[], PetscScalar *fromtdp, PeriodicTimer thetimer);
extern PetscErrorCode TimeDependentTimerIni(const char pre1[], const char pre2pre[], PetscScalar *fromtdt, TimeDependentTimer thetimer);

// extern PetscErrorCode iniStepTimerNew(const char pre1[], const char pre2pre[], PetscInt Iter0, StepTimer *thetimer);
