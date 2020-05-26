
struct timerstruct {
	char	*name;

    PetscBool fixedStep, haveResetStartTimeStep;
    PetscInt count, startTimeStep, startTimeStepResetFreq, numTimeSteps, maxNumIntervals, currInterval;
    PetscInt *timeIntervals;
};

typedef struct timerstruct StepTimer;

struct periodictimestruct {
	char	*name;
    PetscScalar *tdp;    
    PetscScalar cyclePeriod, cycleStep;
    PetscInt numPerPeriod;
};

typedef struct periodictimestruct PeriodicTimer;

struct timedependenttimesstruct {
	char	*name;
    PetscScalar *tdt;
    PetscInt numTimes;
};

typedef struct timedependenttimesstruct TimeDependentTimer;

extern PetscErrorCode iniPeriodicTimer( const char pre[], PeriodicTimer *thetimer );
extern PetscErrorCode iniStepTimer( const char pre[], PetscInt Iter0, StepTimer *thetimer );
extern PetscErrorCode updateStepTimer( const char pre[], PetscInt Iter, StepTimer *thetimer );
extern PetscErrorCode iniTimeDependentTimer( const char pre[], TimeDependentTimer *thetimer );
