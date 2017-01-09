
struct timerstruct {
	char	*name;

    PetscBool fixedStep;
    PetscInt count, startTimeStep, numTimeSteps, maxNumIntervals, currInterval;
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

extern PetscErrorCode iniPeriodicTimer( const char pre[], PeriodicTimer *thetimer );
extern PetscErrorCode iniStepTimer( const char pre[], PetscInt Iter0, StepTimer *thetimer );
extern PetscErrorCode updateStepTimer( const char pre[], PetscInt Iter, StepTimer *thetimer );
