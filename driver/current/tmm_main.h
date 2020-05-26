#define MAXNUMTRACERS 100

extern PetscScalar deltaTClock, time0;
extern PetscInt maxSteps, Iter0;
extern StepTimer writeTimer;
extern PetscInt *gIndices, gLow, gHigh;
extern PetscInt *gBCIndices, lBCSize, gBCSize, gbcLow, gbcHigh;
extern PetscBool rescaleForcing;
