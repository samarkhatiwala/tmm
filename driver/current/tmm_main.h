#define MAXNUMTRACERS 30

extern PetscScalar deltaTClock, time0;
extern PetscInt maxSteps, Iter0, writeSteps;
extern PetscInt *gIndices, gLow, gHigh;
extern PetscInt *gBCIndices, lBCSize, gBCSize, gbcLow, gbcHigh;