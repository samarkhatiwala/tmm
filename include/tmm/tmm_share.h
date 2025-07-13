// #define MAXNUMTRACERS 100

#ifdef DEFINE_VARIABLES
#define EXTERNSHARE
#else
#define EXTERNSHARE extern
#endif

#ifdef DEFINE_PROFILES
#define EXTERNPROF
#else
#define EXTERNPROF extern
#endif

//originally global in main-----------------------------------------------------
EXTERNSHARE PetscScalar deltaTClock, time0;
EXTERNSHARE PetscInt maxSteps, Iter0;
// EXTERNSHARE StepTimer writeTimer;
// EXTERNSHARE PetscBool appendOutput;
EXTERNSHARE PetscInt nb;
EXTERNSHARE PetscInt *gIndices, gLow, gHigh;
EXTERNSHARE PetscInt *gBCIndices, lBCSize, gBCSize, gbcLow, gbcHigh;
// EXTERNSHARE PetscBool doMisfit;
EXTERNSHARE PetscBool rescaleForcing;
//------------------------------------------------------------------------------

EXTERNSHARE PetscBool prescribedBCInUse;
EXTERNSHARE PetscBool calcBCInUse;

//originally in tmm_profile_data.h -----------------------------------------------------
EXTERNPROF PetscInt *gNumProfiles, *gStartIndices, *gEndIndices, *lStartIndices, *lEndIndices, *gProfileLengths, *gSizes;
EXTERNPROF PetscInt *lProfileLength, lNumProfiles, lSize, numPrevProfiles, totalNumProfiles;
EXTERNPROF PetscBool useProfiles;
