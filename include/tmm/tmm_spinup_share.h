#ifdef DEFINE_SPINUP_VARIABLES
#define EXTERNSPINUPSHARE
#else
#define EXTERNSPINUPSHARE extern
#endif

typedef struct {
  PetscInt maxSteps;
  PetscInt Iter0;  
  PetscScalar time0;
  PetscScalar deltaTClock;
  PetscInt itf;
  PetscInt checkpointFreq;
  FILE *logfp;
  Vec Xini;
  PetscScalar *XTovScaleFac;
  PetscScalar *vToXScaleFac;
  PetscBool runExternalModel;
  PetscInt convergedValue;  
  TMMState state;
} TMMSPINUP;
