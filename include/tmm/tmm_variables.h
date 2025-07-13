#ifdef DEFINE_VARIABLES
#define EXTERN
#else
#define EXTERN extern
#endif

//variables made global---------------------------------------------------------

EXTERN Vec templateVec;
EXTERN Vec bcTemplateVec;

/* TM's */
EXTERN Mat Ae, Ai;
EXTERN PeriodicMat Aep, Aip;
EXTERN TimeDependentMat Aetd, Aitd;
EXTERN char mateFile[PETSC_MAX_PATH_LEN], matiFile[PETSC_MAX_PATH_LEN], rfsFile[PETSC_MAX_PATH_LEN];
EXTERN PetscBool periodicMatrix;
EXTERN PetscBool timeDependentMatrix;
EXTERN PetscBool constantMatrix;
EXTERN PeriodicTimer matrixPeriodicTimer;
EXTERN TimeDependentTimer matrixTimeDependentTimer;

/* Rescale forcing */
EXTERN Vec Rfs;
EXTERN PeriodicVec Rfsp;
EXTERN TimeDependentVec Rfstd;

/* BCs */
EXTERN Mat Be, Bi;
EXTERN PeriodicMat Bep, Bip;
EXTERN TimeDependentMat Betd, Bitd; 
EXTERN char matbeFile[PETSC_MAX_PATH_LEN], matbiFile[PETSC_MAX_PATH_LEN];
