#define MAX_MATRIX_NUM_PER_PERIOD 20
// #define MAX_FORCING_NUM_PER_PERIOD 3000

typedef struct _p_PeriodicVec *PeriodicVec;

typedef struct _p_PeriodicVec {
  Vec *qp;
  PetscBool firstTime;
  PetscInt numPerPeriod;
} _p_PeriodicVec;

typedef struct _p_PeriodicArray *PeriodicArray;

typedef struct _p_PeriodicArray {
//   PetscScalar *qp[MAX_FORCING_NUM_PER_PERIOD];
  PetscScalar **qp;
  PetscInt arrayLength;
  PetscInt numValsPerProfile;  
  PetscBool firstTime;
  PetscInt numPerPeriod;
} _p_PeriodicArray;

typedef struct _p_TimeDependentVec *TimeDependentVec;

typedef struct _p_TimeDependentVec {
  Vec utd[2];
  PetscInt vecLength;
  PetscBool firstTime;
  PetscInt itcurr;
} _p_TimeDependentVec;

typedef struct _p_TimeDependentArray *TimeDependentArray;

typedef struct _p_TimeDependentArray {
  PetscScalar *utd[2];
  PetscInt arrayLength;
  PetscInt numValsPerProfile;
  PetscBool firstTime;
  PetscInt itcurr;
} _p_TimeDependentArray;

typedef struct {
  Mat Ap[MAX_MATRIX_NUM_PER_PERIOD];
  PetscBool firstTime;
  PetscInt numPerPeriod;
} PeriodicMat;

typedef struct {
  Mat Atd[2];
  PetscBool firstTime;
  PetscInt itcurr;
} TimeDependentMat;


extern PetscErrorCode calcInterpFactor(PetscInt nmax,PetscScalar t,PetscScalar tarray[],PetscInt *itf,PetscScalar *alpha);
extern PetscErrorCode calcPeriodicInterpFactor(PetscInt n,PetscScalar t,PetscScalar tparr[],PetscInt *itf1,PetscInt *itf2,PetscScalar *al1,PetscScalar *al2);
extern PetscInt findindex(PetscScalar tarr[],PetscInt nmax,PetscScalar t);
extern PetscErrorCode PeriodicMatInterp(PetscScalar tc, Mat *A, PetscScalar cyclePeriod, PetscInt numPerPeriod, 
                                           PetscScalar *tdp, PeriodicMat *user, const char *fileName);


extern PetscErrorCode PeriodicVecCreate(PeriodicVec *c);

extern PetscErrorCode PeriodicVecInterp(PetscScalar tc, Vec *u, PetscScalar cyclePeriod, PetscInt numPerPeriod, 
                                           PetscScalar *tdp, PeriodicVec c, const char *fileName);

extern PetscErrorCode PeriodicVecDestroy(PeriodicVec *c);

extern PetscErrorCode PeriodicArrayCreate(PeriodicArray *arr, PetscInt arrayLength);

extern PetscErrorCode PeriodicArrayDestroy(PeriodicArray *arr);

extern PetscErrorCode TimeDependentVecCreate(TimeDependentVec *c);

extern PetscErrorCode TimeDependentVecInterp(PetscScalar tc, Vec *u, PetscInt N, 
                                                PetscScalar *tdp, TimeDependentVec c, const char *fileName);

extern PetscErrorCode TimeDependentVecDestroy(TimeDependentVec *c);

extern PetscErrorCode TimeDependentArrayCreate(TimeDependentArray *arr, PetscInt arrayLength);

extern PetscErrorCode TimeDependentArrayDestroy(TimeDependentArray *arr);

extern PetscErrorCode PeriodicMatDestroy(PeriodicMat *user);

extern PetscErrorCode TimeDependentMatInterp(PetscScalar tc, Mat *A, PetscInt N, 
                                                PetscScalar *tdp, TimeDependentMat *user, const char *fileName);
extern PetscErrorCode TimeDependentMatDestroy(TimeDependentMat *user);

extern PetscErrorCode writeBinaryScalarData(const char *fileName, PetscScalar *arr, PetscInt N, PetscBool appendToFile);
