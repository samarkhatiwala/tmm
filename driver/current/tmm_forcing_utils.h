#define MAX_MATRIX_NUM_PER_PERIOD 20
#define MAX_FORCING_NUM_PER_PERIOD 2000
typedef struct {
  Vec *up;
  PetscBool firstTime;
  PetscInt numPerPeriod;
} PeriodicVec;

typedef struct {
  Mat Ap[MAX_MATRIX_NUM_PER_PERIOD];
  PetscBool firstTime;
  PetscInt numPerPeriod;
} PeriodicMat;

typedef struct {
  PetscScalar *up[MAX_FORCING_NUM_PER_PERIOD];
  PetscInt arrayLength;
  PetscBool firstTime;
  PetscInt numPerPeriod;
} PeriodicArray;

extern PetscErrorCode calcInterpFactor(PetscInt nmax,PetscScalar t,PetscScalar tarray[],PetscInt *itf,PetscScalar *alpha);
extern PetscErrorCode calcPeriodicInterpFactor(PetscInt n,PetscScalar t,PetscScalar tparr[],PetscInt *itf1,PetscInt *itf2,PetscScalar *al1,PetscScalar *al2);
extern PetscInt findindex(PetscScalar tarr[],PetscInt nmax,PetscScalar t);
extern PetscErrorCode interpPeriodicMatrix(PetscScalar tc, Mat *A, PetscScalar cyclePeriod, PetscInt numPerPeriod, 
                                           PetscScalar *tdp, PeriodicMat *user, char *filename);
extern PetscErrorCode interpPeriodicVector(PetscScalar tc, Vec *u, PetscScalar cyclePeriod, PetscInt numPerPeriod, 
                                           PetscScalar *tdp, PeriodicVec *user, char *filename);

extern PetscErrorCode destroyPeriodicVec(PeriodicVec *user);
extern PetscErrorCode destroyPeriodicMat(PeriodicMat *user);
extern PetscErrorCode destroyPeriodicArray(PeriodicArray *user);

extern PetscErrorCode interpTimeDependentVector(PetscScalar tc, Vec *u, PetscInt numTracers, PetscInt nt, PetscScalar *t, Vec **ut);

extern PetscErrorCode writeBinaryScalarData(char *fileName, PetscScalar *arr, PetscInt N, PetscBool appendToFile);
