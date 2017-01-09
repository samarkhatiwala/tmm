extern PetscErrorCode iniProfileData(PetscInt myId);
extern PetscErrorCode readProfileSurfaceIntData(char *fileName, PetscInt *arr, PetscInt numValsPerProfile);
extern PetscErrorCode readProfileSurfaceScalarData(char *fileName, PetscScalar *arr, PetscInt numValsPerProfile);
extern PetscErrorCode readProfileSurfaceScalarDataRecord(char *fileName, PetscScalar *arr, PetscInt numValsPerProfile, PetscInt iRec);
extern PetscErrorCode writeProfileSurfaceScalarData(char *fileName, PetscScalar *arr, PetscInt numValsPerProfile, PetscBool appendToFile);
extern PetscErrorCode interpPeriodicProfileSurfaceScalarData(PetscScalar tc, PetscScalar *uarr, PetscScalar cyclePeriod,
                                    PetscInt numPerPeriod, PetscScalar *tdp, 
                                    PeriodicArray *user, char *arrFile);
