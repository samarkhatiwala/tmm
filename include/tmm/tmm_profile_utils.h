extern PetscErrorCode iniProfileData(PetscInt myId);
extern PetscErrorCode readProfileSurfaceIntData(const char *fileName, PetscInt *arr, PetscInt numValsPerProfile);
extern PetscErrorCode readProfileSurfaceScalarData(const char *fileName, PetscScalar *arr, PetscInt numValsPerProfile);
extern PetscErrorCode readProfileSurfaceScalarDataRecord(const char *fileName, PetscScalar *arr, PetscInt numValsPerProfile, PetscInt iRec);
extern PetscErrorCode writeProfileSurfaceScalarData(const char *fileName, PetscScalar *arr, PetscInt numValsPerProfile, PetscBool appendToFile);
extern PetscErrorCode interpPeriodicProfileSurfaceScalarData(PetscScalar tc, PetscScalar *uarr, PetscScalar cyclePeriod,
                                    PetscInt numPerPeriod, PetscScalar *tdp, 
                                    PeriodicArray user, const char *fileName);
extern PetscErrorCode interpTimeDependentProfileSurfaceScalarData(PetscScalar tc, PetscScalar *uarr, PetscInt numTimes, PetscScalar *tdt,
                                    TimeDependentArray user, const char *fileName);
extern PetscErrorCode dotProdProfileSurfaceScalarData(PetscScalar *xarr, PetscScalar *yarr, PetscScalar *z);
