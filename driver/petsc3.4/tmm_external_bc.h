extern PetscErrorCode iniCalcBC(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt ntr, Vec *v, Vec *bcc, Vec *bcf);
extern PetscErrorCode calcBC(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, PetscInt ntr, Vec *v, Vec *bcc, Vec *bcf);
extern PetscErrorCode writeBC(PetscScalar tc, PetscInt iLoop,PetscInt numTracers, Vec *v, Vec *bcc, Vec *bcf);
extern PetscErrorCode finalizeCalcBC(PetscScalar tc, PetscInt n, PetscInt ntr);
extern PetscErrorCode reInitializeCalcBC(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt ntr, Vec *v, Vec *bcc, Vec *bcf);
