extern PetscErrorCode iniExternalForcing(PetscScalar tc, PetscInt n, PetscInt ntr, Vec *v, Vec *u);
extern PetscErrorCode calcExternalForcing(PetscScalar tc, PetscInt n, PetscInt it, PetscInt ntr, Vec *v, Vec *u);
extern PetscErrorCode writeExternalForcing(PetscScalar tc, PetscInt it, PetscInt ntr, Vec *v, Vec *u);
extern PetscErrorCode finalizeExternalForcing(PetscScalar tc, PetscInt n, PetscInt ntr);
extern PetscErrorCode reInitializeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt numTracers, Vec *v, Vec *ut);
