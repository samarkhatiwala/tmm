extern PetscErrorCode iniExternalForcing(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
extern PetscErrorCode calcExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *ctx);
extern PetscErrorCode writeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *ctx);
extern PetscErrorCode finalizeExternalForcing(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
extern PetscErrorCode reInitializeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *ctx);
