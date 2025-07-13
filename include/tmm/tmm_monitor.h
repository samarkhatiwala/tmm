extern PetscErrorCode iniMonitor(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
extern PetscErrorCode calcMonitor(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx);
extern PetscErrorCode writeMonitor(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx);
extern PetscErrorCode finalizeMonitor(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
