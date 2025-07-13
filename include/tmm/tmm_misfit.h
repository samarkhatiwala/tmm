extern PetscErrorCode iniMisfit(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
extern PetscErrorCode calcMisfit(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx);
extern PetscErrorCode writeMisfit(PetscScalar tc, PetscInt iLoop, TMMState state, void *ctx);
extern PetscErrorCode finalizeMisfit(PetscScalar tc, PetscInt Iter, TMMState state, void *ctx);
