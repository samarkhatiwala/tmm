extern PetscErrorCode iniCalcBC(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, TMMState state, void *ctx);
extern PetscErrorCode calcCalcBC(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, TMMState state, void *ctx);
extern PetscErrorCode writeCalcBC(PetscScalar tc, PetscInt Iterc, PetscInt iLoop, TMMState state, void *ctx);
extern PetscErrorCode finalizeCalcBC(PetscScalar tc, PetscInt Iterc, TMMState state, void *ctx);
extern PetscErrorCode reInitializeCalcBC(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, TMMState state, void *ctx);
