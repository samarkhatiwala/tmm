extern PetscBool doMisfit;
extern StepTimer misfitTimer;
extern PetscErrorCode iniMisfit(PetscScalar tc, PetscInt n, PetscInt ntr, Vec *v);
extern PetscErrorCode calcMisfit(PetscScalar tc, PetscInt it, PetscInt ntr, Vec *v);
extern PetscErrorCode writeMisfit(PetscScalar tc, PetscInt it, PetscInt ntr, Vec *v);
extern PetscErrorCode finalizeMisfit(PetscScalar tc, PetscInt n, PetscInt ntr);
