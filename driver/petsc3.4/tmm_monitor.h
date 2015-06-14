extern PetscErrorCode iniMonitor(PetscScalar tc, PetscInt n, PetscInt ntr, Vec *v);
extern PetscErrorCode updateMonitor(PetscScalar tc, PetscInt it, PetscInt ntr, Vec *v);
extern PetscErrorCode writeMonitor(PetscScalar tc, PetscInt it, PetscInt ntr, Vec *v);
extern PetscErrorCode finalizeMonitor(PetscScalar tc, PetscInt n, PetscInt ntr);
