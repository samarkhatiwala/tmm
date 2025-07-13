extern PetscInt getExternalSignal(const char cmd[]);
extern PetscInt getExternalSignalMPI(PetscInt isend);
extern PetscInt getExternalSignalFile(PetscInt isend, PetscInt waitTime);
