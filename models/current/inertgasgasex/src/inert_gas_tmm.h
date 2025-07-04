extern void inert_gas_fluxes_(PetscScalar *localTR,PetscScalar *localTs,PetscScalar *localSs, 
                              PetscScalar *localwind, PetscScalar *localfice,PetscScalar *localatmosp, 
                              PetscInt *gasID, PetscScalar *pistonVelocityCoeff,PetscScalar *localVgas,
                              PetscScalar *localFinj,PetscScalar *localFex,
                              PetscScalar *localPTReq);

extern void inert_gas_diagnostics_(PetscInt *Nrloc,PetscScalar localTR[], PetscScalar localTs[], PetscScalar localSs[], 
                                   PetscInt *gasID, PetscScalar localTReq[], PetscScalar localTRsatanom[]);

#if !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) 
#define inert_gas_fluxes_ inert_gas_fluxes
#define inert_gas_diagnostics_ inert_gas_diagnostics
#endif 
