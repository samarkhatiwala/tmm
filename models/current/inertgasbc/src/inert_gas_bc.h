extern void inert_gas_diagnostics_(PetscInt *Nrloc, PetscInt *myIter, PetscScalar *myTime, 
                                 PetscScalar localTR[], PetscScalar localTs[], PetscScalar localSs[], PetscInt *gasID, 
                                 PetscScalar localTReq[], PetscScalar localTRsatanom[]);

extern void inert_gas_bc_(PetscInt *Nrloc, PetscInt *myIter, PetscScalar *myTime, 
                                 PetscScalar localTs[], PetscScalar localSs[], 
                                 PetscScalar localatmosp[], PetscInt *gasID,                                  
                                 PetscScalar localTReq[]);

#if !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) 
#define inert_gas_diagnostics_ inert_gas_diagnostics
#define inert_gas_bc_ inert_gas_bc
#endif 
