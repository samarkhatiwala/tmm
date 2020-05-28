extern void tracegases_model_(PetscInt *myIter, PetscScalar *myTime, 
                                 PetscScalar *localTR,
                                 PetscScalar *localTs, PetscScalar *localSs, PetscScalar *Vgas, PetscScalar *localatmosp, 
                                 PetscInt *gasID, PetscScalar *xTRatm, PetscScalar *mixingRatioScaleFac,
                                 PetscScalar *dzsurf, PetscScalar *localEmP, PetscScalar *TRglobavg,
                                 PetscScalar *localJTR, PetscScalar *localgasexflux, PetscScalar *localtotflux, PetscScalar *localTReq
                                 );

#if !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) 
#define tracegases_model_ tracegases_model
#endif 
