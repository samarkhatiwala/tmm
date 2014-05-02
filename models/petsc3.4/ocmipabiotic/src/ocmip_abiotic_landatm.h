extern void ocmip_abiotic_model_(PetscInt *myIter, PetscScalar *myTime, 
                                 PetscScalar *localDIC,
#ifdef ALLOW_C14                 
                                 PetscScalar *localDIC14,
#endif                                 
                                 PetscScalar *localAlk, PetscScalar *localPO4, PetscScalar *localSiO2,
                                 PetscScalar *localTs, PetscScalar *localSs, PetscScalar *pH, PetscScalar *Vgas,
                                 PetscScalar *localatmosp, PetscScalar *pCO2atm, PetscScalar *dzsurf,
                                 PetscScalar *localEmP, PetscScalar *DICglobavg,
                                 PetscScalar *linearChemistryFactor, PetscScalar *linearChemistryCO2, PetscScalar *linearChemistryDIC,
#ifdef ALLOW_C14                 
                                 PetscScalar *localD14Catm, PetscScalar *DIC14globavg,
#endif                                                                  
                                 PetscScalar *localJDIC, PetscScalar *localgasexflux, PetscScalar *localtotflux
#ifdef ALLOW_C14                 
                                 ,PetscScalar *localJDIC14, PetscScalar *localc14gasexflux, PetscScalar *localc14totflux
#endif                                                                  
                                 );
extern void ini_ocmip_abiotic_model_(PetscInt *myIter, PetscScalar *myTime, 
                                 PetscScalar *localDIC,PetscScalar *localAlk, PetscScalar *localPO4, PetscScalar *localSiO2,
                                 PetscScalar *localTs,PetscScalar *localSs, PetscScalar *pH);
extern void ocmip_abiotic_diagnostics_(PetscInt *Nrloc, PetscInt *myIter, PetscScalar *myTime, PetscScalar *localpco2);
extern void landsource_(PetscScalar *landState, PetscScalar *pCO2atm, PetscScalar *landUseEmission, 
                                 PetscScalar *deltaTsg, PetscScalar *Fland, PetscScalar *landSource);                                                                 

#if !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) 
#define ocmip_abiotic_model_ ocmip_abiotic_model
#define ini_ocmip_abiotic_model_ ini_ocmip_abiotic_model
#define ocmip_abiotic_diagnostics_ ocmip_abiotic_diagnostics
#define landsource_ landsource
#endif 
