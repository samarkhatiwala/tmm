extern void ocmip_abiotic_carbon_model_(PetscInt *myIter, PetscScalar *myTime, 
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
                                 PetscScalar *localJDIC, PetscScalar *localgasexflux, PetscScalar *localtotflux, PetscScalar *localpco2
#ifdef ALLOW_C14                 
                                 ,PetscScalar *localJDIC14, PetscScalar *localc14gasexflux, PetscScalar *localc14totflux
#endif                                                                  
                                 );
extern void ocmip_abiotic_carbon_ini_(PetscInt *myIter, PetscScalar *myTime, 
                                 PetscScalar *localDIC,PetscScalar *localAlk, PetscScalar *localPO4, PetscScalar *localSiO2,
                                 PetscScalar *localTs,PetscScalar *localSs, PetscScalar *pH);
extern void landsource_(PetscScalar *landState, PetscScalar *pCO2atm, PetscScalar *landUseEmission, 
                                 PetscScalar *deltaTsg, PetscScalar *Fland, PetscScalar *landSource);                                                                 

#if !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) 
#define ocmip_abiotic_carbon_model_ ocmip_abiotic_carbon_model
#define ocmip_abiotic_carbon_ini_ ocmip_abiotic_carbon_ini
#define landsource_ landsource
#endif 
