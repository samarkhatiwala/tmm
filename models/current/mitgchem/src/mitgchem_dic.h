extern void mitgchem_copy_data_(PetscInt *nzloc, PetscInt *itr, PetscScalar localTR[], PetscScalar localJTR[], 
                                PetscScalar *DeltaT, PetscInt *direction);

extern void mitgchem_ini_(PetscInt *Nrloc, PetscInt *numTracers, PetscInt *myIter, PetscScalar *myTime, 
						  PetscScalar localTs[],PetscScalar localSs[], PetscScalar *pH, PetscScalar *localsilica,
						  PetscScalar localhFacC[], PetscScalar localrecip_hFacC[], PetscScalar drF[], PetscScalar *deltaT, PetscInt *ip);

extern void mitgchem_model_(PetscInt *Nrloc, PetscInt *myIter, PetscScalar *myTime, 
							PetscScalar localTs[],PetscScalar localSs[], PetscScalar *pH, PetscScalar *localwind,
							PetscScalar *localatmosp, PetscScalar *pco2atmloc,
							PetscScalar *localsilica, PetscScalar *localfice, 
#ifdef READ_PAR
							PetscScalar *localpar, 
#else
							PetscScalar *locallat, 
#endif
#ifdef ALLOW_FE                                    
							PetscScalar *localinputfe, 
#endif
#ifdef LIGHT_CHL
                            PetscScalar *chlloc,
#endif
							PetscScalar *localalpha, PetscScalar *localrain_ratio,
							PetscScalar localhFacC[], PetscScalar localrecip_hFacC[],
							PetscScalar *pco2loc, PetscScalar *co2fluxloc,                                    
						    PetscInt *ip);

extern void mitgchem_diagnostics_(PetscInt *Nrloc, PetscInt *myIter, PetscScalar *myTime, 
                                           PetscScalar localbioac[]);

extern void landsource_(PetscScalar *landState, PetscScalar *pCO2atm, PetscScalar *landUseEmission, 
                                 PetscScalar *deltaTsg, PetscScalar *Fland, PetscScalar *landSource);                                                                 

#if !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) 
#define mitgchem_copy_data_ mitgchem_copy_data
#define mitgchem_ini_ mitgchem_ini
#define mitgchem_model_ mitgchem_model
#define mitgchem_diagnostics_ mitgchem_diagnostics
#define landsource_ landsource
#endif 
                               