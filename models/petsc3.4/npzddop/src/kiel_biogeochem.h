extern void insolation_(PetscInt *N, PetscScalar *myTime, PetscScalar locallatitude[], PetscScalar localswrad[], PetscScalar localstau[]);

extern void kiel_biogeochem_ini_(PetscInt *Nrloc, PetscScalar *DeltaT, 
                                   PetscScalar localTR1[], PetscScalar localTR2[], PetscScalar localTR3[],
                                   PetscScalar localTR4[], PetscScalar localTR5[], PetscScalar localTR6[],                                  
#ifdef CARBON                      
                                   PetscScalar localTR7[], PetscScalar localTR8[],PetscScalar *localph,
#endif
                                   PetscScalar localTs[], PetscScalar localSs[], 
                                   PetscScalar localdz[], PetscScalar drF[], PetscInt *nzmax, PetscInt *nzeuph,
                                   PetscInt *numBiogeochemStepsPerOceanStep,
                                   PetscBool *setDefaults);

extern void kiel_biogeochem_model_(PetscInt *Nrloc, PetscScalar *DeltaT,
                                   PetscScalar localTR1[], PetscScalar localTR2[], PetscScalar localTR3[],
                                   PetscScalar localTR4[], PetscScalar localTR5[], PetscScalar localTR6[],                                   
#ifdef CARBON                      
                                   PetscScalar localTR7[], PetscScalar localTR8[],
				   PetscScalar *DICglobavg, PetscScalar *localEmP, PetscScalar *localpCO2atm,
#endif
                                   PetscScalar localTs[],PetscScalar localSs[], 
                                   PetscScalar *localfice, PetscScalar *localswrad, PetscScalar *localstau,
                                   PetscScalar *localwind, PetscScalar *localatmosp, PetscScalar localdz[], 
                                   PetscScalar localJTR1[],PetscScalar localJTR2[],PetscScalar localJTR3[],
                                   PetscScalar localJTR4[],PetscScalar localJTR5[],PetscScalar localJTR6[],
#ifdef CARBON                      
                                   PetscScalar localJTR7[], PetscScalar localJTR8[], PetscScalar *localph,
				   PetscScalar *localco2airseaflux,
#endif
                                   PetscBool *useSeparateBiogeochemTimeStepping);

extern void kiel_biogeochem_diagnostics_(PetscInt *Nrloc, 
                                         PetscScalar localfbgc1[], PetscScalar localfbgc2[], PetscScalar localfbgc3[], 
					 PetscScalar localfbgc4[], PetscScalar localfbgc5[]);

extern void kiel_biogeochem_set_params_(PetscInt *numbgcparams, PetscScalar bgcparams[]);

#if !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) 
#define insolation_ insolation
#define kiel_biogeochem_ini_ kiel_biogeochem_ini
#define kiel_biogeochem_model_ kiel_biogeochem_model
#define kiel_biogeochem_diagnostics_ kiel_biogeochem_diagnostics
#define kiel_biogeochem_set_params_ kiel_biogeochem_set_params
#endif 

