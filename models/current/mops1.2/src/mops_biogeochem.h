/* $Header: /Users/ikriest/CVS/mops/mops_biogeochem.h,v 1.1.1.1 2015/06/03 17:02:09 ikriest Exp $ */
/* $Name: mops-1_2 $ */

extern void mops_biogeochem_copy_data_(PetscInt *nzloc, PetscInt *itr, PetscScalar localTR[], PetscScalar localJTR[], 
                                PetscScalar *DeltaT, PetscInt *direction);

extern void insolation_(PetscInt *N, PetscScalar *myTime, PetscScalar locallatitude[], PetscScalar localswrad[], PetscScalar localstau[]);

extern void mops_biogeochem_ini_(PetscInt *Nrloc, PetscScalar *DeltaT, 
#ifdef CARBON                      
                                   PetscScalar *localph,
#endif
                                   PetscScalar localTs[], PetscScalar localSs[], 
                                   PetscScalar localdz[], PetscScalar drF[], PetscInt *nzmax, PetscInt *nzeuph,
                                   PetscInt *numBiogeochemStepsPerOceanStep,
                                   PetscBool *setDefaults);

extern void mops_biogeochem_model_(PetscInt *Nrloc, PetscScalar *DeltaT,
#ifdef CARBON                      
				   PetscScalar *DICglobavg, PetscScalar *ALKglobavg, PetscScalar *localEmP, PetscScalar *localpCO2atm,
#endif
                                   PetscScalar localTs[],PetscScalar localSs[], 
                                   PetscScalar *localfice, PetscScalar *localswrad, PetscScalar *localstau,
                                   PetscScalar *localwind, PetscScalar *localatmosp, PetscScalar localdz[], 
#ifdef CARBON                      
                                   PetscScalar *localph,
				   PetscScalar *localco2airseaflux,
#endif
                                   PetscScalar *localburial, PetscScalar *GRunoff, PetscScalar localrunoffvol[],
                                   PetscBool *useSeparateBiogeochemTimeStepping);

extern void mops_biogeochem_diagnostics_(PetscInt *Nrloc, 
                                         PetscScalar localfbgc1[], PetscScalar localfbgc2[], PetscScalar localfbgc3[], 
					 PetscScalar localfbgc4[], PetscScalar localfbgc5[], PetscScalar localfbgc6[], PetscScalar localfbgc7[]);

extern void mops_biogeochem_set_params_(PetscInt *numbgcparams, PetscScalar bgcparams[]);

#if !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) 
#define mops_biogeochem_copy_data_ mops_biogeochem_copy_data
#define insolation_ insolation
#define mops_biogeochem_ini_ mops_biogeochem_ini
#define mops_biogeochem_model_ mops_biogeochem_model
#define mops_biogeochem_diagnostics_ mops_biogeochem_diagnostics
#define mops_biogeochem_set_params_ mops_biogeochem_set_params
#endif 

