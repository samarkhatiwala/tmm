extern void mobi_copy_data_(PetscInt *nzloc, PetscInt *itr, PetscScalar localTR[], PetscInt *direction);

extern void mobi_ini_(PetscScalar zt[], PetscScalar drF[], PetscScalar *DeltaT, 
               PetscScalar *Sglobavg,PetscScalar TRglobavg[], PetscInt *debugFlag);

extern void mobi_calc_(PetscInt *nzloc, PetscScalar *locallatitude, PetscScalar *day, PetscScalar *relyr,
               PetscScalar localTs[], PetscScalar localSs[], PetscScalar TRglobavg[],
# if defined O_carbon
               PetscScalar *pCO2atm, PetscScalar *localwind,
#endif      
#  if defined O_npzd_nitrogen
               PetscScalar localsgbathy[],
#  endif
#  if defined O_npzd_fe_limitation
               PetscScalar localFe[],
#  endif
#  if defined O_embm
               PetscScalar *localswrad,
#  endif
#  if defined O_ice
#   if !defined O_ice_cpts
               PetscScalar *localaice, PetscScalar *localhice, PetscScalar *localhsno,
#   endif
#  endif
               PetscScalar *localEmP, PetscScalar *empglobavg, PetscInt *debugFlag);

#if !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) 
#define mobi_copy_data_ mobi_copy_data
#define mobi_ini_ mobi_ini
#define mobi_calc_ mobi_calc
#endif 

