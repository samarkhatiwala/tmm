extern void mobi_copy_data_(PetscInt *lSize, PetscInt *numLocProfiles, PetscInt *itr, PetscScalar localTR[], PetscInt *direction);

extern void mobi_sed_copy_data_(PetscInt *numLocProfiles, PetscScalar localSedMixTR[], PetscScalar localSedBurTR[], PetscInt *direction);

extern void mobi_ini_(PetscInt *numTracers, PetscInt *lSize, PetscInt *numLocProfiles, PetscInt *nzmax, PetscInt nzarr[],
               PetscScalar zt[], PetscScalar drF[], PetscScalar *DeltaT, PetscScalar locallatitude[], PetscScalar localarea[],
               PetscScalar localsgbathy[],
               PetscScalar *Sglobavg,PetscScalar TRglobavg[],
               PetscScalar *pCO2atm, PetscScalar *dc13atm, PetscScalar *DC14atm,                 
#ifdef O_sed  
               PetscInt *numOceanStepsPerSedStep,
               PetscInt *nzmaxSed, PetscInt *ibmaxSed, PetscInt *numSedMixedTracers,                
               PetscInt *numSedBuriedTracers, PetscScalar *globalweathflx,
               PetscScalar *sedsa, PetscScalar sedmaskloc[],
#endif               
               PetscInt *debugFlag);

extern void mobi_calc_(PetscInt *lSize, PetscInt *numLocProfiles, 
               PetscScalar *day, PetscScalar *relyr,
               PetscScalar localTs[], PetscScalar localSs[], PetscScalar TRglobavg[], 
               PetscScalar localdz[], PetscScalar zt[],
# if defined O_carbon
#if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
#   if defined O_carbon_co2_2d
               PetscScalar pCO2atm[],
#   else
               PetscScalar *pCO2atm,
#   endif               
#endif               
#if defined O_c14ccn_data
               PetscScalar *dc14ccnnhatm, PetscScalar *dc14ccneqatm, PetscScalar *dc14ccnshatm,
#endif
#if defined O_c13ccn_data || defined O_carbon_13_coupled
               PetscScalar *c13o2atm,
#endif
               PetscScalar localwind[],
#endif
#  if defined O_npzd_fe_limitation
               PetscScalar localFe_dissolved[],
#  endif
#ifdef O_npzd_iron
               PetscScalar localFe_adep[], PetscScalar localFe_detr_flux[], PetscScalar localFe_hydr[],
#endif
#  if defined O_embm
               PetscScalar localswrad[],
#  endif
#  if defined O_ice
#   if !defined O_ice_cpts
               PetscScalar localaice[], PetscScalar localhice[], PetscScalar localhsno[],
#   endif
#  endif
               PetscScalar localEmP[], PetscScalar *empglobavg, 
# if defined O_carbon
               PetscScalar gasexfluxloc[], PetscScalar totfluxloc[], 
#if defined O_carbon_13_coupled
               PetscScalar c13gasexfluxloc[],
#endif
# endif
#if defined O_sed
# if defined O_embm
               PetscScalar *globdisch, PetscScalar localdisch[],
# endif
               PetscInt *timeToRunSedModel,
               PetscScalar *globwflx, PetscScalar localwflx[],
#endif
               PetscInt *debugFlag);

extern void mobi_diags_ini_(PetscInt *lNumProfiles, PetscInt *lTotNumPoints, PetscInt *lNum2dDiags,
               PetscInt *lNum3dDiags, PetscInt *debugFlag);

extern void mobi_diags_start_(PetscInt *debugFlag);

extern void mobi_diags_stop_(PetscInt *debugFlag);

extern void mobi_diags_accumulate_(PetscInt *numAvg, PetscInt *avgFlag, PetscInt *debugFlag);

extern void mobi_diags2d_copy_(PetscInt *id, PetscScalar diagArr[], char *fname, PetscInt *debugFlag);

extern void mobi_diags3d_copy_(PetscInt *id, PetscScalar diagArr[], char *fname, PetscInt *debugFlag);

extern void mobi_diags_finalize_(PetscInt *debugFlag);

#if !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) 
#define mobi_copy_data_ mobi_copy_data
#define mobi_sed_copy_data_ mobi_sed_copy_data
#define mobi_ini_ mobi_ini
#define mobi_start_ mobi_start
#define mobi_stop_ mobi_stop
#define mobi_calc_ mobi_calc
#define mobi_diags_ini_ mobi_diags_ini
#define mobi_diags_accumulate_ mobi_diags_accumulate
#define mobi_diags2d_copy_ mobi_diags2d_copy
#define mobi_diags3d_copy_ mobi_diags3d_copy
#define mobi_diags_finalize_ mobi_diags_finalize
#endif 

