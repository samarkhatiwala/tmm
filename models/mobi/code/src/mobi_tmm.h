extern void mobi_copy_data_(PetscInt *lSize, PetscInt *numLocProfiles, PetscInt *itr, PetscScalar localTR[], PetscInt *direction);

extern void mobi_sed_copy_data_(PetscInt *numLocProfiles, PetscScalar localSedMixTR[], PetscScalar localSedBurTR[], PetscInt *direction);

extern void mobi_ini_(PetscInt *numTracers, PetscInt *lSize, PetscInt *numLocProfiles, PetscInt *nzmax, PetscInt nzarr[],
               PetscScalar zt[], PetscScalar drF[], PetscScalar *DeltaT, PetscScalar locallatitude[], PetscScalar localarea[],
               PetscScalar *Sglobavg,PetscScalar TRglobavg[],
#if defined O_carbon               
               PetscScalar *pCO2atm, 
# if defined O_carbon_13               
               PetscScalar *dc13atm, 
# endif     
# if defined O_carbon_14               
               PetscScalar *DC14atm, 
# endif
#endif /* O_carbon */
#if defined O_mobi
               PetscScalar localsgbathy[],
# if defined O_mobi_iron
               PetscScalar localFe_hydr[], 
# endif
# if defined O_mobi_silicon
               PetscScalar localSi_hydr[], 
# endif               
#endif
#if defined O_PaTh
               PetscScalar localPaTh_lith[], PetscScalar wLith[], 
# if !defined O_mobi
               PetscScalar localPaTh_pom[], PetscScalar localPaTh_caco3[], PetscScalar localPaTh_opal[], 
               PetscScalar wPOM[], PetscScalar wCaCO3[], PetscScalar wOpal[], 
# endif               
#endif
#if defined O_sed  
               PetscInt *numOceanStepsPerSedStep,
               PetscInt *nzmaxSed, PetscInt *ibmaxSed, PetscInt *numSedMixedTracers,                
               PetscInt *numSedBuriedTracers, PetscScalar *globalweathflx,
               PetscScalar *sedsa, PetscScalar sedmaskloc[],
#endif /* O_sed */         
               PetscInt *debugFlag);

extern void mobi_calc_(PetscInt *lSize, PetscInt *numLocProfiles, 
               PetscScalar *day, PetscScalar *relyr,
               PetscScalar localTs[], PetscScalar localSs[], PetscScalar TRglobavg[], 
               PetscScalar localEmP[], PetscScalar *empglobavg, 
#if defined BGC
               PetscScalar localwind[], PetscScalar localaice[], PetscScalar localhice[], PetscScalar localhsno[], 
#endif               
#if defined O_carbon
# if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
#  if defined O_carbon_co2_2d
               PetscScalar pCO2atm[],
#  else
               PetscScalar *pCO2atm,
#  endif               
# endif               
# if defined O_c14ccn_data
               PetscScalar *dc14ccnnhatm, PetscScalar *dc14ccneqatm, PetscScalar *dc14ccnshatm,
# endif
# if defined O_c13ccn_data || defined O_carbon_13_coupled
               PetscScalar *c13o2atm,
# endif
#endif /* O_carbon */
#if defined O_mobi
               PetscScalar localswrad[],
               PetscScalar *globdisch, PetscScalar localdisch[],
# if defined O_mobi_iron
               PetscScalar localFe_adep[], 
# endif
# if defined O_mobi_silicon
               PetscScalar localSi_dep[], PetscScalar *globSi_dep, 
# endif
#endif /* O_mobi */
#if defined O_PaTh
               PetscScalar localdust_adep[], 
#endif
#if defined O_carbon
               PetscScalar gasexfluxloc[], PetscScalar totfluxloc[], 
# if defined O_carbon_13_coupled
               PetscScalar c13gasexfluxloc[],
# endif
#endif
#if defined O_sed
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

