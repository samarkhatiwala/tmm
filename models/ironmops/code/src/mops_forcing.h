typedef struct _p_ExternalForcingCtx *ExternalForcingContext;
struct _p_ExternalForcingCtx {
  PETSCHEADER(int);
  PetscInt efctxId;
  PetscInt stateId;
/* Add problem-specific variables below */  
  PetscScalar **localTR, **localJTR;
#ifdef CARBON
  PetscScalar *localph;
#endif

  PetscBool readBGCParams;

#ifdef CARBON
/* atmospheric model variables */
  PetscInt numpCO2atm_hist;
  PetscScalar *TpCO2atm_hist, *pCO2atm_hist;
  PetscBool fixedAtmosCO2;

  PetscBool useAtmModel;
  PetscScalar pCO2atm; /* default initial value */
  PetscScalar atmModelDeltaT;
  PetscScalar Focean;
  PetscScalar Foceanint;

  StepTimer atmWriteTimer;
  PetscBool atmAppendOutput;
  FILE *atmfptime;
  PetscViewer atmfd;
  PetscInt atmfp;
#endif

  PetscScalar GRunoff; /* Global runoff, calculated from burial */
  PetscScalar localFburial;
  PetscScalar Fburial;
  PetscInt burialSumSteps;
  char runoffOutTimeFile[PETSC_MAX_PATH_LEN];  
  char runoffIniOutFile[PETSC_MAX_PATH_LEN];  
  FILE *runofffptime;

  PetscBool calcDiagnostics;
  StepTimer diagTimer;
  PetscBool appendDiagnostics;
/* Add model specific diagnostic variables below */
  Vec fbgc1, fbgc2, fbgc3, fbgc4, fbgc5, fbgc6, fbgc7, fbgc8, fbgc9, fbgc10, fbgc11, fbgc12;
  Vec fbgc1avg, fbgc2avg, fbgc3avg, fbgc4avg, fbgc5avg, fbgc6avg, fbgc7avg, fbgc8avg, fbgc9avg, fbgc10avg, fbgc11avg, fbgc12avg;
  PetscViewer fdfbgc1avg, fdfbgc2avg, fdfbgc3avg, fdfbgc4avg, fdfbgc5avg, fdfbgc6avg, fdfbgc7avg, fdfbgc8avg, fdfbgc9avg, fdfbgc10avg, fdfbgc11avg, fdfbgc12avg;
  PetscScalar *localfbgc1, *localfbgc2, *localfbgc3, *localfbgc4, *localfbgc5, *localfbgc6, *localfbgc7, *localfbgc8, *localfbgc9, *localfbgc10, *localfbgc11, *localfbgc12;
  PetscInt numDiag;

#ifdef CARBON
  PetscScalar *localco2airseafluxdiag, *localco2airseafluxdiagavg;
  char co2airseafluxFile[PETSC_MAX_PATH_LEN];
#endif
};
