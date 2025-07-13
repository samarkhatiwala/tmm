typedef struct _p_ExternalForcingCtx *ExternalForcingContext;
struct _p_ExternalForcingCtx {
  PETSCHEADER(int);
  PetscInt efctxId;
  PetscInt stateId;
/* Add problem-specific variables below */  
  PetscScalar **localTR, **localJTR;
  PetscScalar *TRglobavg;

#if defined O_carbon
PetscScalar *localgasexflux, *localtotflux;
PetscScalar pCO2atm; /* this is initialized to the namelist value */

/* variables for prescribed atmospheric CO2 */
# if defined O_co2ccn_data
char *pCO2atmFiles[2];  
PetscInt numpCO2atm_hist;
PetscScalar *TpCO2atm_hist, *pCO2atm_hist;
# endif

# if defined O_carbon_13
PetscScalar dc13atm; /* this is initialized to the namelist value */
PetscScalar pC13O2atm; /* mole fraction of atmospheric 13CO2; this is only used if O_c13ccn_data or O_carbon_13_coupled are defined */

#  if defined O_c13ccn_data
char *C13atmFiles[2];
PetscScalar *TC13atm_hist, *C13atm_hist;
PetscInt numC13atm_hist;
#  endif
# endif /* O_carbon_13 */

# if defined O_carbon_14
PetscScalar DC14atm; /* this is initialized to the namelist value */
#  if defined O_c14ccn_data
PetscBool timeDependentAtmosphericC14;
PetscScalar dc14ccnnhatm;
PetscScalar dc14ccneqatm;
PetscScalar dc14ccnshatm;
char *C14atmFiles[4];
PetscScalar *TC14atm_hist, *C14atmnh_hist, *C14atmeq_hist, *C14atmsh_hist;
PetscInt numC14atm_hist;
PetscScalar dc14atmvals[3];
#  endif /* O_c14ccn_data */
# endif /* O_carbon_14 */

# if defined O_TMM_interactive_atmosphere
/* Land/Atm model variables */
PetscBool useLandModel;
PetscBool useEmissions;
PetscBool interpEmissions;
char *emFiles[3];  
PetscInt numEmission_hist;
PetscScalar *Tem_hist, *E_hist, *D_hist;
PetscScalar fossilFuelEmission, landUseEmission;
PetscScalar cumulativeEmissions;
char pCO2atmIniFile[PETSC_MAX_PATH_LEN];
PetscScalar landState[3], landSource[3];
PetscScalar deltaTsg;
PetscScalar atmModelDeltaT;
PetscScalar Fland, Focean;

StepTimer atmWriteTimer;
PetscBool atmAppendOutput;
FILE *atmfptime;
char atmOutTimeFile[PETSC_MAX_PATH_LEN];  
PetscScalar pCO2atmavg, Foceanavg, Flandavg, landStateavg[3];

#  if defined O_carbon_13_coupled
char dc13atmIniFile[PETSC_MAX_PATH_LEN];
PetscScalar *localc13gasexflux;
PetscScalar F13ocean;
#  endif /* O_carbon_13_coupled */

# endif /* O_TMM_interactive_atmosphere */

#endif /* O_carbon */

#if defined O_sed
PetscInt numOceanStepsPerSedStep;
PetscInt nzmaxSed, ibmaxSed, numSedMixedTracers, numSedBuriedTracers;
Vec sedMixedTracers, sedBuriedTracers;
PetscScalar *localSedMixedTR, *localSedBuriedTR;
PetscInt *gSedMixedIndices, *gSedBuriedIndices, gSedLow;
PetscInt sedMixedBlockSize, lSedMixedSize, sedBuriedBlockSize, lSedBuriedSize;
PetscScalar globalweathflx; 
PetscScalar localsedsa, globalsedsa;
PetscScalar *localsedmask;
char sedMixedIniFile[PETSC_MAX_PATH_LEN], sedBuriedIniFile[PETSC_MAX_PATH_LEN];
char sedMixedOutFile[PETSC_MAX_PATH_LEN], sedBuriedOutFile[PETSC_MAX_PATH_LEN];
char sedMixedPickupOutFile[PETSC_MAX_PATH_LEN], sedBuriedPickupOutFile[PETSC_MAX_PATH_LEN];
char sedOutTimeFile[PETSC_MAX_PATH_LEN];
PetscFileMode SED_OUTPUT_FILE_MODE;
StepTimer sedWriteTimer;
PetscBool sedAppendOutput;
PetscViewer sedmixedfd, sedburiedfd;
FILE *sedfptime;
#endif /* O_sed */

PetscBool calcDiagnostics;
StepTimer diagTimer;
PetscBool diagFirstTime;
PetscBool diagAppendOutput;
PetscFileMode DIAG_FILE_MODE;
FILE *diagfptime;
PetscViewer diagfd;
int diagfp;
char diagOutTimeFile[PETSC_MAX_PATH_LEN];
/* Add model specific diagnostic variables below */
PetscScalar *localDiag2davg[MAXDIAGS2d];
Vec *Diag3davg;
PetscScalar **localDiag3davg;
PetscInt numDiags2d, numDiags3d;
PetscViewer fddiag3dout[MAXDIAGS3d];
char *Diag2dFile[MAXDIAGS2d], *Diag3dFile[MAXDIAGS3d];
PetscInt doAverage;
};
