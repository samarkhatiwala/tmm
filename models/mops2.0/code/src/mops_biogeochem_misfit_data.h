/* $Header: /Users/ikriest/CVS/mops/mops_biogeochem_misfit_data.h,v 1.2 2016/06/03 09:28:59 ikriest Exp $ */
/* $Name: mops-2_0 $*/
/* Edited by Sophy Oliver Jul 2018 */

/* IK : added for misfit function : start */

#ifdef EXTERNAL_FORCING
  #define ISEXTERN extern
#else 
  #define ISEXTERN
#endif  

ISEXTERN PetscBool averageCost; /* Evaluate misfit by annual averages? */
ISEXTERN StepTimer costTimer;
ISEXTERN PetscScalar *localmbgc1, *localmbgc2, *localmbgc3; /* The pointers to the local model types */
ISEXTERN Vec boxVol; 
ISEXTERN PetscInt *nRegions, maxValsToRead; 
ISEXTERN Vec **mask, **fracMaskVol, *obsbgc, *noise, *obsbgc_norm, *obsbgc_negs, *modbgc, *avmodbgc, *avmodbgc_common, *costbgc;
ISEXTERN PetscScalar **totMaskVol, **fracTotMaskVol, **Gavembgc, **costbgcw1, **costbgcw2, **obsbgcw1, **obsbgcw2, **costnorm; /* global average observed PO4, O2, NO3, ... */
ISEXTERN FILE *misfitf; /* where to write the misfit to - prefer ASCII, for now */
ISEXTERN char misfitFile[PETSC_MAX_PATH_LEN]; 
ISEXTERN PetscViewer fdcbgcavg[3];
ISEXTERN PetscViewer fdbgcavg[3];
ISEXTERN PetscInt iobs, ireg;
ISEXTERN Mat L;
ISEXTERN PetscBool doLoperator, addNoise;
ISEXTERN Vec trTemplateVec;

/* IK : added for misfit function : end */
