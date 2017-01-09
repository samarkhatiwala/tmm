/* $Header: /Users/ikriest/CVS/mops/mops_biogeochem_misfit_data.h,v 1.2 2016/06/03 09:28:59 ikriest Exp $ */
/* $Name: mops-2_0 $*/

/* IK : added for misfit function : start */

#ifdef EXTERNAL_FORCING
  #define ISEXTERN extern
#else 
  #define ISEXTERN
#endif  

ISEXTERN PetscBool averageCost; /* Evaluate misfit by annual averages? */
ISEXTERN StepTimer costTimer;
ISEXTERN PetscScalar *localmbgc1, *localmbgc2, *localmbgc3; /* The pointers to the local model types */
ISEXTERN Vec mbgc1,mbgc2,mbgc3,mbgc1avg,mbgc2avg,mbgc3avg; /* The individual model equivalents/misfits */
ISEXTERN Vec obgc1,obgc2,obgc3; /* The observation arrays */
ISEXTERN Vec cbgc1,cbgc2,cbgc3; /* The individual cost functions (e.g., RMS) */
ISEXTERN Vec globVolFrac; /* volume of each box, as fraction of total volume */
ISEXTERN Vec wbgc1,wbgc2,wbgc3; /* misfit weights for each individual box and tracer type*/
ISEXTERN Vec w1,w2,w3; /* array to store overall weight in */
ISEXTERN PetscScalar Gaveobgc1,Gaveobgc2,Gaveobgc3; /* global average observed PO4, O2, NO3, ... */
ISEXTERN PetscScalar Gavembgc1,Gavembgc2,Gavembgc3; /* global average simulated PO4, O2, NO3, ... */
ISEXTERN PetscScalar Gavecost1,Gavecost2,Gavecost3;  /* global cost function for PO4, O2, NO3 ... */
ISEXTERN PetscScalar RGavecost1,RGavecost2,RGavecost3;  /* square root of global cost function for PO4, O2, NO3 ... */
ISEXTERN PetscScalar Gcost;  /* global overall cost function, to be thrown at optimizer ... */
ISEXTERN FILE *misfitf; /* where to write the misfit to - prefer ASCII, for now */
ISEXTERN char misfitFile[PETSC_MAX_PATH_LEN];  
ISEXTERN PetscViewer fdcbgc1avg, fdcbgc2avg, fdcbgc3avg;

/* IK : added for misfit function : end */
