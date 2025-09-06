/* $Header: /Users/ikriest/CVS/mops/mops_biogeochem_misfit_data.h,v 1.2 2016/06/03 09:28:59 ikriest Exp $ */
/* $Name: mops-2_0 $*/

/* IK : added for misfit function : start */

#ifdef EXTERNAL_FORCING
  #define ISEXTERN extern
#else 
  #define ISEXTERN
#endif  

ISEXTERN PetscBool multiObjective; /* Evaluate misfit by annual averages, boolean operator for second objective */
ISEXTERN StepTimer misfitTimer;

ISEXTERN Vec globVolFrac; /* volume of each box, as fraction of total volume */
ISEXTERN PetscScalar RMSEMetric,OMZMetric,HDMetric;  /* three global overall cost functions, to be thrown at optimizer ... */
ISEXTERN FILE *misfitf; /* where to write the misfit to - prefer ASCII, for now */
ISEXTERN char misfitFile[PETSC_MAX_PATH_LEN];  

/* IK : RMSE for PO4, NO3, O2 */ 

ISEXTERN Vec mbgc1avg,mbgc2avg,mbgc3avg; /* The individual model equivalents/misfits */
ISEXTERN Vec obgc1,obgc2,obgc3; /* The observation arrays */
ISEXTERN Vec cbgc1,cbgc2,cbgc3; /* The individual cost functions (e.g., RMS) */
ISEXTERN Vec wbgc1,wbgc2,wbgc3; /* misfit weights for each individual box and tracer type*/
ISEXTERN Vec w1,w2,w3; /* array to store overall weight in */
ISEXTERN PetscScalar VolSum1,VolSum2,VolSum3; /* global fractional volumes for computation of average */
ISEXTERN PetscScalar NoSum1,NoSum2,NoSum3; /* global fractional volumes for computation of average */
ISEXTERN PetscScalar Gaveobgc1,Gaveobgc2,Gaveobgc3; /* global average observed PO4, O2, NO3, ... */
ISEXTERN PetscScalar Gavembgc1,Gavembgc2,Gavembgc3; /* global average simulated PO4, O2, NO3, ... */
ISEXTERN PetscScalar Gavecost1,Gavecost2,Gavecost3;  /* global cost function for PO4, O2, NO3 ... */
ISEXTERN PetscScalar RGavecost1,RGavecost2,RGavecost3;  /* square root of global cost function for PO4, O2, NO3 ... */
ISEXTERN PetscScalar Ginorgcost;  /* global overall cost function, to be thrown at optimizer ... */

/* IK : RMSE for PHY, ZOO, DET, DOP */

ISEXTERN Vec mbgc4avg,mbgc5avg,mbgc6avg,mbgc7avg; /* The individual model equivalents/misfits */
ISEXTERN Vec obgc4,obgc5,obgc6,obgc7; /* The observation arrays */
ISEXTERN Vec cbgc4,cbgc5,cbgc6,cbgc7; /* The individual cost functions (e.g., RMS) */
ISEXTERN Vec wbgc4,wbgc5,wbgc6,wbgc7; /* misfit weights for each individual box and tracer type*/
ISEXTERN Vec w4,w5,w6,w7; /* array to store overall weight in */
ISEXTERN PetscScalar VolSum4,VolSum5,VolSum6,VolSum7; /* global fractional volumes for computation of average */
ISEXTERN PetscScalar NoSum4,NoSum5,NoSum6,NoSum7; /* global fractional volumes for computation of average */
ISEXTERN PetscScalar Gaveobgc4,Gaveobgc5,Gaveobgc6,Gaveobgc7; /* global average observed PHY, ZOO, DET, DOP */
ISEXTERN PetscScalar Gavembgc4,Gavembgc5,Gavembgc6,Gavembgc7; /* global average simulated PHY, ZOO, DET, DOP */
ISEXTERN PetscScalar Gavecost4,Gavecost5,Gavecost6,Gavecost7;  /* global cost function for PHY, ZOO, DET, DOP */
ISEXTERN PetscScalar RGavecost4,RGavecost5,RGavecost6,RGavecost7;  /* square root of global cost function for PHY, ZOO, DET, DOP */
ISEXTERN PetscScalar Gorgcost;  /* global overall cost function, to be thrown at optimizer */
 
/* IK : OMZ metric */

ISEXTERN Vec onevec, momz, oomz, amomz, aoomz, imomz, ioomz,imoomz;
ISEXTERN PetscScalar GMomz, GOomz, GMOomz;

/* IK : HDMetric for all */ 
ISEXTERN PetscScalar hd1, hd2, hd3, hd4, hd5, hd6, hd7; /* VS HDs for the tracers and the sum of them */  

