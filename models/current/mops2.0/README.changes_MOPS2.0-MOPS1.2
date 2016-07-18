IK, 2016-03-09

This  file documents the changes made to MOPS-Version-1_2, for its 
efficient use within an optimization environment. 
See file README to see how the code works.

******************************************************************************************
******************************************************************************************
Changes to source code files:
******************************************************************************************
******************************************************************************************

NEW FILES:

BGC_MISFIT.h
************
mops_biogeochem_misfit.F
************************
mops_biogeochem_misfit_data.h
*****************************
tmm_misfit.c
************
tmm_misfit.h
************

CHANDGED FILES:

BGC_INI.F
*********
Changed order / position of two parameter definitions to account for adjustment or
----------------------------------------------------------------------------------
non-adjustment during optimization:
-----------------------------------

Moved definition of rhno3ut downwards, to be changed dynamically whenever ro2ut changes during optimization:

         rhno3ut = 0.8d0*ro2ut - rnp ! -HNO3:P ratio for denitrification
 
Moved definition of ACkphy upwards, to stay constant whenever ACmuzoo changes during optimization:
(Note: this varies the initial slope - i.e., the response of zooplankton to changes in food.)

         ACkphy=SQRT(ACmuzoo/1.0d0)/rnp !zoo half-saturation constant [mmol P]


BGC_MODEL.F:
************
Added lines for computation and communication of misfit:
 -------------------------------------------------------
#include "BGC_MISFIT.h"

...
       DO K=1,bgc_kloc
 
         m1_out(k)=0.0d0 
         m2_out(k)=0.0d0 
         m3_out(k)=0.0d0 

       ENDDO
...       
! For now, I am happy to take tracer concentration at the end of each time loop
! for computation of misfit; I might as well sum it up and divide by the number
! of time steps

      DO K=1,bgc_kloc
      
        m1_out(k) = bgc_tracer(k,ipo4)
        m2_out(k) = bgc_tracer(k,ioxy)
        m3_out(k) = bgc_tracer(k,idin)

      ENDDO

external_forcing_mops_biogeochem.c
***********************************
Added/changed lines for computation of misfit function:
-------------------------------------------------------

#define EXTERNAL_FORCING
...
#include "mops_biogeochem_misfit_data.h"
...
#ifdef ASCIIPARAMS
PetscScalar bgcparams[20];
FILE *fpparams;
#else
...
#endif
...
/* IK : added for misfit function : start */
/* IK : added for misfit function : start */
/* IK : added for misfit function : start */	
      if (Iter0+iLoop>=costStartTimeStep) { /* start time averaging (note: costStartTimeStep is ABSOLUTE time step) */	
          mops_biogeochem_misfit_(&nzloc,&localmbgc1[kl],&localmbgc2[kl],&localmbgc3[kl]);
      }
/* IK : added for misfit function : end */
/* IK : added for misfit function : end */
/* IK : added for misfit function : end */
...
/* IK : added for misfit function : start */
/* IK : added for misfit function : start */
/* IK : added for misfit function : start */ 
 	if (Iter0+iLoop>=costStartTimeStep) { /* start time averaging (note: costStartTimeStep is ABSOLUTE time step) */  
 	  ierr = VecSetValues(mbgc1,lSize,gIndices,localmbgc1,INSERT_VALUES);CHKERRQ(ierr);
 	  ierr = VecAssemblyBegin(mbgc1);CHKERRQ(ierr);
 	  ierr = VecAssemblyEnd(mbgc1);CHKERRQ(ierr);    
   
 	  ierr = VecSetValues(mbgc2,lSize,gIndices,localmbgc2,INSERT_VALUES);CHKERRQ(ierr);
 	  ierr = VecAssemblyBegin(mbgc2);CHKERRQ(ierr);
 	  ierr = VecAssemblyEnd(mbgc2);CHKERRQ(ierr);    
   
 	  ierr = VecSetValues(mbgc3,lSize,gIndices,localmbgc3,INSERT_VALUES);CHKERRQ(ierr);
 	  ierr = VecAssemblyBegin(mbgc3);CHKERRQ(ierr);
 	  ierr = VecAssemblyEnd(mbgc3);CHKERRQ(ierr);    
         } 
/* IK : added for misfit function : end */
/* IK : added for misfit function : end */
/* IK : added for misfit function : end */

Added/changed lines for reading of parameters from file:
--------------------------------------------------------
PetscInt it,ipar;
...
ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading %d parameters from file\n",numBGCParams);CHKERRQ(ierr);
#ifdef ASCIIPARAMS
    fpparams=fopen(bgcParamsFile,"r");
    for(ipar=0;ipar<numBGCParams;ipar++) {
      ierr = fscanf(fpparams,"%lf\n",&bgcparams[ipar]);
    }   
    fclose(fpparams);
#else
...
#endif
...
    for (ipar=0; ipar<numBGCParams; ipar++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Parameter no. %d is %f\n",ipar,bgcparams[ipar]);CHKERRQ(ierr);
    }
...
PetscInt ipar;
...
#ifdef ASCIIPARAMS
    fpparams=fopen(bgcParamsFile,"r");
    for(ipar=0;ipar<numBGCParams;ipar++) {
      ierr = fscanf(fpparams,"%lf\n",&bgcparams[ipar]);
    }   
    fclose(fpparams);
#else
...
#endif
... 
    for (ipar=0; ipar<numBGCParams; ipar++) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Parameter no. %d is %f\n",ipar,bgcparams[ipar]);CHKERRQ(ierr);
    }


Added/changed lines for flexible name assignment for runoff pickup:
-------------------------------------------------------------------
char runoffIniOutFile[PETSC_MAX_PATH_LEN];
...
ierr = PetscPrintf(PETSC_COMM_WORLD,"Runoff will be integrated over every %d time steps\n",burialSumSteps);CHKERRQ(ierr);
...
ierr = PetscOptionsGetString(PETSC_NULL,"-pickup_runoff_out",runoffIniOutFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
 if (!flg) {
   strcpy(runoffIniOutFile,"");
   sprintf(runoffIniOutFile,"%s","pickup_runoff.bin");
 }
 ierr = PetscPrintf(PETSC_COMM_WORLD,"Final runoff output will be written to %s\n",runoffIniOutFile);CHKERRQ(ierr);
...
 ierr = writeBinaryScalarData(runoffIniOutFile,&GRunoff,1,PETSC_FALSE);

Added compile option to not by default write out daily runoff:
--------------------------------------------------------------
#ifdef WRITE_RUNOFF
...
#endif
...
#ifdef WRITE_RUNOFF
...
#endif
...
#ifdef WRITE_RUNOFF
...
#endif
...
#ifdef WRITE_RUNOFF
...
#endif

mops_biogeochem.h
*****************
Added lines for computation and communication of misfit:

extern void mops_biogeochem_misfit_(PetscInt *Nrloc,
                                         PetscScalar localmbgc1[], PetscScalar localmbgc2[], PetscScalar localmbgc3[]);
...
#define mops_biogeochem_misfit_ mops_biogeochem_misfit


mops_biogeochem_set_params.F
****************************
Changed lines for reading of parameters from file:

! Third version for Volkmar's optimization tests, 16 October 2015
         ro2ut = parvec(1)
         ACik = parvec(2)
         ACkpo4 = parvec(3)
         ACmuzoo =parvec(4)
         AComniz = parvec(5)
         detmartin = parvec(6)
         
