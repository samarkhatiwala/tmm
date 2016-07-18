/* $Header: /Users/ikriest/CVS/mops/tmm_misfit.c,v 1.1 2015/11/17 14:18:51 ikriest Exp $ */
/* $Name: mops-2_0 $*/

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "petsc_matvec_utils.h"
#include "tmm_main.h"
#include "tmm_forcing_utils.h"
#include "tmm_profile_utils.h"
#include "tmm_profile_data.h"

#include "mops_biogeochem_misfit_data.h"

#define TR v[0]

#undef __FUNCT__
#define __FUNCT__ "iniMisfit"
PetscErrorCode iniMisfit(PetscScalar tc, PetscInt Iter, PetscInt numTracers, Vec *v)
{

  PetscErrorCode ierr;
  PetscBool flg;
  PetscScalar zero = 0.0;
  PetscViewer fd;

/* Add your code here */

   averageCost = PETSC_FALSE;
   
/* IK: added for misfit function - start */
 
/* Type of cost function; so far, only annual averages to PO4, O2 and NO3 available */

  ierr = PetscOptionsGetString(PETSC_NULL,"-misfit_file",misfitFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) {
  strcpy(misfitFile,"");
  sprintf(misfitFile,"%s","misfit.txt");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"misfit will be written to %s\n",misfitFile);CHKERRQ(ierr);
  ierr = PetscFOpen(PETSC_COMM_WORLD,misfitFile,"w",&misfitf);CHKERRQ(ierr);  

  ierr = PetscOptionsHasName(PETSC_NULL,"-average_cost",&averageCost);CHKERRQ(ierr);
	if (averageCost) {
	ierr = PetscOptionsGetInt(PETSC_NULL,"-cost_start_time_step",&costStartTimeStep,&flg);CHKERRQ(ierr);
	if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate (absolute) time step at which to start calculating cost function -cost_start_time_step flag");
	ierr = PetscOptionsGetInt(PETSC_NULL,"-cost_time_steps",&costNumTimeSteps,&flg);CHKERRQ(ierr);
	if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate number of time averaging cost time steps with the -cost_time_steps flag");
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Cost for misfit function will be computed starting at (and including) time step: %d\n", costStartTimeStep);CHKERRQ(ierr);	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Cost for misfit function will be computed over %d time steps\n", costNumTimeSteps);CHKERRQ(ierr);	

/* Read fractional volume of each box for scaling according to box size */

	ierr = VecDuplicate(TR,&globVolFrac);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"volume_fraction.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(globVolFrac,fd);CHKERRQ(ierr);  /* IntoVector */ 
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Read volume file for weighting\n");CHKERRQ(ierr);	

/* Read the observations */ 

  ierr = VecDuplicate(TR,&obgc1);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&obgc2);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&obgc3);CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"po4_obs.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = VecLoad(obgc1,fd);CHKERRQ(ierr);  /* IntoVector */ 
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Read PO4 observations\n");CHKERRQ(ierr);	

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"o2_obs.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = VecLoad(obgc2,fd);CHKERRQ(ierr);  /* IntoVector */ 
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Read O2 observations\n");CHKERRQ(ierr);	

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"no3_obs.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = VecLoad(obgc3,fd);CHKERRQ(ierr);  /* IntoVector */ 
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Read NO3 observations\n");CHKERRQ(ierr);	

/* Initialize the cost function vectors */
	
  ierr = VecDuplicate(TR,&mbgc1);CHKERRQ(ierr);
  ierr = VecSet(mbgc1,zero);CHKERRQ(ierr);
  ierr = VecGetArray(mbgc1,&localmbgc1);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&mbgc1avg);CHKERRQ(ierr);
  ierr = VecSet(mbgc1avg,zero);CHKERRQ(ierr);

  ierr = VecDuplicate(TR,&mbgc2);CHKERRQ(ierr);
  ierr = VecSet(mbgc2,zero);CHKERRQ(ierr);
  ierr = VecGetArray(mbgc2,&localmbgc2);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&mbgc2avg);CHKERRQ(ierr);
  ierr = VecSet(mbgc2avg,zero);CHKERRQ(ierr);

  ierr = VecDuplicate(TR,&mbgc3);CHKERRQ(ierr);
  ierr = VecSet(mbgc3,zero);CHKERRQ(ierr);
  ierr = VecGetArray(mbgc3,&localmbgc3);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&mbgc3avg);CHKERRQ(ierr);
  ierr = VecSet(mbgc3avg,zero);CHKERRQ(ierr);

  ierr = VecDuplicate(TR,&cbgc1);CHKERRQ(ierr);
  ierr = VecSet(cbgc1,zero);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&cbgc2);CHKERRQ(ierr);
  ierr = VecSet(cbgc2,zero);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&cbgc3);CHKERRQ(ierr);
  ierr = VecSet(cbgc3,zero);CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"cbgc1.petsc",FILE_MODE_WRITE,&fdcbgc1avg);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"cbgc2.petsc",FILE_MODE_WRITE,&fdcbgc2avg);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"cbgc3.petsc",FILE_MODE_WRITE,&fdcbgc3avg);CHKERRQ(ierr);

   costCount=0;

        } else {

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Only average cost function available, currently.\n");CHKERRQ(ierr);	        
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Specify -average_cost, and start time step and number of time steps to average.\n");CHKERRQ(ierr);	        
	ierr = PetscPrintf(PETSC_COMM_WORLD,"No misfit will be computed. \n");CHKERRQ(ierr);	        
        costStartTimeStep=10000000;
        costNumTimeSteps=-1;
        
        }
 /* IK: added for misfit function - end */

  return 0;
}

/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "calcMisfit"
PetscErrorCode calcMisfit(PetscScalar tc, PetscInt iLoop, PetscInt numTracers, Vec *v)

{

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  PetscErrorCode ierr;
  PetscScalar zero = 0.0, one = 1.0, minusone = -1.0 ;  

  if (Iter0+iLoop>=costStartTimeStep) { /* note: costStartTimeStep is ABSOLUTE time step */	

/* Add your code here */

   if (averageCost) {   /* only option, currently */

/* IK : Add all tracer snapshots, and increase counter by one */ 
 
    if (costCount<=costNumTimeSteps) {
	  ierr = VecAXPY(mbgc1avg,one,mbgc1);CHKERRQ(ierr);
	  ierr = VecAXPY(mbgc2avg,one,mbgc2);CHKERRQ(ierr);
	  ierr = VecAXPY(mbgc3avg,one,mbgc3);CHKERRQ(ierr);
	  costCount = costCount+1;
	}
   }
 }
 
 return 0;
 }

/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/


#undef __FUNCT__
#define __FUNCT__ "writeMisfit"
PetscErrorCode writeMisfit(PetscScalar tc, PetscInt iLoop, PetscInt numTracers, Vec *v)
{

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  PetscErrorCode ierr;
  PetscScalar zero = 0.0, one = 1.0, minusone = -1.0 ;  

  if (Iter0+iLoop>=costStartTimeStep) { /* note: costStartTimeStep is ABSOLUTE time step */	

/* Add your code here */

  if (averageCost){   /* only option, currently */

/* IK : added for misfit function : start */

	  if (costCount==costNumTimeSteps) { /* time to compute cost function and write to file */
	  
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Computing cost function of time average over %d steps at time %10.5f, step %d\n", costCount, tc, Iter0+iLoop);CHKERRQ(ierr);                      

/* IK : average model over a year; this will be stored again in mbgc1avg, ... */

		ierr = VecScale(mbgc1avg,1.0/costCount);CHKERRQ(ierr);
		ierr = VecScale(mbgc2avg,1.0/costCount);CHKERRQ(ierr);
		ierr = VecScale(mbgc3avg,1.0/costCount);CHKERRQ(ierr);

		ierr = PetscPrintf(PETSC_COMM_WORLD,"Computed annual average concentrations \n");CHKERRQ(ierr);                      

/* IK : (later) write vectors of cost functions for a posteriori analysis */

		ierr = VecView(mbgc1avg,fdcbgc1avg);CHKERRQ(ierr);
		ierr = VecView(mbgc2avg,fdcbgc2avg);CHKERRQ(ierr);
		ierr = VecView(mbgc3avg,fdcbgc3avg);CHKERRQ(ierr);

		ierr = PetscPrintf(PETSC_COMM_WORLD,"Wrote annual average concentrations \n");CHKERRQ(ierr);                      

/* IK : get global average model tracers, e.g., if I ever want to account for the bias  Gavembgc1, ... */

		ierr = VecDot(mbgc1avg,globVolFrac,&Gavembgc1);CHKERRQ(ierr);
		ierr = VecDot(mbgc2avg,globVolFrac,&Gavembgc2);CHKERRQ(ierr);
		ierr = VecDot(mbgc3avg,globVolFrac,&Gavembgc3);CHKERRQ(ierr);

		ierr = PetscPrintf(PETSC_COMM_WORLD,"Computed global mean model tracers \n");CHKERRQ(ierr);                      

/* copy global mean tracer to array cbgc1, ... and perform subsequent operations on this */

		ierr = VecCopy(mbgc1avg,cbgc1);CHKERRQ(ierr);
		ierr = VecCopy(mbgc2avg,cbgc2);CHKERRQ(ierr);
		ierr = VecCopy(mbgc3avg,cbgc3);CHKERRQ(ierr);

/* IK : subtract observations from average model concentrations; this will be stored in cbgc1avg, ....,  again */

		ierr = VecAXPY(cbgc1,minusone,obgc1);CHKERRQ(ierr);
		ierr = VecAXPY(cbgc2,minusone,obgc2);CHKERRQ(ierr);
		ierr = VecAXPY(cbgc3,minusone,obgc3);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Subtracted model-observations  \n");CHKERRQ(ierr);                      
		
/* IK : square the difference between observations and average model concentrations */
/* Note: I use the more VecPow instead of VecNorm, because I want to be as generic as 
possible, e.g., use functions proposed by Evans, 2003; result is again in cbgcavg1, ..., again */

		ierr = VecPow(cbgc1,2.0);CHKERRQ(ierr);
		ierr = VecPow(cbgc2,2.0);CHKERRQ(ierr);
		ierr = VecPow(cbgc3,2.0);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Squared the difference  \n");CHKERRQ(ierr);                      

/* IK : sum all local misfits, and scale with fractional volume at the same time; result is scalar value Gavecost1, ... */

		ierr = VecDot(cbgc1,globVolFrac,&Gavecost1);CHKERRQ(ierr);
		ierr = VecDot(cbgc2,globVolFrac,&Gavecost2);CHKERRQ(ierr);
		ierr = VecDot(cbgc3,globVolFrac,&Gavecost3);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Computed overall sum, weighted by volume  \n");CHKERRQ(ierr);                      

/* IK : get global average observed concentrations Gavebgc1, ... for normalization 
(needed to sum up different quantities; could use any other scaling, e.g., global variance */

		ierr = VecDot(obgc1,globVolFrac,&Gaveobgc1);CHKERRQ(ierr);
		ierr = VecDot(obgc2,globVolFrac,&Gaveobgc2);CHKERRQ(ierr);
		ierr = VecDot(obgc3,globVolFrac,&Gaveobgc3);CHKERRQ(ierr);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Computed global mean observations \n");CHKERRQ(ierr);                      

/* IK : Get the total cost */    

                RGavecost1=sqrt(Gavecost1);
                RGavecost2=sqrt(Gavecost2);
                RGavecost3=sqrt(Gavecost3);
                Gcost = RGavecost1/Gaveobgc1+RGavecost2/Gaveobgc2+RGavecost3/Gaveobgc3;

		ierr = PetscPrintf(PETSC_COMM_WORLD,"Computed overall cost function \n");CHKERRQ(ierr);                      

/* write misfit to file - perhaps better binary? 
                ierr = writeBinaryScalarData("misfit.bin",&Gcost,1,PETSC_TRUE); */

/* some diagnostic output for checking */
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Check translation\n");CHKERRQ(ierr);                      
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing cost function at time %10.5f, step %d:: %10.5f\n", tc, Iter0+iLoop, Gcost);CHKERRQ(ierr);                      
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Components are PO4 %10.5f  O2 %10.5f NO3 %10.5f\n", RGavecost1,RGavecost2,RGavecost3);CHKERRQ(ierr);                      
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Average observations are PO4 %10.5f O2 %10.5f NO3 %10.5f\n", Gaveobgc1,Gaveobgc2,Gaveobgc3);CHKERRQ(ierr);                      
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Average model tracers are PO4 %10.5f O2 %10.5f NO3 %10.5f\n", Gavembgc1,Gavembgc2,Gavembgc3);CHKERRQ(ierr);                      
                ierr = PetscFPrintf(PETSC_COMM_WORLD,misfitf,"%10.5f\n",Gcost);CHKERRQ(ierr);           

                costCount = 0;

      } /* costCount==costNumTimeSteps */ 
      
   } else { /* no other option that average cost, currently */ 
   
     ierr = PetscPrintf(PETSC_COMM_WORLD,"No cost function available \n");CHKERRQ(ierr);                      
   }

} /* Iter0+iLoop>=costStartTimeStep */

/* IK : added for misfit function : end */

  return 0;
}



/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "finalizeMisfit"
PetscErrorCode finalizeMisfit(PetscScalar tc, PetscInt Iter, PetscInt numTracers)
{

  PetscErrorCode ierr;

/* Add your code here */

/* IK : added for misfit function */

    ierr = VecDestroy(&mbgc1);CHKERRQ(ierr);
    ierr = VecDestroy(&mbgc1avg);CHKERRQ(ierr);
    ierr = VecDestroy(&obgc1);CHKERRQ(ierr);
    ierr = VecDestroy(&cbgc1);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fdcbgc1avg);CHKERRQ(ierr);	

    ierr = VecDestroy(&mbgc2);CHKERRQ(ierr);
    ierr = VecDestroy(&mbgc2avg);CHKERRQ(ierr);
    ierr = VecDestroy(&obgc2);CHKERRQ(ierr);
    ierr = VecDestroy(&cbgc2);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fdcbgc2avg);CHKERRQ(ierr);	

    ierr = VecDestroy(&mbgc3);CHKERRQ(ierr);
    ierr = VecDestroy(&mbgc3avg);CHKERRQ(ierr);
    ierr = VecDestroy(&obgc3);CHKERRQ(ierr);
    ierr = VecDestroy(&cbgc3);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fdcbgc3avg);CHKERRQ(ierr);	

    ierr = PetscFClose(PETSC_COMM_WORLD,misfitf);CHKERRQ(ierr);

/* IK : added for misfit function */

  return 0;
}

