/* $Header: /Users/ikriest/CVS/mops/tmm_misfit.c,v 1.1 2015/11/17 14:18:51 ikriest Exp $ */
/* $Name: mops-2_0 $*/

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "tmm_petsc_matvec_utils.h"
#include "tmm_timer.h"
#include "tmm_forcing_utils.h"
#include "tmm_profile_utils.h"
#include "tmm.h"
#include "tmm_misfit.h"
#include "tmm_share.h"
#include "mops_biogeochem_misfit_data.h"

#define TR state->c[0]

static PetscClassId MISFIT_CLASSID;

typedef struct _p_MisfitCtx *MisfitContext;
struct _p_MisfitCtx {
  PETSCHEADER(int);
  PetscInt efctxId;
  PetscInt stateId;
/* Add problem-specific variables below */  
};

// VS: added auxiliary functions for histogram and Hellinger misfit metric:
int calcHist( Vec x, Vec ix, PetscInt bins, PetscScalar delta, PetscInt* hist )
{
/* 
 * calcHist: Histogram w.r.t. a tracer vector,
 * x: tracer vector
 * ix: indicator for considered x components
 * bins: number of bins to be used for the histogram
 * delta: width of each bin of the histogram
 * hist: the histogram 
 */
  PetscErrorCode  ierr;
  Vec             y; 
  PetscInt        i, j, n, nLocal, *histLocal;
  PetscMPIInt     myRank;
  PetscScalar     surplus, lower, upper, *yLocal;

  ierr = MPI_Comm_rank( PETSC_COMM_WORLD, &myRank ); CHKERRQ( ierr );
  // generate vector y from x by setting unconsidered components to zero
  ierr = VecDuplicate( x, &y ); CHKERRQ( ierr );
  ierr = VecPointwiseMult( y, x, ix ); CHKERRQ( ierr );
  ierr = VecGetSize( y, &n ); CHKERRQ( ierr );
  ierr = VecGetLocalSize( y, &nLocal ); CHKERRQ( ierr );
  yLocal = ( PetscScalar* )calloc( nLocal, sizeof( PetscScalar ) );
  histLocal = ( PetscInt* )calloc( bins, sizeof( PetscInt ) );
  // calculate surplus, the number of unconsidered tracer components, i.e.,
  // the vector length n minus the number of 1-components of ix 
  ierr = VecSum( ix, &surplus );
  surplus = ( PetscScalar )n - surplus;
  for ( i = 0; i < bins; i++ ) histLocal[ i ] = 0;

  // generate histogram of local (process) part of y
  ierr = VecGetArray( y, &yLocal ); CHKERRQ( ierr );
  ierr = PetscSortReal( nLocal, yLocal ); CHKERRQ( ierr );
  lower = 0; // start of histogram = lower boundary of class 1
  upper = lower + delta; // upper boundary of class 1
  j = 0; // initialise class 1
  i = 0;
  while ( i < nLocal )
  {
    if ( yLocal[ i ] < upper || j == bins - 1 )
    { // as long as values in vector are lower than the upper boundary, increase counter
      histLocal[ j ] = histLocal[ j ] + 1;
      i++;
    }
    else if ( j < bins - 1 )
    { // if not, shift lower and upper boundary by delta, and increase class number
      lower = upper;
      upper = upper + delta;
      j = j + 1;
    }
  }
  ierr = VecRestoreArray( y, &yLocal ); CHKERRQ( ierr );
  // now add all local histograms to be the global histogram
  i = MPI_Reduce( histLocal, hist, bins, MPI_INT, MPI_SUM, 0, PETSC_COMM_WORLD );
  hist[ 0 ] = hist[ 0 ] - ( int )surplus; // remove surplus of items in the first bin
  for ( i = 0; i < bins; i++ ) hist[ i ]  = hist[ i ];

#ifndef TRANSIENTCOST
  if ( myRank == 0 ) { ierr = PetscPrintf( PETSC_COMM_WORLD, "surplus = %i \n", ( int )surplus ); CHKERRQ( ierr ); }
#endif

  VecDestroy( &y );
  free( yLocal );
  free( histLocal );
  return 0;
}

double calcHD( Vec obs, Vec mod, Vec ix, PetscInt bins, PetscScalar delta, PetscScalar nobs )
{
/* 
 * calcHD: Hellinger distance (HD) w.r.t. two (equal sized) tracer vectors,
 *         e.g., containing observations and model results, respectively
 * obs: first tracer vector, e.g., obgc1
 * mod: second tracer vector, e.g., mbgc1avg
 * ix: indicator of considered components of obs and mod
 * bins: number of bins to be used for histograms the HD is based on
 * delta: width of each bin of the histograms
 * returns HD between obs and mod, acc. to histograms
 */
  PetscErrorCode  ierr;
  PetscInt        i, j, n, *obsHist, *modHist;
  PetscMPIInt     myRank;
  PetscScalar     result;

  ierr = MPI_Comm_rank( PETSC_COMM_WORLD, &myRank ); CHKERRQ( ierr );
  ierr = VecGetSize( obs, &n ); CHKERRQ( ierr );
  obsHist = ( PetscInt* )calloc( bins, sizeof( PetscInt ) );
  modHist = ( PetscInt* )calloc( bins, sizeof( PetscInt ) );

  // generate both histograms
  calcHist( obs, ix, bins, delta, obsHist );
  calcHist( mod, ix, bins, delta, modHist );

#ifndef TRANSIENTCOST
  if ( myRank == 0 )
  { for ( j = 0; j < bins; j++ ) PetscPrintf( PETSC_COMM_WORLD, "obsHist[ %i ] = %12.6f\n", j, obsHist[ j ]/nobs );
    for ( j = 0; j < bins; j++ ) PetscPrintf( PETSC_COMM_WORLD, "modHist[ %i ] = %12.6f\n", j, modHist[ j ]/nobs );
  }
#endif

  // calculate Hellinger distance from both global histograms
  result = 0;
  for ( j = 0; j < bins; j++ ) result = result + pow( pow( modHist[ j ], 0.5 ) - pow( obsHist[ j ], 0.5 ), 2.0 );
  result = 0.5 * result / nobs;
  result = pow( result, 0.5 );

  free( obsHist );
  free( modHist );
  return result;
}

#undef __FUNCT__
#define __FUNCT__ "iniMisfit"
PetscErrorCode iniMisfit(PetscScalar tc, PetscInt Iter, TMMState state, void *userctx)
{

  PetscInt numTracers;
  const char *prefix;
  void *ctx;

  PetscErrorCode ierr;
  PetscBool flg;
  PetscScalar zero = 0.0, one = 1.0;
  PetscViewer fd;

  static PetscBool registered = PETSC_FALSE;
  static PetscInt efctxId = 0;

/* Add your code here */
  MisfitContext ef;

  MPI_Comm comm = PETSC_COMM_WORLD;
  
  if (!registered) {
    PetscClassIdRegister("Misfit context", &MISFIT_CLASSID);
    registered = PETSC_TRUE;
  }
  PetscHeaderCreate(ef, MISFIT_CLASSID, "Misfit", "Misfit context", "Misfit", comm, 0, 0); 

  efctxId++;
  ef->efctxId=efctxId;
  ef->stateId=state->stateId;

  PetscContainer ctxcontainer;
  PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
  PetscCall(PetscContainerSetPointer(ctxcontainer, (void*)ef));
  PetscCall(PetscObjectCompose((PetscObject)state, "misfit ctx", (PetscObject)ctxcontainer));
  state->misfitctxcontainer = ctxcontainer;
  PetscCall(PetscContainerDestroy(&ctxcontainer));

  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

// Now set problem specific data
// Common data (only initialize/read once)
  if (ef->efctxId==1) {
// etc etc
  }
  
  ierr = StepTimerCreate(&misfitTimer);CHKERRQ(ierr);
  ierr = StepTimerIni("misfit_", prefix, Iter0+1, misfitTimer);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Misfit will be computed starting at (and including) time step: %d\n", misfitTimer->startTimeStep);CHKERRQ(ierr);	
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Misfit will be computed over %d time steps\n", misfitTimer->numTimeSteps);CHKERRQ(ierr);
   
  multiObjective = PETSC_FALSE;
   
/* IK: added for misfit function - start */
 

  ierr = PetscOptionsGetString(NULL,prefix,"-misfit_file",misfitFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) {
  strcpy(misfitFile,"");
  sprintf(misfitFile,"%s","misfit.txt");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"misfit will be written to %s\n",misfitFile);CHKERRQ(ierr);
  ierr = PetscFOpen(PETSC_COMM_WORLD,misfitFile,"w",&misfitf);CHKERRQ(ierr);  

  ierr = PetscOptionsHasName(NULL,prefix,"-multi_objective",&multiObjective);CHKERRQ(ierr);
  if(multiObjective) {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing second objective for misfit function\n");CHKERRQ(ierr);	
   }

/* Read fractional volume of each box for scaling according to box size */

	ierr = VecDuplicate(TR,&globVolFrac);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"volume_fraction.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(globVolFrac,fd);CHKERRQ(ierr);  /* IntoVector */ 
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Read volume file for weighting\n");CHKERRQ(ierr);	

/* Read vector of weights for individual tracer types and gridboxes */
	ierr = VecDuplicate(TR,&wbgc1);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"weights.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(wbgc1,fd);CHKERRQ(ierr);  /* IntoVector */ 
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
	
        ierr = VecDuplicate(TR,&wbgc2);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"weights.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(wbgc2,fd);CHKERRQ(ierr);  /* IntoVector */ 
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
	
	ierr = VecDuplicate(TR,&wbgc3);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"weights.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(wbgc3,fd);CHKERRQ(ierr);  /* IntoVector */ 
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Read misfit file(s) for weighting\n");CHKERRQ(ierr);	

/* construct total weight: for now, multiply weight with fractional volume */

      ierr = VecDuplicate(TR,&w1);CHKERRQ(ierr);
      ierr = VecDuplicate(TR,&w2);CHKERRQ(ierr);
      ierr = VecDuplicate(TR,&w3);CHKERRQ(ierr);
      ierr = VecSet(w1,zero);CHKERRQ(ierr);
      ierr = VecSet(w2,zero);CHKERRQ(ierr);
      ierr = VecSet(w3,zero);CHKERRQ(ierr);
      ierr = VecPointwiseMult(w1,wbgc1,globVolFrac);
      ierr = VecPointwiseMult(w2,wbgc2,globVolFrac);
      ierr = VecPointwiseMult(w3,wbgc3,globVolFrac);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Set up weight arrays\n");CHKERRQ(ierr);	

/* IK : Check number of observations */

      ierr = VecSum(wbgc1,&NoSum1);CHKERRQ(ierr);
      ierr = VecSum(wbgc2,&NoSum2);CHKERRQ(ierr);
      ierr = VecSum(wbgc3,&NoSum3);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of observations: PO4: %12.0f O2: %12.0f NO3: %12.0f \n",NoSum1, NoSum2, NoSum3);CHKERRQ(ierr);	

/* IK : Compute global fractional volume of all locations where observations exist for weighting */

      ierr = VecSum(w1,&VolSum1);CHKERRQ(ierr);
      ierr = VecSum(w2,&VolSum2);CHKERRQ(ierr);
      ierr = VecSum(w3,&VolSum3);CHKERRQ(ierr);

/* Scale to get volume fraction of total observed volume */
             
      ierr = VecScale(w1,1.0/VolSum1);CHKERRQ(ierr);
      ierr = VecScale(w2,1.0/VolSum2);CHKERRQ(ierr);
      ierr = VecScale(w3,1.0/VolSum3);CHKERRQ(ierr);

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
	
  ierr = VecDuplicate(TR,&mbgc1avg);CHKERRQ(ierr);
  ierr = VecSet(mbgc1avg,zero);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&mbgc2avg);CHKERRQ(ierr);
  ierr = VecSet(mbgc2avg,zero);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&mbgc3avg);CHKERRQ(ierr);
  ierr = VecSet(mbgc3avg,zero);CHKERRQ(ierr);

  ierr = VecDuplicate(TR,&cbgc1);CHKERRQ(ierr);
  ierr = VecSet(cbgc1,zero);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&cbgc2);CHKERRQ(ierr);
  ierr = VecSet(cbgc2,zero);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&cbgc3);CHKERRQ(ierr);
  ierr = VecSet(cbgc3,zero);CHKERRQ(ierr);

/* IK Start OMZ misfit */
/* IK Start OMZ misfit */

  ierr = VecDuplicate(TR,&onevec);CHKERRQ(ierr);
  ierr = VecSet(onevec,one);CHKERRQ(ierr);

  ierr = VecDuplicate(TR,&momz);CHKERRQ(ierr);
  ierr = VecSet(momz,zero);CHKERRQ(ierr);
        
  ierr = VecDuplicate(TR,&oomz);CHKERRQ(ierr);
  ierr = VecSet(oomz,zero);CHKERRQ(ierr);

  ierr = VecDuplicate(TR,&amomz);CHKERRQ(ierr);
  ierr = VecSet(amomz,zero);CHKERRQ(ierr);

  ierr = VecDuplicate(TR,&aoomz);CHKERRQ(ierr);
  ierr = VecSet(aoomz,zero);CHKERRQ(ierr);

  ierr = VecDuplicate(TR,&imomz);CHKERRQ(ierr);
  ierr = VecSet(imomz,zero);CHKERRQ(ierr);

  ierr = VecDuplicate(TR,&ioomz);CHKERRQ(ierr);
  ierr = VecSet(ioomz,zero);CHKERRQ(ierr);

  ierr = VecDuplicate(TR,&imoomz);CHKERRQ(ierr);
  ierr = VecSet(imoomz,zero);CHKERRQ(ierr);

/* IK End OMZ misfit */
/* IK End OMZ misfit */

/* Start organic matter misfit */
/* Start organic matter misfit */

/* Read vector of weights for individual tracer types and gridboxes */
	ierr = VecDuplicate(TR,&wbgc4);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"phy_weight.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(wbgc4,fd);CHKERRQ(ierr);  /* IntoVector */ 
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
	
	ierr = VecDuplicate(TR,&wbgc5);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"zoo_weight.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(wbgc5,fd);CHKERRQ(ierr);  /* IntoVector */ 
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
	
	ierr = VecDuplicate(TR,&wbgc6);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"pop_weight.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(wbgc6,fd);CHKERRQ(ierr);  /* IntoVector */ 
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
	
	ierr = VecDuplicate(TR,&wbgc7);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"dop_weight.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(wbgc7,fd);CHKERRQ(ierr);  /* IntoVector */ 
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Read organic matter misfit file(s) for weighting\n");CHKERRQ(ierr);	

/* construct total weight: for now, multiply weight with fractional volume */

      ierr = VecDuplicate(TR,&w4);CHKERRQ(ierr);
      ierr = VecDuplicate(TR,&w5);CHKERRQ(ierr);
      ierr = VecDuplicate(TR,&w6);CHKERRQ(ierr);
      ierr = VecDuplicate(TR,&w7);CHKERRQ(ierr);
      ierr = VecSet(w4,zero);CHKERRQ(ierr);
      ierr = VecSet(w5,zero);CHKERRQ(ierr);
      ierr = VecSet(w6,zero);CHKERRQ(ierr);
      ierr = VecSet(w7,zero);CHKERRQ(ierr);
      ierr = VecPointwiseMult(w4,wbgc4,globVolFrac); /* IK So far, weigh these with volume fraction */
      ierr = VecPointwiseMult(w5,wbgc5,globVolFrac); /* IK So far, weigh these with volume fraction */
      ierr = VecPointwiseMult(w6,wbgc6,globVolFrac); /* IK So far, weigh these with volume fraction */
      ierr = VecPointwiseMult(w7,wbgc7,globVolFrac); /* IK So far, weigh these with volume fraction */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Set up weight arrays for organic matter\n");CHKERRQ(ierr);	

/* IK : Check number of observations */

      ierr = VecSum(wbgc4,&NoSum4);CHKERRQ(ierr);
      ierr = VecSum(wbgc5,&NoSum5);CHKERRQ(ierr);
      ierr = VecSum(wbgc6,&NoSum6);CHKERRQ(ierr);
      ierr = VecSum(wbgc7,&NoSum7);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of observations: PHY: %12.0f ZOO: %10.0f POP: %12.0f DOP: %12.0f \n",NoSum4, NoSum5, NoSum6, NoSum7);CHKERRQ(ierr);

/* IK : Compute global fractional volume of all locations where observations exist for weighting */

      ierr = VecSum(w4,&VolSum4);CHKERRQ(ierr);
      ierr = VecSum(w5,&VolSum5);CHKERRQ(ierr);
      ierr = VecSum(w6,&VolSum6);CHKERRQ(ierr);
      ierr = VecSum(w7,&VolSum7);CHKERRQ(ierr);

/* Scale to get volume fraction of total observed volume */
             
      ierr = VecScale(w4,1.0/VolSum4);CHKERRQ(ierr);
      ierr = VecScale(w5,1.0/VolSum5);CHKERRQ(ierr);
      ierr = VecScale(w6,1.0/VolSum6);CHKERRQ(ierr);
      ierr = VecScale(w7,1.0/VolSum7);CHKERRQ(ierr);

/* Read the observations */ 

  ierr = VecDuplicate(TR,&obgc4);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&obgc5);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&obgc6);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&obgc7);CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"phy_obs.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = VecLoad(obgc4,fd);CHKERRQ(ierr);  /* IntoVector */ 
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Read phytoplankton observations\n");CHKERRQ(ierr);	

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"zoo_obs.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = VecLoad(obgc5,fd);CHKERRQ(ierr);  /* IntoVector */ 
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Read zooplankton observations\n");CHKERRQ(ierr);	

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"pop_obs.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = VecLoad(obgc6,fd);CHKERRQ(ierr);  /* IntoVector */ 
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Read POP observations\n");CHKERRQ(ierr);	

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"dop_obs.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = VecLoad(obgc7,fd);CHKERRQ(ierr);  /* IntoVector */ 
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Read DOP observations\n");CHKERRQ(ierr);	

/* Initialize the cost function vectors */
	
  ierr = VecDuplicate(TR,&mbgc4avg);CHKERRQ(ierr);
  ierr = VecSet(mbgc4avg,zero);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&mbgc5avg);CHKERRQ(ierr);
  ierr = VecSet(mbgc5avg,zero);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&mbgc6avg);CHKERRQ(ierr);
  ierr = VecSet(mbgc6avg,zero);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&mbgc7avg);CHKERRQ(ierr);
  ierr = VecSet(mbgc7avg,zero);CHKERRQ(ierr);

  ierr = VecDuplicate(TR,&cbgc4);CHKERRQ(ierr);
  ierr = VecSet(cbgc4,zero);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&cbgc5);CHKERRQ(ierr);
  ierr = VecSet(cbgc5,zero);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&cbgc6);CHKERRQ(ierr);
  ierr = VecSet(cbgc6,zero);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&cbgc7);CHKERRQ(ierr);
  ierr = VecSet(cbgc7,zero);CHKERRQ(ierr);

/* End organic matter misfit */
/* End organic matter misfit */
 
 /* IK: added for misfit function - end */

  return 0;
}

/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "calcMisfit"
PetscErrorCode calcMisfit(PetscScalar tc, PetscInt iLoop, TMMState state, void *userctx)
{

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  PetscInt numTracers;
  const char *prefix;
  void *ctx;

  PetscErrorCode ierr;
  PetscScalar zero = 0.0, one = 1.0, minusone = -1.0 ;  
  PetscScalar zoodetfrac = 0.5;
  PetscInt ipo4=0,idop=1,ioxy=2,iphy=3,izoo=4,idet=5,ino3=6 ;  

  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

  if (Iter0+iLoop>=misfitTimer->startTimeStep) { /* start computing misfit (note: startTimeStep is ABSOLUTE time step) */	

/* IK : Add all tracer snapshots, and increase counter by one */ 
 
//     if (misfitTimer->count<misfitTimer->numTimeSteps) {
	  ierr = VecAXPY(mbgc1avg,one,state->c[ipo4]);CHKERRQ(ierr); /* PO4 */
	  ierr = VecAXPY(mbgc2avg,one,state->c[ioxy]);CHKERRQ(ierr); /* O2 */
	  ierr = VecAXPY(mbgc3avg,one,state->c[ino3]);CHKERRQ(ierr); /* NO3 */
	  ierr = VecAXPY(mbgc4avg,one,state->c[iphy]);CHKERRQ(ierr); /* phytoplankton */
	  ierr = VecAXPY(mbgc5avg,one,state->c[izoo]);CHKERRQ(ierr); /* zooplankton */
	  ierr = VecAXPY(mbgc6avg,one,state->c[iphy]);CHKERRQ(ierr); /* POP: phytoplankton component */
	  ierr = VecAXPY(mbgc6avg,one,state->c[idet]);CHKERRQ(ierr); /* POP: detritus component */
	  ierr = VecAXPY(mbgc6avg,zoodetfrac,state->c[izoo]);CHKERRQ(ierr); /* POP: zooplankton component */
	  ierr = VecAXPY(mbgc7avg,one,state->c[idop]);CHKERRQ(ierr); /* DOP: phytoplankton component */
	  misfitTimer->count++;
// 	}
 }
 
 return 0;
 }

/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/


#undef __FUNCT__
#define __FUNCT__ "writeMisfit"
PetscErrorCode writeMisfit(PetscScalar tc, PetscInt iLoop, TMMState state, void *userctx)
{

/* Note: tc and iLoop are the time and step at the end of the current time step. */
  PetscInt numTracers;
  const char *prefix;
  void *ctx;

  PetscErrorCode ierr;
  PetscScalar zero = 0.0, one = 1.0, minusone = -1.0 ;  
  PetscScalar omzcrit = 50.0; /* IK Edit this if necessary */

  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

/* IK: for HD */

  PetscInt nbins = 50; /* VS number of bins per histogram */
  PetscScalar delta1 = 0.08, delta2 = 8.0, delta3 = 1.0, delta4 = 0.002, delta5 = 0.002, delta6 = 0.01, delta7 = 0.02; /* VS bin widths for histograms of NO3, O2, PO4, Phy, Zoo, Det, DOP */

  if (Iter0+iLoop>=misfitTimer->startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */	

/* IK : added for misfit function : start */

	  if (misfitTimer->count==misfitTimer->numTimeSteps) { /* time to compute cost function and write to file */

/* IK : average model over a year; this will be stored again in mbgc1avg, ... */

		ierr = VecScale(mbgc1avg,1.0/misfitTimer->count);CHKERRQ(ierr);
		ierr = VecScale(mbgc2avg,1.0/misfitTimer->count);CHKERRQ(ierr);
		ierr = VecScale(mbgc3avg,1.0/misfitTimer->count);CHKERRQ(ierr);
		ierr = VecScale(mbgc4avg,1.0/misfitTimer->count);CHKERRQ(ierr);
		ierr = VecScale(mbgc5avg,1.0/misfitTimer->count);CHKERRQ(ierr);
		ierr = VecScale(mbgc6avg,1.0/misfitTimer->count);CHKERRQ(ierr);
		ierr = VecScale(mbgc7avg,1.0/misfitTimer->count);CHKERRQ(ierr);

/* copy global mean tracer to array cbgc1, ... and perform subsequent operations on this */

		ierr = VecCopy(mbgc1avg,cbgc1);CHKERRQ(ierr);
		ierr = VecCopy(mbgc2avg,cbgc2);CHKERRQ(ierr);
		ierr = VecCopy(mbgc3avg,cbgc3);CHKERRQ(ierr);
		ierr = VecCopy(mbgc4avg,cbgc4);CHKERRQ(ierr);
		ierr = VecCopy(mbgc5avg,cbgc5);CHKERRQ(ierr);
		ierr = VecCopy(mbgc6avg,cbgc6);CHKERRQ(ierr);
		ierr = VecCopy(mbgc7avg,cbgc7);CHKERRQ(ierr);

/* IK Begin standard RMSE misfit */
/* IK Begin standard RMSE misfit */

/* IK : Compute global fractional volume of all locations where observations exist for calculation of global average concentrations */

/* IK : get global average model tracers, e.g., if I ever want to account for the bias  Gavembgc1, ... */

		ierr = VecDot(mbgc1avg,w1,&Gavembgc1);CHKERRQ(ierr);
		ierr = VecDot(mbgc2avg,w2,&Gavembgc2);CHKERRQ(ierr);
		ierr = VecDot(mbgc3avg,w3,&Gavembgc3);CHKERRQ(ierr);

/* IK : get global average observed concentrations Gavebgc1, ... for normalization 
(needed to sum up different quantities; could use any other scaling, e.g., global variance */

		ierr = VecDot(obgc1,w1,&Gaveobgc1);CHKERRQ(ierr);
		ierr = VecDot(obgc2,w2,&Gaveobgc2);CHKERRQ(ierr);
		ierr = VecDot(obgc3,w3,&Gaveobgc3);CHKERRQ(ierr);

/* IK : subtract observations from average model concentrations; this will be stored in cbgc1avg, ....,  again */

		ierr = VecAXPY(cbgc1,minusone,obgc1);CHKERRQ(ierr);
		ierr = VecAXPY(cbgc2,minusone,obgc2);CHKERRQ(ierr);
		ierr = VecAXPY(cbgc3,minusone,obgc3);CHKERRQ(ierr);
		
/* IK : square the difference between observations and average model concentrations */
/* Note: I use the more VecPow instead of VecNorm, because I want to be as generic as 
possible, e.g., use functions proposed by Evans, 2003; result is again in cbgcavg1, ..., again */

		ierr = VecPow(cbgc1,2.0);CHKERRQ(ierr);
		ierr = VecPow(cbgc2,2.0);CHKERRQ(ierr);
		ierr = VecPow(cbgc3,2.0);CHKERRQ(ierr);

/* IK : sum all local misfits, and scale with fractional volume at the same time; result is scalar value Gavecost1, ... */

		ierr = VecDot(cbgc1,w1,&Gavecost1);CHKERRQ(ierr);
		ierr = VecDot(cbgc2,w2,&Gavecost2);CHKERRQ(ierr);
		ierr = VecDot(cbgc3,w3,&Gavecost3);CHKERRQ(ierr);             

/* IK : Get the total cost */    

                RGavecost1=sqrt(Gavecost1);
                RGavecost2=sqrt(Gavecost2);
                RGavecost3=sqrt(Gavecost3);

            Ginorgcost = RGavecost1/Gaveobgc1+RGavecost2/Gaveobgc2+RGavecost3/Gaveobgc3;

/* IK End standard RMSE misfit */
/* IK End standard RMSE misfit */

/* IK Start OMZ misfit */
/* IK Start OMZ misfit */

/* IK Copy model and observed concentrations to new vector for OMZs */ 
                ierr = VecCopy(mbgc2avg,momz);CHKERRQ(ierr);
                ierr = VecCopy(obgc2,oomz);CHKERRQ(ierr);

/* IK Multiply concentrations by -1:  [x]= -1 * [x] */ 

                ierr = VecScale(momz,minusone);CHKERRQ(ierr);
                ierr = VecScale(oomz,minusone);CHKERRQ(ierr);

/* IK Add criterion for for OMZ. This will be negative outside OMZ.  [x]= crit - [x] */ 

                ierr = VecShift(momz,omzcrit);CHKERRQ(ierr);
                ierr = VecShift(oomz,omzcrit);CHKERRQ(ierr);

/* IK Create a new vector that contains the absolute values of the difference. */ 

                ierr = VecCopy(momz,amomz);CHKERRQ(ierr);
                ierr = VecCopy(oomz,aoomz);CHKERRQ(ierr);

                ierr = VecAbs(amomz);CHKERRQ(ierr);
                ierr = VecAbs(aoomz);CHKERRQ(ierr);

/* IK Add the two vectors. This will contain twice the difference for waters within OMZ, and 0 outside */ 

                ierr = VecAXPY(momz,one,amomz);CHKERRQ(ierr);
                ierr = VecAXPY(oomz,one,aoomz);CHKERRQ(ierr);

/* IK Set all regions within OMZ to 1. Outside will be 0. */ 

                ierr = VecPointwiseMin(imomz,momz,onevec);CHKERRQ(ierr);
                ierr = VecPointwiseMin(ioomz,oomz,onevec);CHKERRQ(ierr);

/* IK Get the region of overlap and sum. This will be 1 in regions of overlap, and 0 outside. Sum the volume. */

                ierr = VecPointwiseMult(imoomz,imomz,ioomz);;CHKERRQ(ierr);
                ierr = VecDot(imoomz,globVolFrac,&GMOomz);CHKERRQ(ierr);

/* IK Get simulated and observed total OMZ volume */

                ierr = VecDot(imomz,globVolFrac,&GMomz);CHKERRQ(ierr);
                ierr = VecDot(ioomz,globVolFrac,&GOomz);CHKERRQ(ierr);

/* IK Compute the global metric */

                OMZMetric = one - GMOomz/(GMomz+GOomz-GMOomz);

                 
/* IK End OMZ misfit */
/* IK End OMZ misfit */

/* IK Start organic matter misfit */
/* IK Start organic matter misfit */

/* IK : Calculate global average model tracers */

		ierr = VecDot(mbgc4avg,w4,&Gavembgc4);CHKERRQ(ierr);
		ierr = VecDot(mbgc5avg,w5,&Gavembgc5);CHKERRQ(ierr);
		ierr = VecDot(mbgc6avg,w6,&Gavembgc6);CHKERRQ(ierr);
		ierr = VecDot(mbgc7avg,w7,&Gavembgc7);CHKERRQ(ierr);

/* IK : Calculate global average concentrations */

		ierr = VecDot(obgc4,w4,&Gaveobgc4);CHKERRQ(ierr); 
		ierr = VecDot(obgc5,w5,&Gaveobgc5);CHKERRQ(ierr); 
		ierr = VecDot(obgc6,w6,&Gaveobgc6);CHKERRQ(ierr); 
		ierr = VecDot(obgc7,w7,&Gaveobgc7);CHKERRQ(ierr); 

/* IK : subtract observations from average model concentrations; this will be stored in cbgc4avg, ....,  again */

		ierr = VecAXPY(cbgc4,minusone,obgc4);CHKERRQ(ierr);
		ierr = VecAXPY(cbgc5,minusone,obgc5);CHKERRQ(ierr);
		ierr = VecAXPY(cbgc6,minusone,obgc6);CHKERRQ(ierr);
		ierr = VecAXPY(cbgc7,minusone,obgc7);CHKERRQ(ierr);
		
/* IK : square the difference between observations and average model concentrations */

		ierr = VecPow(cbgc4,2.0);CHKERRQ(ierr);
		ierr = VecPow(cbgc5,2.0);CHKERRQ(ierr);
		ierr = VecPow(cbgc6,2.0);CHKERRQ(ierr);
		ierr = VecPow(cbgc7,2.0);CHKERRQ(ierr);

/* IK : sum all local misfits, and scale with weight volume at the same time; result is scalar value Gavecost4, ... */

		ierr = VecDot(cbgc4,w4,&Gavecost4);CHKERRQ(ierr);
		ierr = VecDot(cbgc5,w5,&Gavecost5);CHKERRQ(ierr);
		ierr = VecDot(cbgc6,w6,&Gavecost6);CHKERRQ(ierr);
		ierr = VecDot(cbgc7,w7,&Gavecost7);CHKERRQ(ierr);

/* IK : Get the total cost: normalise by specific total volume and take square root */    

            RGavecost4=sqrt(Gavecost4);
            RGavecost5=sqrt(Gavecost5);
            RGavecost6=sqrt(Gavecost6);
            RGavecost7=sqrt(Gavecost7);
            Gorgcost = RGavecost4/Gaveobgc4+RGavecost5/Gaveobgc5+RGavecost6/Gaveobgc6;

/* IK End organic matter misfit */
/* IK End organic matter misfit */

/*IK Total RMSE Cost */ 
            RMSEMetric = Ginorgcost + Gorgcost;

/* VS Start Hellinger misfit */
/* VS Start Hellinger misfit */

            hd1 = calcHD( obgc1, mbgc1avg, wbgc1, nbins, delta1, NoSum1 );                
            hd2 = calcHD( obgc2, mbgc2avg, wbgc2, nbins, delta2, NoSum2 );                
            hd3 = calcHD( obgc3, mbgc3avg, wbgc3, nbins, delta3, NoSum3 );  
            hd4 = calcHD( obgc4, mbgc4avg, wbgc4, nbins, delta4, NoSum4 );  
            hd5 = calcHD( obgc5, mbgc5avg, wbgc5, nbins, delta5, NoSum5 );  
            hd6 = calcHD( obgc6, mbgc6avg, wbgc6, nbins, delta6, NoSum6 );  
            hd7 = calcHD( obgc7, mbgc7avg, wbgc7, nbins, delta7, NoSum7 );  
            HDMetric = hd1 + hd2 + hd3 + hd4 + hd5 + hd6;              

/* VS End Hellinger misfit */
/* VS End Hellinger misfit */

/* IK Write single objective and, if desired second objective for multiobjective */

    ierr = PetscFPrintf(PETSC_COMM_WORLD,misfitf,"%10.5f\n",RMSEMetric);CHKERRQ(ierr);           
    if (multiObjective) {
                ierr = PetscFPrintf(PETSC_COMM_WORLD,misfitf,"%10.5f\n",HDMetric);CHKERRQ(ierr);                     
    }
    
/* some diagnostic output for checking - only do this if not in transient cost function mode */
#ifndef TRANSIENTCOST
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing  cost function at time %10.5f, step %d:: %10.5f %10.5f %10.5f\n", tc, Iter0+iLoop,RMSEMetric,Ginorgcost,Gorgcost);CHKERRQ(ierr);    
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Components are PO4 %10.5f  O2 %10.5f NO3 %10.5f OMZ %10.5f\n", RGavecost1,RGavecost2,RGavecost3,GMOomz);CHKERRQ(ierr);                      
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Average observations are PO4 %10.5f O2 %10.5f NO3 %10.5f OMZ %10.5f\n", Gaveobgc1,Gaveobgc2,Gaveobgc3,GOomz);CHKERRQ(ierr);                      
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Average model tracers are PO4 %10.5f O2 %10.5f NO3 %10.5f OMZ %10.5f\n", Gavembgc1,Gavembgc2,Gavembgc3,GMomz);CHKERRQ(ierr);                      
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Organic components are PHY %10.5f  ZOO %10.5f DET %10.5f DOP %10.5f\n", RGavecost4,RGavecost5,RGavecost6,RGavecost7);CHKERRQ(ierr);                      
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Average organic observations are PHY %10.5f ZOO %10.5f DET %10.5f DOP %10.5f\n", Gaveobgc4,Gaveobgc5,Gaveobgc6,Gaveobgc7);CHKERRQ(ierr);                      
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Average organic model tracers are PHY %10.5f ZOO %10.5f DET %10.5f DOP %10.5f\n", Gavembgc4,Gavembgc5,Gavembgc6,Gavembgc7);CHKERRQ(ierr);                      
		ierr = PetscPrintf(PETSC_COMM_WORLD,"HD components model tracers are PO4 %10.5f  O2 %10.5f NO3 %10.5f PHY %10.5f ZOO %10.5f DET %10.5f DOP %10.5f Total %10.5f\n", hd1,hd2,hd3,hd4,hd5,hd6,hd7,HDMetric);CHKERRQ(ierr);                      
#else
                ierr = PetscPrintf(PETSC_COMM_WORLD,"%10.5f", tc);CHKERRQ(ierr);             
		ierr = PetscPrintf(PETSC_COMM_WORLD,"%10.5f %10.5f %10.5f", GMomz,GOomz,GMOomz);CHKERRQ(ierr); 
		ierr = PetscPrintf(PETSC_COMM_WORLD,"%10.5f %10.5f %10.5f", Ginorgcost,Gorgcost,RMSEMetric);CHKERRQ(ierr); 
		ierr = PetscPrintf(PETSC_COMM_WORLD,"%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f",RGavecost1,RGavecost2,RGavecost3,RGavecost4,RGavecost5,RGavecost6,RGavecost7);
		ierr = PetscPrintf(PETSC_COMM_WORLD,"%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f",Gavembgc1,Gavembgc2,Gavembgc3,Gavembgc4,Gavembgc5,Gavembgc6,Gavembgc7);CHKERRQ(ierr); 
		ierr = PetscPrintf(PETSC_COMM_WORLD,"%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f",Gaveobgc1,Gaveobgc2,Gaveobgc3,Gaveobgc4,Gaveobgc5,Gaveobgc6,Gaveobgc7);CHKERRQ(ierr); 
		ierr = PetscPrintf(PETSC_COMM_WORLD,"%10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n", hd1,hd2,hd3,hd4,hd5,hd6,hd7,HDMetric);CHKERRQ(ierr);   
#endif

/* IK: set simulated tracer array to zero again */

		ierr = VecSet(mbgc1avg,zero);CHKERRQ(ierr);
		ierr = VecSet(mbgc2avg,zero);CHKERRQ(ierr);
		ierr = VecSet(mbgc3avg,zero);CHKERRQ(ierr);
                ierr = VecSet(mbgc4avg,zero);CHKERRQ(ierr);
                ierr = VecSet(mbgc5avg,zero);CHKERRQ(ierr);
                ierr = VecSet(mbgc6avg,zero);CHKERRQ(ierr);
                ierr = VecSet(mbgc7avg,zero);CHKERRQ(ierr);


        ierr = StepTimerUpdate(Iter0+iLoop, misfitTimer);CHKERRQ(ierr);
    
      } /* misfitTimer.Count==misfitTimer.numTimeSteps */ 
   
   }

/* IK : added for misfit function : end */

  return 0;
}



/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/
/* -----------------------------------------------------------------------------------------------------------*/

#undef __FUNCT__
#define __FUNCT__ "finalizeMisfit"
PetscErrorCode finalizeMisfit(PetscScalar tc, PetscInt Iter, TMMState state, void *userctx)
{
  PetscInt numTracers;
  const char *prefix;
  void *ctx;

  PetscErrorCode ierr;

/* Add your code here */
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

/* IK : added for misfit function */

    ierr = VecDestroy(&mbgc1avg);CHKERRQ(ierr);
    ierr = VecDestroy(&obgc1);CHKERRQ(ierr);
    ierr = VecDestroy(&cbgc1);CHKERRQ(ierr);

    ierr = VecDestroy(&mbgc2avg);CHKERRQ(ierr);
    ierr = VecDestroy(&obgc2);CHKERRQ(ierr);
    ierr = VecDestroy(&cbgc2);CHKERRQ(ierr);

    ierr = VecDestroy(&mbgc3avg);CHKERRQ(ierr);
    ierr = VecDestroy(&obgc3);CHKERRQ(ierr);
    ierr = VecDestroy(&cbgc3);CHKERRQ(ierr);

/* IK Start OMZ misfit */
/* IK Start OMZ misfit */

    ierr = VecDestroy(&onevec);CHKERRQ(ierr);
    ierr = VecDestroy(&momz);CHKERRQ(ierr);
    ierr = VecDestroy(&oomz);CHKERRQ(ierr);
    ierr = VecDestroy(&amomz);CHKERRQ(ierr);
    ierr = VecDestroy(&aoomz);CHKERRQ(ierr);
    ierr = VecDestroy(&imomz);CHKERRQ(ierr);
    ierr = VecDestroy(&ioomz);CHKERRQ(ierr);
    ierr = VecDestroy(&imoomz);CHKERRQ(ierr);

/* IK End OMZ misfit */
/* IK End OMZ misfit */

/* IK Start organic matter misfit */
/* IK Start organic matter misfit */

    ierr = VecDestroy(&mbgc4avg);CHKERRQ(ierr);
    ierr = VecDestroy(&obgc4);CHKERRQ(ierr);
    ierr = VecDestroy(&cbgc4);CHKERRQ(ierr);

    ierr = VecDestroy(&mbgc5avg);CHKERRQ(ierr);
    ierr = VecDestroy(&obgc5);CHKERRQ(ierr);
    ierr = VecDestroy(&cbgc5);CHKERRQ(ierr);

    ierr = VecDestroy(&mbgc6avg);CHKERRQ(ierr);
    ierr = VecDestroy(&obgc6);CHKERRQ(ierr);
    ierr = VecDestroy(&cbgc6);CHKERRQ(ierr);

    ierr = VecDestroy(&mbgc7avg);CHKERRQ(ierr);
    ierr = VecDestroy(&obgc7);CHKERRQ(ierr);
    ierr = VecDestroy(&cbgc7);CHKERRQ(ierr);

/* IK End organic matter misfit */
/* IK End organic matter misfit */

    ierr = PetscFClose(PETSC_COMM_WORLD,misfitf);CHKERRQ(ierr);

/* IK : added for misfit function */

  return 0;
}

