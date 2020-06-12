#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "petscmat.h"
#include "petsc_matvec_utils.h"
#include "tmm_timer.h"
#include "tmm_forcing_utils.h"
#include "tmm_profile_utils.h"
#include "tmm_profile_data.h"
#include "tmm_misfit.h"
#include "tmm_main.h"
#include "inert_gas_bc_tmm.h"

/* Macros to map tracer names to vectors */
#define TR v[0]
#define BCc bcc[0]
#define BCf bcf[0]

PetscScalar *localTs,*localSs;
PetscScalar *localTR;
PetscScalar *localBCc,*localBCf;
PetscScalar DeltaT;
PetscInt gasID = 0;

Vec atmosp,Tss,Sss;
PeriodicVec atmospp, Tssp, Sssp;
PetscScalar *localatmosp,*localTss,*localSss;

PetscBool periodicBiogeochemForcing = PETSC_FALSE;
PeriodicTimer biogeochemTimer;

PetscBool calcDiagnostics = PETSC_FALSE;
StepTimer diagTimer;
PetscBool appendDiagnostics = PETSC_FALSE;
/* Add model specific diagnostic variables below */
Vec Ts,Ss;
PeriodicVec Tsp, Ssp;
Vec TReqdiag, TRsatanomdiag, TRsatanomdiagavg, TReqdiagavg, BCdiagavg;
PetscScalar *localTReqdiag, *localTRsatanomdiag;
PetscViewer fdTRsatanomdiagavg, fdTReqdiagavg, fdBCdiagavg;

#undef __FUNCT__
#define __FUNCT__ "iniCalcBC"
PetscErrorCode iniCalcBC(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt numTracers, Vec *v, Vec *bcc, Vec *bcf)
{

  PetscErrorCode ierr;
  PetscInt itr;
  PetscViewer fd;
  PetscBool flg;
  PetscInt it;
  PetscScalar myTime;
  PetscScalar zero = 0.0;

/*   v[0]=TR1,v[1]=TR2,... */
  ierr = VecGetArray(TR,&localTR);CHKERRQ(ierr);

  for (itr=0; itr<numTracers; itr++) {
    ierr = VecSet(bcc[itr],zero); CHKERRQ(ierr);
    ierr = VecSet(bcf[itr],zero); CHKERRQ(ierr);    
  }
  ierr = VecGetArray(BCc,&localBCc);CHKERRQ(ierr);
  ierr = VecGetArray(BCf,&localBCf);CHKERRQ(ierr);

  ierr = PetscOptionsGetInt(NULL,NULL,"-inert_gas_id",&gasID,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate inert gas tracer ID with the -inert_gas_id option");  

  ierr = PetscOptionsGetReal(NULL,NULL,"-biogeochem_deltat",&DeltaT,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate biogeochemical time step in seconds with the -biogeochem_deltat option");  

  ierr = PetscOptionsHasName(NULL,NULL,"-periodic_biogeochem_forcing",&periodicBiogeochemForcing);CHKERRQ(ierr);

  if (periodicBiogeochemForcing) {    
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic biogeochemical forcing specified\n");CHKERRQ(ierr);
    ierr = iniPeriodicTimer("periodic_biogeochem_", &biogeochemTimer);CHKERRQ(ierr);
  }
  
/* Grid arrays */

/* Forcing fields */  
  ierr = VecDuplicate(BCc,&Tss);CHKERRQ(ierr);
  ierr = VecDuplicate(BCc,&Sss);CHKERRQ(ierr);  
  ierr = VecDuplicate(BCc,&atmosp);CHKERRQ(ierr);

  if (periodicBiogeochemForcing) {    
    Tssp.firstTime = PETSC_TRUE;
    Sssp.firstTime = PETSC_TRUE;
    atmospp.firstTime = PETSC_TRUE;
  } else {
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Tss.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr = VecLoad(Tss,fd);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Sss.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr = VecLoad(Sss,fd);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"atmosp.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr = VecLoad(atmosp,fd);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
  }  
  
  ierr = VecGetArray(Tss,&localTss);CHKERRQ(ierr);
  ierr = VecGetArray(Sss,&localSss);CHKERRQ(ierr);
  ierr = VecGetArray(atmosp,&localatmosp);CHKERRQ(ierr);
/*   ierr = PetscPrintf(PETSC_COMM_WORLD,"Done reading surface T, S, and atmospheric pressure\n");CHKERRQ(ierr); */

  ierr = PetscOptionsHasName(NULL,NULL,"-calc_diagnostics",&calcDiagnostics);CHKERRQ(ierr);
  if (calcDiagnostics) {    
/*   Read T and S */
    ierr = VecDuplicate(TR,&Ts);CHKERRQ(ierr);
    ierr = VecDuplicate(TR,&Ss);CHKERRQ(ierr);  
    if (periodicBiogeochemForcing) {    
      Tsp.firstTime = PETSC_TRUE;
      Ssp.firstTime = PETSC_TRUE;

      ierr = interpPeriodicVector(tc,&Ts,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Tsp,"Ts_");
      ierr = interpPeriodicVector(tc,&Ss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Ssp,"Ss_");	
      
    } else {
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Ts.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = VecLoad(Ts,fd);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Ss.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = VecLoad(Ss,fd);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
    }  
    
    ierr = VecGetArray(Ts,&localTs);CHKERRQ(ierr);
    ierr = VecGetArray(Ss,&localSs);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Done reading T/S\n");CHKERRQ(ierr);

    ierr = VecDuplicate(TR,&TReqdiag);CHKERRQ(ierr);
    ierr = VecSet(TReqdiag,zero);CHKERRQ(ierr);
    ierr = VecGetArray(TReqdiag,&localTReqdiag);CHKERRQ(ierr);
  
/*Data for diagnostics */
    ierr = iniStepTimer("diag_", Iter0, &diagTimer);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed starting at (and including) time step: %d\n", diagTimer.startTimeStep);CHKERRQ(ierr);	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed over %d time steps\n", diagTimer.numTimeSteps);CHKERRQ(ierr);	
    
	ierr = VecDuplicate(TR,&TRsatanomdiag);CHKERRQ(ierr);
	ierr = VecSet(TRsatanomdiag,zero);CHKERRQ(ierr);
	ierr = VecGetArray(TRsatanomdiag,&localTRsatanomdiag);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&TRsatanomdiagavg);CHKERRQ(ierr);
	ierr = VecSet(TRsatanomdiagavg,zero);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&TReqdiagavg);CHKERRQ(ierr);
	ierr = VecSet(TReqdiagavg,zero);CHKERRQ(ierr);
	ierr = VecDuplicate(BCc,&BCdiagavg);CHKERRQ(ierr);
	ierr = VecSet(BCdiagavg,zero);CHKERRQ(ierr);
    
  }

/* Initialize biogeochem model */  
  if (periodicBiogeochemForcing) {   
    ierr = interpPeriodicVector(tc,&Tss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Tssp,"Tss_");
    ierr = interpPeriodicVector(tc,&Sss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Sssp,"Sss_");	
    ierr = interpPeriodicVector(tc,&atmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&atmospp,"atmosp_");	        
  }  

  myTime = DeltaT*Iterc; /* Iter should start at 0 */
  inert_gas_bc_(&lBCSize,&Iterc,&myTime,&localTss[0],&localSss[0],
                        &localatmosp[0],&gasID,
                        &localBCc[0]);
    
  if (periodicBiogeochemForcing) {   
    ierr = interpPeriodicVector(tf,&Tss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Tssp,"Tss_");
    ierr = interpPeriodicVector(tf,&Sss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Sssp,"Sss_");	
    ierr = interpPeriodicVector(tf,&atmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&atmospp,"atmosp_");	        
  }  

  myTime = DeltaT*Iterf; /* Iter should start at 0 */
  inert_gas_bc_(&lBCSize,&Iterf,&myTime,&localTss[0],&localSss[0],
                        &localatmosp[0],&gasID,
                        &localBCf[0]);      

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "calcBC"
PetscErrorCode calcBC(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt iLoop, PetscInt numTracers, Vec *v, Vec *bcc, Vec *bcf)
{

  PetscErrorCode ierr;
  PetscScalar myTime;
  PetscScalar zero = 0.0, one = 1.0;

/*  Recompute BC */
  if (periodicBiogeochemForcing) {   
    ierr = interpPeriodicVector(tc,&Tss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Tssp,"Tss_");
    ierr = interpPeriodicVector(tc,&Sss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Sssp,"Sss_");	
    ierr = interpPeriodicVector(tc,&atmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&atmospp,"atmosp_");	        

    myTime = DeltaT*Iterc; /* Iter should start at 0 */
    inert_gas_bc_(&lBCSize,&Iterc,&myTime,&localTss[0],&localSss[0],
                          &localatmosp[0],&gasID,
                          &localBCc[0]);
    
    ierr = interpPeriodicVector(tf,&Tss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Tssp,"Tss_");
    ierr = interpPeriodicVector(tf,&Sss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Sssp,"Sss_");	
    ierr = interpPeriodicVector(tf,&atmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&atmospp,"atmosp_");	        

    myTime = DeltaT*Iterf; /* Iter should start at 0 */
    inert_gas_bc_(&lBCSize,&Iterf,&myTime,&localTss[0],&localSss[0],
                          &localatmosp[0],&gasID,
                          &localBCf[0]);
    
    ierr = VecSetValues(BCc,lBCSize,gBCIndices,localBCc,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(BCc);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(BCc);CHKERRQ(ierr);    
  
    ierr = VecSetValues(BCf,lBCSize,gBCIndices,localBCf,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(BCf);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(BCf);CHKERRQ(ierr);        
  }

  myTime = DeltaT*Iterc; /* Iter should start at 0 */
  if (calcDiagnostics) {  
	if (Iter0+iLoop>=diagTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */	
/*      Note: this will recompute localTReqdiag with an atmospheric pressure of 1 because now we  */
/*            are computing localTReqdiag for the full 3-d grid and not just the surface points as  */
/*            done in S/R inert_gas_fluxes. */

      if (periodicBiogeochemForcing) {    
        ierr = interpPeriodicVector(tc,&Ts,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Tsp,"Ts_");
        ierr = interpPeriodicVector(tc,&Ss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Ssp,"Ss_");	
      }

      inert_gas_diagnostics_(&lSize,&Iterc,&myTime,&localTR[0],&localTs[0],&localSs[0],&gasID,
      &localTReqdiag[0],&localTRsatanomdiag[0]);
      
      ierr = VecSetValues(TReqdiag,lSize,gIndices,localTReqdiag,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecAssemblyBegin(TReqdiag);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(TReqdiag);CHKERRQ(ierr);    

      ierr = VecSetValues(TRsatanomdiag,lSize,gIndices,localTRsatanomdiag,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecAssemblyBegin(TRsatanomdiag);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(TRsatanomdiag);CHKERRQ(ierr);    
	}  
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "writeBC"
PetscErrorCode writeBC(PetscScalar tc, PetscInt iLoop, PetscInt numTracers, Vec *v, Vec *bcc, Vec *bcf)
{

  PetscErrorCode ierr;
  PetscInt ip;
  PetscScalar zero = 0.0, one = 1.0;  

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  if (calcDiagnostics) {  
	if (Iter0+iLoop>=diagTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */  
  
	  if (diagTimer.count<=diagTimer.numTimeSteps) { /* still within same averaging block so accumulate */

		ierr = VecAXPY(TRsatanomdiagavg,one,TRsatanomdiag);CHKERRQ(ierr);	  
		ierr = VecAXPY(TReqdiagavg,one,TReqdiag);CHKERRQ(ierr);	  
		ierr = VecAXPY(BCdiagavg,one,BCc);CHKERRQ(ierr);	  

		diagTimer.count++;
	  }

	  if (diagTimer.count==diagTimer.numTimeSteps) { /* time to write averages to file */
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing diagnostics time average at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);                      

		ierr = VecScale(TRsatanomdiagavg,1.0/diagTimer.count);CHKERRQ(ierr);
		ierr = VecScale(TReqdiagavg,1.0/diagTimer.count);CHKERRQ(ierr);
		ierr = VecScale(BCdiagavg,1.0/diagTimer.count);CHKERRQ(ierr);

        if (!appendDiagnostics) {
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"TRsatanomdiagavg.petsc",FILE_MODE_WRITE,&fdTRsatanomdiagavg);CHKERRQ(ierr);
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"TReqdiagavg.petsc",FILE_MODE_WRITE,&fdTReqdiagavg);CHKERRQ(ierr);
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"BCdiagavg.petsc",FILE_MODE_WRITE,&fdBCdiagavg);CHKERRQ(ierr);        
        }
        
		ierr = VecView(TRsatanomdiagavg,fdTRsatanomdiagavg);CHKERRQ(ierr);
		ierr = VecView(TReqdiagavg,fdTReqdiagavg);CHKERRQ(ierr);
		ierr = VecView(BCdiagavg,fdBCdiagavg);CHKERRQ(ierr);

        appendDiagnostics=PETSC_TRUE;

/*      reset diagnostic arrays */
		ierr = VecSet(TRsatanomdiagavg,zero); CHKERRQ(ierr);
		ierr = VecSet(TReqdiagavg,zero); CHKERRQ(ierr);
		ierr = VecSet(BCdiagavg,zero); CHKERRQ(ierr);

        ierr = updateStepTimer("diag_", Iter0+iLoop, &diagTimer);CHKERRQ(ierr);

	  }
	}  
  }


  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "finalizeCalcBC"
PetscErrorCode finalizeCalcBC(PetscScalar tc, PetscInt Iter, PetscInt numTracers)
{

  PetscErrorCode ierr;


  ierr = VecDestroy(&Tss);CHKERRQ(ierr);  
  ierr = VecDestroy(&Sss);CHKERRQ(ierr);  
  ierr = VecDestroy(&atmosp);CHKERRQ(ierr);  

  if (periodicBiogeochemForcing) {    
    ierr = destroyPeriodicVec(&Tssp);CHKERRQ(ierr);
    ierr = destroyPeriodicVec(&Sssp);CHKERRQ(ierr);
    ierr = destroyPeriodicVec(&atmospp);CHKERRQ(ierr);
  }    
  
  if (calcDiagnostics) {    
    ierr = VecDestroy(&Ts);CHKERRQ(ierr);
    ierr = VecDestroy(&Ss);CHKERRQ(ierr);
    if (periodicBiogeochemForcing) {    
      ierr = destroyPeriodicVec(&Tsp);CHKERRQ(ierr);
      ierr = destroyPeriodicVec(&Ssp);CHKERRQ(ierr);
    }    
	ierr = VecDestroy(&TRsatanomdiag);CHKERRQ(ierr);
	ierr = VecDestroy(&TRsatanomdiagavg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fdTRsatanomdiagavg);CHKERRQ(ierr);	  
	ierr = VecDestroy(&TReqdiagavg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fdTReqdiagavg);CHKERRQ(ierr);	  
	ierr = VecDestroy(&BCdiagavg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fdBCdiagavg);CHKERRQ(ierr);	  
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "reInitializeCalcBC"
PetscErrorCode reInitializeCalcBC(PetscScalar tc, PetscInt Iterc, PetscScalar tf, PetscInt Iterf, PetscInt numTracers, Vec *v, Vec *bcc, Vec *bcf)
{

  PetscErrorCode ierr;
  PetscScalar myTime;

/*  Recompute BC */
  if (periodicBiogeochemForcing) {   
    ierr = interpPeriodicVector(tc,&Tss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Tssp,"Tss_");
    ierr = interpPeriodicVector(tc,&Sss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Sssp,"Sss_");	
    ierr = interpPeriodicVector(tc,&atmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&atmospp,"atmosp_");	        

    myTime = DeltaT*Iterc; /* Iter should start at 0 */
    inert_gas_bc_(&lBCSize,&Iterc,&myTime,&localTss[0],&localSss[0],
                          &localatmosp[0],&gasID,
                          &localBCc[0]);
    
    ierr = interpPeriodicVector(tf,&Tss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Tssp,"Tss_");
    ierr = interpPeriodicVector(tf,&Sss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Sssp,"Sss_");	
    ierr = interpPeriodicVector(tf,&atmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&atmospp,"atmosp_");	        

    myTime = DeltaT*Iterf; /* Iter should start at 0 */
    inert_gas_bc_(&lBCSize,&Iterf,&myTime,&localTss[0],&localSss[0],
                          &localatmosp[0],&gasID,
                          &localBCf[0]);
    
    ierr = VecSetValues(BCc,lBCSize,gBCIndices,localBCc,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(BCc);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(BCc);CHKERRQ(ierr);    
  
    ierr = VecSetValues(BCf,lBCSize,gBCIndices,localBCf,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(BCf);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(BCf);CHKERRQ(ierr);        
  }

  return 0;
}

