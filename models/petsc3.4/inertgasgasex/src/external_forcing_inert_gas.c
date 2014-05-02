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
#include "tmm_main.h"
#include "inert_gas.h"

/* Macros to map tracer names to vectors */
#define TR v[0]
#define JTR ut[0]

Vec Ts,Ss;
PetscInt *gIndices;
PetscScalar *localTs,*localSs;
PetscScalar *localTR;
PetscScalar *localJTR;
PetscScalar *localfice,*localwind, *localatmosp;
Vec TReq;
PetscScalar *localVgas,*localFgas,*localFinj,*localFex,*localTReq;
PetscScalar DeltaT, dzsurf;
PetscInt gasID = 0;
PetscScalar fluxScaling[3];

PeriodicVec Tsp, Ssp;
PeriodicArray localficep, localwindp, localatmospp;

PetscBool useSeparateWinds = PETSC_FALSE;
PetscScalar *localuwind,*localvwind;
PeriodicArray localuwindp, localvwindp;
PetscInt numWindsPeriods;
PetscScalar *tdpWinds; /* arrays for periodic forcing */
PetscScalar windsCyclePeriod, windsCycleStep;

PetscInt numBiogeochemPeriods;
PetscScalar *tdpBiogeochem; /* arrays for periodic forcing */
PetscBool periodicBiogeochemForcing = PETSC_FALSE;
PetscScalar biogeochemCyclePeriod, biogeochemCycleStep;

PetscBool calcDiagnostics = PETSC_FALSE;
PetscInt diagNumTimeSteps, diagStartTimeStep, diagCount;
PetscBool appendDiagnostics = PETSC_FALSE;
/* Add model specific diagnostic variables below */
PetscScalar *localgasfluxdiagavg, *localinjfluxdiagavg, *localexfluxdiagavg;
Vec TRsatanomdiag, TRsatanomdiagavg, TReqdiagavg;
PetscScalar *localTRsatanomdiag;
PetscViewer fdTRsatanomdiagavg, fdTReqdiagavg;

#if defined (FORSPINUP) || defined (FORJACOBIAN)
PetscScalar relaxTau[50], relaxLambda[50], relaxValue[50];
PetscBool relaxTracer = PETSC_FALSE;
#endif

#undef __FUNCT__
#define __FUNCT__ "iniExternalForcing"
PetscErrorCode iniExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt numTracers, Vec *v, Vec *ut)
{
  PetscErrorCode ierr;
  PetscInt gLow, gHigh, il;
  PetscInt ip, kl, nzloc;
  PetscInt itr;
  PetscViewer fd;
  PetscBool flg;
  PetscInt it;
  PetscScalar myTime;
  PetscScalar zero = 0.0;
  PetscInt maxValsToRead;
    
#if defined (FORSPINUP) || defined (FORJACOBIAN)
  ierr = PetscOptionsHasName(PETSC_NULL,"-relax_tracer",&relaxTracer);CHKERRQ(ierr);
  if (relaxTracer) {  

    maxValsToRead = numTracers;
    ierr = PetscOptionsGetRealArray(PETSC_NULL,"-relax_tau",relaxTau,&maxValsToRead,&flg);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate tracer relaxation tau with the -relax_tau option");
    if (maxValsToRead != numTracers) {
      SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of relaxation tau values specified");
    }

    maxValsToRead = numTracers;
    ierr = PetscOptionsGetRealArray(PETSC_NULL,"-relax_value",relaxValue,&maxValsToRead,&flg);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate relaxation values with the -relax_value option");
    if (maxValsToRead != numTracers) {
      SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of relaxation values specified");
    }
    
    for (itr=0; itr<numTracers; itr++) {
      if (relaxTau[itr]>0.0) {
        relaxLambda[itr]=1.0/relaxTau[itr];
      } else {
        relaxLambda[itr]=0.0;
        relaxValue[itr]=0.0;
      }      
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Tracer %d relaxation lambda=%15.11f, relaxation value=%10.8f\n",itr,relaxLambda[itr],relaxValue[itr]);CHKERRQ(ierr);
    }      
    
  }
#endif
    
/*   v[0]=TR1,v[1]=TR2,... */
  ierr = VecGetArray(TR,&localTR);CHKERRQ(ierr);

  for (itr=0; itr<numTracers; itr++) {
    ierr = VecSet(ut[itr],zero); CHKERRQ(ierr);
  }
  ierr = VecGetArray(JTR,&localJTR);CHKERRQ(ierr);

  ierr = PetscOptionsGetInt(PETSC_NULL,"-inert_gas_id",&gasID,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate inert gas tracer ID with the -inert_gas_id option");  

  ierr = PetscOptionsGetReal(PETSC_NULL,"-biogeochem_deltat",&DeltaT,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate biogeochemical time step in seconds with the -biogeochem_deltat option");  

  ierr = PetscOptionsGetReal(PETSC_NULL,"-dzsurf",&dzsurf,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate surface layer thickness with the -dzsurf option");  

  fluxScaling[0]=1.0;
  fluxScaling[1]=1.0;
  fluxScaling[2]=1.0;
  maxValsToRead = 3;
  ierr = PetscOptionsGetRealArray(PETSC_NULL,"-flux_scaling",fluxScaling,&maxValsToRead,&flg);
  if (flg) {
    if (maxValsToRead != 3) {
      SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of scaling values specified");
    }  
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Fluxes will be scaled as follows:\n");CHKERRQ(ierr);
    ierr=PetscPrintf(PETSC_COMM_WORLD,"  Finj: %15.11f\n",fluxScaling[0]);CHKERRQ(ierr);
    ierr=PetscPrintf(PETSC_COMM_WORLD,"  Fex : %15.11f\n",fluxScaling[1]);CHKERRQ(ierr);
    ierr=PetscPrintf(PETSC_COMM_WORLD,"  Fgas: %15.11f\n",fluxScaling[2]);CHKERRQ(ierr);
  }

  ierr = PetscOptionsHasName(PETSC_NULL,"-periodic_biogeochem_forcing",&periodicBiogeochemForcing);CHKERRQ(ierr);

  if (periodicBiogeochemForcing) {    
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic biogeochemical forcing specified\n");CHKERRQ(ierr);

/*  read time data */
/*  IMPORTANT: time units must be the same as that used by the toplevel driver */
    ierr = PetscOptionsGetReal(PETSC_NULL,"-periodic_biogeochem_cycle_period",&biogeochemCyclePeriod,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate biogeochemical forcing cycling time with the -periodic_biogeochem_cycle_period option");
    ierr = PetscOptionsGetReal(PETSC_NULL,"-periodic_biogeochem_cycle_step",&biogeochemCycleStep,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate biogeochemical forcing cycling step with the -periodic_biogeochem_cycle_step option");
    numBiogeochemPeriods=biogeochemCyclePeriod/biogeochemCycleStep;
/*  array for holding extended time array */
    PetscMalloc((numBiogeochemPeriods+2)*sizeof(PetscScalar), &tdpBiogeochem); 
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Periodic biogeochemical forcing specified at times:\n");CHKERRQ(ierr);            
    for (it=0; it<=numBiogeochemPeriods+1; it++) {
      tdpBiogeochem[it]=(-biogeochemCycleStep/2.0) + it*biogeochemCycleStep;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"tdpBiogeochem=%10.5f\n", tdpBiogeochem[it]);CHKERRQ(ierr);        
    }    
  }
  
/*   Read T and S */
  ierr = VecDuplicate(TR,&Ts);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&Ss);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    Tsp.firstTime = PETSC_TRUE;
    Ssp.firstTime = PETSC_TRUE;
  } else {
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Ts.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(Ts,fd);CHKERRQ(ierr); /* IntoVector */
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Ss.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(Ss,fd);CHKERRQ(ierr); /* IntoVector */
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
  }  
  ierr = VecGetArray(Ts,&localTs);CHKERRQ(ierr);
  ierr = VecGetArray(Ss,&localSs);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Done reading T/S\n");CHKERRQ(ierr);

/*   Compute global indices for local piece of vectors */
  ierr = VecGetOwnershipRange(Ts,&gLow,&gHigh);CHKERRQ(ierr);
  gHigh = gHigh - 1; /* Note: gHigh is one more than the last local element */
  ierr = PetscMalloc(lSize*sizeof(PetscInt),&gIndices);CHKERRQ(ierr);  
  for (il=0; il<lSize; il++) {
    gIndices[il] = il + gLow;
  }  

/* Grid arrays */

/* Forcing fields */  
  ierr = PetscOptionsHasName(PETSC_NULL,"-use_separate_winds",&useSeparateWinds);CHKERRQ(ierr);
  if (useSeparateWinds) {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Separate winds specified: gas transfer velocity will be computed using winds\n");CHKERRQ(ierr);        
	if (!periodicBiogeochemForcing) SETERRQ(PETSC_COMM_WORLD,1,"Separate winds can only be used with periodic biogeochemical forcing!");          

	ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localwind);CHKERRQ(ierr);  
	
	ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localuwind);CHKERRQ(ierr);    
	localuwindp.firstTime = PETSC_TRUE;
	localuwindp.arrayLength = lNumProfiles;    
	ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localvwind);CHKERRQ(ierr);    
	localvwindp.firstTime = PETSC_TRUE;
	localvwindp.arrayLength = lNumProfiles;    
	
	ierr = PetscOptionsGetReal(PETSC_NULL,"-periodic_winds_cycle_period",&windsCyclePeriod,&flg);CHKERRQ(ierr);
	if (flg) {
	  ierr = PetscOptionsGetReal(PETSC_NULL,"-periodic_winds_cycle_step",&windsCycleStep,&flg);CHKERRQ(ierr);
	  numWindsPeriods=windsCyclePeriod/windsCycleStep;
  /*  array for holding extended time array */
	  PetscMalloc((numWindsPeriods+2)*sizeof(PetscScalar), &tdpWinds); 
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Periodic winds forcing specified at times:\n");CHKERRQ(ierr);            
	  for (it=0; it<=numWindsPeriods+1; it++) {
		tdpWinds[it]=(-windsCycleStep/2.0) + it*windsCycleStep;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"tdpWinds=%10.5f\n", tdpWinds[it]);CHKERRQ(ierr);        
	  }    
/*         ierr = PetscPrintf(PETSC_COMM_WORLD,"Winds: %d, %g, %g\n",numWindsPeriods,windsCyclePeriod,windsCycleStep);CHKERRQ(ierr);             */
	  
	} else {
	  windsCyclePeriod=biogeochemCyclePeriod;
	  windsCycleStep=biogeochemCycleStep;
	  numWindsPeriods=windsCyclePeriod/windsCycleStep;
  /*  array for holding extended time array */
	  PetscMalloc((numWindsPeriods+2)*sizeof(PetscScalar), &tdpWinds); 
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Periodic winds forcing specified at times:\n");CHKERRQ(ierr);            
	  for (it=0; it<=numWindsPeriods+1; it++) {
		tdpWinds[it]=(-windsCycleStep/2.0) + it*windsCycleStep;
		ierr = PetscPrintf(PETSC_COMM_WORLD,"tdpWinds=%10.5f\n", tdpWinds[it]);CHKERRQ(ierr);        
	  }    
	}
      
  } else {
	ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localwind);CHKERRQ(ierr);  
	if (periodicBiogeochemForcing) {    
	  localwindp.firstTime = PETSC_TRUE;
	  localwindp.arrayLength = lNumProfiles;
	} else {  
	  ierr = readProfileSurfaceScalarData("wind.bin",localwind,1);  
	}
  }
  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localfice);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localficep.firstTime = PETSC_TRUE;
    localficep.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("fice.bin",localfice,1);  
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localVgas);CHKERRQ(ierr);  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localFgas);CHKERRQ(ierr);  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localFinj);CHKERRQ(ierr);  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localFex);CHKERRQ(ierr);    

  ierr = VecDuplicate(TR,&TReq);CHKERRQ(ierr);
  ierr = VecSet(TReq,zero);CHKERRQ(ierr);
  ierr = VecGetArray(TReq,&localTReq);CHKERRQ(ierr);

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localatmosp);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localatmospp.firstTime = PETSC_TRUE;
    localatmospp.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("atmosp.bin",localatmosp,1);  
  }

  
  if (periodicBiogeochemForcing) {   
	ierr = interpPeriodicVector(tc,&Ts,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Ssp,"Ss_");	
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localficep,"fice_");
    if (useSeparateWinds) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localuwind,windsCyclePeriod,numWindsPeriods,
                                                    tdpWinds,&localuwindp,"uwind_");                                                    
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localvwind,windsCyclePeriod,numWindsPeriods,
                                                    tdpWinds,&localvwindp,"vwind_");                                                        
    } else {                                                  
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemCyclePeriod,numBiogeochemPeriods,
													tdpBiogeochem,&localwindp,"wind_");
	}												
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemCyclePeriod,numBiogeochemPeriods,
									              tdpBiogeochem,&localatmospp,"atmosp_");
  }  

/* Initialize biogeochem model */
  myTime = DeltaT*Iter; /* Iter should start at 0 */
  for (ip=0; ip<lNumProfiles; ip++) {
	kl=lStartIndices[ip];  

    if (useSeparateWinds) localwind[ip] = sqrt(pow(localuwind[ip],2) + pow(localvwind[ip],2));      
   
    inert_gas_fluxes_(&Iter,&myTime,&localTR[kl],&localTs[kl],&localSs[kl],
                          &localwind[ip],&localfice[ip],&localatmosp[ip],
                          &gasID,&localVgas[ip],&localFinj[ip],&localFex[ip],
                          &localTReq[kl]);

    localFinj[ip]=fluxScaling[0]*localFinj[ip]*(1.0-localfice[ip]);
    localFex[ip]=fluxScaling[1]*localFex[ip]*(1.0-localfice[ip]);                              
  }

  ierr = PetscOptionsHasName(PETSC_NULL,"-calc_diagnostics",&calcDiagnostics);CHKERRQ(ierr);
  if (calcDiagnostics) {    
/*Data for diagnostics */
	ierr = PetscOptionsGetInt(PETSC_NULL,"-diag_start_time_step",&diagStartTimeStep,&flg);CHKERRQ(ierr);
	if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate (absolute) time step at which to start storing diagnostics with the -diag_start_time_step flag");
	ierr = PetscOptionsGetInt(PETSC_NULL,"-diag_time_steps",&diagNumTimeSteps,&flg);CHKERRQ(ierr);
	if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate number of time averaging diagnostics time steps with the -diag_time_step flag");
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed starting at (and including) time step: %d\n", diagStartTimeStep);CHKERRQ(ierr);	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed over %d time steps\n", diagNumTimeSteps);CHKERRQ(ierr);	

    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localgasfluxdiagavg);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localinjfluxdiagavg);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localexfluxdiagavg);CHKERRQ(ierr);  
    
	ierr = VecDuplicate(TR,&TRsatanomdiag);CHKERRQ(ierr);
	ierr = VecSet(TRsatanomdiag,zero);CHKERRQ(ierr);
	ierr = VecGetArray(TRsatanomdiag,&localTRsatanomdiag);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&TRsatanomdiagavg);CHKERRQ(ierr);
	ierr = VecSet(TRsatanomdiagavg,zero);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&TReqdiagavg);CHKERRQ(ierr);
	ierr = VecSet(TReqdiagavg,zero);CHKERRQ(ierr);
    
    for (ip=0; ip<lNumProfiles; ip++) {
      localgasfluxdiagavg[ip]=0.0;
      localinjfluxdiagavg[ip]=0.0;
      localexfluxdiagavg[ip]=0.0;
    }

	diagCount=0;
  }

  return 0;
}


#undef __FUNCT__
#define __FUNCT__ "calcExternalForcing"
PetscErrorCode calcExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop, PetscInt numTracers, Vec *v, Vec *ut)
{

  PetscErrorCode ierr;
  PetscInt itr, ip, nzloc, kl;
  PetscScalar myTime;
  PetscScalar Fas;  
  PetscScalar zero = 0.0, one = 1.0;

  myTime = DeltaT*Iter; /* Iter should start at 0 */

  if (periodicBiogeochemForcing) {   
	ierr = interpPeriodicVector(tc,&Ts,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Ssp,"Ss_");	
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localficep,"fice_");
    if (useSeparateWinds) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localuwind,windsCyclePeriod,numWindsPeriods,
                                                    tdpWinds,&localuwindp,"uwind_");                                                    
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localvwind,windsCyclePeriod,numWindsPeriods,
                                                    tdpWinds,&localvwindp,"vwind_");                                                        
    } else {                                                                                                    
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemCyclePeriod,numBiogeochemPeriods,
													tdpBiogeochem,&localwindp,"wind_");
	}
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemCyclePeriod,numBiogeochemPeriods,
									              tdpBiogeochem,&localatmospp,"atmosp_");

/*  Recompute fluxes */
    for (ip=0; ip<lNumProfiles; ip++) {
      nzloc=lProfileLength[ip];    
      kl=lStartIndices[ip];  
  
      if (useSeparateWinds) localwind[ip] = sqrt(pow(localuwind[ip],2) + pow(localvwind[ip],2));      

      inert_gas_fluxes_(&Iter,&myTime,&localTR[kl],&localTs[kl],&localSs[kl],
                            &localwind[ip],&localfice[ip],&localatmosp[ip],
                            &gasID,&localVgas[ip],&localFinj[ip],&localFex[ip],
                            &localTReq[kl]);

      localFinj[ip]=fluxScaling[0]*localFinj[ip]*(1.0-localfice[ip]);
      localFex[ip]=fluxScaling[1]*localFex[ip]*(1.0-localfice[ip]);                            
    }
  }

/* Compute air-sea gas exchange term */
  for (ip=0; ip<lNumProfiles; ip++) {
    nzloc=lProfileLength[ip];  
    kl=lStartIndices[ip];
    
    localFgas[ip]=-fluxScaling[2]*localVgas[ip]*(1.0-localfice[ip])*(localTR[kl] - localTReq[kl]);    
    Fas = localFgas[ip] + localFinj[ip] + localFex[ip];
    
    localJTR[kl] = Fas/dzsurf;

	if (calcDiagnostics) {  
	  if (Iter0+iLoop>=diagStartTimeStep) { /* start time averaging (note: diagStartTimeStep is ABSOLUTE time step) */	
/*      Note: this will recompute localTReq with an atmospheric pressure of 1 because now we  */
/*            are computing localTReq for the full 3-d grid and not just the surface points as  */
/*            done in S/R inert_gas_fluxes. */
        inert_gas_diagnostics_(&nzloc,&Iter,&myTime,&localTR[kl],&localTs[kl],&localSs[kl],&gasID,
        &localTReq[kl],&localTRsatanomdiag[kl]);
      }
	}
                             
  } /* end loop over profiles */

  if (calcDiagnostics) {  
	if (Iter0+iLoop>=diagStartTimeStep) { /* start time averaging (note: diagStartTimeStep is ABSOLUTE time step) */  
      ierr = VecSetValues(TReq,lSize,gIndices,localTReq,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecAssemblyBegin(TReq);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(TReq);CHKERRQ(ierr);    

      ierr = VecSetValues(TRsatanomdiag,lSize,gIndices,localTRsatanomdiag,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecAssemblyBegin(TRsatanomdiag);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(TRsatanomdiag);CHKERRQ(ierr);    
    }
  }

  ierr = VecSetValues(JTR,lSize,gIndices,localJTR,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(JTR);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(JTR);CHKERRQ(ierr);    

#if defined (FORSPINUP) || defined (FORJACOBIAN)
/* add relaxation term: ut = ut - lambda*(v-vr) = ut -lambda*v + lambda*vr */
    if (relaxTracer) {
      for (itr=0; itr<numTracers; itr++) {
        ierr = VecAXPY(ut[itr],-relaxLambda[itr],v[itr]);CHKERRQ(ierr); /* ut = ut - lambda*v */
        ierr = VecShift(ut[itr],relaxLambda[itr]*relaxValue[itr]);CHKERRQ(ierr); /* ut = ut + lambda*vr */
      }
    }
#endif

/* Convert to discrete tendency */
  ierr = VecScale(JTR,DeltaT);CHKERRQ(ierr);  
    
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "writeExternalForcing"
PetscErrorCode writeExternalForcing(PetscScalar tc, PetscInt iLoop, PetscInt numTracers, Vec *v, Vec *ut)
{

  PetscErrorCode ierr;
  PetscInt ip;
  PetscScalar zero = 0.0, one = 1.0;  

/* Note: tc and iLoop are the time and step at the end of the current time step. */

  if (calcDiagnostics) {  
	if (Iter0+iLoop>=diagStartTimeStep) { /* start time averaging (note: diagStartTimeStep is ABSOLUTE time step) */  
  
	  if (diagCount<=diagNumTimeSteps) { /* still within same averaging block so accumulate */

		ierr = VecAXPY(TRsatanomdiagavg,one,TRsatanomdiag);CHKERRQ(ierr);	  
		ierr = VecAXPY(TReqdiagavg,one,TReq);CHKERRQ(ierr);	  

        for (ip=0; ip<lNumProfiles; ip++) {
          localgasfluxdiagavg[ip]=localFgas[ip]+localgasfluxdiagavg[ip];
          localinjfluxdiagavg[ip]=localFinj[ip]+localinjfluxdiagavg[ip];
          localexfluxdiagavg[ip]=localFex[ip]+localexfluxdiagavg[ip];
        }	  

		diagCount = diagCount+1;
	  }

	  if (diagCount==diagNumTimeSteps) { /* time to write averages to file */
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing diagnostics time average at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);                      

		ierr = VecScale(TRsatanomdiagavg,1.0/diagCount);CHKERRQ(ierr);
		ierr = VecScale(TReqdiagavg,1.0/diagCount);CHKERRQ(ierr);

        for (ip=0; ip<lNumProfiles; ip++) {
          localgasfluxdiagavg[ip]=localgasfluxdiagavg[ip]/diagCount;
          localinjfluxdiagavg[ip]=localinjfluxdiagavg[ip]/diagCount;
          localexfluxdiagavg[ip]=localexfluxdiagavg[ip]/diagCount;
        }	          

        if (!appendDiagnostics) {
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"TRsatanomdiagavg.petsc",FILE_MODE_WRITE,&fdTRsatanomdiagavg);CHKERRQ(ierr);
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"TReqdiagavg.petsc",FILE_MODE_WRITE,&fdTReqdiagavg);CHKERRQ(ierr);
        }
		ierr = VecView(TRsatanomdiagavg,fdTRsatanomdiagavg);CHKERRQ(ierr);
		ierr = VecView(TReqdiagavg,fdTReqdiagavg);CHKERRQ(ierr);

		ierr = writeProfileSurfaceScalarData("gasfluxavg.bin",localgasfluxdiagavg,1,appendDiagnostics);  		
		ierr = writeProfileSurfaceScalarData("injfluxavg.bin",localinjfluxdiagavg,1,appendDiagnostics);  		
		ierr = writeProfileSurfaceScalarData("exfluxavg.bin",localexfluxdiagavg,1,appendDiagnostics);  		

        appendDiagnostics=PETSC_TRUE;

/*      reset diagnostic arrays */
		ierr = VecSet(TRsatanomdiagavg,zero); CHKERRQ(ierr);
		ierr = VecSet(TReqdiagavg,zero); CHKERRQ(ierr);
        for (ip=0; ip<lNumProfiles; ip++) {
          localgasfluxdiagavg[ip]=0.0;
          localinjfluxdiagavg[ip]=0.0;
          localexfluxdiagavg[ip]=0.0;
        }	  

		diagCount = 0;        
	  }
	}  
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "finalizeExternalForcing"
PetscErrorCode finalizeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt numTracers)
{

  PetscErrorCode ierr;
  
  ierr = VecDestroy(&Ts);CHKERRQ(ierr);
  ierr = VecDestroy(&Ss);CHKERRQ(ierr);
  ierr = PetscFree(gIndices);CHKERRQ(ierr);  

  if (periodicBiogeochemForcing) {    
    ierr = destroyPeriodicVec(&Tsp);CHKERRQ(ierr);
    ierr = destroyPeriodicVec(&Ssp);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localficep);CHKERRQ(ierr);
    if (useSeparateWinds) {
      ierr = destroyPeriodicArray(&localuwindp);CHKERRQ(ierr);
      ierr = destroyPeriodicArray(&localvwindp);CHKERRQ(ierr);        
    } else {
	  ierr = destroyPeriodicArray(&localwindp);CHKERRQ(ierr);
	}
    ierr = destroyPeriodicArray(&localatmospp);CHKERRQ(ierr);
  }    
  
  if (calcDiagnostics) {    
	ierr = VecDestroy(&TRsatanomdiag);CHKERRQ(ierr);
	ierr = VecDestroy(&TRsatanomdiagavg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fdTRsatanomdiagavg);CHKERRQ(ierr);	  
	ierr = VecDestroy(&TReqdiagavg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fdTReqdiagavg);CHKERRQ(ierr);	  
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "reInitializeExternalForcing"
PetscErrorCode reInitializeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt numTracers, Vec *v, Vec *ut)
{
  PetscErrorCode ierr;
  PetscInt ip, kl;
  PetscScalar myTime;

  if (periodicBiogeochemForcing) {   
	ierr = interpPeriodicVector(tc,&Ts,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Ssp,"Ss_");	
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localficep,"fice_");
    if (useSeparateWinds) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localuwind,windsCyclePeriod,numWindsPeriods,
                                                    tdpWinds,&localuwindp,"uwind_");                                                    
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localvwind,windsCyclePeriod,numWindsPeriods,
                                                    tdpWinds,&localvwindp,"vwind_");                                                        
    } else {                                                  
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemCyclePeriod,numBiogeochemPeriods,
													tdpBiogeochem,&localwindp,"wind_");
	}												                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemCyclePeriod,numBiogeochemPeriods,
									              tdpBiogeochem,&localatmospp,"atmosp_");
  }  

/* Initialize biogeochem model */
  myTime = DeltaT*Iter; /* Iter should start at 0 */
  for (ip=0; ip<lNumProfiles; ip++) {
	kl=lStartIndices[ip];  

    if (useSeparateWinds) localwind[ip] = sqrt(pow(localuwind[ip],2) + pow(localvwind[ip],2));      

    inert_gas_fluxes_(&Iter,&myTime,&localTR[kl],&localTs[kl],&localSs[kl],
                          &localwind[ip],&localfice[ip],&localatmosp[ip],
                          &gasID,&localVgas[ip],&localFinj[ip],&localFex[ip],
                          &localTReq[kl]);

    localFinj[ip]=fluxScaling[0]*localFinj[ip]*(1.0-localfice[ip]);
    localFex[ip]=fluxScaling[1]*localFex[ip]*(1.0-localfice[ip]);                          
  }
    
  return 0;
}
