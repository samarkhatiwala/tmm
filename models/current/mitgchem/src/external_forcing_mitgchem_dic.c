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
#include "mitgchem_dic_tmm.h"

/* Macros to map tracer names to vectors */
/* v[0]=DIC,v[1]=Alk,v[2]=PO4,v[3]=DOP,v[4]=O2,v[5]=Fe */
#define TR v[0]
#define DIC v[0]
#define Alk v[1]
#define JDIC ut[0]
#define JAlk ut[1]
#define localDIC localTR[0]
#define localAlk localTR[1]
#define localJDIC localJTR[0]
#define localJAlk localJTR[1]

Vec Ts,Ss;
PetscScalar *localTs,*localSs;
PetscScalar **localTR, **localJTR;
PetscScalar *pH;
PetscScalar *localwind,*localatmosp,*localsilica,*localfice,*localEmP;
#ifdef READ_PAR
PetscScalar *localpar;
#else
PetscScalar *locallat;
#endif
#ifdef ALLOW_FE
PetscScalar *localinputfe;
#endif
#ifdef LIGHT_CHL
PetscScalar *localchl;
#endif
PetscScalar *localhFacC,*localrecip_hFacC;
PetscInt nzmax;
PetscScalar DeltaT;
PetscScalar *drF;

PetscBool useSeparateBiogeochemTimeStepping = PETSC_FALSE;

PetscInt toModel = 1; 
PetscInt fromModel = 2;

char fileName[PETSC_MAX_PATH_LEN];  
PetscScalar *localalpha, *localrain_ratio;

PeriodicVec Tsp, Ssp;
PeriodicArray localwindp, localatmospp, localsilicap,localficep,localEmPp;
#ifdef READ_PAR
PeriodicArray localparp;
#endif
#ifdef ALLOW_FE
PeriodicArray localinputfep;
#endif
#ifdef LIGHT_CHL
PeriodicArray localchlp;
#endif
PetscBool periodicBiogeochemForcing = PETSC_FALSE;
PeriodicTimer biogeochemTimer;

PetscInt maxValsToRead;
PetscInt dummyInt = 0;
char *pCO2atmFiles[2];  
PetscInt numpCO2atm_hist = 0;
PetscScalar *TpCO2atm_hist, *pCO2atm_hist;
PetscBool fixedAtmosCO2 = PETSC_TRUE;
char pCO2atmIniFile[PETSC_MAX_PATH_LEN];  

/* Land/Atm model variables */
PetscBool useAtmModel = PETSC_FALSE;
PetscBool useLandModel = PETSC_FALSE;
PetscBool useEmissions = PETSC_FALSE;
PetscBool interpEmissions = PETSC_FALSE;
char *emFiles[3];  
PetscInt numEmission_hist = 0;
PetscScalar *Tem_hist, *E_hist, *D_hist;
PetscScalar fossilFuelEmission = 0.0, landUseEmission = 0.0;
PetscScalar cumulativeEmissions = 0.0;
PetscScalar pCO2atm_ini = 280.0; /* default initial value */
PetscScalar pCO2atm = 280.0; /* default initial value */
PetscScalar *localdA;
PetscScalar landState[3], landSource[3];
PetscScalar deltaTsg = 0.0;
PetscScalar ppmToPgC=2.1324;
PetscScalar atmModelDeltaT;
PetscScalar daysPerYear, secondsPerYear;
PetscScalar Fland = 0.0, Focean=0.0;

StepTimer atmWriteTimer;
PetscBool atmAppendOutput;
FILE *atmfptime;
PetscViewer atmfd;
PetscInt atmfp;
char atmOutTimeFile[PETSC_MAX_PATH_LEN];  
PetscScalar pCO2atmavg, Foceanavg, Flandavg, landStateavg[3];

PetscBool calcDiagnostics = PETSC_FALSE;
StepTimer diagTimer;
PetscBool appendDiagnostics = PETSC_FALSE;
Vec bioac, bioacavg;
PetscViewer fdbioacavg;
PetscScalar *localbioac;
PetscScalar *localpco2diag, *localpco2diagavg;
PetscScalar *localgasexfluxdiag, *localgasexfluxdiagavg, *localempfluxdiag, *localempfluxdiagavg;

#if defined (FORSPINUP) || defined (FORJACOBIAN)
PetscScalar relaxTau[50], relaxLambda[50], relaxValue[50];
PetscBool relaxTracer = PETSC_FALSE;
#endif

#undef __FUNCT__
#define __FUNCT__ "iniExternalForcing"
PetscErrorCode iniExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt numTracers, Vec *v, Vec *ut)
{
  PetscErrorCode ierr;
  PetscInt ip, kl, nzloc;
  PetscInt itr;  
  PetscViewer fd;
  PetscInt fp;
  PetscBool flg;
  PetscInt it;
  PetscScalar myTime;
  PetscScalar zero = 0.0;

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

  for (itr=0; itr<numTracers; itr++) {
    ierr = VecSet(ut[itr],zero); CHKERRQ(ierr);
  }
  
/*   v[0]=DIC,v[1]=Alk,v[2]=PO4,v[3]=DOP,v[4]=O2,v[5]=Fe */
  ierr = VecGetArrays(v,numTracers,&localTR);CHKERRQ(ierr);
  ierr = VecGetArrays(ut,numTracers,&localJTR);CHKERRQ(ierr);

#ifdef ALLOW_OLD_VIRTUALFLUX
  SETERRQ(PETSC_COMM_WORLD,1,"Virtual flux option not supported!");  
#endif

  ierr = PetscOptionsGetReal(PETSC_NULL,"-biogeochem_deltat",&DeltaT,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate biogeochemical time step in seconds with the -biogeochem_deltat option");  

  ierr = PetscOptionsGetReal(PETSC_NULL,"-days_per_year",&daysPerYear,&flg);CHKERRQ(ierr);
  if (!flg) {
    daysPerYear = 360.0;
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of days per year is %12.7f\n",daysPerYear);CHKERRQ(ierr);
  secondsPerYear = 86400.0*daysPerYear;

/* Note: The mitgchem package always internally time-steps tracers and returns the updated tracer field at the next time step. */
/* To use this internal time-stepping, switch on the useSeparateBiogeochemTimeStepping flag with -separate_biogeochem_time_stepping. */
/* Alternatively, we can compute a tendency from the new and old tracer field (done in mitgchem_copy_data) and use that within the */
/* TMM driver code to time-step tracers. This is the default. */
  ierr = PetscOptionsGetBool(PETSC_NULL,"-separate_biogeochem_time_stepping",&useSeparateBiogeochemTimeStepping,0);CHKERRQ(ierr);
#if defined (FORSPINUP) || defined (FORJACOBIAN)
  if (useSeparateBiogeochemTimeStepping) {
    SETERRQ(PETSC_COMM_WORLD,1,"Cannot use the -separate_biogeochem_time_stepping option with SPINUP or JACOBIAN ");    
  }
#endif
  if (useSeparateBiogeochemTimeStepping) {
    fromModel = 3;  
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Biogeochem model will be time-stepped independently\n");CHKERRQ(ierr);
  }  

  ierr = PetscOptionsHasName(PETSC_NULL,"-periodic_biogeochem_forcing",&periodicBiogeochemForcing);CHKERRQ(ierr);

  if (periodicBiogeochemForcing) {    
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic biogeochemical forcing specified\n");CHKERRQ(ierr);
    ierr = iniPeriodicTimer("periodic_biogeochem_", &biogeochemTimer);CHKERRQ(ierr);
  }

/*   Read T and S */
  ierr = VecDuplicate(TR,&Ts);CHKERRQ(ierr);
  ierr = VecDuplicate(TR,&Ss);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    Tsp.firstTime = PETSC_TRUE;
    Ssp.firstTime = PETSC_TRUE;
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

/* Land/Atm model data */
  ierr = PetscOptionsHasName(PETSC_NULL,"-use_atm_model",&useAtmModel);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(PETSC_NULL,"-use_land_model",&useLandModel);CHKERRQ(ierr);

  if ((useLandModel) && (!useAtmModel)) SETERRQ(PETSC_COMM_WORLD,1,"ERROR: Land model cannot be used without the atmospheric model");

  if (useAtmModel) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using interactive atmospheric model\n");CHKERRQ(ierr);  

    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localdA);CHKERRQ(ierr);
    ierr = readProfileSurfaceScalarData("dA.bin",localdA,1);  

/* overwrite default value */
	ierr = PetscOptionsGetReal(PETSC_NULL,"-pco2atm_ini",&pCO2atm_ini,&flg);CHKERRQ(ierr); /* read from command line */
    if (!flg) {
      ierr = PetscOptionsGetString(PETSC_NULL,"-pco2atm_ini_file",pCO2atmIniFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
      if (flg) { /* read from binary file */
        ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,pCO2atmIniFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
        ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
        ierr = PetscBinaryRead(fp,&pCO2atm_ini,1,PETSC_SCALAR);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      }
    }
    pCO2atm = pCO2atm_ini;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using initial atmospheric pCO2 of %g ppm\n",pCO2atm);CHKERRQ(ierr);
      
    atmModelDeltaT = DeltaT/secondsPerYear; /* time step in years */

	ierr = iniStepTimer("atm_write_", Iter0, &atmWriteTimer);CHKERRQ(ierr);

    ierr = PetscOptionsHasName(PETSC_NULL,"-atm_append",&atmAppendOutput);CHKERRQ(ierr);
    if (atmAppendOutput) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output will be appended\n");CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output will overwrite existing file(s)\n");CHKERRQ(ierr);
    }    

/* Output times */
    ierr = PetscOptionsGetString(PETSC_NULL,"-atm_time_file",atmOutTimeFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
    if (!flg) {
      strcpy(atmOutTimeFile,"");
      sprintf(atmOutTimeFile,"%s","atm_output_time.txt");
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output times will be written to %s\n",atmOutTimeFile);CHKERRQ(ierr);

    if (!atmAppendOutput) {
      ierr = PetscFOpen(PETSC_COMM_WORLD,atmOutTimeFile,"w",&atmfptime);CHKERRQ(ierr);  
      ierr = PetscFPrintf(PETSC_COMM_WORLD,atmfptime,"%d   %10.5f\n",Iter0,time0);CHKERRQ(ierr);     
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing atmospheric output at time %10.5f, step %d\n", tc,Iter);CHKERRQ(ierr);  
      ierr = writeBinaryScalarData("pCO2atm_output.bin",&pCO2atm,1,PETSC_FALSE);
    } else {
      ierr = PetscFOpen(PETSC_COMM_WORLD,atmOutTimeFile,"a",&atmfptime);CHKERRQ(ierr);  
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output will be appended. Initial condition will NOT be written\n");CHKERRQ(ierr);      
    }

    if (useLandModel) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Using interactive land model\n");CHKERRQ(ierr);      
      ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"land_ini.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
      ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
      ierr = PetscBinaryRead(fp,landState,3,PETSC_SCALAR);CHKERRQ(ierr);
      ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    

      if (!atmAppendOutput) {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing land output at time %10.5f, step %d\n", tc,Iter);CHKERRQ(ierr);  
        ierr = writeBinaryScalarData("land_state_output.bin",landState,3,PETSC_FALSE);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Land model output will be appended. Initial condition will NOT be written\n");CHKERRQ(ierr);      
      }

    }

    /* CO2 emissions */
	maxValsToRead = 3;
	emFiles[0] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* time file */
	emFiles[1] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* fossil fuel emissions file */
	emFiles[2] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* land use emissions file */
    ierr = PetscOptionsGetStringArray(PETSC_NULL,"-emissions_history",emFiles,&maxValsToRead,&useEmissions);CHKERRQ(ierr);
    if (useEmissions) { /* Read emissions history */
      if (maxValsToRead != 3) {
        SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of file names specified for emissions");
      }      
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Using prescribed emissions\n");CHKERRQ(ierr);     
      /* read time data */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,emFiles[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	  ierr = PetscBinaryRead(fp,&numEmission_hist,1,PETSC_INT);CHKERRQ(ierr);  
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of points in emission files is %d \n",numEmission_hist);CHKERRQ(ierr);  
      ierr = PetscMalloc(numEmission_hist*sizeof(PetscScalar),&Tem_hist);CHKERRQ(ierr); 
      ierr = PetscBinaryRead(fp,Tem_hist,numEmission_hist,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      /* read fossil fuel emissions */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,emFiles[1],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
      ierr = PetscMalloc(numEmission_hist*sizeof(PetscScalar),&E_hist);CHKERRQ(ierr); 
      ierr = PetscBinaryRead(fp,E_hist,numEmission_hist,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      /* read land use emissions */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,emFiles[2],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
      ierr = PetscMalloc(numEmission_hist*sizeof(PetscScalar),&D_hist);CHKERRQ(ierr); 
      ierr = PetscBinaryRead(fp,D_hist,numEmission_hist,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      ierr = PetscOptionsHasName(PETSC_NULL,"-interp_emissions",&interpEmissions);CHKERRQ(ierr);
      if (interpEmissions) {      
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: Emissions will be interpolated in time. If you're prescribing annual emissions\n");CHKERRQ(ierr);     
        ierr = PetscPrintf(PETSC_COMM_WORLD,"         interpolation may lead to a different net emission input than what is prescribed\n");CHKERRQ(ierr);     
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Emissions will NOT be interpolated in time. It is assumed that you're prescribing annual emissions and that\n");CHKERRQ(ierr);     
        ierr = PetscPrintf(PETSC_COMM_WORLD,"the time data in file %s are the beginning of the year for which the corresponding emission is prescribed.\n",emFiles[0]);CHKERRQ(ierr);     
      }
    }  
        
  } else {  /* not using atm model */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using prescribed atmospheric pCO2\n");CHKERRQ(ierr);

    /* prescribed atmospheric CO2 */
  	maxValsToRead = 2;
	pCO2atmFiles[0] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* time file */
	pCO2atmFiles[1] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* atmospheric pCO2 history file */
    ierr = PetscOptionsGetStringArray(PETSC_NULL,"-pco2atm_history",pCO2atmFiles,&maxValsToRead,&flg);CHKERRQ(ierr);
    if (flg) { /* Read atmospheric pCO2 history */
      if (maxValsToRead != 2) {
        SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of file names specified for atmospheric pCO2 history");
      }      
      fixedAtmosCO2 = PETSC_FALSE;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading time-dependent atmospheric pCO2 history\n");CHKERRQ(ierr);      
      /* read time data */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,pCO2atmFiles[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	  ierr = PetscBinaryRead(fp,&numpCO2atm_hist,1,PETSC_INT);CHKERRQ(ierr);  
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of points in atmospheric history file is %d \n",numpCO2atm_hist);CHKERRQ(ierr);  
      ierr = PetscMalloc(numpCO2atm_hist*sizeof(PetscScalar),&TpCO2atm_hist);CHKERRQ(ierr); 
      ierr = PetscBinaryRead(fp,TpCO2atm_hist,numpCO2atm_hist,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      /* read atmospheric pCO2 data */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,pCO2atmFiles[1],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
      ierr = PetscMalloc(numpCO2atm_hist*sizeof(PetscScalar),&pCO2atm_hist);CHKERRQ(ierr); 
      ierr = PetscBinaryRead(fp,pCO2atm_hist,numpCO2atm_hist,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      
      pCO2atm = pCO2atm_hist[0];

    }	else {
      ierr = PetscOptionsGetReal(PETSC_NULL,"-pco2atm",&pCO2atm,&flg);CHKERRQ(ierr); /* overwrite default value */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Using fixed atmospheric pCO2 of %g ppm\n",pCO2atm);CHKERRQ(ierr);
      
    }    

  }
/* Grid arrays */
  ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localhFacC);CHKERRQ(ierr);  
  ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localrecip_hFacC);CHKERRQ(ierr);  
  
  ierr = VecLoadVecIntoArray(TR,"hFacC.petsc",localhFacC);CHKERRQ(ierr);
  ierr = VecLoadVecIntoArray(TR,"recip_hFacC.petsc",localrecip_hFacC);CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"drF.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fp,&nzmax,1,PETSC_INT);CHKERRQ(ierr);  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of vertical layers is %d \n",nzmax);CHKERRQ(ierr);  
  ierr = PetscMalloc(nzmax*sizeof(PetscScalar),&drF);CHKERRQ(ierr); 
  ierr = PetscBinaryRead(fp,drF,nzmax,PETSC_SCALAR);CHKERRQ(ierr);  
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

/* Forcing fields */  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&pH);CHKERRQ(ierr);

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localwind);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localwindp.firstTime = PETSC_TRUE;
    localwindp.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("wind.bin",localwind,1);  
  }
  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localatmosp);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localatmospp.firstTime = PETSC_TRUE;
    localatmospp.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("atmosp.bin",localatmosp,1);  
  }
  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localsilica);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localsilicap.firstTime = PETSC_TRUE;
    localsilicap.arrayLength = lNumProfiles;  
  } else {  
    ierr = readProfileSurfaceScalarData("silica.bin",localsilica,1);  
  }
  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localfice);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localficep.firstTime = PETSC_TRUE;
    localficep.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("fice.bin",localfice,1);  
  }

#ifdef ALLOW_FE
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localinputfe);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localinputfep.firstTime = PETSC_TRUE;
    localinputfep.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("inputfe.bin",localinputfe,1);  
  }
#endif

#ifdef LIGHT_CHL
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localchl);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localchlp.firstTime = PETSC_TRUE;
    localchlp.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("chl.bin",localchl,1);  
  }
#endif

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localEmP);CHKERRQ(ierr);
  if (periodicBiogeochemForcing) {    
	localEmPp.firstTime = PETSC_TRUE;
	localEmPp.arrayLength = lNumProfiles;
  } else {  
	ierr = readProfileSurfaceScalarData("EmP.bin",localEmP,1);  
  }

#ifdef READ_PAR
/* Read PAR from file */
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localpar);CHKERRQ(ierr);
  if (periodicBiogeochemForcing) {    
    localparp.firstTime = PETSC_TRUE;
    localparp.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("par.bin",localpar,1);  
  }  
#else
/* Read latitude from file */
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&locallat);CHKERRQ(ierr);  
  ierr = readProfileSurfaceScalarData("latitude.bin",locallat,1);  
#endif	  
  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localalpha);CHKERRQ(ierr);  
  for (ip=0; ip<lNumProfiles; ip++) { /* default is to let the model set alpha internally */
	localalpha[ip]=-1.0;
  }
  ierr = PetscOptionsGetString(PETSC_NULL,"-alpha_file",fileName,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (flg) { /* read from binary file and overwrite */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"alpha field is being read from file: %s\n",fileName);CHKERRQ(ierr);  
	ierr = readProfileSurfaceScalarData(fileName,localalpha,1);  
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localrain_ratio);CHKERRQ(ierr);  
  for (ip=0; ip<lNumProfiles; ip++) { /* default is to let the model set rain_ratio internally */
	localrain_ratio[ip]=-1.0;
  }
  ierr = PetscOptionsGetString(PETSC_NULL,"-rain_ratio_file",fileName,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (flg) { /* read from binary file and overwrite */
    ierr = PetscPrintf(PETSC_COMM_WORLD,"rain ratio field is being read from file: %s\n",fileName);CHKERRQ(ierr);    
	ierr = readProfileSurfaceScalarData(fileName,localrain_ratio,1);  
  }
  
  if (periodicBiogeochemForcing) {   
	ierr = interpPeriodicVector(tc,&Ts,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Ssp,"Ss_");	
	ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localEmPp,"EmP_");                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localwindp,"wind_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localatmospp,"atmosp_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localsilica,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localsilicap,"silica_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localficep,"fice_");
#ifdef READ_PAR                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localpar,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localparp,"par_");
#endif                                                                                                   
#ifdef ALLOW_FE                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localinputfe,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localinputfep,"inputfe_");
#endif                                                 
#ifdef LIGHT_CHL                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localchl,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localchlp,"chl_");
#endif                                                 
  }
    
/* Initialize biogeochem model */
  myTime = DeltaT*Iter; /* Iter should start at 0 */  
  for (ip=0; ip<lNumProfiles; ip++) {
    nzloc=lProfileLength[ip];
    kl=lStartIndices[ip];
/* 
    if (ip==0) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"ip=%d,kl=%d,nzloc=%d,th=%g,s=%g\n",ip,kl,nzloc,localTs[kl],localSs[kl]);CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"dic=%g,fe=%g\n",localDIC[kl],localFe[kl]);CHKERRQ(ierr);      
    }  

 */
	for (itr=0; itr<numTracers; itr++) {    	
	  mitgchem_copy_data_(&nzloc,&itr,&localTR[itr][kl],&localJTR[itr][kl],&DeltaT,&toModel);
	}  
 
    mitgchem_ini_(&nzloc,&numTracers,&Iter,&myTime,
                           &localTs[kl],&localSs[kl],&pH[ip],&localsilica[ip],
                           &localhFacC[kl],&localrecip_hFacC[kl],&drF[0],&DeltaT,&ip);
                           
  }

  ierr = PetscOptionsHasName(PETSC_NULL,"-calc_diagnostics",&calcDiagnostics);CHKERRQ(ierr);
  if (calcDiagnostics) {    
/*Data for diagnostics */
    ierr = iniStepTimer("diag_", Iter0, &diagTimer);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed starting at (and including) time step: %d\n", diagTimer.startTimeStep);CHKERRQ(ierr);	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed over %d time steps\n", diagTimer.numTimeSteps);CHKERRQ(ierr);	

	ierr = VecDuplicate(TR,&bioac);CHKERRQ(ierr);
	ierr = VecSet(bioac,zero);CHKERRQ(ierr);
	ierr = VecGetArray(bioac,&localbioac);CHKERRQ(ierr);
	ierr = VecDuplicate(TR,&bioacavg);CHKERRQ(ierr);
	ierr = VecSet(bioacavg,zero);CHKERRQ(ierr);
	
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localpco2diag);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localpco2diagavg);CHKERRQ(ierr);  

    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localgasexfluxdiag);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localgasexfluxdiagavg);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localempfluxdiag);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localempfluxdiagavg);CHKERRQ(ierr);  	

    for (ip=0; ip<lNumProfiles; ip++) {
      localpco2diag[ip]=0.0;
      localpco2diagavg[ip]=0.0;      
      localgasexfluxdiag[ip]=0.0;
      localgasexfluxdiagavg[ip]=0.0;      
      localempfluxdiag[ip]=0.0;
      localempfluxdiagavg[ip]=0.0;      
    }

    if (useAtmModel) {
      pCO2atmavg=0.0;
      Foceanavg=0.0;

	  if (useLandModel) {
		Flandavg=0.0;
		landStateavg[0]=0.0;
		landStateavg[1]=0.0;
		landStateavg[2]=0.0;
	  }      
    }
    
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
  PetscScalar zero = 0.0, one = 1.0;  
  PetscInt k;
  PetscScalar interpfac;  
  PetscInt itf;
  PetscScalar localFocean;
  PetscScalar recipdzsurfloc;
  PetscScalar localpco2 = 0.0, localgasexflux = 0.0;
  PetscScalar localempflux = 0.0, localempalkflux = 0.0;

  myTime = DeltaT*Iter; /* Iter should start at 0 */

  if (periodicBiogeochemForcing) {   
	ierr = interpPeriodicVector(tc,&Ts,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Ssp,"Ss_");	
	ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localEmPp,"EmP_");                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localwindp,"wind_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localatmospp,"atmosp_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localsilica,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localsilicap,"silica_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localficep,"fice_");
#ifdef ALLOW_FE                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localinputfe,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localinputfep,"inputfe_");
#endif                                                 
#ifdef READ_PAR                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localpar,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localparp,"par_");
#endif                                                                                                   
#ifdef LIGHT_CHL                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localchl,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localchlp,"chl_");
#endif                                                 
  }

  if (useAtmModel) {
    if (useEmissions) {  
/*   Interpolate emissions */
      if (tc>=Tem_hist[0]) {
        if (interpEmissions) {
		  ierr = calcInterpFactor(numEmission_hist,tc,Tem_hist,&itf,&interpfac); CHKERRQ(ierr);
		  fossilFuelEmission = interpfac*E_hist[itf] + (1.0-interpfac)*E_hist[itf+1];	  
		  landUseEmission = interpfac*D_hist[itf] + (1.0-interpfac)*D_hist[itf+1];	          
        } else {
          itf=findindex(Tem_hist,numEmission_hist,floor(tc));
		  fossilFuelEmission = E_hist[itf];
		  landUseEmission = D_hist[itf];          
        }
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: time < %10.5f. Setting emissions to 0\n",Tem_hist[0]);CHKERRQ(ierr);
        fossilFuelEmission = 0.0;
        landUseEmission = 0.0;
      }
      cumulativeEmissions = cumulativeEmissions + atmModelDeltaT*(fossilFuelEmission + landUseEmission); /* PgC */
    }
  } else {  
/* Interpolate atmospheric pCO2   */
    if (!fixedAtmosCO2) { 
      if (tc>=TpCO2atm_hist[0]) {
        ierr = calcInterpFactor(numpCO2atm_hist,tc,TpCO2atm_hist,&itf,&interpfac); CHKERRQ(ierr);
        pCO2atm = interpfac*pCO2atm_hist[itf] + (1.0-interpfac)*pCO2atm_hist[itf+1];	  
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: time < %10.5f. Assuming pCO2atm=%g\n",TpCO2atm_hist[0],pCO2atm);CHKERRQ(ierr);
      }
    }  
  }

  localFocean = 0.0;  
  Focean = 0.0;
  Fland = 0.0;

  for (ip=0; ip<lNumProfiles; ip++) {
    nzloc=lProfileLength[ip];
    kl=lStartIndices[ip];

	recipdzsurfloc = localrecip_hFacC[kl]/drF[0];
	localempflux = localEmP[ip]*localDIC[kl]*recipdzsurfloc;
	localempalkflux = localEmP[ip]*localAlk[kl]*recipdzsurfloc;

    if (useSeparateBiogeochemTimeStepping) { /* update DIC and Alk first */
      localDIC[kl] = localDIC[kl] + localempflux*DeltaT;
      localAlk[kl] = localAlk[kl] + localempalkflux*DeltaT;      
    }

	for (itr=0; itr<numTracers; itr++) {    	
	  mitgchem_copy_data_(&nzloc,&itr,&localTR[itr][kl],&localJTR[itr][kl],&DeltaT,&toModel);
	}  
 
    mitgchem_model_(&nzloc,&Iter,&myTime,
                            &localTs[kl],&localSs[kl],&pH[ip],&localwind[ip],
                            &localatmosp[ip],&pCO2atm,
                            &localsilica[ip],&localfice[ip],
#ifdef READ_PAR                            
                            &localpar[ip],
#else
                            &locallat[ip],
#endif
#ifdef ALLOW_FE                            
                            &localinputfe[ip],
#endif
#ifdef LIGHT_CHL                            
                            &localchl[ip],
#endif
                            &localalpha[ip],&localrain_ratio[ip],
                            &localhFacC[kl],&localrecip_hFacC[kl],
                            &localpco2,&localgasexflux,
                            &ip);

	for (itr=0; itr<numTracers; itr++) {    
	  mitgchem_copy_data_(&nzloc,&itr,&localTR[itr][kl],&localJTR[itr][kl],&DeltaT,&fromModel);
	}  

    if (!useSeparateBiogeochemTimeStepping) {
	  localJDIC[kl] = localJDIC[kl] + localempflux;
	  localJAlk[kl] = localJAlk[kl] + localempalkflux;
    }
    
    if (useAtmModel) {                 
      localFocean = localFocean + (localgasexflux+localempflux)*localdA[ip]*(12.0/1.e15)*secondsPerYear; /* PgC/y */
    }

	if (calcDiagnostics) {  
	  if (Iter0+iLoop>=diagTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */	
        mitgchem_diagnostics_(&nzloc,&Iter,&myTime,&localbioac[kl]);
        localpco2diag[ip]=localpco2;
        localgasexfluxdiag[ip]=localgasexflux;
        localempfluxdiag[ip]=localempflux;
      }
	}
  } /* end loop over profiles */

  if (calcDiagnostics) {  
	if (Iter0+iLoop>=diagTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */  
	  ierr = VecSetValues(bioac,lSize,gIndices,localbioac,INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecAssemblyBegin(bioac);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(bioac);CHKERRQ(ierr);    	  
    }
  }
  
  if (useAtmModel) {
    MPI_Allreduce(&localFocean, &Focean, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);    

	if (useLandModel) {
	  landsource_(&landState[0],&pCO2atm,&landUseEmission,&deltaTsg,&Fland,&landSource[0]); /* returns S and Fl in PgC/y */
/*    time step land */
	  for (k=0;k<=2;k++) {
		landState[k] = landState[k] + atmModelDeltaT*landSource[k];
	  }
	}

/*  time step atmosphere */
    pCO2atm = pCO2atm + atmModelDeltaT*(fossilFuelEmission + landUseEmission - Focean - Fland)/ppmToPgC;
  }  

  if (useSeparateBiogeochemTimeStepping) {  /* return updated tracer field */
	for (itr=0; itr<numTracers; itr++) {  
	  ierr = VecSetValues(v[itr],lSize,gIndices,localTR[itr],INSERT_VALUES);CHKERRQ(ierr);
	  ierr = VecAssemblyBegin(v[itr]);CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(v[itr]);CHKERRQ(ierr);    
	}
  } else {
	for (itr=0; itr<numTracers; itr++) {  
	  ierr = VecSetValues(ut[itr],lSize,gIndices,localJTR[itr],INSERT_VALUES);CHKERRQ(ierr);
	  ierr = VecAssemblyBegin(ut[itr]);CHKERRQ(ierr);
	  ierr = VecAssemblyEnd(ut[itr]);CHKERRQ(ierr);    
	}
  
#if defined (FORSPINUP) || defined (FORJACOBIAN)
/* add relaxation term: ut = ut - lambda*(v-vr) = ut -lambda*v + lambda*vr */
    if (relaxTracer) {
      for (itr=0; itr<numTracers; itr++) {
        ierr = VecAXPY(ut[itr],-relaxLambda[itr],v[itr]);CHKERRQ(ierr); /* ut = ut - lambda*v */
        ierr = VecShift(ut[itr],relaxLambda[itr]*relaxValue[itr]);CHKERRQ(ierr); /* ut = ut + lambda*vr */
      }
    }
#endif
  
/*  Convert to discrete tendency */
	for (itr=0; itr<numTracers; itr++) {
	  ierr = VecScale(ut[itr],DeltaT);CHKERRQ(ierr);
	}
  }
  
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

  if (useAtmModel) {
/* write instantaneous atmos model state */
	if (Iter0+iLoop>=atmWriteTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
	  if (atmWriteTimer.count<=atmWriteTimer.numTimeSteps) {
		atmWriteTimer.count++;
	  }
	  if (atmWriteTimer.count==atmWriteTimer.numTimeSteps) { /* time to write out */
    
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing atmospheric model output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
		ierr = PetscFPrintf(PETSC_COMM_WORLD,atmfptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);           
		ierr = writeBinaryScalarData("pCO2atm_output.bin",&pCO2atm,1,PETSC_TRUE);

		if (useEmissions) {
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Cumulative emissions at time %10.5f, step %d = %10.6f PgC\n", tc, Iter0+iLoop, cumulativeEmissions);CHKERRQ(ierr);
		}  

		if (useLandModel) {
/*        write instantaneous land model state */
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing land model output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
		  ierr = writeBinaryScalarData("land_state_output.bin",landState,3,PETSC_TRUE);
		}

		ierr = updateStepTimer("atm_write_", Iter0+iLoop, &atmWriteTimer);CHKERRQ(ierr);

      }
    }
  }

  if (calcDiagnostics) {  
	if (Iter0+iLoop>=diagTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */  

	  if (diagTimer.count<=diagTimer.numTimeSteps) { /* still within same averaging block so accumulate */

		ierr = VecAXPY(bioacavg,one,bioac);CHKERRQ(ierr);
        for (ip=0; ip<lNumProfiles; ip++) {
          localpco2diagavg[ip]=localpco2diag[ip]+localpco2diagavg[ip];
          localgasexfluxdiagavg[ip]=localgasexfluxdiag[ip]+localgasexfluxdiagavg[ip];
          localempfluxdiagavg[ip]=localempfluxdiag[ip]+localempfluxdiagavg[ip];
        }

        if (useAtmModel) {
          pCO2atmavg=pCO2atm+pCO2atmavg;
          Foceanavg=Focean+Foceanavg;

		  if (useLandModel) {
			Flandavg=Fland+Flandavg;
			landStateavg[0]=landState[0]+landStateavg[0];
			landStateavg[1]=landState[1]+landStateavg[1];
			landStateavg[2]=landState[2]+landStateavg[2];          
		  }                  
        }        
        
		diagTimer.count++;
	  }

	  if (diagTimer.count==diagTimer.numTimeSteps) { /* time to write averages to file */
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing diagnostics time average at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);                      
		ierr = VecScale(bioacavg,1.0/diagTimer.count);CHKERRQ(ierr);
        for (ip=0; ip<lNumProfiles; ip++) {
          localpco2diagavg[ip]=localpco2diagavg[ip]/diagTimer.count;
          localgasexfluxdiagavg[ip]=localgasexfluxdiagavg[ip]/diagTimer.count;
          localempfluxdiagavg[ip]=localempfluxdiagavg[ip]/diagTimer.count;
        }	  

        if (!appendDiagnostics) {
  		  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"bioac.petsc",FILE_MODE_WRITE,&fdbioacavg);CHKERRQ(ierr);		
        }
		ierr = VecView(bioacavg,fdbioacavg);CHKERRQ(ierr);

		ierr = writeProfileSurfaceScalarData("pco2_surf.bin",localpco2diagavg,1,appendDiagnostics);  		
		ierr = writeProfileSurfaceScalarData("gasexflux_surf.bin",localgasexfluxdiagavg,1,appendDiagnostics);  		
		ierr = writeProfileSurfaceScalarData("empflux_surf.bin",localempfluxdiagavg,1,appendDiagnostics);  		

        if (useAtmModel) {
          pCO2atmavg=pCO2atmavg/diagTimer.count;
          Foceanavg=Foceanavg/diagTimer.count;   
          ierr = writeBinaryScalarData("pCO2atm_avg.bin",&pCO2atmavg,1,appendDiagnostics);  		
          ierr = writeBinaryScalarData("Focean_avg.bin",&Foceanavg,1,appendDiagnostics);  		

		  if (useLandModel) {
			Flandavg=Flandavg/diagTimer.count;        
			landStateavg[0]=landStateavg[0]/diagTimer.count;
			landStateavg[1]=landStateavg[1]/diagTimer.count;
			landStateavg[2]=landStateavg[2]/diagTimer.count;           
			ierr = writeBinaryScalarData("Fland_avg.bin",&Flandavg,1,appendDiagnostics);  		
			ierr = writeBinaryScalarData("land_state_avg.bin",landStateavg,3,appendDiagnostics);
		  }          
        }

        appendDiagnostics=PETSC_TRUE;

/*      reset diagnostic arrays */
		ierr = VecSet(bioacavg,zero); CHKERRQ(ierr);
        for (ip=0; ip<lNumProfiles; ip++) {
          localpco2diagavg[ip]=0.0;
          localgasexfluxdiagavg[ip]=0.0;
          localempfluxdiagavg[ip]=0.0;
        }	  

        if (useAtmModel) {
          pCO2atmavg=0.0;
          Foceanavg=0.0;

		  if (useLandModel) {
			Flandavg=0.0;
			landStateavg[0]=0.0;
			landStateavg[1]=0.0;
			landStateavg[2]=0.0;          
		  }          
        }
        ierr = updateStepTimer("diag_", Iter0+iLoop, &diagTimer);CHKERRQ(ierr);

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

/* write final pickup */
  if (useAtmModel) {
/* write instantaneous atmos model state */
    ierr = writeBinaryScalarData("pickup_pCO2atm.bin",&pCO2atm,1,PETSC_FALSE);

	if (useLandModel) {
/*    write instantaneous land model state */
	  ierr = writeBinaryScalarData("pickup_land_state.bin",landState,3,PETSC_FALSE);
	}    
  }
  
  ierr = VecDestroy(&Ts);CHKERRQ(ierr);
  ierr = VecDestroy(&Ss);CHKERRQ(ierr);

  if (periodicBiogeochemForcing) {    
    ierr = destroyPeriodicVec(&Tsp);CHKERRQ(ierr);
    ierr = destroyPeriodicVec(&Ssp);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localwindp);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localatmospp);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localsilicap);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localficep);CHKERRQ(ierr);
#ifdef READ_PAR
    ierr = destroyPeriodicArray(&localparp);CHKERRQ(ierr);
#endif
#ifdef ALLOW_FE    
    ierr = destroyPeriodicArray(&localinputfep);CHKERRQ(ierr);
#endif
#ifdef LIGHT_CHL
    ierr = destroyPeriodicArray(&localchlp);CHKERRQ(ierr);
#endif
  }    

  if (useAtmModel) {
    ierr = PetscFClose(PETSC_COMM_WORLD,atmfptime);CHKERRQ(ierr);
  }

  if (calcDiagnostics) {  
	ierr = VecDestroy(&bioac);CHKERRQ(ierr);
	ierr = VecDestroy(&bioacavg);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fdbioacavg);CHKERRQ(ierr);	
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "reInitializeExternalForcing"
PetscErrorCode reInitializeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt numTracers, Vec *v, Vec *ut)
{
  PetscErrorCode ierr;
  PetscInt ip, kl, nzloc;
  PetscScalar myTime;

  if (periodicBiogeochemForcing) {   
	ierr = interpPeriodicVector(tc,&Ts,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Ssp,"Ss_");	
	ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localEmPp,"EmP_");                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localwindp,"wind_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localatmospp,"atmosp_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localsilica,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localsilicap,"silica_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localficep,"fice_");
#ifdef READ_PAR                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localpar,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localparp,"par_");
#endif                                                                                                   
#ifdef ALLOW_FE                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localinputfe,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localinputfep,"inputfe_");
#endif                                                 
#ifdef LIGHT_CHL                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localchl,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localchlp,"chl_");
#endif                                                 
  }
    
/* Initialize biogeochem model */
  myTime = DeltaT*Iter; /* Iter should start at 0 */
  for (ip=0; ip<lNumProfiles; ip++) {
    nzloc=lProfileLength[ip];
    kl=lStartIndices[ip];
    mitgchem_ini_(&nzloc,&numTracers,&Iter,&myTime,
                           &localTs[kl],&localSs[kl],&pH[ip],&localsilica[ip],
                           &localhFacC[kl],&localrecip_hFacC[kl],&drF[0],&DeltaT,&ip);
  }
  
  return 0;
}
