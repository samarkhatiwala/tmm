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
#include "TRACEGASES_OPTIONS.h"
#include "tracegases.h"

/* Macros to map tracer names to vectors */
#define TR v[0]
#define JTR ut[0]

Vec surfVolFrac;
Vec Ts,Ss;
PetscScalar *localTs,*localSs;
PetscScalar *localTR;
PetscScalar *localJTR;
PetscScalar *localfice,*localxkw, *localatmosp;
PetscScalar *localVgas660;
PetscScalar DeltaT, *localdzsurf;
PetscBool useVirtualFlux = PETSC_FALSE;
PetscScalar *localEmP;
PeriodicArray localEmPp;
PetscBool useWinds = PETSC_FALSE;
PetscScalar *localuwind,*localvwind;
PeriodicArray localuwindp, localvwindp;
PetscInt numWindsPeriods;
PetscScalar *tdpWinds; /* arrays for periodic forcing */
PetscScalar windsCyclePeriod, windsCycleStep;
PetscScalar pistonVelocity=0.31; /* default piston velocity when using winds [cm/hr] */

PeriodicVec Tsp, Ssp;
PeriodicArray localficep, localxkwp, localatmospp;
PetscInt numBiogeochemPeriods;
PetscScalar *tdpBiogeochem; /* arrays for periodic forcing */
PetscBool periodicBiogeochemForcing = PETSC_FALSE;
PetscScalar biogeochemCyclePeriod, biogeochemCycleStep;

PetscInt maxValsToRead;
PetscInt dummyInt = 0;

PetscInt gasID=0;
PetscScalar mixingRatioScaleFactor=0.0; /* scale factor for mixing ratio, e.g., if mixing ratio is in ppm, mixingRatioScaleFactor=1.e-6 */
char *xTRatmFiles[2];  
PetscInt numxTRatm_hist = 0;
PetscScalar *TxTRatm_hist, *xTRatm_hist, *xTRatm_hist0, *xTRatm_hist1;
PetscInt ithist = -1;
PetscBool fixedAtmosMixingRatio = PETSC_TRUE;
PetscBool spatiallyVariableMixingRatio = PETSC_FALSE;
char xTRatmIniFile[PETSC_MAX_PATH_LEN];  

/* Atm model variables */
PetscBool useAtmModel = PETSC_FALSE;
PetscBool useEmissions = PETSC_FALSE;
PetscBool interpEmissions = PETSC_FALSE;
char *emFiles[2];  
PetscInt numEmission_hist = 0;
PetscScalar *Tem_hist, *E_hist;
PetscScalar annualEmission = 0.0;
PetscScalar cumulativeEmissions = 0.0;
PetscScalar xTRatm_ini = 0.0; /* default initial value */
PetscScalar xTRatm = 0.0; /* default initial value */
PetscScalar *localdA;
PetscScalar atmModelDeltaT;
PetscScalar secPerYear=86400.0*360.0;
PetscScalar Focean=0.0;
PetscScalar mixingRatioToMass=0.0; /* This converts the mixing ratio of the gas to mass in the units used to specify emissions */
PetscScalar moleToMass=0.0; /* This converts 1 mole of the gas to mass in the units used to specify emissions */

PetscInt atmWriteSteps;
PetscBool atmAppendOutput;
FILE *atmfptime;
PetscViewer atmfd;
PetscInt atmfp;
char atmOutTimeFile[PETSC_MAX_PATH_LEN];  
PetscScalar xTRatmavg, Foceanavg;

PetscBool calcDiagnostics = PETSC_FALSE;
PetscInt diagNumTimeSteps, diagStartTimeStep, diagCount;
PetscBool appendDiagnostics = PETSC_FALSE;
PetscScalar *localTReqdiag, *localTReqdiagavg;
PetscScalar *localgasexfluxdiag, *localgasexfluxdiagavg, *localtotfluxdiag, *localtotfluxdiagavg;

#undef __FUNCT__
#define __FUNCT__ "iniExternalForcing"
PetscErrorCode iniExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt numTracers, Vec *v, Vec *ut)
{
  PetscErrorCode ierr;
  PetscInt ip, kl;
  PetscViewer fd;
  PetscInt fp;
  PetscBool flg, flg1;
  PetscInt it;
  PetscScalar myTime;
  PetscScalar zero = 0.0;

  ierr = VecGetArray(TR,&localTR);CHKERRQ(ierr);

  ierr = VecSet(JTR,zero); CHKERRQ(ierr);    
  ierr = VecGetArray(JTR,&localJTR);CHKERRQ(ierr);

  ierr = PetscOptionsGetReal(PETSC_NULL,"-mixing_ratio_scale_factor",&mixingRatioScaleFactor,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate the scale factor to convert atmospheric mixing ratio to atm with the -mixing_ratio_scale_factor option");  

  ierr = PetscOptionsGetInt(PETSC_NULL,"-gas_id",&gasID,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate the gas ID as an integer with the -gas_id option: 1=N2O, 2=CFC11, 3=CFC12, 4=SF6");

  ierr = PetscOptionsGetReal(PETSC_NULL,"-biogeochem_deltat",&DeltaT,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate biogeochemical time step in seconds with the -biogeochem_deltat option");  

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
	ierr = VecLoad(Ts,fd);CHKERRQ(ierr);  
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Ss.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(Ss,fd);CHKERRQ(ierr);    
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
  }  
  ierr = VecGetArray(Ts,&localTs);CHKERRQ(ierr);
  ierr = VecGetArray(Ss,&localSs);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Done reading T/S\n");CHKERRQ(ierr);

/* Atm model data */
  ierr = PetscOptionsHasName(PETSC_NULL,"-use_atm_model",&useAtmModel);CHKERRQ(ierr);

  if (useAtmModel) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using interactive atmospheric model\n");CHKERRQ(ierr);  

    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localdA);CHKERRQ(ierr);
    ierr = readProfileSurfaceScalarData("dA.bin",localdA,1);  

    ierr = PetscOptionsGetReal(PETSC_NULL,"-mixing_ratio_to_mass",&mixingRatioToMass,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate the mixing ratio to mass value for this tracer with the -mixing_ratio_to_mass option");  

    ierr = PetscOptionsGetReal(PETSC_NULL,"-mole_to_mass",&moleToMass,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate the mole to mass ratio for this tracer with the -mole_to_mass option");  

/* overwrite default value */
	ierr = PetscOptionsGetReal(PETSC_NULL,"-xTRatm_ini",&xTRatm_ini,&flg);CHKERRQ(ierr); /* read from command line */
    if (!flg) {
      ierr = PetscOptionsGetString(PETSC_NULL,"-xTRatm_ini_file",xTRatmIniFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
      if (flg) { /* read from binary file */
        ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,xTRatmIniFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
        ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
        ierr = PetscBinaryRead(fp,&xTRatm_ini,1,PETSC_SCALAR);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      }
    }
    xTRatm = xTRatm_ini;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using initial atmospheric mixing ratio of %g\n",xTRatm);CHKERRQ(ierr);
      
    atmModelDeltaT = DeltaT/secPerYear; /* time step in years */

    ierr = PetscOptionsGetInt(PETSC_NULL,"-atm_write_steps",&atmWriteSteps,&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate atmospheric model output step with the -atm_write_steps option");

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
      ierr = writeBinaryScalarData("xTRatm_output.bin",&xTRatm,1,PETSC_FALSE);
    } else {
      ierr = PetscFOpen(PETSC_COMM_WORLD,atmOutTimeFile,"a",&atmfptime);CHKERRQ(ierr);  
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output will be appended. Initial condition will NOT be written\n");CHKERRQ(ierr);      
    }

    /* TR emissions */
	maxValsToRead = 3;
	emFiles[0] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* time file */
	emFiles[1] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* emissions file */
    ierr = PetscOptionsGetStringArray(PETSC_NULL,"-emissions_history",emFiles,&maxValsToRead,&useEmissions);CHKERRQ(ierr);
    if (useEmissions) { /* Read emissions history */
      if (maxValsToRead != 2) {
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
      /* read emissions */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,emFiles[1],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
      ierr = PetscMalloc(numEmission_hist*sizeof(PetscScalar),&E_hist);CHKERRQ(ierr); 
      ierr = PetscBinaryRead(fp,E_hist,numEmission_hist,PETSC_SCALAR);CHKERRQ(ierr);
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
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using prescribed atmospheric mixing ratio\n");CHKERRQ(ierr);

    /* prescribed atmospheric mixing ratio */
  	maxValsToRead = 2;
	xTRatmFiles[0] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* time file */
	xTRatmFiles[1] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* atmospheric mixing ratio history file */
    ierr = PetscOptionsGetStringArray(PETSC_NULL,"-xTRatm_history",xTRatmFiles,&maxValsToRead,&flg);CHKERRQ(ierr);
    if (flg) { /* Read atmospheric xTR history */
      if (maxValsToRead != 2) {
        SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of file names specified for atmospheric mixing ratio history");
      }
      fixedAtmosMixingRatio = PETSC_FALSE;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading time-dependent atmospheric mixing ratio history\n");CHKERRQ(ierr);
      ierr = PetscOptionsHasName(PETSC_NULL,"-spatially_variable_mixing_ratio",&spatiallyVariableMixingRatio);CHKERRQ(ierr);
      if (spatiallyVariableMixingRatio) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"   Spatially-uniform atmospheric mixing ratio history has been specified\n");CHKERRQ(ierr);
      } else {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"   Spatially-variable atmospheric mixing ratio history has been specified\n");CHKERRQ(ierr);      
      }
      /* read time data */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,xTRatmFiles[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	  ierr = PetscBinaryRead(fp,&numxTRatm_hist,1,PETSC_INT);CHKERRQ(ierr);  
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of points in atmospheric mixing ratio history file is %d \n",numxTRatm_hist);CHKERRQ(ierr);  
      ierr = PetscMalloc(numxTRatm_hist*sizeof(PetscScalar),&TxTRatm_hist);CHKERRQ(ierr); 
      ierr = PetscBinaryRead(fp,TxTRatm_hist,numxTRatm_hist,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      /* read atmospheric xTR data */
      if (spatiallyVariableMixingRatio) {
		ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&xTRatm_hist);CHKERRQ(ierr);
		ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&xTRatm_hist0);CHKERRQ(ierr);
		ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&xTRatm_hist1);CHKERRQ(ierr);
		ierr = readProfileSurfaceScalarData(xTRatmFiles[1],xTRatm_hist,1); /* read initial slice */ 
		ithist=-1;
      } else {      
		ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,xTRatmFiles[1],FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
		ierr = PetscMalloc(numxTRatm_hist*sizeof(PetscScalar),&xTRatm_hist);CHKERRQ(ierr); 
		ierr = PetscBinaryRead(fp,xTRatm_hist,numxTRatm_hist,PETSC_SCALAR);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);            
        xTRatm = xTRatm_hist[0];
      }  

    } else {
      ierr = PetscOptionsGetReal(PETSC_NULL,"-xTRatm",&xTRatm,&flg);CHKERRQ(ierr); /* overwrite default value */
	  if (flg) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Using fixed, spatially-constant atmospheric mixing ratio of %g\n",xTRatm);CHKERRQ(ierr);
	  } else {	
	    SETERRQ(PETSC_COMM_WORLD,1,"No atmospheric mixing ratio has been specified!");      
	  }  
    }  

  } /* useAtmModel */

/* Grid arrays */
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localdzsurf);CHKERRQ(ierr);
  ierr = readProfileSurfaceScalarData("dzsurf.bin",localdzsurf,1);  

/* Forcing fields */  

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localVgas660);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) { /* read monthly mean ice and xkw or winds */
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localfice);CHKERRQ(ierr);    
    localficep.firstTime = PETSC_TRUE;
    localficep.arrayLength = lNumProfiles;

    ierr = PetscOptionsHasName(PETSC_NULL,"-use_winds",&useWinds);CHKERRQ(ierr);
    if (useWinds) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Periodic winds specified: gas transfer velocity will be computed using winds\n");CHKERRQ(ierr);      
      ierr = PetscOptionsGetReal(PETSC_NULL,"-piston_velocity",&pistonVelocity,&flg);CHKERRQ(ierr); /* overwrite default value */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Piston velocity of %10.5f cm/hr will be used\n", pistonVelocity);CHKERRQ(ierr);      
      pistonVelocity = pistonVelocity*0.01/3600.0; /* convert cm/hr to m/s */
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
      ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localxkw);CHKERRQ(ierr);    
      localxkwp.firstTime = PETSC_TRUE;
      localxkwp.arrayLength = lNumProfiles;
    }   
  } else { /* read annual mean Vgas660 */
    ierr = readProfileSurfaceScalarData("Vgas660.bin",localVgas660,1);  
  }  

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localatmosp);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localatmospp.firstTime = PETSC_TRUE;
    localatmospp.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("atmosp.bin",localatmosp,1);  
  }

  ierr = PetscOptionsHasName(PETSC_NULL,"-use_virtual_flux",&useVirtualFlux);CHKERRQ(ierr);

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localEmP);CHKERRQ(ierr);
  if (periodicBiogeochemForcing) {    
	localEmPp.firstTime = PETSC_TRUE;
	localEmPp.arrayLength = lNumProfiles;
  } else {  
	ierr = readProfileSurfaceScalarData("EmP.bin",localEmP,1);  
  }

  if (useVirtualFlux) {
	ierr = VecDuplicate(TR,&surfVolFrac);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"surface_volume_fraction.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(surfVolFrac,fd);CHKERRQ(ierr);  
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
  }
    
  if (periodicBiogeochemForcing) {   
	ierr = interpPeriodicVector(tc,&Ts,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Ssp,"Ss_");	
	ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemCyclePeriod,numBiogeochemPeriods,
												  tdpBiogeochem,&localEmPp,"EmP_");                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localficep,"fice_");
    if (useWinds) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localuwind,windsCyclePeriod,numWindsPeriods,
                                                    tdpWinds,&localuwindp,"uwind_");                                                    
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localvwind,windsCyclePeriod,numWindsPeriods,
                                                    tdpWinds,&localvwindp,"vwind_");                                                        
    } else {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localxkw,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                    tdpBiogeochem,&localxkwp,"xkw_");
    }
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemCyclePeriod,numBiogeochemPeriods,
									              tdpBiogeochem,&localatmospp,"atmosp_");
  }  

/* Initialize biogeochem model */
  myTime = DeltaT*Iter; /* Iter should start at 0 */  

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localTReqdiag);CHKERRQ(ierr);  /* always need to pass this */

  ierr = PetscOptionsHasName(PETSC_NULL,"-calc_diagnostics",&calcDiagnostics);CHKERRQ(ierr);
  if (calcDiagnostics) {    
/*Data for diagnostics */
	ierr = PetscOptionsGetInt(PETSC_NULL,"-diag_start_time_step",&diagStartTimeStep,&flg);CHKERRQ(ierr);
	if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate (absolute) time step at which to start storing diagnostics with the -diag_start_time_step flag");
	ierr = PetscOptionsGetInt(PETSC_NULL,"-diag_time_steps",&diagNumTimeSteps,&flg);CHKERRQ(ierr);
	if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate number of time averaging diagnostics time steps with the -diag_time_step flag");
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed starting at (and including) time step: %d\n", diagStartTimeStep);CHKERRQ(ierr);	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed over %d time steps\n", diagNumTimeSteps);CHKERRQ(ierr);	

    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localTReqdiagavg);CHKERRQ(ierr);  

    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localgasexfluxdiag);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localgasexfluxdiagavg);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localtotfluxdiag);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localtotfluxdiagavg);CHKERRQ(ierr);  

    for (ip=0; ip<lNumProfiles; ip++) {
      localTReqdiag[ip]=0.0;
      localTReqdiagavg[ip]=0.0;      
      localgasexfluxdiag[ip]=0.0;
      localgasexfluxdiagavg[ip]=0.0;      
      localtotfluxdiag[ip]=0.0;
      localtotfluxdiagavg[ip]=0.0;      
    }

    if (useAtmModel) {
      xTRatmavg=0.0;
      Foceanavg=0.0;
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
  static PetscScalar TRemp = 0.0;
  PetscInt itr, ip, nzloc, kl;
  PetscScalar myTime;
  PetscScalar alpha;
  PetscInt itf;
  PetscScalar zero = 0.0, one = 1.0;
  PetscInt k;
  PetscScalar localFocean;
  PetscScalar localgasexflux = 0.0, localtotflux = 0.0;
  PetscScalar w2;
  
  myTime = DeltaT*Iter; /* Iter should start at 0 */

  if (periodicBiogeochemForcing) {   
	ierr = interpPeriodicVector(tc,&Ts,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemCyclePeriod,numBiogeochemPeriods,tdpBiogeochem,&Ssp,"Ss_");	
	ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemCyclePeriod,numBiogeochemPeriods,
												  tdpBiogeochem,&localEmPp,"EmP_");                                                      
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localficep,"fice_");
    if (useWinds) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localuwind,windsCyclePeriod,numWindsPeriods,
                                                    tdpWinds,&localuwindp,"uwind_");
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localvwind,windsCyclePeriod,numWindsPeriods,
                                                    tdpWinds,&localvwindp,"vwind_");    
    } else {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localxkw,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                    tdpBiogeochem,&localxkwp,"xkw_");
    }
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemCyclePeriod,numBiogeochemPeriods,
									              tdpBiogeochem,&localatmospp,"atmosp_");

/*  Recompute gas exchange coeff */
    if (useWinds) {
      for (ip=0; ip<lNumProfiles; ip++) {
        kl=lStartIndices[ip];
        w2 = pow(localuwind[ip],2) + pow(localvwind[ip],2);
        localVgas660[ip] = (1.0-localfice[ip])*pistonVelocity*w2;
      }                                                      
    } else {
      for (ip=0; ip<lNumProfiles; ip++) {
        kl=lStartIndices[ip];    
        localVgas660[ip] = (1.0-localfice[ip])*localxkw[ip];
      }
    }
  }

  if (useAtmModel) {
    if (useEmissions) {  
/*   Interpolate emissions */
      if (tc>=Tem_hist[0]) {
        if (interpEmissions) {
		  ierr = calcInterpFactor(numEmission_hist,tc,Tem_hist,&itf,&alpha); CHKERRQ(ierr);
		  annualEmission = alpha*E_hist[itf] + (1.0-alpha)*E_hist[itf+1];	  
        } else {
          itf=findindex(Tem_hist,numEmission_hist,floor(tc));
		  annualEmission = E_hist[itf];
        }
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: time < %10.5f. Setting emissions to 0\n",Tem_hist[0]);CHKERRQ(ierr);
        annualEmission = 0.0;
      }
      cumulativeEmissions = cumulativeEmissions + atmModelDeltaT*annualEmission; /* mass */
    }
  } else {  
/* Interpolate atmospheric xTR   */
    if (!fixedAtmosMixingRatio) { 
      if (tc>=TxTRatm_hist[0]) {
        ierr = calcInterpFactor(numxTRatm_hist,tc,TxTRatm_hist,&itf,&alpha); CHKERRQ(ierr);
        if (spatiallyVariableMixingRatio) {
		  if (itf != ithist) { /* time to read new bracketing slices: itf uses 0-based, while readProfileSurfaceScalarDataRecord uses 1-based indexing*/
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading new bracketing slices for spatially-variable mixing ratio at time = %g: %d and %d\n",tc,itf+1,itf+2);CHKERRQ(ierr);        
			ierr = readProfileSurfaceScalarDataRecord(xTRatmFiles[1],xTRatm_hist0,1,itf+1);
			ierr = readProfileSurfaceScalarDataRecord(xTRatmFiles[1],xTRatm_hist1,1,itf+2);
			ithist=itf;
		  }	
		  /* interpolate in time */
		  for (ip=0; ip<lNumProfiles; ip++) {
			xTRatm_hist[ip] = alpha*xTRatm_hist0[ip] + (1.0-alpha)*xTRatm_hist1[ip];	          
		  }			
        } else {
          /* interpolate in time */
          xTRatm = alpha*xTRatm_hist[itf] + (1.0-alpha)*xTRatm_hist[itf+1];
        }  
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: time < %10.5f. Using xTRatm at initial time.\n",TxTRatm_hist[0]);CHKERRQ(ierr);
      }
    }  
  }

  if (useVirtualFlux) { /* use the global surface mean value to calculate E-P contribution */
    ierr = VecDot(surfVolFrac,TR,&TRemp);CHKERRQ(ierr); /* volume weighted mean surface TR */									              
  }

  localFocean = 0.0;  
  Focean = 0.0;
  
/* Compute air-sea gas exchange term */
  for (ip=0; ip<lNumProfiles; ip++) {
    kl=lStartIndices[ip];
    nzloc=lProfileLength[ip];
	if (spatiallyVariableMixingRatio) xTRatm = xTRatm_hist[ip];
	if (!useVirtualFlux) TRemp=localTR[kl]; /* use the local surface value to calculate E-P contribution */
    tracegases_model_(&Iter,&myTime,
                         &localTR[kl],
                         &localTs[kl],&localSs[kl],&localVgas660[ip],
                         &localatmosp[ip],&gasID,&xTRatm,&mixingRatioScaleFactor,&localdzsurf[ip],
                         &localEmP[ip],&TRemp,
                         &localJTR[kl],&localgasexflux,&localtotflux,&localTReqdiag[ip]
                         );

    if (useAtmModel) {                 
      localFocean = localFocean + localtotflux*localdA[ip]*moleToMass*secPerYear; /* mass/y */
    }
    
	if (calcDiagnostics) {  
	  if (Iter0+iLoop>=diagStartTimeStep) { /* start time averaging (note: diagStartTimeStep is ABSOLUTE time step) */	
        localgasexfluxdiag[ip]=localgasexflux;
        localtotfluxdiag[ip]=localtotflux;
      }
	}                         
  } /* end loop over profiles */
  
  if (useAtmModel) {
	MPI_Allreduce(&localFocean, &Focean, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);    

/*  time step atmosphere */
    xTRatm = xTRatm + atmModelDeltaT*(annualEmission - Focean)/mixingRatioToMass;
  }  

  ierr = VecSetValues(JTR,lSize,gIndices,localJTR,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(JTR);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(JTR);CHKERRQ(ierr);    

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

  if (useAtmModel) {
/* write instantaneous atmos model state */
    if ((iLoop % atmWriteSteps)==0) {  /*  time to write out */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing atmospheric model output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
      ierr = PetscFPrintf(PETSC_COMM_WORLD,atmfptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);           
      ierr = writeBinaryScalarData("xTRatm_output.bin",&xTRatm,1,PETSC_TRUE);

	  if (useEmissions) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Cumulative emissions at time %10.5f, step %d = %10.6f\n", tc, Iter0+iLoop, cumulativeEmissions);CHKERRQ(ierr);
      }  
    }

  }

  if (calcDiagnostics) {  
	if (Iter0+iLoop>=diagStartTimeStep) { /* start time averaging (note: diagStartTimeStep is ABSOLUTE time step) */  

	  if (diagCount<=diagNumTimeSteps) { /* still within same averaging block so accumulate */
	  
        for (ip=0; ip<lNumProfiles; ip++) {
          localTReqdiagavg[ip]=localTReqdiag[ip]+localTReqdiagavg[ip];
          localgasexfluxdiagavg[ip]=localgasexfluxdiag[ip]+localgasexfluxdiagavg[ip];
          localtotfluxdiagavg[ip]=localtotfluxdiag[ip]+localtotfluxdiagavg[ip];      
        }	  

        if (useAtmModel) {
          xTRatmavg=xTRatm+xTRatmavg;
          Foceanavg=Focean+Foceanavg;
        }        
        
		diagCount = diagCount+1;
	  }

	  if (diagCount==diagNumTimeSteps) { /* time to write averages to file */
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing diagnostics time average at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);                      

        for (ip=0; ip<lNumProfiles; ip++) {
          localTReqdiagavg[ip]=localTReqdiagavg[ip]/diagCount;
          localgasexfluxdiagavg[ip]=localgasexfluxdiagavg[ip]/diagCount;
          localtotfluxdiagavg[ip]=localtotfluxdiagavg[ip]/diagCount;
        }	  

        ierr = writeProfileSurfaceScalarData("TReq_surf.bin",localTReqdiagavg,1,appendDiagnostics);  		
        ierr = writeProfileSurfaceScalarData("gasexflux_surf.bin",localgasexfluxdiagavg,1,appendDiagnostics);  		
        ierr = writeProfileSurfaceScalarData("totalflux_surf.bin",localtotfluxdiagavg,1,appendDiagnostics);  		
        if (useAtmModel) {
          xTRatmavg=xTRatmavg/diagCount;
          Foceanavg=Foceanavg/diagCount;   
          ierr = writeBinaryScalarData("xTRatm_avg.bin",&xTRatmavg,1,appendDiagnostics);  		
          ierr = writeBinaryScalarData("Focean_avg.bin",&Foceanavg,1,appendDiagnostics);  		
        }

        appendDiagnostics=PETSC_TRUE;

/*      reset diagnostic arrays */
        for (ip=0; ip<lNumProfiles; ip++) {
          localTReqdiagavg[ip]=0.0;
          localgasexfluxdiagavg[ip]=0.0;
          localtotfluxdiagavg[ip]=0.0;
        }	  


        if (useAtmModel) {
          xTRatmavg=0.0;
          Foceanavg=0.0;
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

/* write final pickup */
  if (useAtmModel) {
/* write instantaneous atmos model state */
    ierr = writeBinaryScalarData("pickup_xTRatm.bin",&xTRatm,1,PETSC_FALSE);
  }
  
  ierr = VecDestroy(&Ts);CHKERRQ(ierr);
  ierr = VecDestroy(&Ss);CHKERRQ(ierr);

  if (periodicBiogeochemForcing) {    
    ierr = destroyPeriodicVec(&Tsp);CHKERRQ(ierr);
    ierr = destroyPeriodicVec(&Ssp);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localficep);CHKERRQ(ierr);
    if (useWinds) {
      ierr = destroyPeriodicArray(&localuwindp);CHKERRQ(ierr);
      ierr = destroyPeriodicArray(&localvwindp);CHKERRQ(ierr);
    } else {   
      ierr = destroyPeriodicArray(&localxkwp);CHKERRQ(ierr);
    }
    ierr = destroyPeriodicArray(&localatmospp);CHKERRQ(ierr);
	ierr = destroyPeriodicArray(&localEmPp);CHKERRQ(ierr);
  }    

  if (useAtmModel) {
    ierr = PetscFClose(PETSC_COMM_WORLD,atmfptime);CHKERRQ(ierr);
  }
  
  if (useVirtualFlux) {
    ierr = VecDestroy(&surfVolFrac);CHKERRQ(ierr);  
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
	ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemCyclePeriod,numBiogeochemPeriods,
												  tdpBiogeochem,&localEmPp,"EmP_");                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                  tdpBiogeochem,&localficep,"fice_");
    if (useWinds) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localuwind,windsCyclePeriod,numWindsPeriods,
                                                    tdpWinds,&localuwindp,"uwind_");
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localvwind,windsCyclePeriod,numWindsPeriods,
                                                    tdpWinds,&localvwindp,"vwind_");    
    } else {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localxkw,biogeochemCyclePeriod,numBiogeochemPeriods,
                                                    tdpBiogeochem,&localxkwp,"xkw_");
    }                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemCyclePeriod,numBiogeochemPeriods,
									              tdpBiogeochem,&localatmospp,"atmosp_");
  }  

/* Initialize biogeochem model */
  myTime = DeltaT*Iter; /* Iter should start at 0 */  

  return 0;
}
