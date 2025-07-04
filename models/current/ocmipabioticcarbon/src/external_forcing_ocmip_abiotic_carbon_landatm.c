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
#include "ocmip_abiotic_carbon_landatm_tmm.h"

/* Macros to map tracer names to vectors */
/* v[0]=DIC */
/* v[1]=DIC14 */
#define DIC v[0]
#define JDIC ut[0]
#ifdef ALLOW_C14
#define DIC14 v[1]
#define JDIC14 ut[1]
#endif

Vec surfVolFrac;
PetscScalar *localDIC;
PetscScalar *localJDIC;
#ifdef ALLOW_C14
PetscScalar *localDIC14;
PetscScalar *localJDIC14;
#endif
PetscScalar *localTs,*localSs;
PetscScalar *localAlk, *localfice, *localEmP;
PetscScalar *localPO4,*localSiO2;
PetscScalar *localxkw, *localatmosp;
PetscScalar *localVgas;
PetscScalar DeltaT, *localdzsurf;
PetscBool useVirtualFlux = PETSC_FALSE;
PetscScalar *pH;
PetscBool useWinds = PETSC_FALSE;
PeriodicTimer windsTimer;
PetscScalar *localuwind,*localvwind;
PeriodicArray localuwindp, localvwindp;
PetscScalar pistonVelocityCoeff=0.31; /* default piston velocity coefficient when using winds [cm/hr]*[s^2/m^2] */

PetscBool useLinearChemistry = PETSC_FALSE;
PetscScalar *linearChemistryFactor, *linearChemistryCO2, *linearChemistryDIC;

PeriodicArray localTsp, localSsp;
PeriodicArray localAlkp, localficep, localEmPp;
PeriodicArray localPO4p, localSiO2p;
PeriodicArray localxkwp, localatmospp;
TimeDependentArray localTstd, localSstd;
TimeDependentArray localAlktd, localficetd, localEmPtd;

PetscBool periodicBiogeochemForcing = PETSC_FALSE;
PetscBool timeDependentBiogeochemForcing = PETSC_FALSE;
PeriodicTimer biogeochemTimer;
TimeDependentTimer timeDependentBiogeochemTimer;

PetscInt maxValsToRead;
PetscInt dummyInt = 0;
char *pCO2atmFiles[2];  
PetscInt numpCO2atm_hist = 0;
PetscScalar *TpCO2atm_hist, *pCO2atm_hist;
PetscBool fixedAtmosCO2 = PETSC_TRUE;
char pCO2atmIniFile[PETSC_MAX_PATH_LEN];  
PetscScalar *localpCO2atm, *localpCO2atm0, *localpCO2atm1;
PetscBool spatiallyVariableCO2 = PETSC_FALSE;
PetscInt itfCO2 = -1;

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
char atmOutTimeFile[PETSC_MAX_PATH_LEN];  
PetscScalar pCO2atmavg, Foceanavg, Flandavg, landStateavg[3];

#ifdef ALLOW_C14
/* PetscBool useC14 = PETSC_FALSE; */
PetscScalar lambdaDecayC14 = 1.2097e-4; /* y^-1 we convert units below */
PetscScalar DC14atm = 0.0; /* default initial value */
PetscScalar *localDC14atm, *localDC14atm0, *localDC14atm1;
PetscBool fixedAtmosC14 = PETSC_TRUE;
char *C14atmFiles[2];
PetscScalar *TC14atm_hist;
PetscInt numC14atm_hist = 0;
PetscInt itfC14 = -1;
#endif

PetscBool calcDiagnostics = PETSC_FALSE;
StepTimer diagTimer;
PetscBool diagAppendOutput = PETSC_FALSE;
PetscFileMode DIAG_FILE_MODE;
FILE *diagfptime;
char diagOutTimeFile[PETSC_MAX_PATH_LEN];
PetscScalar *localpco2diag, *localpco2diagavg;
PetscScalar *localgasexfluxdiag, *localgasexfluxdiagavg, *localtotfluxdiag, *localtotfluxdiagavg;
PetscScalar *localc14gasexfluxdiag, *localc14gasexfluxdiagavg, *localc14totfluxdiag, *localc14totfluxdiagavg;

/* Vec JDICavg; */
/* PetscViewer fdjdic; */
#ifdef ALLOW_C14
Vec Rocn, Rocnavg;
PetscViewer fdrocn;
#endif

#undef __FUNCT__
#define __FUNCT__ "iniExternalForcing"
PetscErrorCode iniExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt numTracers, Vec *v, Vec *ut)
{
  PetscErrorCode ierr;
  PetscInt ip, kl;
  PetscViewer fd;
  int fp;
  PetscBool flg, flg1;
  PetscInt it;
  PetscScalar myTime;
  PetscScalar zero = 0.0;

#ifdef ALLOW_C14
  if (numTracers!=2) SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of tracers specified! Must specify 2 tracers for C14");
#else
  if (numTracers==2) SETERRQ(PETSC_COMM_WORLD,1,"Too many tracers selected!");
#endif

  ierr = VecGetArray(DIC,&localDIC);CHKERRQ(ierr);

  ierr = VecSet(JDIC,zero); CHKERRQ(ierr);    
  ierr = VecGetArray(JDIC,&localJDIC);CHKERRQ(ierr);

  ierr = PetscOptionsGetReal(NULL,NULL,"-biogeochem_deltat",&DeltaT,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate biogeochemical time step in seconds with the -biogeochem_deltat option");  

  ierr = PetscOptionsGetReal(NULL,NULL,"-days_per_year",&daysPerYear,&flg);CHKERRQ(ierr);
  if (!flg) {
    daysPerYear = 360.0;
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of days per year is %12.7f\n",daysPerYear);CHKERRQ(ierr);
  secondsPerYear = 86400.0*daysPerYear;

#ifdef ALLOW_C14
  ierr = VecGetArray(DIC14,&localDIC14);CHKERRQ(ierr);

  ierr = VecSet(JDIC14,zero); CHKERRQ(ierr);    
  ierr = VecGetArray(JDIC14,&localJDIC14);CHKERRQ(ierr);

  ierr=PetscPrintf(PETSC_COMM_WORLD,"C14 will also be simulated\n");CHKERRQ(ierr);

  lambdaDecayC14 = lambdaDecayC14/secondsPerYear; /* convert y^-1 to s^-1 */
  
#endif

  ierr = PetscOptionsHasName(NULL,NULL,"-periodic_biogeochem_forcing",&periodicBiogeochemForcing);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-time_dependent_biogeochem_forcing",&timeDependentBiogeochemForcing);CHKERRQ(ierr);

  if (periodicBiogeochemForcing) {    
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic biogeochemical forcing specified\n");CHKERRQ(ierr);
    ierr = iniPeriodicTimer("periodic_biogeochem_", &biogeochemTimer);CHKERRQ(ierr);
  }

  if (timeDependentBiogeochemForcing) {    
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Time-dependent biogeochemical forcing specified\n");CHKERRQ(ierr);
    ierr = iniTimeDependentTimer("time_dependent_biogeochem_", &timeDependentBiogeochemTimer);CHKERRQ(ierr);
  }

/* Land/Atm model data */
  ierr = PetscOptionsHasName(NULL,NULL,"-use_atm_model",&useAtmModel);CHKERRQ(ierr);
  ierr = PetscOptionsHasName(NULL,NULL,"-use_land_model",&useLandModel);CHKERRQ(ierr);

  if ((useLandModel) && (!useAtmModel)) SETERRQ(PETSC_COMM_WORLD,1,"ERROR: Land model cannot be used without the atmospheric model");

#ifdef ALLOW_C14
  if (useAtmModel) SETERRQ(PETSC_COMM_WORLD,1,"ERROR: C14 is not supported with atmospheric model!");
#endif

  if (useAtmModel) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using interactive atmospheric model\n");CHKERRQ(ierr);  

    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localdA);CHKERRQ(ierr);
    ierr = readProfileSurfaceScalarData("dA.bin",localdA,1);  

/* overwrite default value */
	ierr = PetscOptionsGetReal(NULL,NULL,"-pco2atm_ini",&pCO2atm_ini,&flg);CHKERRQ(ierr); /* read from command line */
    if (!flg) {
      ierr = PetscOptionsGetString(NULL,NULL,"-pco2atm_ini_file",pCO2atmIniFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
      if (flg) { /* read from binary file */
        ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,pCO2atmIniFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
        ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
        ierr = PetscBinaryRead(fp,&pCO2atm_ini,1,NULL,PETSC_SCALAR);CHKERRQ(ierr);
        ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      }
    }
    pCO2atm = pCO2atm_ini;
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using initial atmospheric pCO2 of %g ppm\n",pCO2atm);CHKERRQ(ierr);
      
    atmModelDeltaT = DeltaT/secondsPerYear; /* time step in years */

	ierr = iniStepTimer("atm_write_", Iter0, &atmWriteTimer);CHKERRQ(ierr);

    ierr = PetscOptionsHasName(NULL,NULL,"-atm_append",&atmAppendOutput);CHKERRQ(ierr);
    if (atmAppendOutput) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output will be appended\n");CHKERRQ(ierr);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output will overwrite existing file(s)\n");CHKERRQ(ierr);
    }    

/* Output times */
    ierr = PetscOptionsGetString(NULL,NULL,"-atm_time_file",atmOutTimeFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
    if (!flg) {
      strcpy(atmOutTimeFile,"");
      sprintf(atmOutTimeFile,"%s","atm_output_time.txt");
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output times will be written to %s\n",atmOutTimeFile);CHKERRQ(ierr);

    if (!atmAppendOutput) {
      ierr = PetscFOpen(PETSC_COMM_WORLD,atmOutTimeFile,"w",&atmfptime);CHKERRQ(ierr);  
	  if (Iter0==atmWriteTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
		ierr = PetscFPrintf(PETSC_COMM_WORLD,atmfptime,"%d   %10.5f\n",Iter0,time0);CHKERRQ(ierr);     
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing atmospheric output at time %10.5f, step %d\n", tc,Iter);CHKERRQ(ierr);  
		ierr = writeBinaryScalarData("pCO2atm_output.bin",&pCO2atm,1,PETSC_FALSE);
	  }	
    } else {
      ierr = PetscFOpen(PETSC_COMM_WORLD,atmOutTimeFile,"a",&atmfptime);CHKERRQ(ierr);  
	  if (Iter0==atmWriteTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */            
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output will be appended. Initial condition will NOT be written\n");CHKERRQ(ierr);      
	  }	
    }

    if (useLandModel) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Using interactive land model\n");CHKERRQ(ierr);      
      ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"land_ini.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
      ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
      ierr = PetscBinaryRead(fp,landState,3,NULL,PETSC_SCALAR);CHKERRQ(ierr);
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
    ierr = PetscOptionsGetStringArray(NULL,NULL,"-emissions_history",emFiles,&maxValsToRead,&useEmissions);CHKERRQ(ierr);
    if (useEmissions) { /* Read emissions history */
      if (maxValsToRead != 3) {
        SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of file names specified for emissions");
      }      
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Using prescribed emissions\n");CHKERRQ(ierr);     
      /* read time data */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,emFiles[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	  ierr = PetscBinaryRead(fp,&numEmission_hist,1,NULL,PETSC_INT);CHKERRQ(ierr);  
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of points in emission files is %d \n",numEmission_hist);CHKERRQ(ierr);  
      ierr = PetscMalloc(numEmission_hist*sizeof(PetscScalar),&Tem_hist);CHKERRQ(ierr); 
      ierr = PetscBinaryRead(fp,Tem_hist,numEmission_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      /* read fossil fuel emissions */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,emFiles[1],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
      ierr = PetscMalloc(numEmission_hist*sizeof(PetscScalar),&E_hist);CHKERRQ(ierr); 
      ierr = PetscBinaryRead(fp,E_hist,numEmission_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      /* read land use emissions */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,emFiles[2],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
      ierr = PetscMalloc(numEmission_hist*sizeof(PetscScalar),&D_hist);CHKERRQ(ierr); 
      ierr = PetscBinaryRead(fp,D_hist,numEmission_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      ierr = PetscOptionsHasName(NULL,NULL,"-interp_emissions",&interpEmissions);CHKERRQ(ierr);
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
	ierr = PetscOptionsHasName(NULL,NULL,"-spatially_variable_atmospheric_co2",&spatiallyVariableCO2);CHKERRQ(ierr);
    
    if (spatiallyVariableCO2) {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Spatially-variable atmospheric pCO2 has been specified\n");CHKERRQ(ierr);          
	  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localpCO2atm);CHKERRQ(ierr);    
    }
    
	pCO2atmFiles[0] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* time file */
	pCO2atmFiles[1] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* atmospheric pCO2 history file */
  	maxValsToRead = 2;
    ierr = PetscOptionsGetStringArray(NULL,NULL,"-pco2atm_history",pCO2atmFiles,&maxValsToRead,&flg);CHKERRQ(ierr);
    if (flg) { /* Read atmospheric pCO2 history */
      if (maxValsToRead != 2) {
        SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of file names specified for atmospheric pCO2 history");
      }      
      fixedAtmosCO2 = PETSC_FALSE;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading time-dependent atmospheric pCO2 history\n");CHKERRQ(ierr);      
      /* read time data */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,pCO2atmFiles[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	  ierr = PetscBinaryRead(fp,&numpCO2atm_hist,1,NULL,PETSC_INT);CHKERRQ(ierr);  
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of points in atmospheric history file is %d \n",numpCO2atm_hist);CHKERRQ(ierr);  
      ierr = PetscMalloc(numpCO2atm_hist*sizeof(PetscScalar),&TpCO2atm_hist);CHKERRQ(ierr); 
      ierr = PetscBinaryRead(fp,TpCO2atm_hist,numpCO2atm_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      /* read atmospheric pCO2 data */
      if (spatiallyVariableCO2) {
		ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localpCO2atm0);CHKERRQ(ierr);  
		ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localpCO2atm1);CHKERRQ(ierr);  
		ierr = readProfileSurfaceScalarData(pCO2atmFiles[1],localpCO2atm,1); /* read initial slice */ 
		itfCO2=-1;      
      } else {
		ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,pCO2atmFiles[1],FILE_MODE_READ,&fd);CHKERRQ(ierr);
		ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
		ierr = PetscMalloc(numpCO2atm_hist*sizeof(PetscScalar),&pCO2atm_hist);CHKERRQ(ierr); 
		ierr = PetscBinaryRead(fp,pCO2atm_hist,numpCO2atm_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
		ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	  
		pCO2atm = pCO2atm_hist[0];
      }
    } else {
      if (spatiallyVariableCO2) {
		ierr = PetscOptionsGetString(NULL,NULL,"-pco2atm",pCO2atmFiles[1],PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
		if (flg1) {
  		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading spatially-variable fixed atmospheric pCO2\n");CHKERRQ(ierr);		
		  ierr = readProfileSurfaceScalarData(pCO2atmFiles[1],localpCO2atm,1);  
		} else {
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: no atmospheric pCO2 file given! Using a fixed, uniform value of %g ppm\n",pCO2atm);CHKERRQ(ierr);
		  for (ip=0; ip<lNumProfiles; ip++) {
			localpCO2atm[ip]=pCO2atm;
		  }
		}      
      } else {
		ierr = PetscOptionsGetReal(NULL,NULL,"-pco2atm",&pCO2atm,&flg);CHKERRQ(ierr); /* overwrite default value */
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Using fixed atmospheric pCO2 of %g ppm\n",pCO2atm);CHKERRQ(ierr);
      }    
    }    

#ifdef ALLOW_C14
    /* prescribed atmospheric DeltaC14 */
	ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localDC14atm);CHKERRQ(ierr);
	C14atmFiles[0] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* time file */
	C14atmFiles[1] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* atmospheric DeltaC14 history file */
  	maxValsToRead = 2;
    ierr = PetscOptionsGetStringArray(NULL,NULL,"-c14atm_history",C14atmFiles,&maxValsToRead,&flg);CHKERRQ(ierr);
    if (flg) { /* Read atmospheric C14 history */
      if (maxValsToRead != 2) {
        SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of file names specified for atmospheric DeltaC14 history");
      }      
      fixedAtmosC14 = PETSC_FALSE;
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading time-dependent atmospheric DeltaC14 history\n");CHKERRQ(ierr);
      /* read time data */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,C14atmFiles[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	  ierr = PetscBinaryRead(fp,&numC14atm_hist,1,NULL,PETSC_INT);CHKERRQ(ierr);  
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of points in atmospheric history file is %d \n",numC14atm_hist);CHKERRQ(ierr);  
      ierr = PetscMalloc(numC14atm_hist*sizeof(PetscScalar),&TC14atm_hist);CHKERRQ(ierr); 
      ierr = PetscBinaryRead(fp,TC14atm_hist,numC14atm_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
      /* read atmospheric DeltaC14 data */
      ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localDC14atm0);CHKERRQ(ierr);  
      ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localDC14atm1);CHKERRQ(ierr);  
      ierr = readProfileSurfaceScalarData(C14atmFiles[1],localDC14atm,1); /* read initial slice */ 
      itfC14=-1;

    }	else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Using fixed atmospheric DeltaC14\n");CHKERRQ(ierr);
	  ierr = PetscOptionsGetString(NULL,NULL,"-c14atm",C14atmFiles[1],PETSC_MAX_PATH_LEN-1,&flg1);CHKERRQ(ierr);
	  if (flg1) {
	    ierr = readProfileSurfaceScalarData(C14atmFiles[1],localDC14atm,1);  
	  } else {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: no atmospheric DeltaC14 file given! Using a fixed, uniform value of %g per mil\n",DC14atm);CHKERRQ(ierr);
		for (ip=0; ip<lNumProfiles; ip++) {
		  localDC14atm[ip]=DC14atm;
		}
	  }
    }    
#endif    
  
  }

/* Grid arrays */
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localdzsurf);CHKERRQ(ierr);
  ierr = readProfileSurfaceScalarData("dzsurf.bin",localdzsurf,1);  

/* Forcing fields */  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&pH);CHKERRQ(ierr);

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localVgas);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) { /* read winds or transfer coefficient (xkw) */
    ierr = PetscOptionsHasName(NULL,NULL,"-use_winds",&useWinds);CHKERRQ(ierr);
    if (useWinds) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Periodic winds specified: gas transfer velocity will be computed using winds\n");CHKERRQ(ierr);      
      ierr = PetscOptionsGetReal(NULL,NULL,"-piston_velocity_coeff",&pistonVelocityCoeff,&flg);CHKERRQ(ierr); /* overwrite default value */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Piston velocity coefficient of %10.5f [cm/hr]*[s^2/m^2] will be used\n", pistonVelocityCoeff);CHKERRQ(ierr);      
      pistonVelocityCoeff = pistonVelocityCoeff*0.01/3600.0; /* convert [cm/hr]*[s^2/m^2] to [m/s]*[s^2/m^2] */
      ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localuwind);CHKERRQ(ierr);    
      localuwindp.firstTime = PETSC_TRUE;
      localuwindp.arrayLength = lNumProfiles;    
      ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localvwind);CHKERRQ(ierr);    
      localvwindp.firstTime = PETSC_TRUE;
      localvwindp.arrayLength = lNumProfiles;    

      ierr = PetscOptionsHasName(NULL,NULL,"-periodic_winds_cycle_period",&flg);CHKERRQ(ierr);
      if (flg) {
	    ierr = iniPeriodicTimer("periodic_winds_", &windsTimer);CHKERRQ(ierr);
      } else {
	    ierr = iniPeriodicTimer("periodic_biogeochem_", &windsTimer);CHKERRQ(ierr);
      }
      
    } else {
      ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localxkw);CHKERRQ(ierr);    
      localxkwp.firstTime = PETSC_TRUE;
      localxkwp.arrayLength = lNumProfiles;
    }   
  } else { /* read annual mean Vgas */
    ierr = readProfileSurfaceScalarData("Vgas.bin",localVgas,1);  
  }  

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localfice);CHKERRQ(ierr); /* NOTE: fice is only use for periodic or time-dependent biogeochem forcing */ 
  if (timeDependentBiogeochemForcing) {
    localficetd.firstTime = PETSC_TRUE;
    localficetd.arrayLength = lNumProfiles;
  } else if (periodicBiogeochemForcing) {    
    localficep.firstTime = PETSC_TRUE;
    localficep.arrayLength = lNumProfiles;
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localTs);CHKERRQ(ierr);  
  if (timeDependentBiogeochemForcing) {
    localTstd.firstTime = PETSC_TRUE;
    localTstd.arrayLength = lNumProfiles;  
  } else if (periodicBiogeochemForcing) {    
    localTsp.firstTime = PETSC_TRUE;
    localTsp.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("Ts.bin",localTs,1);  
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localSs);CHKERRQ(ierr);
  if (timeDependentBiogeochemForcing) {
    localSstd.firstTime = PETSC_TRUE;
    localSstd.arrayLength = lNumProfiles;  
  } else if (periodicBiogeochemForcing) {    
    localSsp.firstTime = PETSC_TRUE;
    localSsp.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("Ss.bin",localSs,1);  
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localAlk);CHKERRQ(ierr);
  if (timeDependentBiogeochemForcing) {
    localAlktd.firstTime = PETSC_TRUE;
    localAlktd.arrayLength = lNumProfiles;
  } else if (periodicBiogeochemForcing) {    
    localAlkp.firstTime = PETSC_TRUE;
    localAlkp.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("Alk.bin",localAlk,1);  
  }

  ierr = PetscOptionsHasName(NULL,NULL,"-use_virtual_flux",&useVirtualFlux);CHKERRQ(ierr);

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localEmP);CHKERRQ(ierr);
  if (timeDependentBiogeochemForcing) {
	localEmPtd.firstTime = PETSC_TRUE;
	localEmPtd.arrayLength = lNumProfiles;
  } else if (periodicBiogeochemForcing) {    
	localEmPp.firstTime = PETSC_TRUE;
	localEmPp.arrayLength = lNumProfiles;
  } else {  
	ierr = readProfileSurfaceScalarData("EmP.bin",localEmP,1);  
  }

  if (useVirtualFlux) {
	ierr = VecDuplicate(DIC,&surfVolFrac);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"surface_volume_fraction.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(surfVolFrac,fd);CHKERRQ(ierr);  
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      
  }
  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localatmosp);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localatmospp.firstTime = PETSC_TRUE;
    localatmospp.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("atmosp.bin",localatmosp,1);  
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localPO4);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localPO4p.firstTime = PETSC_TRUE;
    localPO4p.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("PO4.bin",localPO4,1);  
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localSiO2);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localSiO2p.firstTime = PETSC_TRUE;
    localSiO2p.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("SiO2.bin",localSiO2,1);  
  }

  ierr = PetscOptionsHasName(NULL,NULL,"-use_linear_chemistry",&useLinearChemistry);CHKERRQ(ierr);
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&linearChemistryFactor);CHKERRQ(ierr);  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&linearChemistryCO2);CHKERRQ(ierr);  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&linearChemistryDIC);CHKERRQ(ierr);  
  if (useLinearChemistry) {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Linear chemistry will be used when computing CO2aq from DIC\n");CHKERRQ(ierr);              
    ierr = readProfileSurfaceScalarData("linear_chemistry_factor.bin",linearChemistryFactor,1);  
    ierr = readProfileSurfaceScalarData("linear_chemistry_co2.bin",linearChemistryCO2,1);  
    ierr = readProfileSurfaceScalarData("linear_chemistry_dic.bin",linearChemistryDIC,1);  
  } else {
	for (ip=0; ip<lNumProfiles; ip++) {
	  linearChemistryFactor[ip]=-999.0;
	  linearChemistryCO2[ip]=-999.0;
	  linearChemistryDIC[ip]=-999.0;	  
    }
  }

  if (timeDependentBiogeochemForcing) {
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localTs,timeDependentBiogeochemTimer.numTimes,timeDependentBiogeochemTimer.tdt,&localTstd,"Ts.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localSs,timeDependentBiogeochemTimer.numTimes,timeDependentBiogeochemTimer.tdt,&localSstd,"Ss.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localEmP,timeDependentBiogeochemTimer.numTimes,timeDependentBiogeochemTimer.tdt,&localEmPtd,"EmP.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localfice,timeDependentBiogeochemTimer.numTimes,timeDependentBiogeochemTimer.tdt,&localficetd,"fice.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localAlk,timeDependentBiogeochemTimer.numTimes,timeDependentBiogeochemTimer.tdt,&localAlktd,"Alk.bin");    
  } else if (periodicBiogeochemForcing) {    
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localTs,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localTsp,"Ts_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localSs,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localSsp,"Ss_");
	ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localEmPp,"EmP_");                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localficep,"fice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localAlk,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localAlkp,"Alk_");    
  }    
  if (periodicBiogeochemForcing) {   
    if (useWinds) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localuwind,windsTimer.cyclePeriod,windsTimer.numPerPeriod,windsTimer.tdp,&localuwindp,"uwind_");                                                    
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localvwind,windsTimer.cyclePeriod,windsTimer.numPerPeriod,windsTimer.tdp,&localvwindp,"vwind_");                                                        
    } else {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localxkw,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localxkwp,"xkw_");
    }
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localatmospp,"atmosp_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localPO4,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localPO4p,"PO4_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localSiO2,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localSiO2p,"SiO2_");
  }  

/* Initialize biogeochem model */
  myTime = DeltaT*Iter; /* Iter should start at 0 */  
  for (ip=0; ip<lNumProfiles; ip++) {
	kl=lStartIndices[ip];  

    ocmip_abiotic_carbon_ini_(&Iter,&myTime,&localDIC[kl],&localAlk[ip],&localPO4[ip],&localSiO2[ip],                 
                           &localTs[ip],&localSs[ip],&pH[ip]);
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localpco2diag);CHKERRQ(ierr);  /* always need to pass this */

  ierr = PetscOptionsHasName(NULL,NULL,"-calc_diagnostics",&calcDiagnostics);CHKERRQ(ierr);
  if (calcDiagnostics) {    
/*Data for diagnostics */
    ierr = iniStepTimer("diag_", Iter0+1, &diagTimer);CHKERRQ(ierr);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed starting at and including (absolute) time step: %d\n", diagTimer.startTimeStep);CHKERRQ(ierr);	
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed over %d time steps\n", diagTimer.numTimeSteps);CHKERRQ(ierr);	

	ierr = PetscOptionsHasName(NULL,NULL,"-diag_append",&diagAppendOutput);CHKERRQ(ierr);
	if (diagAppendOutput) {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostic output will be appended\n");CHKERRQ(ierr);
	  DIAG_FILE_MODE=FILE_MODE_APPEND;
	} else {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostic output will overwrite existing file(s)\n");CHKERRQ(ierr);
	  DIAG_FILE_MODE=FILE_MODE_WRITE;
	}

/* Output times */
	ierr = PetscOptionsGetString(NULL,NULL,"-diag_time_file",diagOutTimeFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
	if (!flg) {
	  strcpy(diagOutTimeFile,"");
	  sprintf(diagOutTimeFile,"%s","diagnostic_output_time.txt");
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostic output times will be written to %s\n",diagOutTimeFile);CHKERRQ(ierr);

	if (!diagAppendOutput) {
	  ierr = PetscFOpen(PETSC_COMM_WORLD,diagOutTimeFile,"w",&diagfptime);CHKERRQ(ierr);  
	} else {
	  ierr = PetscFOpen(PETSC_COMM_WORLD,diagOutTimeFile,"a",&diagfptime);CHKERRQ(ierr);  
	}

    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localpco2diagavg);CHKERRQ(ierr);  

    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localgasexfluxdiag);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localgasexfluxdiagavg);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localtotfluxdiag);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localtotfluxdiagavg);CHKERRQ(ierr);  

#ifdef ALLOW_C14
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localc14gasexfluxdiag);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localc14gasexfluxdiagavg);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localc14totfluxdiag);CHKERRQ(ierr);  
    ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localc14totfluxdiagavg);CHKERRQ(ierr);
#endif

    for (ip=0; ip<lNumProfiles; ip++) {
      localpco2diag[ip]=0.0;
      localpco2diagavg[ip]=0.0;      
      localgasexfluxdiag[ip]=0.0;
      localgasexfluxdiagavg[ip]=0.0;      
      localtotfluxdiag[ip]=0.0;
      localtotfluxdiagavg[ip]=0.0;      
#ifdef ALLOW_C14
      localc14gasexfluxdiag[ip]=0.0;
      localc14gasexfluxdiagavg[ip]=0.0;      
      localc14totfluxdiag[ip]=0.0;
      localc14totfluxdiagavg[ip]=0.0;      
#endif
    }

#ifdef ALLOW_C14
    ierr = VecDuplicate(DIC,&Rocn);CHKERRQ(ierr);
    ierr = VecDuplicate(DIC,&Rocnavg);CHKERRQ(ierr);
    ierr = VecSet(Rocnavg,zero); CHKERRQ(ierr);        
    /* Open file here for subsequent output */
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Rocn.petsc",DIAG_FILE_MODE,&fdrocn);CHKERRQ(ierr);    
#endif      

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
  PetscScalar DICemp = 0.0, DIC14emp = 0.0;
  PetscScalar Vgas660, Sc;
  PetscInt itr, ip, nzloc, kl;
  PetscScalar myTime;
  PetscScalar alpha;
  PetscInt itf;
  PetscScalar zero = 0.0, one = 1.0;
  PetscInt k;
  PetscScalar localFocean;
  PetscScalar localgasexflux = 0.0, localtotflux = 0.0, localc14gasexflux = 0.0, localc14totflux = 0.0;
  PetscScalar w2;
  
  myTime = DeltaT*Iter; /* Iter should start at 0 */

  if (timeDependentBiogeochemForcing) {
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localTs,timeDependentBiogeochemTimer.numTimes,timeDependentBiogeochemTimer.tdt,&localTstd,"Ts.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localSs,timeDependentBiogeochemTimer.numTimes,timeDependentBiogeochemTimer.tdt,&localSstd,"Ss.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localEmP,timeDependentBiogeochemTimer.numTimes,timeDependentBiogeochemTimer.tdt,&localEmPtd,"EmP.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localfice,timeDependentBiogeochemTimer.numTimes,timeDependentBiogeochemTimer.tdt,&localficetd,"fice.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localAlk,timeDependentBiogeochemTimer.numTimes,timeDependentBiogeochemTimer.tdt,&localAlktd,"Alk.bin");    
  } else if (periodicBiogeochemForcing) {    
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localTs,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localTsp,"Ts_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localSs,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localSsp,"Ss_");
	ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localEmPp,"EmP_");                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localficep,"fice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localAlk,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localAlkp,"Alk_");    
  }    
  if (periodicBiogeochemForcing) {   
    if (useWinds) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localuwind,windsTimer.cyclePeriod,windsTimer.numPerPeriod,windsTimer.tdp,&localuwindp,"uwind_");                                                    
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localvwind,windsTimer.cyclePeriod,windsTimer.numPerPeriod,windsTimer.tdp,&localvwindp,"vwind_");                                                        
    } else {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localxkw,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localxkwp,"xkw_");
    }
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localatmospp,"atmosp_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localPO4,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localPO4p,"PO4_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localSiO2,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localSiO2p,"SiO2_");
  }  

  if (periodicBiogeochemForcing) {   
/*  Recompute gas exchange coeff */
    if (useWinds) {
      for (ip=0; ip<lNumProfiles; ip++) {
        kl=lStartIndices[ip];
        w2 = pow(localuwind[ip],2) + pow(localvwind[ip],2);
        Vgas660 = (1.0-localfice[ip])*pistonVelocityCoeff*w2;
        Sc = 2073.1 - 125.62*localTs[ip] + 3.6276*pow(localTs[ip],2) - 0.043219*pow(localTs[ip],3);
        localVgas[ip]=Vgas660/sqrt(Sc/660.0);
      }                                                      
    } else {
      for (ip=0; ip<lNumProfiles; ip++) {
        kl=lStartIndices[ip];
    
        Vgas660 = (1.0-localfice[ip])*localxkw[ip];
        Sc = 2073.1 - 125.62*localTs[ip] + 3.6276*pow(localTs[ip],2) - 0.043219*pow(localTs[ip],3);
        localVgas[ip]=Vgas660/sqrt(Sc/660.0);
      }
    }
  }

  if (useAtmModel) {
    if (useEmissions) {  
/*   Interpolate emissions */
      if (tc>=Tem_hist[0]) {
        if (interpEmissions) {
		  ierr = calcInterpFactor(numEmission_hist,tc,Tem_hist,&itf,&alpha); CHKERRQ(ierr);
		  fossilFuelEmission = alpha*E_hist[itf] + (1.0-alpha)*E_hist[itf+1];	  
		  landUseEmission = alpha*D_hist[itf] + (1.0-alpha)*D_hist[itf+1];	          
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
      if (spatiallyVariableCO2) {
		if (tc>=TpCO2atm_hist[0]) {
		  ierr = calcInterpFactor(numpCO2atm_hist,tc,TpCO2atm_hist,&itf,&alpha); CHKERRQ(ierr);
		  if (itf != itfCO2) { /* time to read new bracketing slices: itf uses 0-based, while readProfileSurfaceScalarDataRecord uses 1-based indexing*/
			ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading new bracketing slices for pCO2atm at time = %g: %d and %d\n",tc,itf+1,itf+2);CHKERRQ(ierr);        
			ierr = readProfileSurfaceScalarDataRecord(pCO2atmFiles[1],localpCO2atm0,1,itf+1);
			ierr = readProfileSurfaceScalarDataRecord(pCO2atmFiles[1],localpCO2atm1,1,itf+2);
			itfCO2=itf;
		  }
/*        interpolate in time */
		  for (ip=0; ip<lNumProfiles; ip++) {
			localpCO2atm[ip] = alpha*localpCO2atm0[ip] + (1.0-alpha)*localpCO2atm1[ip];	          
		  }        
		} else {
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: time < %10.5f. Assuming pCO2atm fixed at initial value\n",TpCO2atm_hist[0]);CHKERRQ(ierr);
		}      
      } else {
		if (tc>=TpCO2atm_hist[0]) {
		  ierr = calcInterpFactor(numpCO2atm_hist,tc,TpCO2atm_hist,&itf,&alpha); CHKERRQ(ierr);
		  pCO2atm = alpha*pCO2atm_hist[itf] + (1.0-alpha)*pCO2atm_hist[itf+1];	  
		} else {
		  ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: time < %10.5f. Assuming pCO2atm=%g\n",TpCO2atm_hist[0],pCO2atm);CHKERRQ(ierr);
		}
	  }	
    }  
    
#ifdef ALLOW_C14
    if (!fixedAtmosC14) { 
      if (tc>=TC14atm_hist[0]) {
        ierr = calcInterpFactor(numC14atm_hist,tc,TC14atm_hist,&itf,&alpha); CHKERRQ(ierr);
        if (itf != itfC14) { /* time to read new bracketing slices: itf uses 0-based, while readProfileSurfaceScalarDataRecord uses 1-based indexing*/
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading new bracketing slices for DeltaC14atm at time = %g: %d and %d\n",tc,itf+1,itf+2);CHKERRQ(ierr);        
          ierr = readProfileSurfaceScalarDataRecord(C14atmFiles[1],localDC14atm0,1,itf+1);
          ierr = readProfileSurfaceScalarDataRecord(C14atmFiles[1],localDC14atm1,1,itf+2);
          itfC14=itf;
        }
/*         interpolate in time */
        for (ip=0; ip<lNumProfiles; ip++) {
          localDC14atm[ip] = alpha*localDC14atm0[ip] + (1.0-alpha)*localDC14atm1[ip];	          
        }        
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: time < %10.5f. Assuming DeltaC14atm fixed at initial value\n",TC14atm_hist[0]);CHKERRQ(ierr);
      }
    }  
#endif
  }

#ifdef ALLOW_C14
  ierr = VecSet(JDIC14,zero); CHKERRQ(ierr); /* This is needed for C14 because we are adding an interior source (decay) term */
#endif

  if (useVirtualFlux) { /* use the global surface mean value to calculate E-P contribution */
    ierr = VecDot(surfVolFrac,DIC,&DICemp);CHKERRQ(ierr); /* volume weighted mean surface DIC */									              
#ifdef ALLOW_C14
    ierr = VecDot(surfVolFrac,DIC14,&DIC14emp);CHKERRQ(ierr); /* volume weighted mean surface DIC14 */									              
#endif    
  }

  localFocean = 0.0;  
  Focean = 0.0;
  Fland = 0.0;
  
/* Compute air-sea gas exchange term */
  for (ip=0; ip<lNumProfiles; ip++) {
    kl=lStartIndices[ip];
    nzloc=lProfileLength[ip];
    
    if (spatiallyVariableCO2) pCO2atm = localpCO2atm[ip];
    
	if (!useVirtualFlux) { /* use the local surface value to calculate E-P contribution */
      DICemp=localDIC[kl];
#ifdef ALLOW_C14      
      DIC14emp=localDIC14[kl];
#endif
    }
    ocmip_abiotic_carbon_model_(&Iter,&myTime,
                         &localDIC[kl],
#ifdef ALLOW_C14
                         &localDIC14[kl],
#endif                         
                         &localAlk[ip],&localPO4[ip],&localSiO2[ip],
                         &localTs[ip],&localSs[ip],&pH[ip],&localVgas[ip],
                         &localatmosp[ip],&pCO2atm,&localdzsurf[ip],
                         &localEmP[ip],&DICemp,
                         &linearChemistryFactor[ip],&linearChemistryCO2[ip],&linearChemistryDIC[ip],
#ifdef ALLOW_C14
                         &localDC14atm[ip],&DIC14emp,
#endif                                                  
                         &localJDIC[kl],&localgasexflux,&localtotflux,&localpco2diag[ip]
#ifdef ALLOW_C14
                         ,&localJDIC14[kl],&localc14gasexflux,&localc14totflux
#endif                                                                           
                         );

    if (useAtmModel) {                 
      localFocean = localFocean + localtotflux*localdA[ip]*(12.0/1.e15)*secondsPerYear; /* PgC/y */
    }
    
	if (calcDiagnostics) {  
	  if (Iter0+iLoop>=diagTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */
        localgasexfluxdiag[ip]=localgasexflux;
        localtotfluxdiag[ip]=localtotflux;
#ifdef ALLOW_C14
        localc14gasexfluxdiag[ip]=localc14gasexflux;
        localc14totfluxdiag[ip]=localc14totflux;
#endif
      }
	}                         
  } /* end loop over profiles */
  
  if (useAtmModel) {
	MPI_Allreduce(&localFocean, &Focean, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);    

	if (useLandModel) {
	  landsource_(&landState[0],&pCO2atm,&landUseEmission,&deltaTsg,&Fland,&landSource[0]); /* returns S and Fl in PgC/y */
/*    time step land */
	  for (k=0;k<=2;k++) {
		landState[k] = landState[k] + atmModelDeltaT*landSource[k];
	  }
	}

/* time step atmosphere */
    pCO2atm = pCO2atm + atmModelDeltaT*(fossilFuelEmission + landUseEmission - Focean - Fland)/ppmToPgC;
  }  

  ierr = VecSetValues(JDIC,lSize,gIndices,localJDIC,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(JDIC);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(JDIC);CHKERRQ(ierr);    
#ifdef ALLOW_C14
  ierr = VecSetValues(JDIC14,lSize,gIndices,localJDIC14,INSERT_VALUES);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(JDIC14);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(JDIC14);CHKERRQ(ierr);    
#endif

/* Radioactive decay term */
#ifdef ALLOW_C14
    ierr = VecAXPY(JDIC14,-1.0*lambdaDecayC14,DIC14);CHKERRQ(ierr); /* JDIC14 <- JDIC14 - lambda*DIC14 */  
#endif

/* Convert to discrete tendency */
  ierr = VecScale(JDIC,DeltaT);CHKERRQ(ierr);
#ifdef ALLOW_C14
  ierr = VecScale(JDIC14,DeltaT);CHKERRQ(ierr);
#endif
  
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
	if (Iter0+iLoop>atmWriteTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
	  if (atmWriteTimer.count<atmWriteTimer.numTimeSteps) {
		atmWriteTimer.count++;
	  }
	}

	if (Iter0+iLoop>=atmWriteTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
      if ((atmWriteTimer.count==atmWriteTimer.numTimeSteps) || (Iter0+iLoop==atmWriteTimer.startTimeStep)) { /* time to write out */

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
		
      }

      if (atmWriteTimer.count==atmWriteTimer.numTimeSteps) {
		ierr = updateStepTimer("atm_write_", Iter0+iLoop, &atmWriteTimer);CHKERRQ(ierr);
	  }
    }

  }

  if (calcDiagnostics) {  
	if (Iter0+iLoop>=diagTimer.startTimeStep) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */

#ifdef ALLOW_C14
      ierr = VecPointwiseDivide(Rocn,DIC14,DIC);CHKERRQ(ierr);
#endif

	  if (diagTimer.count<diagTimer.numTimeSteps) { /* still within same averaging block so accumulate */
	  
        for (ip=0; ip<lNumProfiles; ip++) {
          localpco2diagavg[ip]=localpco2diag[ip]+localpco2diagavg[ip];
          localgasexfluxdiagavg[ip]=localgasexfluxdiag[ip]+localgasexfluxdiagavg[ip];
          localtotfluxdiagavg[ip]=localtotfluxdiag[ip]+localtotfluxdiagavg[ip];      
#ifdef ALLOW_C14
          localc14gasexfluxdiagavg[ip]=localc14gasexfluxdiag[ip]+localc14gasexfluxdiagavg[ip];
          localc14totfluxdiagavg[ip]=localc14totfluxdiag[ip]+localc14totfluxdiagavg[ip];
#endif              
        }	  

#ifdef ALLOW_C14
        ierr = VecAXPY(Rocnavg,one,Rocn);CHKERRQ(ierr);
#endif              

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
		ierr = PetscFPrintf(PETSC_COMM_WORLD,diagfptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);           

        for (ip=0; ip<lNumProfiles; ip++) {
          localpco2diagavg[ip]=localpco2diagavg[ip]/diagTimer.count;
          localgasexfluxdiagavg[ip]=localgasexfluxdiagavg[ip]/diagTimer.count;
          localtotfluxdiagavg[ip]=localtotfluxdiagavg[ip]/diagTimer.count;
#ifdef ALLOW_C14
          localc14gasexfluxdiagavg[ip]=localc14gasexfluxdiagavg[ip]/diagTimer.count;
          localc14totfluxdiagavg[ip]=localc14totfluxdiagavg[ip]/diagTimer.count;
#endif              
        }	  

#ifdef ALLOW_C14
        ierr = VecScale(Rocnavg,1.0/diagTimer.count);CHKERRQ(ierr);
#endif

        ierr = writeProfileSurfaceScalarData("pCO2_surf.bin",localpco2diagavg,1,diagAppendOutput);  		
        ierr = writeProfileSurfaceScalarData("gasexflux_surf.bin",localgasexfluxdiagavg,1,diagAppendOutput);  		
        ierr = writeProfileSurfaceScalarData("totalflux_surf.bin",localtotfluxdiagavg,1,diagAppendOutput);  		
#ifdef ALLOW_C14
        ierr = writeProfileSurfaceScalarData("c14gasexflux_surf.bin",localc14gasexfluxdiagavg,1,diagAppendOutput);  		
        ierr = writeProfileSurfaceScalarData("c14totalflux_surf.bin",localc14totfluxdiagavg,1,diagAppendOutput);  		
        ierr = VecView(Rocnavg,fdrocn);CHKERRQ(ierr);
#endif		  
        if (useAtmModel) {
          pCO2atmavg=pCO2atmavg/diagTimer.count;
          Foceanavg=Foceanavg/diagTimer.count;   
          ierr = writeBinaryScalarData("pCO2atm_avg.bin",&pCO2atmavg,1,diagAppendOutput);  		
          ierr = writeBinaryScalarData("Focean_avg.bin",&Foceanavg,1,diagAppendOutput);  		

		  if (useLandModel) {
			Flandavg=Flandavg/diagTimer.count;        
			landStateavg[0]=landStateavg[0]/diagTimer.count;
			landStateavg[1]=landStateavg[1]/diagTimer.count;
			landStateavg[2]=landStateavg[2]/diagTimer.count;           
			ierr = writeBinaryScalarData("Fland_avg.bin",&Flandavg,1,diagAppendOutput);  		
			ierr = writeBinaryScalarData("land_state_avg.bin",landStateavg,3,diagAppendOutput);
		  }          
        }

		diagAppendOutput=PETSC_TRUE;
		DIAG_FILE_MODE=FILE_MODE_APPEND;

/*      reset diagnostic arrays */
        for (ip=0; ip<lNumProfiles; ip++) {
          localpco2diagavg[ip]=0.0;
          localgasexfluxdiagavg[ip]=0.0;
          localtotfluxdiagavg[ip]=0.0;
#ifdef ALLOW_C14
          localc14gasexfluxdiagavg[ip]=0.0;
          localc14totfluxdiagavg[ip]=0.0;
#endif              
        }	  

#ifdef ALLOW_C14
       ierr = VecSet(Rocnavg,zero);CHKERRQ(ierr);        
#endif

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
/*   write instantaneous land model state */
	  ierr = writeBinaryScalarData("pickup_land_state.bin",landState,3,PETSC_FALSE);
	}    
  }

  if (timeDependentBiogeochemForcing) {
    ierr = destroyTimeDependentArray(&localTstd);CHKERRQ(ierr);
    ierr = destroyTimeDependentArray(&localSstd);CHKERRQ(ierr);
    ierr = destroyTimeDependentArray(&localEmPtd);CHKERRQ(ierr);
    ierr = destroyTimeDependentArray(&localficetd);CHKERRQ(ierr);
    ierr = destroyTimeDependentArray(&localAlktd);CHKERRQ(ierr);
  } else if (periodicBiogeochemForcing) {    
    ierr = destroyPeriodicArray(&localTsp);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localSsp);CHKERRQ(ierr);
	ierr = destroyPeriodicArray(&localEmPp);CHKERRQ(ierr);    
    ierr = destroyPeriodicArray(&localficep);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localAlkp);CHKERRQ(ierr);   
  }    
  
  if (periodicBiogeochemForcing) {    
    if (useWinds) {
      ierr = destroyPeriodicArray(&localuwindp);CHKERRQ(ierr);
      ierr = destroyPeriodicArray(&localvwindp);CHKERRQ(ierr);
    } else {   
      ierr = destroyPeriodicArray(&localxkwp);CHKERRQ(ierr);
    }
    ierr = destroyPeriodicArray(&localatmospp);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localPO4p);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localSiO2p);CHKERRQ(ierr);
  }    

  if (useAtmModel) {
    ierr = PetscFClose(PETSC_COMM_WORLD,atmfptime);CHKERRQ(ierr);
  }
  
  if (useVirtualFlux) {
    ierr = VecDestroy(&surfVolFrac);CHKERRQ(ierr);  
  }

  if (calcDiagnostics) {      
	ierr = PetscFClose(PETSC_COMM_WORLD,diagfptime);CHKERRQ(ierr);  	    
#ifdef ALLOW_C14  
    ierr = VecDestroy(&Rocn);CHKERRQ(ierr);
    ierr = VecDestroy(&Rocnavg);CHKERRQ(ierr);    
    ierr = PetscViewerDestroy(&fdrocn);CHKERRQ(ierr);
#endif  
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

  if (timeDependentBiogeochemForcing) {
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localTs,timeDependentBiogeochemTimer.numTimes,timeDependentBiogeochemTimer.tdt,&localTstd,"Ts.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localSs,timeDependentBiogeochemTimer.numTimes,timeDependentBiogeochemTimer.tdt,&localSstd,"Ss.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localEmP,timeDependentBiogeochemTimer.numTimes,timeDependentBiogeochemTimer.tdt,&localEmPtd,"EmP.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localfice,timeDependentBiogeochemTimer.numTimes,timeDependentBiogeochemTimer.tdt,&localficetd,"fice.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localAlk,timeDependentBiogeochemTimer.numTimes,timeDependentBiogeochemTimer.tdt,&localAlktd,"Alk.bin");    
  } else if (periodicBiogeochemForcing) {    
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localTs,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localTsp,"Ts_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localSs,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localSsp,"Ss_");
	ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localEmPp,"EmP_");                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localfice,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localficep,"fice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localAlk,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localAlkp,"Alk_");    
  }    
  if (periodicBiogeochemForcing) {   
    if (useWinds) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localuwind,windsTimer.cyclePeriod,windsTimer.numPerPeriod,windsTimer.tdp,&localuwindp,"uwind_");                                                    
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localvwind,windsTimer.cyclePeriod,windsTimer.numPerPeriod,windsTimer.tdp,&localvwindp,"vwind_");                                                        
    } else {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localxkw,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localxkwp,"xkw_");
    }
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localatmospp,"atmosp_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localPO4,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localPO4p,"PO4_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localSiO2,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&localSiO2p,"SiO2_");
  }  

/* Initialize biogeochem model */
  myTime = DeltaT*Iter; /* Iter should start at 0 */  
  for (ip=0; ip<lNumProfiles; ip++) {
	kl=lStartIndices[ip];  

    ocmip_abiotic_carbon_ini_(&Iter,&myTime,&localDIC[kl],&localAlk[ip],&localPO4[ip],&localSiO2[ip],                 
                           &localTs[ip],&localSs[ip],&pH[ip]);
  }

  return 0;
}
