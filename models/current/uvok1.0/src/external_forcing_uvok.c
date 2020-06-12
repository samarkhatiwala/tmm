#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define READ_SWRAD
#include "petscmat.h"
#include "petsc_matvec_utils.h"
#include "tmm_timer.h"
#include "tmm_forcing_utils.h"
#include "tmm_profile_utils.h"
#include "tmm_profile_data.h"
#include "tmm_misfit.h"
#include "tmm_main.h"
#include "uvok_tmm.h"

#define TR1 v[0]

Vec Ts,Ss;
PetscScalar *localTs,*localSs;
PetscScalar **localTR, **localJTR;
PetscBool useEmP = PETSC_TRUE;
PetscScalar *localEmP, EmPglobavg;
PeriodicArray localEmPp;
Vec surfVolFrac;
PetscScalar Sglobavg = 0.0, *TRglobavg;
PetscScalar *localwind,*localaice,*localhice,*localhsno,*localdz,*localatmosp;
PetscScalar *localswrad;
PetscScalar *locallatitude;
/* PetscScalar *localsgbathy; */

#ifdef O_npzd_fe_limitation
Vec Fe_dissolved;
PetscScalar *localFe_dissolved;
PeriodicVec Fe_dissolvedp;
#endif
#ifdef O_npzd_iron
PetscScalar *localFe_adep, *localFe_detr_flux;
PeriodicArray localFe_adepp, localFe_detr_fluxp;
Vec Fe_hydr;
PetscScalar *localFe_hydr;
#endif

PetscScalar daysPerYear, secondsPerYear;

PetscInt nzmax;
PetscScalar DeltaT;
PetscScalar *drF, *zt;

PeriodicVec Tsp, Ssp;
PeriodicArray localwindp,localaicep,localhicep,localhsnop,localatmospp;
#ifdef READ_SWRAD
PeriodicArray localswradp;
#endif

PetscBool periodicBiogeochemForcing = PETSC_FALSE;
PeriodicTimer biogeochemTimer;

PetscInt maxValsToRead;

/* variables for prescribed atmospheric CO2 */
#if defined O_carbon
#if defined O_co2ccn_data
char *pCO2atmFiles[2];  
PetscInt numpCO2atm_hist = 0;
PetscScalar *TpCO2atm_hist, *pCO2atm_hist;
PetscScalar pCO2atm = 277.0; /* default initial value */
#endif

#if defined O_TMM_interactive_atmosphere
/* Land/Atm model variables */
PetscBool useLandModel = PETSC_FALSE;
PetscBool useEmissions = PETSC_FALSE;
PetscBool interpEmissions = PETSC_FALSE;
char *emFiles[3];  
PetscInt numEmission_hist = 0;
PetscScalar *Tem_hist, *E_hist, *D_hist;
PetscScalar fossilFuelEmission = 0.0, landUseEmission = 0.0;
PetscScalar cumulativeEmissions = 0.0;
PetscScalar pCO2atm_ini = 277.0; /* default initial value */
PetscScalar pCO2atm = 277.0; /* default initial value */
char pCO2atmIniFile[PETSC_MAX_PATH_LEN];
PetscScalar *localdA;
PetscScalar landState[3], landSource[3];
PetscScalar deltaTsg = 0.0;
PetscScalar ppmToPgC=2.1324;
PetscScalar atmModelDeltaT;
PetscScalar Fland = 0.0, Focean=0.0;

StepTimer atmWriteTimer;
PetscBool atmAppendOutput;
FILE *atmfptime;
char atmOutTimeFile[PETSC_MAX_PATH_LEN];  
PetscScalar pCO2atmavg, Foceanavg, Flandavg, landStateavg[3];
#endif

# if defined O_c14ccn_data
PetscScalar dc14ccnnatm = 0.0;
PetscScalar dc14ccnsatm = 0.0;
PetscScalar dc14ccneatm = 0.0;
PetscScalar DC14atm = 0.0; /* default initial value */
char *C14atmFiles[2];
PetscScalar *TC14atm_hist, *C14atm_hist;
PetscInt numC14atm_hist = 0;
#endif                 
#endif /* O_carbon */

PetscBool calcDiagnostics = PETSC_FALSE;
StepTimer diagTimer;
PetscBool diagFirstTime = PETSC_TRUE;
PetscBool diagAppendOutput = PETSC_FALSE;
PetscFileMode DIAG_FILE_MODE;
FILE *diagfptime;
PetscViewer diagfd;
PetscInt diagfp;
char diagOutTimeFile[PETSC_MAX_PATH_LEN];
/* Add model specific diagnostic variables below */
#define MAXDIAGS2d 40
#define MAXDIAGS3d 40
PetscScalar *localDiag2d[MAXDIAGS2d], *localDiag2davg[MAXDIAGS2d];
Vec *Diag3d, *Diag3davg;
PetscScalar **localDiag3d, **localDiag3davg;
PetscInt numDiags2d = 0, numDiags3d = 0, id2d, id3d;
PetscViewer fddiag3dout[MAXDIAGS3d];
char *Diag2dFile[MAXDIAGS2d], *Diag3dFile[MAXDIAGS3d];
PetscInt doAverage=0;

PetscMPIInt myId;
PetscInt debugFlag = 0;

PetscBool MYTRUE = PETSC_TRUE, MYFALSE = PETSC_FALSE;

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
  int fp;
  PetscBool flg;
  PetscInt it, m;
  PetscScalar myTime;
  PetscScalar zero = 0.0;
  
#if defined (FORSPINUP) || defined (FORJACOBIAN)
  ierr = PetscOptionsHasName(NULL,NULL,"-relax_tracer",&relaxTracer);CHKERRQ(ierr);
  if (relaxTracer) {  

    maxValsToRead = numTracers;
    ierr = PetscOptionsGetRealArray(NULL,NULL,"-relax_tau",relaxTau,&maxValsToRead,&flg);
    if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate tracer relaxation tau with the -relax_tau option");
    if (maxValsToRead != numTracers) {
      SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of relaxation tau values specified");
    }

    maxValsToRead = numTracers;
    ierr = PetscOptionsGetRealArray(NULL,NULL,"-relax_value",relaxValue,&maxValsToRead,&flg);
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

/* This sequence should be identical to that in S/R tracer_init in UVic_ESCM.F */
  m = 0;
# if defined O_carbon
#define dic v[0]
#define sdic ut[0]
  m++;
#  if defined O_carbon_13
#define dic13 v[1]
#define sdic13 ut[1]
  m++;
#  endif
#  if defined O_carbon_14
#define c14 v[1]
#define sc14 ut[1]
  m++;
#  endif
# endif
# if defined O_npzd_alk
#define alk v[2]
#define salk ut[2]
  m++;
# endif
# if defined O_npzd_o2
#define o2 v[3]
#define so2 ut[3]
  m++;
# endif
# if defined O_npzd
#define po4 v[4]
#define spo4 ut[4]
  m++;
  /*#define dop v[6]
#define sdop ut[6]
m++; */
#define phyt v[5]
#define sphyt ut[5]
  m++;
#define zoop v[6]
#define szoop ut[6]
  m++;
#define detr v[7]
#define sdetr ut[7]
  m++;
#  if defined O_npzd_nitrogen
#define no3 v[8]
#define sno3 ut[8]
  m++;
  /*#define don v[11]
#define sdon ut[11]
m++; */
#define diaz v[9]
#define sdiaz ut[9]
  m++;
#   if defined O_npzd_nitrogen_15
#define din15 v[13]
#define sdin15 ut[13]
  m++;
#define don15 v[14]
#define sdon15 ut[14]
  m++;
#define phytn15 v[15]
#define sphytn15 ut[15]
  m++;
#define zoopn15 v[16]
#define szoopn15 ut[16]
  m++;
#define detrn15 v[17]
#define sdetrn15 ut[17]
  m++;
#define diazn15 v[18]
#define sdiazn15 ut[18]
  m++;
#   endif
#  endif
#  if defined O_carbon_13
#define doc13 v[19]
#define sdoc13 ut[19]
  m++;
#define phytc13 v[20]
#define sphytc13 ut[20]
  m++;
#define zoopc13 v[21]
#define szoopc13 ut[21]
  m++;
#define detrc13 v[22]
#define sdetrc13 ut[22]
  m++;
#   if defined O_npzd_nitrogen
#define diazc13 v[23]
#define sdiazc13 ut[23]
  m++;
#   endif
#  if defined O_npzd_iron
#define dfe v[24]
#define sdfe ut[24]
  m++;
#define detrfe v[25]
#define sdetrfe ut[25]
  m++;
#  endif
#  endif
# endif
  
  if (m != numTracers) {
    SETERRQ(PETSC_COMM_WORLD,1,"Error: numTracers does not match expected number of tracers!");    
  }
  
  for (itr=0; itr<numTracers; itr++) {
    ierr = VecSet(ut[itr],zero); CHKERRQ(ierr);
  }

  ierr = VecGetArrays(v,numTracers,&localTR);CHKERRQ(ierr);
  ierr = VecGetArrays(ut,numTracers,&localJTR);CHKERRQ(ierr);
  
  ierr = PetscOptionsGetReal(NULL,NULL,"-biogeochem_deltat",&DeltaT,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate biogeochemical time step in seconds with the -biogeochem_deltat option");  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Ocean time step for BGC length is  %12.7f seconds\n",DeltaT);CHKERRQ(ierr);

  ierr = PetscOptionsGetReal(NULL,NULL,"-days_per_year",&daysPerYear,&flg);CHKERRQ(ierr);
  if (!flg) {
    daysPerYear = 360.0;
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of days per year is %12.7f\n",daysPerYear);CHKERRQ(ierr);
  secondsPerYear = 86400.0*daysPerYear;

  ierr = PetscOptionsHasName(NULL,NULL,"-periodic_biogeochem_forcing",&periodicBiogeochemForcing);CHKERRQ(ierr);

  if (periodicBiogeochemForcing) {    
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic biogeochemical forcing specified\n");CHKERRQ(ierr);

    ierr = iniPeriodicTimer("periodic_biogeochem_", &biogeochemTimer);CHKERRQ(ierr);

/*  IMPORTANT: time units must be the same as that used by the toplevel driver */
  }
  
/*   Read T and S */
  ierr = VecDuplicate(TR1,&Ts);CHKERRQ(ierr);
  ierr = VecDuplicate(TR1,&Ss);CHKERRQ(ierr);  
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

#if defined O_carbon
#if defined O_co2ccn_data && defined O_TMM_interactive_atmosphere
  SETERRQ(PETSC_COMM_WORLD,1,"Cannot use both O_co2ccn_data and O_TMM_interactive_atmosphere options at the same time!");
#endif

#if defined O_co2ccn_data
  /* prescribed atmospheric CO2 */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using prescribed atmospheric CO2\n");CHKERRQ(ierr);
  maxValsToRead = 2;
  pCO2atmFiles[0] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* time file */
  pCO2atmFiles[1] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* atmospheric pCO2 history file */
  ierr = PetscOptionsGetStringArray(NULL,NULL,"-pco2atm_history",pCO2atmFiles,&maxValsToRead,&flg);CHKERRQ(ierr);
  if (flg) { /* Read atmospheric pCO2 history */
	if (maxValsToRead != 2) {
	  SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of file names specified for atmospheric pCO2 history");
	}      
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
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,pCO2atmFiles[1],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscMalloc(numpCO2atm_hist*sizeof(PetscScalar),&pCO2atm_hist);CHKERRQ(ierr); 
	ierr = PetscBinaryRead(fp,pCO2atm_hist,numpCO2atm_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	
	pCO2atm = pCO2atm_hist[0];

  }	else {
	SETERRQ(PETSC_COMM_WORLD,1,"Must specify atmospheric CO2 history with -pco2atm_history when using the O_co2ccn_data option");
  }    
#endif

#if defined O_TMM_interactive_atmosphere
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
	ierr = PetscFPrintf(PETSC_COMM_WORLD,atmfptime,"%d   %10.5f\n",Iter0,time0);CHKERRQ(ierr);     
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing atmospheric output at time %10.5f, step %d\n", tc,Iter);CHKERRQ(ierr);  
	ierr = writeBinaryScalarData("pCO2atm_output.bin",&pCO2atm,1,PETSC_FALSE);
  } else {
	ierr = PetscFOpen(PETSC_COMM_WORLD,atmOutTimeFile,"a",&atmfptime);CHKERRQ(ierr);  
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output will be appended. Initial condition will NOT be written\n");CHKERRQ(ierr);      
  }

  ierr = PetscOptionsHasName(NULL,NULL,"-use_land_model",&useLandModel);CHKERRQ(ierr);
  if (useLandModel) {
	SETERRQ(PETSC_COMM_WORLD,1,"Land model not yet supported!");
  
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
#endif

# if defined O_c14ccn_data
  /* prescribed atmospheric DeltaC14 */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using prescribed atmospheric Delta C14\n");CHKERRQ(ierr);
  maxValsToRead = 2;
  C14atmFiles[0] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* time file */
  C14atmFiles[1] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* atmospheric DeltaC14 history file */
  ierr = PetscOptionsGetStringArray(NULL,NULL,"-c14atm_history",C14atmFiles,&maxValsToRead,&flg);CHKERRQ(ierr);
  if (flg) { /* Read atmospheric Delta C14 history */
	if (maxValsToRead != 2) {
	  SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of file names specified for atmospheric Delta C14 history");
	}      
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading time-dependent atmospheric Delta C14 history\n");CHKERRQ(ierr);      
	/* read time data */
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,C14atmFiles[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscBinaryRead(fp,&numC14atm_hist,1,NULL,PETSC_INT);CHKERRQ(ierr);  
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of points in atmospheric history file is %d \n",numC14atm_hist);CHKERRQ(ierr);  
	ierr = PetscMalloc(numC14atm_hist*sizeof(PetscScalar),&TC14atm_hist);CHKERRQ(ierr); 
	ierr = PetscBinaryRead(fp,TC14atm_hist,numC14atm_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	/* read atmospheric Delta C14 data */
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,C14atmFiles[1],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscMalloc(numC14atm_hist*sizeof(PetscScalar),&C14atm_hist);CHKERRQ(ierr); 
	ierr = PetscBinaryRead(fp,C14atm_hist,numC14atm_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	
	DC14atm = C14atm_hist[0];
	dc14ccnnatm = DC14atm;
	dc14ccnsatm = DC14atm;
	dc14ccneatm = DC14atm;
	
  }	else {
	SETERRQ(PETSC_COMM_WORLD,1,"Must specify atmospheric Delta C14 history with -c14atm_history when using the O_c14ccn_data option");		
  }  
#endif
#endif /* O_carbon */

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localEmP);CHKERRQ(ierr);
  for (ip=0; ip<lNumProfiles; ip++) { /* initialize to zero to be safe */
    localEmP[ip]=0.0;
  }
  EmPglobavg = 0.0; /* set this to zero for now */  

  ierr = PetscOptionsHasName(NULL,NULL,"-no_use_emp",&flg);CHKERRQ(ierr);
  if (flg) {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: -no_use_emp has been specified. E-P will be set to zero.\n");CHKERRQ(ierr);
	useEmP = PETSC_FALSE;
  }  
  
  if (useEmP) {
   if (periodicBiogeochemForcing) {    
	 localEmPp.firstTime = PETSC_TRUE;
	 localEmPp.arrayLength = lNumProfiles;
   } else {  
	 ierr = readProfileSurfaceScalarData("EmP.bin",localEmP,1);  
   }
  }
  
 /* always need this */  
  ierr = VecDuplicate(TR1,&surfVolFrac);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"surface_volume_fraction.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = VecLoad(surfVolFrac,fd);CHKERRQ(ierr);  
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      

  ierr = PetscMalloc(numTracers*sizeof(PetscScalar),&TRglobavg);CHKERRQ(ierr);	  

/* Grid arrays */
/*  ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localsgbathy);CHKERRQ(ierr);    
    ierr = VecLoadVecIntoArray(TR1,"sgbathy.petsc",localsgbathy);CHKERRQ(ierr); */

  ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localdz);CHKERRQ(ierr);    
  ierr = VecLoadVecIntoArray(TR1,"dz.petsc",localdz);CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"zt.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fp,&nzmax,1,NULL,PETSC_INT);CHKERRQ(ierr);  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of vertical layers is %d \n",nzmax);CHKERRQ(ierr);
  ierr = PetscMalloc(nzmax*sizeof(PetscScalar),&zt);CHKERRQ(ierr); 
  ierr = PetscBinaryRead(fp,zt,nzmax,NULL,PETSC_SCALAR);CHKERRQ(ierr);  
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"drF.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
  ierr = PetscBinaryRead(fp,&nzmax,1,NULL,PETSC_INT);CHKERRQ(ierr);  
  ierr = PetscMalloc(nzmax*sizeof(PetscScalar),&drF);CHKERRQ(ierr); 
  ierr = PetscBinaryRead(fp,drF,nzmax,NULL,PETSC_SCALAR);CHKERRQ(ierr);  
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

/* Forcing fields */
#ifdef O_npzd_fe_limitation
  ierr = VecDuplicate(TR1,&Fe_dissolved);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    Fe_dissolvedp.firstTime = PETSC_TRUE;
  } else {
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Fe_dissolved.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = VecLoad(Fe_dissolved,fd);CHKERRQ(ierr);    
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
  }  
  ierr = VecGetArray(Fe_dissolved,&localFe_dissolved);CHKERRQ(ierr);
#endif
#ifdef O_npzd_iron
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localFe_adep);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localFe_adepp.firstTime = PETSC_TRUE;
    localFe_adepp.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("Fe_adep.bin",localFe_adep,1);  
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localFe_detr_flux);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localFe_detr_fluxp.firstTime = PETSC_TRUE;
    localFe_detr_fluxp.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("Fe_detr_flux.bin",localFe_detr_flux,1);  
  }

  ierr = VecDuplicate(TR1,&Fe_hydr);CHKERRQ(ierr);  
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Fe_hydr.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = VecLoad(Fe_hydr,fd);CHKERRQ(ierr);    
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
  ierr = VecGetArray(Fe_hydr,&localFe_hydr);CHKERRQ(ierr);
#endif

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localswrad);CHKERRQ(ierr);  
#ifdef READ_SWRAD
  if (periodicBiogeochemForcing) {    
    localswradp.firstTime = PETSC_TRUE;
    localswradp.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("swrad.bin",localswrad,1);  
  }
#endif	  

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&locallatitude);CHKERRQ(ierr);  
  ierr = readProfileSurfaceScalarData("latitude.bin",locallatitude,1);  
  
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localaice);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localaicep.firstTime = PETSC_TRUE;
    localaicep.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("aice.bin",localaice,1);  
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localhice);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localhicep.firstTime = PETSC_TRUE;
    localhicep.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("hice.bin",localhice,1);  
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localhsno);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    localhsnop.firstTime = PETSC_TRUE;
    localhsnop.arrayLength = lNumProfiles;
  } else {  
    ierr = readProfileSurfaceScalarData("hsno.bin",localhsno,1);  
  }

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

/* Initialize biogeochem model */
  myTime = DeltaT*Iter; /* Iter should start at 0 */

  if (periodicBiogeochemForcing) {   
	ierr = interpPeriodicVector(tc,&Ts,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Ssp,"Ss_");
#ifdef O_npzd_fe_limitation
	ierr = interpPeriodicVector(tc,&Fe_dissolved,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Fe_dissolvedp,"Fe_dissolved_");		
#endif
#ifdef O_npzd_iron
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localFe_adep,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localFe_adepp,"Fe_adep_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localFe_detr_flux,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localFe_detr_fluxp,"Fe_detr_flux_");
#endif
#ifdef READ_SWRAD
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localswradp,"swrad_");
#else
   insolation_(&lNumProfiles,&myTime,&locallatitude[0],&localswrad[0],&localtau[0]);
#endif                                                 
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localaice,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localaicep,"aice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhice,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localhicep,"hice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhsno,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localhsnop,"hsno_");                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localwindp,"wind_");   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localatmospp,"atmosp_");					                              
    if (useEmP) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                    biogeochemTimer.tdp,&localEmPp,"EmP_");                                                  
    }
  } else {
#ifndef READ_SWRAD
    insolation_(&lNumProfiles,&myTime,&locallatitude[0],&localswrad[0],&localtau[0]);
#endif    
  }

/* compute global means */
  ierr = VecDot(surfVolFrac,Ss,&Sglobavg);CHKERRQ(ierr); /* volume weighted mean surface salinity */   
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Global average of salinity =%10.5f\n", Sglobavg);CHKERRQ(ierr); 
  for (itr=0; itr<numTracers; itr++) {    
	ierr = VecDot(surfVolFrac,v[itr],&TRglobavg[itr]);CHKERRQ(ierr); /* volume weighted mean surface TR */
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Global average of tracer %d = %10.5f\n", itr, TRglobavg[itr]);CHKERRQ(ierr);									                    
  }  

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myId);CHKERRQ(ierr);
  if (myId == 0) debugFlag = 1;
  
  uvok_ini_(&nzmax,zt,drF,&DeltaT,&Sglobavg,TRglobavg,&debugFlag);  

  ierr = PetscOptionsHasName(NULL,NULL,"-calc_diagnostics",&calcDiagnostics);CHKERRQ(ierr);
  if (calcDiagnostics) {    
/*Data for diagnostics */
    ierr = iniStepTimer("diag_", Iter0, &diagTimer);CHKERRQ(ierr);

	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed starting at (and including) time step: %d\n", diagTimer.startTimeStep);CHKERRQ(ierr);	
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
	
    diagFirstTime=PETSC_TRUE;
    doAverage=0;

    uvok_diags_ini_(&lNumProfiles, &lSize, &numDiags2d, &numDiags3d, &debugFlag);

    if (numDiags2d > MAXDIAGS2d) {
      SETERRQ(PETSC_COMM_WORLD,1,"Number of 2-d diagnostics requested exceeds maximum. You must increase MAXDIAGS2d!");
    }    
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of 2-d diagnostics requested: %d\n",numDiags2d);CHKERRQ(ierr);	                            	  
    
    if (numDiags2d>0) {
	  for (id2d=0; id2d<numDiags2d; id2d++) {
		ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localDiag2davg[id2d]);CHKERRQ(ierr);
		Diag2dFile[id2d] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
		strcpy(Diag2dFile[id2d],"");
		for (ip=0; ip<lNumProfiles; ip++) {
		  localDiag2davg[id2d][ip]=0.0;      
		}
	  }
    }
   
/* Set the number of 3-d diagnostics */
/*	uvok_diags3d_ini_(&numDiags3d, &minusone, &debugFlag); */
	if (numDiags3d > MAXDIAGS3d) {
	  SETERRQ(PETSC_COMM_WORLD,1,"Number of 3-d diagnostics requested exceeds maximum. You must increase MAXDIAGS3d!");
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of 3-d diagnostics requested: %d\n",numDiags3d);CHKERRQ(ierr);	                            	  
	
	if (numDiags2d == 0 & numDiags3d == 0) {
	  SETERRQ(PETSC_COMM_WORLD,1,"You have specified the -calc_diagnostics flag but no diagnostics are being computed!");
    }

	if (numDiags3d>0) {
	  ierr = VecDuplicateVecs(TR1,numDiags3d,&Diag3davg);CHKERRQ(ierr);

	  for (id3d=0; id3d<numDiags3d; id3d++) {
		ierr = VecSet(Diag3davg[id3d],zero);CHKERRQ(ierr);
		Diag3dFile[id3d] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
		strcpy(Diag3dFile[id3d],"");
	  }
	  ierr = VecGetArrays(Diag3davg,numDiags3d,&localDiag3davg);CHKERRQ(ierr);        
    }
    
    doAverage=0;
    
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
  PetscScalar relyr, day;
  PetscInt toModel = 1; 
  PetscInt fromModel = 2;
  PetscInt itf;
  PetscScalar alpha;
  PetscInt k;
  PetscScalar localFocean;  
  PetscScalar localgasexflux = 0.0, localtotflux = 0.0;

  myTime = DeltaT*Iter; /* Iter should start at 0 */

  if (periodicBiogeochemForcing) {   
	ierr = interpPeriodicVector(tc,&Ts,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Ssp,"Ss_");	
#ifdef O_npzd_fe_limitation	
	ierr = interpPeriodicVector(tc,&Fe_dissolved,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Fe_dissolvedp,"Fe_dissolved_");	
#endif	
#ifdef O_npzd_iron
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localFe_adep,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localFe_adepp,"Fe_adep_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localFe_detr_flux,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localFe_detr_fluxp,"Fe_detr_flux_");
#endif
#ifdef READ_SWRAD
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localswradp,"swrad_");
#else
    insolation_(&lNumProfiles,&myTime,&locallatitude[0],&localswrad[0],&localtau[0]);
#endif                                                 
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localaice,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localaicep,"aice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhice,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localhicep,"hice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhsno,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localhsnop,"hsno_");                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localwindp,"wind_");  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localatmospp,"atmosp_");					                              
    if (useEmP) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                    biogeochemTimer.tdp,&localEmPp,"EmP_");                                                      
    }
  }

#if defined O_carbon
#if defined O_co2ccn_data
/* Interpolate atmospheric pCO2   */
  if (tc>=TpCO2atm_hist[0]) {
	ierr = calcInterpFactor(numpCO2atm_hist,tc,TpCO2atm_hist,&itf,&alpha); CHKERRQ(ierr);
	pCO2atm = alpha*pCO2atm_hist[itf] + (1.0-alpha)*pCO2atm_hist[itf+1];	  
  } else {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: time < %10.5f. Assuming pCO2atm=%g\n",TpCO2atm_hist[0],pCO2atm);CHKERRQ(ierr);
  }
#endif

#if defined O_TMM_interactive_atmosphere
  if (useEmissions) {  
/* Interpolate emissions */
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
#endif

# if defined O_c14ccn_data
  if (tc>=TC14atm_hist[0]) {
	ierr = calcInterpFactor(numC14atm_hist,tc,TC14atm_hist,&itf,&alpha); CHKERRQ(ierr);
	DC14atm = alpha*C14atm_hist[itf] + (1.0-alpha)*C14atm_hist[itf+1];	  
	dc14ccnnatm = DC14atm;
	dc14ccnsatm = DC14atm;
	dc14ccneatm = DC14atm;
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: time < %10.5f. Assuming DC14atm=%g\n",TC14atm_hist[0],DC14atm);CHKERRQ(ierr);
  }
#endif
#endif /* O_carbon */

/* compute global means */
  if (useEmP) {
#  if !defined O_constant_flux_reference 
	for (itr=0; itr<numTracers; itr++) {    
      ierr = VecDot(surfVolFrac,v[itr],&TRglobavg[itr]);CHKERRQ(ierr); /* volume weighted mean surface TR */									              
    }
# endif    
/*    EmPglobavg = 0.0; */ /* this is set to zero above */
  }

  relyr = myTime/secondsPerYear; /* number of years (and fractional years) of model */
  day = myTime/86400.0 - floor(relyr)*daysPerYear; /* relative day number referenced to the beginning of the current year */

#if defined O_carbon
#if defined O_TMM_interactive_atmosphere
  localFocean = 0.0;  
  Focean = 0.0;
  Fland = 0.0;
#endif
#endif /* O_carbon */

	if (calcDiagnostics) {  
	  if (Iter0+iLoop==diagTimer.startTimeStep) { /* start time averaging (note: diagStartTimeStep is ABSOLUTE time step) */	
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Switching on diagnostics accumulation at step %d\n", Iter0+iLoop);CHKERRQ(ierr);
		uvok_diags_start_(&debugFlag); /* Switch on diagnostics accumulation at start of cycle */
	  }	
	}                         

  for (ip=0; ip<lNumProfiles; ip++) {
    nzloc=lProfileLength[ip];
    kl=lStartIndices[ip];

	for (itr=0; itr<numTracers; itr++) {    	
	  uvok_copy_data_(&nzloc,&itr,&localTR[itr][kl],&toModel);
	}  

    uvok_calc_(&nzloc,&locallatitude[ip],&day,&relyr,
               &localTs[kl],&localSs[kl],&TRglobavg[0],&localdz[kl],zt,
# if defined O_carbon
#if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
               &pCO2atm,
#endif               
               &localwind[ip],
#endif
# if defined O_c14ccn_data
               &dc14ccnnatm,&dc14ccnsatm,&dc14ccneatm,
#endif                 
#  if defined O_npzd_nitrogen_15
               &localsgbathy[kl],
#  endif
#  if defined O_npzd_fe_limitation
               &localFe_dissolved[kl],
#  endif
#ifdef O_npzd_iron
               &localFe_adep[ip], &localFe_detr_flux[ip], &localFe_hydr[kl],
#endif
#  if defined O_embm
               &localswrad[ip],
#  endif
#  if defined O_ice
#   if !defined O_ice_cpts
               &localaice[ip], &localhice[ip], &localhsno[ip],
#   endif
#  endif
               &localEmP[ip], &EmPglobavg, 
# if defined O_carbon
               &localgasexflux, &localtotflux, 
# endif                    
               &debugFlag);

	for (itr=0; itr<numTracers; itr++) {    
	  uvok_copy_data_(&nzloc,&itr,&localJTR[itr][kl],&fromModel);
	}  

	if (calcDiagnostics) {  
	  if (Iter0+iLoop>=diagTimer.startTimeStep) { /* start time averaging (note: diagStartTimeStep is ABSOLUTE time step) */	
	    uvok_diags_accumulate_(&ip, &kl, &nzloc, &diagTimer.count, &doAverage, &debugFlag);
	  }	
	}                         

#if defined O_carbon
#if defined O_TMM_interactive_atmosphere
/*  Note: following UVic (gasbc.F) we use the gas exchange flux to compute atmospheric CO2 evolution. */
/*  The virtual flux term included in localtotflux should go be zero when integrated over the global ocean surface */
	localFocean = localFocean + localgasexflux*localdA[ip]*(12.0/1.e15)*secondsPerYear; /* PgC/y */
#endif
#endif /* O_carbon */

  } /* end loop over profiles */

  if (calcDiagnostics) {  
	if (Iter0+iLoop>=diagTimer.startTimeStep) { /* start time averaging (note: diagStartTimeStep is ABSOLUTE time step) */	
/*    still within same averaging block so increment accumulate count */	
	  diagTimer.count++;
	}	
  }
    
#if defined O_carbon
#if defined O_TMM_interactive_atmosphere
  MPI_Allreduce(&localFocean, &Focean, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  if (Iter0+iLoop>=atmWriteTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
	if (atmWriteTimer.count<=atmWriteTimer.numTimeSteps) {
	  atmWriteTimer.count++;
	}
	if (atmWriteTimer.count==atmWriteTimer.numTimeSteps) { /* time to write out */
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Focean: %d = %10.5f\n", Iter, Focean);CHKERRQ(ierr);
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"pCO2atm: %d = %10.5f\n", Iter, pCO2atm);CHKERRQ(ierr);
//    note: we update the StepTimer later in writeExternalForcing	  
	}
  }
      
  if (useLandModel) {
/*	landsource_(&landState[0],&pCO2atm,&landUseEmission,&deltaTsg,&Fland,&landSource[0]); */ /* returns S and Fl in PgC/y */
/*    time step land */
	for (k=0;k<=2;k++) {
	  landState[k] = landState[k] + atmModelDeltaT*landSource[k];
	}
  }

/* time step atmosphere */
  pCO2atm = pCO2atm + atmModelDeltaT*(fossilFuelEmission + landUseEmission - Focean - Fland)/ppmToPgC;
#endif
#endif /* O_carbon */

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
  
/* Convert to discrete tendency */
  for (itr=0; itr<numTracers; itr++) {
	ierr = VecScale(ut[itr],DeltaT);CHKERRQ(ierr);
  }
    
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "writeExternalForcing"
PetscErrorCode writeExternalForcing(PetscScalar tc, PetscInt iLoop, PetscInt numTracers, Vec *v, Vec *ut)
{

  PetscErrorCode ierr;
  PetscInt ip, nzloc, kl;
  PetscScalar zero = 0.0, one = 1.0;  

/* Note: tc and iLoop are the time and step at the end of the current time step. */

#if defined O_carbon
#if defined O_TMM_interactive_atmosphere
/* write instantaneous atmos model state */
  if (Iter0+iLoop>=atmWriteTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
	if (atmWriteTimer.count==atmWriteTimer.numTimeSteps) { /* time to write out */
//    note: we've already incremented the count above	
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing atmospheric model output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
	  ierr = PetscFPrintf(PETSC_COMM_WORLD,atmfptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);           
	  ierr = writeBinaryScalarData("pCO2atm_output.bin",&pCO2atm,1,PETSC_TRUE);

	  if (useEmissions) {
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Cumulative emissions at time %10.5f, step %d = %10.6f PgC\n", tc, Iter0+iLoop, cumulativeEmissions);CHKERRQ(ierr);
	  }  

	  if (useLandModel) {
/*      write instantaneous land model state */
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing land model output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
		ierr = writeBinaryScalarData("land_state_output.bin",landState,3,PETSC_TRUE);
	  }

	  ierr = updateStepTimer("atm_write_", Iter0+iLoop, &atmWriteTimer);CHKERRQ(ierr);
	  
	}
  }
#endif
#endif /* O_carbon */

  if (calcDiagnostics) {  
	if (Iter0+iLoop>=diagTimer.startTimeStep) { /* start time averaging (note: diagStartTimeStep is ABSOLUTE time step) */  

	  if (diagTimer.count==diagTimer.numTimeSteps) { /* time to write averages to file */
        doAverage=1;
		for (ip=0; ip<lNumProfiles; ip++) {
		  nzloc=lProfileLength[ip];
		  kl=lStartIndices[ip];
		  uvok_diags_accumulate_(&ip, &kl, &nzloc, &diagTimer.count, &doAverage, &debugFlag);
        }
        
		ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing diagnostics time average at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
		ierr = PetscFPrintf(PETSC_COMM_WORLD,diagfptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);           

		if (numDiags2d>0) {
		  for (id2d=0; id2d<numDiags2d; id2d++) {
			uvok_diags2d_copy_(&id2d,&localDiag2davg[id2d][0],Diag2dFile[id2d],&debugFlag);          		  
			ierr = writeProfileSurfaceScalarData(Diag2dFile[id2d],localDiag2davg[id2d],1,diagAppendOutput);		   
		  }  
        }
        
		if (numDiags3d>0) {        
		  for (id3d=0; id3d<numDiags3d; id3d++) {
			uvok_diags3d_copy_(&id3d,&localDiag3davg[id3d][0],Diag3dFile[id3d],&debugFlag);
			ierr = VecSetValues(Diag3davg[id3d],lSize,gIndices,localDiag3davg[id3d],INSERT_VALUES);CHKERRQ(ierr);
			ierr = VecAssemblyBegin(Diag3davg[id3d]);CHKERRQ(ierr);
			ierr = VecAssemblyEnd(Diag3davg[id3d]);CHKERRQ(ierr);
            /* Open file here if first time. The reason for doing it here is that we don't know the file name */
            /* until the first call to uvok_diags3d_copy above */
			if (diagFirstTime) {
		      ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,Diag3dFile[id3d],DIAG_FILE_MODE,&fddiag3dout[id3d]);CHKERRQ(ierr);
			}
			ierr = VecView(Diag3davg[id3d],fddiag3dout[id3d]);CHKERRQ(ierr);	      
		  }  
        }
        
        if (diagFirstTime) {
          diagFirstTime=PETSC_FALSE;
          diagAppendOutput=PETSC_TRUE;
          DIAG_FILE_MODE=FILE_MODE_APPEND;
        }

        ierr = updateStepTimer("diag_", Iter0+iLoop, &diagTimer);CHKERRQ(ierr);          
		doAverage = 0;

		if (diagTimer.haveResetStartTimeStep) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Switching off diagnostics accumulation at step %d\n", Iter0+iLoop);CHKERRQ(ierr);		
  	      uvok_diags_stop_(&debugFlag); /* Switch off diagnostics accumulation until next cycle */
		}
		
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
#if defined O_carbon
#if defined O_TMM_interactive_atmosphere
/* write instantaneous atmos model state */
  ierr = writeBinaryScalarData("pickup_pCO2atm.bin",&pCO2atm,1,PETSC_FALSE);

  if (useLandModel) {
/*   write instantaneous land model state */
	ierr = writeBinaryScalarData("pickup_land_state.bin",landState,3,PETSC_FALSE);
  }

  ierr = PetscFClose(PETSC_COMM_WORLD,atmfptime);CHKERRQ(ierr);  
#endif
#endif /* O_carbon */

  ierr = VecDestroy(&Ts);CHKERRQ(ierr);
  ierr = VecDestroy(&Ss);CHKERRQ(ierr);
#ifdef O_npzd_fe_limitation  
  ierr = VecDestroy(&Fe_dissolved);CHKERRQ(ierr);
  if (periodicBiogeochemForcing) {    
    ierr = destroyPeriodicVec(&Fe_dissolvedp);CHKERRQ(ierr);    
  }
#endif
#ifdef O_npzd_iron
  if (periodicBiogeochemForcing) {    
    ierr = destroyPeriodicArray(&localFe_adepp);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localFe_detr_fluxp);CHKERRQ(ierr);    
  }
  ierr = VecDestroy(&Fe_hydr);CHKERRQ(ierr);
#endif
  
  if (periodicBiogeochemForcing) {    
    ierr = destroyPeriodicVec(&Tsp);CHKERRQ(ierr);
    ierr = destroyPeriodicVec(&Ssp);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localaicep);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localhicep);CHKERRQ(ierr);
    ierr = destroyPeriodicArray(&localhsnop);CHKERRQ(ierr);    
    ierr = destroyPeriodicArray(&localwindp);CHKERRQ(ierr);    
    ierr = destroyPeriodicArray(&localatmospp);CHKERRQ(ierr);    
#ifdef READ_SWRAD    
    ierr = destroyPeriodicArray(&localswradp);CHKERRQ(ierr);
#endif    
  }    

  if (calcDiagnostics) {      
    uvok_diags_finalize_(&debugFlag);
	if (numDiags3d>0) {    
	  ierr = VecDestroyVecs(numDiags3d,&Diag3davg);CHKERRQ(ierr);  
	  for (id3d=0; id3d<numDiags3d; id3d++) {
		ierr = PetscViewerDestroy(&fddiag3dout[id3d]);CHKERRQ(ierr);
	  }
	} 
    ierr = PetscFClose(PETSC_COMM_WORLD,diagfptime);CHKERRQ(ierr);  	   
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
  PetscViewer fd;

  myTime = DeltaT*Iter; /* Iter should start at 0 */

  if (periodicBiogeochemForcing) {   
	ierr = interpPeriodicVector(tc,&Ts,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Tsp,"Ts_");
	ierr = interpPeriodicVector(tc,&Ss,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Ssp,"Ss_");	
#ifdef O_npzd_fe_limitation	
	ierr = interpPeriodicVector(tc,&Fe_dissolved,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,biogeochemTimer.tdp,&Fe_dissolvedp,"Fe_dissolved_");	
#endif	
#ifdef O_npzd_iron
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localFe_adep,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localFe_adepp,"Fe_adep_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localFe_detr_flux,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localFe_detr_fluxp,"Fe_detr_flux_");
#endif
#ifdef READ_SWRAD
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localswradp,"swrad_");
#else
    insolation_(&lNumProfiles,&myTime,&locallatitude[0],&localswrad[0],&localtau[0]);
#endif                                                 
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localaice,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localaicep,"aice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhice,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localhicep,"hice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhsno,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localhsnop,"hsno_");                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localwindp,"wind_");   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                  biogeochemTimer.tdp,&localatmospp,"atmosp_");	
    if (useEmP) {
      ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer.cyclePeriod,biogeochemTimer.numPerPeriod,
                                                    biogeochemTimer.tdp,&localEmPp,"EmP_");                                                      
    }
  }
    
  return 0;
}
