#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "CPP_OPTIONS_CHECK.h"

#if defined O_carbon || defined O_mobi_alk || defined O_mobi_o2 || defined O_mobi
#define BGC
#endif

#if defined BGC || (defined O_PaTh && defined O_PaTh_vflux)
#define NEEDEMP
#endif 

#define MAXDIAGS2d 40
#define MAXDIAGS3d 60

#define rc13std 0.0112372
#define r13Func(dc13) ((dc13*1e-3 + 1)*rc13std)
#define dc13Func(c13,c) (1.e3*(c13/(c-c13)/rc13std - 1.))

// This driver code interfaces the TMM to MOBI. Depending on which options (tracers) 
// are switched on in MOBI_TMM_OPTIONS.h, different forcing data are required. Unless 
// stated otherwise, these can be constant in time, periodic or fully time-dependent.
// ifdef NEEDEMP (see above): EmP
// ifdef BGC (see above): T, S, aice, hice, hsno, wind, atmosp
// ifdef O_mobi: swrad, disch, Fe_adep
//               periodic or constant only: Si_dep
//               constant only: Fe_hydr, Si_hydr
//               in addition to a periodic or time-dependent disch, annual mean 
//               disch and Si_dep are also always read
// ifdef O_PaTh: periodic or constant only: dust_dep
//               constant only: PaTh_pom, PaTh_caco3, PaTh_opal, PaTh_lith

#define READ_SWRAD
#include "petscmat.h"
#include "tmm_petsc_matvec_utils.h"
#include "tmm_timer.h"
#include "tmm_forcing_utils.h"
#include "tmm_profile_utils.h"
#include "tmm.h"
#include "tmm_share.h"
#include "mobi_tmm.h"
#include "mobi_forcing.h"

#define TR1 state->c[0]

static PetscClassId EXTERNALFORCING_CLASSID;

// Common variables
PetscScalar zero = 0.0, one = 1.0;  
PetscScalar ppmToPgC=2.1324;

PetscInt toModel = 1; 
PetscInt fromModel = 2;

PetscScalar daysPerYear, secondsPerYear;

PetscBool periodicBiogeochemForcing = PETSC_FALSE;
PetscBool timeDependentBiogeochemForcing = PETSC_FALSE;
PeriodicTimer biogeochemTimer;
TimeDependentTimer timeDependentBiogeochemTimer;

PetscInt nzmax;
PetscScalar DeltaT;
PetscScalar *drF, *zt;
PetscScalar *localdA;
PetscScalar *locallatitude;
Vec surfVolFrac;
PetscScalar totalOceanSurfaceArea;

// Always need these
Vec Ts,Ss;
PetscScalar *localTs,*localSs;
PetscScalar Sglobavg = 0.0;
PetscScalar *localEmP, EmPglobavg;

#if defined NEEDEMP
PetscBool useEmP = PETSC_TRUE;
PeriodicArray localEmPp;
TimeDependentArray localEmPtd;
#endif

#if defined BGC
PeriodicVec Tsp, Ssp;
TimeDependentVec Tstd, Sstd;
PetscScalar *localwind, *localaice, *localhice, *localhsno, *localatmosp;
PeriodicArray localwindp, localaicep, localhicep, localhsnop, localatmospp;
TimeDependentArray localwindtd, localaicetd, localhicetd, localhsnotd, localatmosptd;
#endif

#if defined O_mobi
PetscScalar *localswrad;
# if defined READ_SWRAD
PeriodicArray localswradp;
TimeDependentArray localswradtd;
# endif
PetscScalar *localdisch, totlocaldisch, globaldisch;
PeriodicArray localdischp;
TimeDependentArray localdischtd;
PetscScalar *localsgbathy;

# if defined O_mobi_iron
PetscScalar *localFe_adep;
PeriodicArray localFe_adepp;
TimeDependentArray localFe_adeptd;
PetscScalar *localFe_hydr;
# endif /* O_mobi_iron */

# if defined O_mobi_silicon
PetscScalar *localSi_dep;
PeriodicArray localSi_depp;
PetscScalar globalSi_dep, totlocalSi_dep;
PetscScalar *localSi_hydr;
# endif /* O_mobi_silicon */

#endif /* O_mobi */

#if defined O_PaTh
PetscScalar *localdust_adep;
PeriodicArray localdust_adepp;
PetscScalar *localPaTh_lith, *wLith;
# if !defined O_mobi
PetscScalar *localPaTh_pom, *localPaTh_caco3, *localPaTh_opal;
PetscScalar *wPOM, *wCaCO3, *wOpal;
# endif
#endif /* O_PaTh */

PetscInt id2d, id3d;

PetscMPIInt myId;
PetscInt debugFlag = 0;

#undef __FUNCT__
#define __FUNCT__ "iniExternalForcing"
PetscErrorCode iniExternalForcing(PetscScalar tc, PetscInt Iter, TMMState state, void *userctx)
{
  PetscErrorCode ierr;
  PetscInt numTracers;
  const char *prefix;
  PetscInt ip, il;
  PetscInt itr;
  PetscViewer fd;
  int fp;
  PetscInt maxValsToRead;  
  PetscBool flg;
  PetscScalar myTime;
  PetscScalar localtotdA;

  static PetscBool registered = PETSC_FALSE;
  static PetscInt efctxId = 0;

  ExternalForcingContext ef;

  MPI_Comm comm = PETSC_COMM_WORLD;
  
  if (!registered) {
    PetscClassIdRegister("ExternalForcing context", &EXTERNALFORCING_CLASSID);
    registered = PETSC_TRUE;
  }
  PetscHeaderCreate(ef, EXTERNALFORCING_CLASSID, "ExternalForcing", "ExternalForcing context", "ExternalForcing", comm, 0, 0); 

  efctxId++;
  ef->efctxId=efctxId;
  ef->stateId=state->stateId;

  PetscContainer ctxcontainer;
  PetscCall(PetscContainerCreate(PetscObjectComm((PetscObject)state), &ctxcontainer));
  PetscCall(PetscContainerSetPointer(ctxcontainer, (void*)ef));
  PetscCall(PetscObjectCompose((PetscObject)state, "external forcing ctx", (PetscObject)ctxcontainer));
  state->extforcctxcontainer = ctxcontainer;
  PetscCall(PetscContainerDestroy(&ctxcontainer));

//   if (ef->stateId != state->id) {
//     SETERRQ(PETSC_COMM_WORLD,1,"ERROR: Mismatch between user context and state!");
//   }
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

#if defined O_carbon
# if !defined O_co2ccn_user
  SETERRQ(PETSC_COMM_WORLD,1,"You must define O_co2ccn_user in MOBI_TMM_OPTIONS.h!");
# endif

# if defined O_co2ccn_data && defined O_TMM_interactive_atmosphere
  SETERRQ(PETSC_COMM_WORLD,1,"Cannot use both O_co2ccn_data and O_TMM_interactive_atmosphere options at the same time!");
# endif

# if defined O_carbon_13_coupled && !defined O_TMM_interactive_atmosphere
  SETERRQ(PETSC_COMM_WORLD,1,"Must define O_TMM_interactive_atmosphere when O_carbon_13_coupled is defined!");
# endif

# if defined O_c13ccn_data && defined O_carbon_13_coupled
  SETERRQ(PETSC_COMM_WORLD,1,"Cannot use both O_c13ccn_data and O_carbon_13_coupled options at the same time!");
# endif

#endif /* O_carbon */

// Now set problem specific data
// Common data (only initialize/read once)
  if (ef->efctxId==1) {

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
  ierr = PetscOptionsHasName(NULL,NULL,"-time_dependent_biogeochem_forcing",&timeDependentBiogeochemForcing);CHKERRQ(ierr);


  if (periodicBiogeochemForcing) {    
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Periodic biogeochemical forcing specified\n");CHKERRQ(ierr);
    ierr = PeriodicTimerCreate(&biogeochemTimer);CHKERRQ(ierr);
    ierr = PeriodicTimerIni("periodic_biogeochem_", NULL, NULL, biogeochemTimer);CHKERRQ(ierr);
  }

  if (timeDependentBiogeochemForcing) {    
    ierr=PetscPrintf(PETSC_COMM_WORLD,"Time-dependent biogeochemical forcing specified\n");CHKERRQ(ierr);
    ierr = TimeDependentTimerCreate(&timeDependentBiogeochemTimer);CHKERRQ(ierr);
    ierr = TimeDependentTimerIni("time_dependent_biogeochem_", NULL, NULL, timeDependentBiogeochemTimer);CHKERRQ(ierr);
  }

#if defined O_sed
  ierr = PetscOptionsGetInt(NULL,NULL,"-num_ocean_steps_per_sed_step",&numOceanStepsPerSedStep,&flg);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of ocean steps per sediment step = %d\n",numOceanStepsPerSedStep);CHKERRQ(ierr);
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localsedmask);CHKERRQ(ierr);   
#endif
  
/* Grid arrays */
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

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localdA);CHKERRQ(ierr);
  ierr = readProfileSurfaceScalarData("dA.bin",localdA,1);  /* cm^2 */
  localtotdA = 0.0;  
  totalOceanSurfaceArea = 0.0;
  for (ip=0; ip<lNumProfiles; ip++) {
	localtotdA = localtotdA + localdA[ip];
  }
  MPI_Allreduce(&localtotdA, &totalOceanSurfaceArea, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);    
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Total ocean surface area = %g [cm^2]\n",totalOceanSurfaceArea);CHKERRQ(ierr);

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&locallatitude);CHKERRQ(ierr);  
  ierr = readProfileSurfaceScalarData("latitude.bin",locallatitude,1);  

  ierr = VecDuplicate(TR1,&surfVolFrac);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"surface_volume_fraction.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = VecLoad(surfVolFrac,fd);CHKERRQ(ierr);  
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);      

/* Forcing fields */
/* Initialize T and S */
  ierr = VecDuplicate(TR1,&Ts);CHKERRQ(ierr);
  ierr = VecDuplicate(TR1,&Ss);CHKERRQ(ierr);
  ierr = VecSet(Ts,zero); CHKERRQ(ierr);
  ierr = VecSet(Ss,zero); CHKERRQ(ierr);
  ierr = VecGetArray(Ts,&localTs);CHKERRQ(ierr);
  ierr = VecGetArray(Ss,&localSs);CHKERRQ(ierr);

/* Initialize EmP */
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localEmP);CHKERRQ(ierr);
  for (ip=0; ip<lNumProfiles; ip++) { /* initialize to zero to be safe */
    localEmP[ip]=0.0;
  }
  EmPglobavg = 0.0; /* set this to zero for now */  

#if defined NEEDEMP
  ierr = PetscOptionsHasName(NULL,NULL,"-no_use_emp",&flg);CHKERRQ(ierr);
  if (flg) {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: -no_use_emp has been specified. E-P will be set to zero.\n");CHKERRQ(ierr);
	useEmP = PETSC_FALSE;
  }  
  
  if (useEmP) {
	if (timeDependentBiogeochemForcing) {
    ierr = TimeDependentArrayCreate(&localEmPtd, lNumProfiles);  
	} else if (periodicBiogeochemForcing) {    
    ierr = PeriodicArrayCreate(&localEmPp, lNumProfiles);  
	} else {  
	  ierr = readProfileSurfaceScalarData("EmP.bin",localEmP,1);  
	  ierr = dotProdProfileSurfaceScalarData(localEmP,localdA,&EmPglobavg);
	  EmPglobavg=EmPglobavg/totalOceanSurfaceArea;
	}
  }
#endif

#if defined O_mobi
  ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localsgbathy);CHKERRQ(ierr);    
  ierr = VecLoadVecIntoArray(TR1,"sgbathy.petsc",localsgbathy);CHKERRQ(ierr);

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localswrad);CHKERRQ(ierr);  
# if defined READ_SWRAD
  if (timeDependentBiogeochemForcing) {
    ierr = TimeDependentArrayCreate(&localswradtd, lNumProfiles);  
  } else if (periodicBiogeochemForcing) {   
    ierr = PeriodicArrayCreate(&localswradp, lNumProfiles);  
  } else {  
    ierr = readProfileSurfaceScalarData("swrad.bin",localswrad,1);  
  }
# endif	/* READ_SWRAD */

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localdisch);CHKERRQ(ierr);  /* g/cm^2/s */
  if (timeDependentBiogeochemForcing) {
    ierr = TimeDependentArrayCreate(&localdischtd, lNumProfiles);  
  } else if (periodicBiogeochemForcing) {
    ierr = PeriodicArrayCreate(&localdischp, lNumProfiles);  
  }
// We always read the annual mean discharge to calculate the global annual mean discharge. 
// This avoids having to do an MPI_Allreduce at each time step.
  ierr = readProfileSurfaceScalarData("disch.bin",localdisch,1); 
  totlocaldisch = 0.0;
  for (ip=0; ip<lNumProfiles; ip++) {
	totlocaldisch = totlocaldisch + localdisch[ip]*localdA[ip]; /* g/s */
  }    
  MPI_Allreduce(&totlocaldisch, &globaldisch, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Global annual mean discharge = %g [g/s]\n", globaldisch);CHKERRQ(ierr);    

# if defined O_mobi_iron
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localFe_adep);CHKERRQ(ierr);  
  if (timeDependentBiogeochemForcing) {    
    ierr = TimeDependentArrayCreate(&localFe_adeptd, lNumProfiles);  
  } else if (periodicBiogeochemForcing) {    
    ierr = PeriodicArrayCreate(&localFe_adepp, lNumProfiles);  
  } else {  
    ierr = readProfileSurfaceScalarData("Fe_adep.bin",localFe_adep,1);  
  }

  ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localFe_hydr);CHKERRQ(ierr);
  ierr = VecLoadVecIntoArray(TR1,"Fe_hydr.petsc",localFe_hydr);CHKERRQ(ierr);
# endif /* O_mobi_iron */

# if defined O_mobi_silicon
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localSi_dep);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    ierr = PeriodicArrayCreate(&localSi_depp, lNumProfiles);  
  }
// We always read the annual mean deposition to calculate the global annual mean deposition. 
// This avoids having to do an MPI_Allreduce at each time step.
  ierr = readProfileSurfaceScalarData("Si_dep.bin",localSi_dep,1);  
  totlocalSi_dep = 0.0;
  for (ip=0; ip<lNumProfiles; ip++) {
	totlocalSi_dep = totlocalSi_dep + localSi_dep[ip]*localdA[ip]; /* units ??? */
  }    
  MPI_Allreduce(&totlocalSi_dep, &globalSi_dep, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Global annual silicon deposition = %g [units ???]\n", globalSi_dep);CHKERRQ(ierr);    

  ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localSi_hydr);CHKERRQ(ierr);
  ierr = VecLoadVecIntoArray(TR1,"Si_hydr.petsc",localSi_hydr);CHKERRQ(ierr);
# endif /* O_mobi_silicon */

#endif /* O_mobi */

#if defined O_PaTh
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localdust_adep);CHKERRQ(ierr);  
  if (periodicBiogeochemForcing) {    
    ierr = PeriodicArrayCreate(&localdust_adepp, lNumProfiles);  
  } else {  
    ierr = readProfileSurfaceScalarData("dust_dep.bin",localdust_adep,1);
  }

  ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localPaTh_lith);CHKERRQ(ierr);
  ierr = VecLoadVecIntoArray(TR1,"PaTh_lith.petsc",localPaTh_lith);CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"wLith.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscMalloc(nzmax*sizeof(PetscScalar),&wLith);CHKERRQ(ierr); 
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
//   ierr = PetscBinaryRead(fp,&nzmax,1,NULL,PETSC_INT);CHKERRQ(ierr);  
  ierr = PetscBinaryRead(fp,wLith,nzmax,NULL,PETSC_SCALAR);CHKERRQ(ierr);  
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

# if !defined O_mobi
  ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localPaTh_pom);CHKERRQ(ierr);
  ierr = VecLoadVecIntoArray(TR1,"PaTh_pom.petsc",localPaTh_pom);CHKERRQ(ierr);

  ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localPaTh_caco3);CHKERRQ(ierr);
  ierr = VecLoadVecIntoArray(TR1,"PaTh_caco3.petsc",localPaTh_caco3);CHKERRQ(ierr);

  ierr = PetscMalloc(lSize*sizeof(PetscScalar),&localPaTh_opal);CHKERRQ(ierr);
  ierr = VecLoadVecIntoArray(TR1,"PaTh_opal.petsc",localPaTh_opal);CHKERRQ(ierr);
  
  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"wPOM.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscMalloc(nzmax*sizeof(PetscScalar),&wPOM);CHKERRQ(ierr); 
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
//   ierr = PetscBinaryRead(fp,&nzmax,1,NULL,PETSC_INT);CHKERRQ(ierr);  
  ierr = PetscBinaryRead(fp,wPOM,nzmax,NULL,PETSC_SCALAR);CHKERRQ(ierr);  
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"wCaCO3.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscMalloc(nzmax*sizeof(PetscScalar),&wCaCO3);CHKERRQ(ierr); 
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
//   ierr = PetscBinaryRead(fp,&nzmax,1,NULL,PETSC_INT);CHKERRQ(ierr);  
  ierr = PetscBinaryRead(fp,wCaCO3,nzmax,NULL,PETSC_SCALAR);CHKERRQ(ierr);  
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"wOpal.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
  ierr = PetscMalloc(nzmax*sizeof(PetscScalar),&wOpal);CHKERRQ(ierr); 
  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
//   ierr = PetscBinaryRead(fp,&nzmax,1,NULL,PETSC_INT);CHKERRQ(ierr);  
  ierr = PetscBinaryRead(fp,wOpal,nzmax,NULL,PETSC_SCALAR);CHKERRQ(ierr);  
  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
  
# endif  
#endif /* O_PaTh */

#if defined BGC
  if (timeDependentBiogeochemForcing) {
    TimeDependentVecCreate(&Tstd);
    TimeDependentVecCreate(&Sstd);
  } else if (periodicBiogeochemForcing) {
	   PeriodicVecCreate(&Tsp);
 	  PeriodicVecCreate(&Ssp);
  } else {
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Ts.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr = VecLoad(Ts,fd);CHKERRQ(ierr);  
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"Ss.petsc",FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr = VecLoad(Ss,fd);CHKERRQ(ierr);    
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Done reading T/S\n");CHKERRQ(ierr);	
  }  

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localaice);CHKERRQ(ierr);  
  if (timeDependentBiogeochemForcing) {
    ierr = TimeDependentArrayCreate(&localaicetd, lNumProfiles);  
  } else if (periodicBiogeochemForcing) {
    ierr = PeriodicArrayCreate(&localaicep, lNumProfiles);  
  } else {  
    ierr = readProfileSurfaceScalarData("aice.bin",localaice,1);  
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localhice);CHKERRQ(ierr);  
  if (timeDependentBiogeochemForcing) {    
    ierr = TimeDependentArrayCreate(&localhicetd, lNumProfiles);  
  } else if (periodicBiogeochemForcing) {   
    ierr = PeriodicArrayCreate(&localhicep, lNumProfiles);  
  } else {  
    ierr = readProfileSurfaceScalarData("hice.bin",localhice,1);  
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localhsno);CHKERRQ(ierr);  
  if (timeDependentBiogeochemForcing) {  
    ierr = TimeDependentArrayCreate(&localhsnotd, lNumProfiles);  
  } else if (periodicBiogeochemForcing) {   
    ierr = PeriodicArrayCreate(&localhsnop, lNumProfiles);  
  } else {  
    ierr = readProfileSurfaceScalarData("hsno.bin",localhsno,1);  
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localwind);CHKERRQ(ierr);  
  if (timeDependentBiogeochemForcing) {    
    ierr = TimeDependentArrayCreate(&localwindtd, lNumProfiles);  
  } else if (periodicBiogeochemForcing) { 
    ierr = PeriodicArrayCreate(&localwindp, lNumProfiles);  
  } else {  
    ierr = readProfileSurfaceScalarData("wind.bin",localwind,1);  
  }

  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&localatmosp);CHKERRQ(ierr);  
  if (timeDependentBiogeochemForcing) { 
    ierr = TimeDependentArrayCreate(&localatmosptd, lNumProfiles);  
  } else if (periodicBiogeochemForcing) {  
    ierr = PeriodicArrayCreate(&localatmospp, lNumProfiles);  
  } else {  
    ierr = readProfileSurfaceScalarData("atmosp.bin",localatmosp,1);  
  }

#endif /* BGC */

  } /* end common data */
  
  for (itr=0; itr<numTracers; itr++) {
    ierr = VecSet(state->qef[itr],zero); CHKERRQ(ierr);
  }

  ierr = VecGetArrays(state->c,numTracers,&ef->localTR);CHKERRQ(ierr);
  ierr = VecGetArrays(state->qef,numTracers,&ef->localJTR);CHKERRQ(ierr);
  
/* Initialize biogeochem model */
  myTime = DeltaT*Iter; /* Iter should start at 0 */

#if defined NEEDEMP
  if (useEmP) {
	if (timeDependentBiogeochemForcing) {    
	  ierr = interpTimeDependentProfileSurfaceScalarData(tc,localEmP,timeDependentBiogeochemTimer->numTimes,
														 timeDependentBiogeochemTimer->tdt,localEmPtd,"EmPtd.bin");
	  ierr = dotProdProfileSurfaceScalarData(localEmP,localdA,&EmPglobavg);
	  EmPglobavg=EmPglobavg/totalOceanSurfaceArea;														 
	} else if (periodicBiogeochemForcing) {   
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
													biogeochemTimer->tdp,localEmPp,"EmP_");
	  ierr = dotProdProfileSurfaceScalarData(localEmP,localdA,&EmPglobavg);
	  EmPglobavg=EmPglobavg/totalOceanSurfaceArea;													
	}
  }	
#endif
  
#if defined O_mobi
# if defined READ_SWRAD
  if (timeDependentBiogeochemForcing) {    
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localswrad,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localswradtd,"swradtd.bin");
  } else if (periodicBiogeochemForcing) {   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localswradp,"swrad_");
  }
# else
   insolation_(&lNumProfiles,&myTime,&locallatitude[0],&localswrad[0],&localtau[0]);
# endif /* READ_SWRAD */

  if (timeDependentBiogeochemForcing) {    
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localdisch,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localdischtd,"dischtd.bin");
  } else if (periodicBiogeochemForcing) {   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localdisch,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localdischp,"disch_");
  }
  	
# if defined O_mobi_iron
  if (timeDependentBiogeochemForcing) {    
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localFe_adep,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localFe_adeptd,"Fe_adeptd.bin");
  } else if (periodicBiogeochemForcing) {   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localFe_adep,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localFe_adepp,"Fe_adep_");
  }                                                
# endif /* O_mobi_iron */

# if defined O_mobi_silicon
  if (periodicBiogeochemForcing) {
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localSi_dep,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localSi_depp,"Si_dep_");
  }
# endif /* O_mobi_silicon */
#endif /* O_mobi */

#if defined O_PaTh
  if (periodicBiogeochemForcing) {
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localdust_adep,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localdust_adepp,"dust_dep_");
  }
#endif /* O_PaTh */

#if defined BGC
  if (timeDependentBiogeochemForcing) {    
    ierr = TimeDependentVecInterp(tc,&Ts,timeDependentBiogeochemTimer->numTimes,timeDependentBiogeochemTimer->tdt,Tstd,"Tstd.petsc");
    ierr = TimeDependentVecInterp(tc,&Ss,timeDependentBiogeochemTimer->numTimes,timeDependentBiogeochemTimer->tdt,Sstd,"Sstd.petsc");
  
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localaice,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localaicetd,"aicetd.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localhice,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localhicetd,"hicetd.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localhsno,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localhsnotd,"hsnotd.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localwind,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localwindtd,"windtd.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localatmosp,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localatmosptd,"atmosptd.bin");
  } else if (periodicBiogeochemForcing) {   
    ierr = PeriodicVecInterp(tc,&Ts,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,Tsp,"Ts_");
    ierr = PeriodicVecInterp(tc,&Ss,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,Ssp,"Ss_");	
  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localaice,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localaicep,"aice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhice,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localhicep,"hice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhsno,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localhsnop,"hsno_");                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localwindp,"wind_");   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localatmospp,"atmosp_");	
  }
#endif /* BGC */  
                                                			                                                                                				                              
/* compute global means */
  ierr = PetscMalloc(numTracers*sizeof(PetscScalar),&ef->TRglobavg);CHKERRQ(ierr);	  
  ierr = VecDot(surfVolFrac,Ss,&Sglobavg);CHKERRQ(ierr); /* volume weighted mean surface salinity */   
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Global average of salinity =%10.5f\n", Sglobavg);CHKERRQ(ierr); 
  for (itr=0; itr<numTracers; itr++) {    
	   ierr = VecDot(surfVolFrac,state->c[itr],&ef->TRglobavg[itr]);CHKERRQ(ierr); /* volume weighted mean surface TR */
   	ierr = PetscPrintf(PETSC_COMM_WORLD,"Global average of tracer %d = %10.5f\n", itr, ef->TRglobavg[itr]);CHKERRQ(ierr);									                    
  }  

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&myId);CHKERRQ(ierr);
  if (myId == 0) debugFlag = 1;

  mobi_ini_(&numTracers,&lSize,&lNumProfiles,&nzmax,lProfileLength,
            zt,drF,&DeltaT,locallatitude,localdA,
            &Sglobavg,ef->TRglobavg,
#if defined O_carbon
            &ef->pCO2atm,
# if defined O_carbon_13
            &ef->dc13atm,
# endif
# if defined O_carbon_14
            &ef->DC14atm,
# endif
#endif /* O_carbon */
#if defined O_mobi
            localsgbathy,
# if defined O_mobi_iron
			localFe_hydr,
# endif
# if defined O_mobi_silicon
			localSi_hydr,
# endif
#endif
#if defined O_PaTh
			localPaTh_lith, wLith, 
# if !defined O_mobi
			localPaTh_pom, localPaTh_caco3, localPaTh_opal, 
			wPOM, wCaCO3, wOpal, 
# endif
#endif
#if defined O_sed
            &numOceanStepsPerSedStep,
            &nzmaxSed,&ibmaxSed,&numSedMixedTracers,
            &numSedBuriedTracers,&globalweathflx,
            &localsedsa,localsedmask,
#endif /* O_sed */
            &debugFlag);

#if defined O_carbon
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&ef->localgasexflux);CHKERRQ(ierr);
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&ef->localtotflux);CHKERRQ(ierr);
  for (ip=0; ip<lNumProfiles; ip++) {
     ef->localgasexflux[ip]=0.0;
     ef->localtotflux[ip]=0.0;
  }

# if defined O_co2ccn_data
  /* prescribed atmospheric CO2 */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using prescribed atmospheric CO2\n");CHKERRQ(ierr);
  maxValsToRead = 2;
  ef->pCO2atmFiles[0] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* time file */
  ef->pCO2atmFiles[1] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* atmospheric pCO2 history file */
  ierr = PetscOptionsGetStringArray(NULL,NULL,"-pco2atm_history",ef->pCO2atmFiles,&maxValsToRead,&flg);CHKERRQ(ierr);
  if (flg) { /* Read atmospheric pCO2 history */
	if (maxValsToRead != 2) {
	  SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of file names specified for atmospheric pCO2 history");
	}      
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading time-dependent atmospheric pCO2 history\n");CHKERRQ(ierr);      
	/* read time data */
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,ef->pCO2atmFiles[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscBinaryRead(fp,&ef->numpCO2atm_hist,1,NULL,PETSC_INT);CHKERRQ(ierr);  
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of points in atmospheric pCO2 history file is %d \n",ef->numpCO2atm_hist);CHKERRQ(ierr);  
	ierr = PetscMalloc(ef->numpCO2atm_hist*sizeof(PetscScalar),&ef->TpCO2atm_hist);CHKERRQ(ierr); 
	ierr = PetscBinaryRead(fp,ef->TpCO2atm_hist,ef->numpCO2atm_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	/* read atmospheric pCO2 data */
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,ef->pCO2atmFiles[1],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscMalloc(ef->numpCO2atm_hist*sizeof(PetscScalar),&ef->pCO2atm_hist);CHKERRQ(ierr); 
	ierr = PetscBinaryRead(fp,ef->pCO2atm_hist,ef->numpCO2atm_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	
	ef->pCO2atm = ef->pCO2atm_hist[0];

  }	else {
	SETERRQ(PETSC_COMM_WORLD,1,"Must specify atmospheric CO2 history with -pco2atm_history when using the O_co2ccn_data option");
  }    
# endif /* O_co2ccn_data */

# if defined O_TMM_interactive_atmosphere
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using interactive atmospheric model\n");CHKERRQ(ierr);  

/* overwrite namelist value */
  ierr = PetscOptionsGetReal(NULL,NULL,"-pco2atm_ini",&ef->pCO2atm,&flg);CHKERRQ(ierr); /* read from command line */
  if (!flg) {
	ierr = PetscOptionsGetString(NULL,NULL,"-pco2atm_ini_file",ef->pCO2atmIniFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
	if (flg) { /* read from binary file */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,ef->pCO2atmIniFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	  ierr = PetscBinaryRead(fp,&ef->pCO2atm,1,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	}
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using initial atmospheric pCO2 of %g ppm\n",ef->pCO2atm);CHKERRQ(ierr);

#  if defined O_carbon_13_coupled
/* overwrite namelist value */
  ierr = PetscOptionsGetReal(NULL,NULL,"-dc13atm_ini",&ef->dc13atm,&flg);CHKERRQ(ierr); /* read from command line */
  if (!flg) {
	ierr = PetscOptionsGetString(NULL,NULL,"-dc13atm_ini_file",ef->dc13atmIniFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
	if (flg) { /* read from binary file */
	  ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,ef->dc13atmIniFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
	  ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	  ierr = PetscBinaryRead(fp,&ef->dc13atm,1,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	  ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	}
  }
  pC13O2atm = (r13Func(ef->dc13atm)/(1+r13Func(ef->dc13atm)))*ef->pCO2atm; /* mole fraction of atmospheric 13CO2; this is the prognostic (time-stepped) variable */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using initial atmospheric pC13O2 (d13C) of %g ppm (%g permille)\n",ef->pC13O2atm,ef->dc13atm);CHKERRQ(ierr);
  ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&ef->localc13gasexflux);CHKERRQ(ierr);
  for (ip=0; ip<lNumProfiles; ip++) {
     ef->localc13gasexflux[ip]=0.0;
  }  
#  endif /* O_carbon_13_coupled */

  atmModelDeltaT = DeltaT/secondsPerYear; /* time step in years */

  ierr = StepTimerCreate(&ef->atmWriteTimer);CHKERRQ(ierr);
 	ierr = StepTimerIni("atm_write_", prefix, Iter0, ef->atmWriteTimer);CHKERRQ(ierr);

  ierr = PetscOptionsHasName(NULL,NULL,"-atm_append",&ef->atmAppendOutput);CHKERRQ(ierr);
  if (ef->atmAppendOutput) {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output will be appended\n");CHKERRQ(ierr);
  } else {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output will overwrite existing file(s)\n");CHKERRQ(ierr);
  }    

/* Output times */
  ierr = PetscOptionsGetString(NULL,NULL,"-atm_time_file",ef->atmOutTimeFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) {
	strcpy(ef->atmOutTimeFile,"");
	sprintf(ef->atmOutTimeFile,"%s","atm_output_time.txt");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output times will be written to %s\n",ef->atmOutTimeFile);CHKERRQ(ierr);

  if (!ef->atmAppendOutput) {
    ierr = PetscFOpen(PETSC_COMM_WORLD,ef->atmOutTimeFile,"w",&ef->atmfptime);CHKERRQ(ierr);  
    if (Iter0==(ef->atmWriteTimer->startTimeStep)) { /* note: startTimeStep is ABSOLUTE time step */
      ierr = PetscFPrintf(PETSC_COMM_WORLD,atmfptime,"%d   %10.5f\n",Iter0,time0);CHKERRQ(ierr);     
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing atmospheric output at time %10.5f, step %d\n", tc,Iter);CHKERRQ(ierr);  
      ierr = writeBinaryScalarData("pCO2atm_output.bin",&ef->pCO2atm,1,PETSC_FALSE);
#  if defined O_carbon_13_coupled
      ierr = writeBinaryScalarData("pC13O2atm_output.bin",&ef->pC13O2atm,1,PETSC_FALSE);
      ierr = writeBinaryScalarData("dc13atm_output.bin",&ef->dc13atm,1,PETSC_FALSE);
#  endif	
    }
  } else {
    ierr = PetscFOpen(PETSC_COMM_WORLD,ef->atmOutTimeFile,"a",&atmfptime);CHKERRQ(ierr);  
    if (Iter0==(ef->atmWriteTimer->startTimeStep)) { /* note: startTimeStep is ABSOLUTE time step */      	
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Atmospheric model output will be appended. Initial condition will NOT be written\n");CHKERRQ(ierr);      
    }  
  }

  ierr = PetscOptionsHasName(NULL,NULL,"-use_land_model",&ef->useLandModel);CHKERRQ(ierr);
  if (ef->useLandModel) {
    SETERRQ(PETSC_COMM_WORLD,1,"Land model not yet supported!");
  
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Using interactive land model\n");CHKERRQ(ierr);      
    ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,"land_ini.bin",FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
    ierr = PetscBinaryRead(fp,ef->landState,3,NULL,PETSC_SCALAR);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);    

    if (!ef->atmAppendOutput) {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing land output at time %10.5f, step %d\n", tc,Iter);CHKERRQ(ierr);  
      ierr = writeBinaryScalarData("land_state_output.bin",ef->landState,3,PETSC_FALSE);
    } else {
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Land model output will be appended. Initial condition will NOT be written\n");CHKERRQ(ierr);      
    }
  }

  /* CO2 emissions */
  maxValsToRead = 3;
  ef->emFiles[0] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* time file */
  ef->emFiles[1] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* fossil fuel emissions file */
  ef->emFiles[2] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* land use emissions file */
  ierr = PetscOptionsGetStringArray(NULL,NULL,"-emissions_history",ef->emFiles,&maxValsToRead,&ef->useEmissions);CHKERRQ(ierr);
  if (ef->useEmissions) { /* Read emissions history */
	if (maxValsToRead != 3) {
	  SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of file names specified for emissions");
	}      
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Using prescribed emissions\n");CHKERRQ(ierr);     
	/* read time data */
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,ef->emFiles[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscBinaryRead(fp,&ef->numEmission_hist,1,NULL,PETSC_INT);CHKERRQ(ierr);  
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of points in CO2 emission files is %d \n",ef->numEmission_hist);CHKERRQ(ierr);  
	ierr = PetscMalloc(ef->numEmission_hist*sizeof(PetscScalar),&ef->Tem_hist);CHKERRQ(ierr); 
	ierr = PetscBinaryRead(fp,ef->Tem_hist,ef->numEmission_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	/* read fossil fuel emissions */
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,ef->emFiles[1],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscMalloc(ef->numEmission_hist*sizeof(PetscScalar),&ef->E_hist);CHKERRQ(ierr); 
	ierr = PetscBinaryRead(fp,ef->E_hist,ef->numEmission_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	/* read land use emissions */
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,ef->emFiles[2],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscMalloc(ef->numEmission_hist*sizeof(PetscScalar),&ef->D_hist);CHKERRQ(ierr); 
	ierr = PetscBinaryRead(fp,ef->D_hist,numEmission_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	ierr = PetscOptionsHasName(NULL,NULL,"-interp_emissions",&ef->interpEmissions);CHKERRQ(ierr);
	if (ef->interpEmissions) {      
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: Emissions will be interpolated in time. If you're prescribing annual emissions\n");CHKERRQ(ierr);     
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"         interpolation may lead to a different net emission input than what is prescribed\n");CHKERRQ(ierr);     
	} else {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Emissions will NOT be interpolated in time. It is assumed that you're prescribing annual emissions and that\n");CHKERRQ(ierr);     
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"the time data in file %s are the beginning of the year for which the corresponding emission is prescribed.\n",emFiles[0]);CHKERRQ(ierr);     
	}
  }  
# endif /* O_TMM_interactive_atmosphere */

# if defined O_c14ccn_data
  /* prescribed atmospheric DeltaC14 */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using prescribed atmospheric Delta C14\n");CHKERRQ(ierr);
  ef->timeDependentAtmosphericC14 = PETSC_FALSE;  
  maxValsToRead = 4;
  ef->C14atmFiles[0] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* time file */
  ef->C14atmFiles[1] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* NH atmospheric DeltaC14 history file */
  ef->C14atmFiles[2] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* EQ atmospheric DeltaC14 history file */
  ef->C14atmFiles[3] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* SH atmospheric DeltaC14 history file */
  ierr = PetscOptionsGetStringArray(NULL,NULL,"-c14atm_history",ef->C14atmFiles,&maxValsToRead,&ef->timeDependentAtmosphericC14);CHKERRQ(ierr);
  if (ef->timeDependentAtmosphericC14) { /* Read atmospheric Delta C14 history */
	if (maxValsToRead != 4) {
	  SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of file names specified for atmospheric Delta C14 history");
	}      
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading time-dependent atmospheric Delta C14 history\n");CHKERRQ(ierr);      
	/* read time data */
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,ef->C14atmFiles[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscBinaryRead(fp,&ef->numC14atm_hist,1,NULL,PETSC_INT);CHKERRQ(ierr);  
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of points in atmospheric Delta C14 history files is %d \n",ef->numC14atm_hist);CHKERRQ(ierr);  
	ierr = PetscMalloc(ef->numC14atm_hist*sizeof(PetscScalar),&ef->TC14atm_hist);CHKERRQ(ierr); 
	ierr = PetscBinaryRead(fp,ef->TC14atm_hist,ef->numC14atm_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	/* read atmospheric Delta C14 data */
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading NH atmospheric Delta C14 history from %s\n",ef->C14atmFiles[1]);CHKERRQ(ierr);      
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,ef->C14atmFiles[1],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscMalloc(ef->numC14atm_hist*sizeof(PetscScalar),&ef->C14atmnh_hist);CHKERRQ(ierr); 
	ierr = PetscBinaryRead(fp,ef->C14atmnh_hist,ef->numC14atm_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	dc14ccnnhatm = C14atmnh_hist[0];
//
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading EQ atmospheric Delta C14 history from %s\n",ef->C14atmFiles[2]);CHKERRQ(ierr);
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,ef->C14atmFiles[2],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscMalloc(ef->numC14atm_hist*sizeof(PetscScalar),&ef->C14atmeq_hist);CHKERRQ(ierr); 
	ierr = PetscBinaryRead(fp,ef->C14atmeq_hist,ef->numC14atm_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	ef->dc14ccneqatm = ef->C14atmeq_hist[0];
//
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading SH atmospheric Delta C14 history from %s\n",ef->C14atmFiles[3]);CHKERRQ(ierr);      
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,ef->C14atmFiles[3],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscMalloc(ef->numC14atm_hist*sizeof(PetscScalar),&ef->C14atmsh_hist);CHKERRQ(ierr); 
	ierr = PetscBinaryRead(fp,ef->C14atmsh_hist,ef->numC14atm_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	ef->dc14ccnshatm = ef->C14atmsh_hist[0];	
  }	else {
    /* No atmospheric DeltaC14 history specified */  
    ef->dc14atmvals[0]=0.0;
    ef->dc14atmvals[1]=0.0;
    ef->dc14atmvals[2]=0.0;    
    maxValsToRead = 3;
    ierr = PetscOptionsGetRealArray(NULL,NULL,"-c14atm_vals",ef->dc14atmvals,&maxValsToRead,&flg);
    ef->dc14ccnnhatm=ef->dc14atmvals[0];
    ef->dc14ccneqatm=ef->dc14atmvals[1];
    ef->dc14ccnshatm=ef->dc14atmvals[2];
    if (flg) {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using fixed atmospheric Delta C14 values of %g (NH), %g (EQ) and %g (SH)\n",ef->dc14ccnnhatm,ef->dc14ccneqatm,ef->dc14ccnshatm);CHKERRQ(ierr);
    } else {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Option O_c14ccn_data specified but no atmospheric Delta C14 history or values given\n");CHKERRQ(ierr);      
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using default fixed atmospheric Delta C14 values of %g (NH), %g (EQ) and %g (SH)\n",ef->dc14ccnnhatm,ef->dc14ccneqatm,ef->dc14ccnshatm);CHKERRQ(ierr);
    }
  }  
# endif /* O_c14ccn_data */

# if defined O_c13ccn_data
  /* prescribed atmospheric deltaC13 */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using prescribed atmospheric delta C13\n");CHKERRQ(ierr);
  maxValsToRead = 2;
  ef->C13atmFiles[0] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* time file */
  ef->C13atmFiles[1] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char)); /* atmospheric deltaC13 history file */
  ierr = PetscOptionsGetStringArray(NULL,NULL,"-c13atm_history",ef->C13atmFiles,&maxValsToRead,&flg);CHKERRQ(ierr);
  if (flg) { /* Read atmospheric delta C13 history */
	if (maxValsToRead != 2) {
	  SETERRQ(PETSC_COMM_WORLD,1,"Insufficient number of file names specified for atmospheric delta C13 history");
	}      
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading time-dependent atmospheric delta C13 history\n");CHKERRQ(ierr);      
	/* read time data */
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,ef->C13atmFiles[0],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscBinaryRead(fp,&ef->numC13atm_hist,1,NULL,PETSC_INT);CHKERRQ(ierr);  
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of points in atmospheric delta C13 history file is %d \n",ef->numC13atm_hist);CHKERRQ(ierr);  
	ierr = PetscMalloc(ef->numC13atm_hist*sizeof(PetscScalar),&ef->TC13atm_hist);CHKERRQ(ierr); 
	ierr = PetscBinaryRead(fp,ef->TC13atm_hist,ef->numC13atm_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	/* read atmospheric delta C13 data */
	ierr = PetscViewerBinaryOpen(PETSC_COMM_SELF,ef->C13atmFiles[1],FILE_MODE_READ,&fd);CHKERRQ(ierr);
	ierr = PetscViewerBinaryGetDescriptor(fd,&fp);CHKERRQ(ierr);
	ierr = PetscMalloc(ef->numC13atm_hist*sizeof(PetscScalar),&ef->C13atm_hist);CHKERRQ(ierr); 
	ierr = PetscBinaryRead(fp,ef->C13atm_hist,ef->numC13atm_hist,NULL,PETSC_SCALAR);CHKERRQ(ierr);
	ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);
	
	ef->dc13atm = ef->C13atm_hist[0];
	ef->pC13O2atm = (r13Func(ef->dc13atm)/(1+r13Func(ef->dc13atm)))*ef->pCO2atm; /* mole fraction of atmospheric 13CO2 */
	
  }	else {
	SETERRQ(PETSC_COMM_WORLD,1,"Must specify atmospheric delta C13 history with -c13atm_history when using the O_c13ccn_data option");		
  }  
# endif /* O_c13ccn_data */

# if !defined O_co2ccn_data && !defined O_TMM_interactive_atmosphere
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using a fixed, namelist value of atmospheric CO2 of %g ppm\n",ef->pCO2atm);CHKERRQ(ierr);
# endif

# if defined O_carbon_13
# if !defined O_c13ccn_data && !defined O_carbon_13_coupled
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Using a fixed, namelist value of delta C13 of %g permille\n",ef->dc13atm);CHKERRQ(ierr);
# endif
# endif

#endif /* O_carbon */

#if defined O_sed
  MPI_Allreduce(&localsedsa, &globalsedsa, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Global sediment area = %g [cm^2]\n", globalsedsa);CHKERRQ(ierr); 

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Global weathering flux is %g umol C/s\n",globalweathflx);CHKERRQ(ierr);  

  ierr = writeProfileSurfaceScalarData("sedmask.bin",localsedmask,1,PETSC_FALSE);		   

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Sediment grid for mixed tracers is (nzmax x numSedMixedTracers): %d x %d\n",nzmaxSed,numSedMixedTracers);CHKERRQ(ierr);  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Sediment grid for buried tracers is (ibmax x numSedBuriedTracers): %d x %d\n",ibmaxSed,numSedBuriedTracers);CHKERRQ(ierr);

  sedMixedBlockSize = nzmaxSed*numSedMixedTracers;
  lSedMixedSize = lNumProfiles*sedMixedBlockSize;
  ierr = VecCreate(PETSC_COMM_WORLD,&sedMixedTracers);CHKERRQ(ierr);
  ierr = VecSetSizes(sedMixedTracers,lSedMixedSize,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(sedMixedTracers);CHKERRQ(ierr);
/*   Compute global indices for local piece of sediment vector */
  ierr = VecGetOwnershipRange(sedMixedTracers,&gSedLow,NULL);CHKERRQ(ierr);
  ierr = PetscMalloc(lSedMixedSize*sizeof(PetscInt),&gSedMixedIndices);CHKERRQ(ierr);  
  for (il=0; il<lSedMixedSize; il++) {
    gSedMixedIndices[il] = il + gSedLow;
  }  

  ierr = VecGetArray(sedMixedTracers,&localSedMixedTR);CHKERRQ(ierr);

  sedBuriedBlockSize = ibmaxSed*numSedBuriedTracers;
  lSedBuriedSize = lNumProfiles*sedBuriedBlockSize;
  ierr = VecCreate(PETSC_COMM_WORLD,&sedBuriedTracers);CHKERRQ(ierr);
  ierr = VecSetSizes(sedBuriedTracers,lSedBuriedSize,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = VecSetFromOptions(sedBuriedTracers);CHKERRQ(ierr);
/*   Compute global indices for local piece of sediment vector */
  ierr = VecGetOwnershipRange(sedBuriedTracers,&gSedLow,NULL);CHKERRQ(ierr);
  ierr = PetscMalloc(lSedBuriedSize*sizeof(PetscInt),&gSedBuriedIndices);CHKERRQ(ierr);  
  for (il=0; il<lSedBuriedSize; il++) {
    gSedBuriedIndices[il] = il + gSedLow;
  }  

  ierr = VecGetArray(sedBuriedTracers,&localSedBuriedTR);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished initializing sediment arrays\n");CHKERRQ(ierr);  

  mobi_sed_copy_data_(&lNumProfiles,&localSedMixedTR[0],&localSedBuriedTR[0],&fromModel);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished copying sediment arrays\n");CHKERRQ(ierr);  
  
//   ierr = VecSetValues(sedMixedTracers,lSedMixedSize,gSedMixedIndices,localSedMixedTR,INSERT_VALUES);CHKERRQ(ierr);
//   ierr = VecAssemblyBegin(sedMixedTracers);CHKERRQ(ierr);
//   ierr = VecAssemblyEnd(sedMixedTracers);CHKERRQ(ierr);    

//   ierr = VecSetValues(sedBuriedTracers,lSedBuriedSize,gSedBuriedIndices,localSedBuriedTR,INSERT_VALUES);CHKERRQ(ierr);
//   ierr = VecAssemblyBegin(sedBuriedTracers);CHKERRQ(ierr);
//   ierr = VecAssemblyEnd(sedBuriedTracers);CHKERRQ(ierr);    

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished pushing sediment arrays\n");CHKERRQ(ierr);  
  
/* Initial conditions     */
  ierr = PetscOptionsGetString(NULL,NULL,"-ised_mixed",sedMixedIniFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading sediment mixed tracers initial conditions from %s\n", sedMixedIniFile);CHKERRQ(ierr);    
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,sedMixedIniFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr = VecLoad(sedMixedTracers,fd);CHKERRQ(ierr); /* IntoVector */
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);          
  } else {  /* set to zero */
// 	 ierr = VecSet(sedMixedTracers,zero);CHKERRQ(ierr);
  }

  ierr = PetscOptionsGetString(NULL,NULL,"-ised_Buried",sedBuriedIniFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (flg) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Reading sediment buried tracers initial conditions from %s\n", sedBuriedIniFile);CHKERRQ(ierr);    
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,sedBuriedIniFile,FILE_MODE_READ,&fd);CHKERRQ(ierr);
    ierr = VecLoad(sedBuriedTracers,fd);CHKERRQ(ierr); /* IntoVector */
    ierr = PetscViewerDestroy(&fd);CHKERRQ(ierr);          
  } else {  /* set to zero */
// 	 ierr = VecSet(sedBuriedTracers,zero);CHKERRQ(ierr);
  }

  mobi_sed_copy_data_(&lNumProfiles,&localSedMixedTR[0],&localSedBuriedTR[0],&toModel);

  ierr = StepTimerCreate(&ef->sedWriteTimer);CHKERRQ(ierr);
 	ierr = StepTimerIni("sed_write_", prefix, Iter0, ef->sedWriteTimer);CHKERRQ(ierr);
    
/* Output file */
  ierr = PetscOptionsGetString(NULL,NULL,"-osed_mixed",sedMixedOutFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate sediment mixed tracers output file name with the -osed_mixed option");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Sediment mixed tracers will be written to %s\n", sedMixedOutFile);CHKERRQ(ierr);    

  ierr = PetscOptionsGetString(NULL,NULL,"-osed_buried",sedBuriedOutFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_WORLD,1,"Must indicate sediment buried tracers output file name with the -osed_buried option");
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Sediment buried tracers will be written to %s\n", sedBuriedOutFile);CHKERRQ(ierr);    
  
  ierr = PetscOptionsHasName(NULL,NULL,"-sed_append",&sedAppendOutput);CHKERRQ(ierr);
  if (sedAppendOutput) {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Sediment model output will be appended\n");CHKERRQ(ierr);
    SED_OUTPUT_FILE_MODE=FILE_MODE_APPEND;
  } else {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Sediment model output will overwrite existing file(s)\n");CHKERRQ(ierr);
    SED_OUTPUT_FILE_MODE=FILE_MODE_WRITE;
  }    

/* Output times */
  ierr = PetscOptionsGetString(NULL,NULL,"-sed_time_file",sedOutTimeFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) {
	strcpy(sedOutTimeFile,"");
	sprintf(sedOutTimeFile,"%s","sed_output_time.txt");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Sediment model output times will be written to %s\n",sedOutTimeFile);CHKERRQ(ierr);

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,sedMixedOutFile,SED_OUTPUT_FILE_MODE,&sedmixedfd);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,sedBuriedOutFile,SED_OUTPUT_FILE_MODE,&sedburiedfd);CHKERRQ(ierr);

  if (!sedAppendOutput) {
   	ierr = PetscFOpen(PETSC_COMM_WORLD,sedOutTimeFile,"w",&sedfptime);CHKERRQ(ierr);  
    if (Iter0==sedWriteTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
      ierr = PetscFPrintf(PETSC_COMM_WORLD,sedfptime,"%d   %10.5f\n",Iter0,time0);CHKERRQ(ierr);     
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing sediment output at time %10.5f, step %d\n", tc,Iter);CHKERRQ(ierr);
      ierr = VecView(sedMixedTracers,sedmixedfd);CHKERRQ(ierr);
      ierr = VecView(sedBuriedTracers,sedburiedfd);CHKERRQ(ierr);
    }  
  } else {
    ierr = PetscFOpen(PETSC_COMM_WORLD,sedOutTimeFile,"a",&sedfptime);CHKERRQ(ierr);  
    if (Iter0==sedWriteTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */      		
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Sediment model output will be appended. Initial condition will NOT be written\n");CHKERRQ(ierr);      
    }  
  }

/* File name for final pickup */
  ierr = PetscOptionsGetString(NULL,NULL,"-sed_mixed_pickup_out",sedMixedPickupOutFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) {
	strcpy(sedMixedPickupOutFile,"");
    sprintf(sedMixedPickupOutFile,"%s","sed_mixed_pickup.petsc");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Final sediment mixed tracers pickup will be written to %s\n",sedMixedPickupOutFile);CHKERRQ(ierr);

  ierr = PetscOptionsGetString(NULL,NULL,"-sed_buried_pickup_out",sedBuriedPickupOutFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
  if (!flg) {
	strcpy(sedBuriedPickupOutFile,"");
    sprintf(sedBuriedPickupOutFile,"%s","sed_buried_pickup.petsc");
  }
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Final sediment buried tracers pickup will be written to %s\n",sedBuriedPickupOutFile);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Finished initializing sediment arrays\n");CHKERRQ(ierr);  
#endif /* O_sed */

  ierr = PetscOptionsHasName(NULL,NULL,"-calc_diagnostics",&ef->calcDiagnostics);CHKERRQ(ierr);
  if (ef->calcDiagnostics) {    
/*Data for diagnostics */
  ierr = StepTimerCreate(&ef->diagTimer);CHKERRQ(ierr);
 	ierr = StepTimerIni("diag_", prefix, Iter0+1, ef->diagTimer);CHKERRQ(ierr);
 	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed starting at and including (absolute) time step: %d\n", ef->diagTimer->startTimeStep);CHKERRQ(ierr);	
	 ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostics will be computed over %d time steps\n", ef->diagTimer->numTimeSteps);CHKERRQ(ierr);	

	ierr = PetscOptionsHasName(NULL,NULL,"-diag_append",&ef->diagAppendOutput);CHKERRQ(ierr);
	if (ef->diagAppendOutput) {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostic output will be appended\n");CHKERRQ(ierr);
	  ef->DIAG_FILE_MODE=FILE_MODE_APPEND;
	} else {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostic output will overwrite existing file(s)\n");CHKERRQ(ierr);
	  ef->DIAG_FILE_MODE=FILE_MODE_WRITE;
	}

/* Output times */
	ierr = PetscOptionsGetString(NULL,NULL,"-diag_time_file",ef->diagOutTimeFile,PETSC_MAX_PATH_LEN-1,&flg);CHKERRQ(ierr);
	if (!flg) {
	  strcpy(ef->diagOutTimeFile,"");
	  sprintf(ef->diagOutTimeFile,"%s","diagnostic_output_time.txt");
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Diagnostic output times will be written to %s\n",ef->diagOutTimeFile);CHKERRQ(ierr);

	if (!ef->diagAppendOutput) {
	  ierr = PetscFOpen(PETSC_COMM_WORLD,ef->diagOutTimeFile,"w",&ef->diagfptime);CHKERRQ(ierr);  
	} else {
	  ierr = PetscFOpen(PETSC_COMM_WORLD,ef->diagOutTimeFile,"a",&ef->diagfptime);CHKERRQ(ierr);  
	}
	
    ef->diagFirstTime=PETSC_TRUE;
    ef->doAverage=0;

    mobi_diags_ini_(&lNumProfiles, &lSize, &ef->numDiags2d, &ef->numDiags3d, &debugFlag);

    if (ef->numDiags2d > MAXDIAGS2d) {
      SETERRQ(PETSC_COMM_WORLD,1,"Number of 2-d diagnostics requested exceeds maximum. You must increase MAXDIAGS2d!");
    }    
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of 2-d diagnostics requested: %d\n",ef->numDiags2d);CHKERRQ(ierr);	                            	  
    
    if (ef->numDiags2d>0) {
	  for (id2d=0; id2d<ef->numDiags2d; id2d++) {
		ierr = PetscMalloc(lNumProfiles*sizeof(PetscScalar),&ef->localDiag2davg[id2d]);CHKERRQ(ierr);
		ef->Diag2dFile[id2d] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
		strcpy(ef->Diag2dFile[id2d],"");
// 		sprintf(Diag2dFile[id2d],"%s%03d%s","diag2d_",id2d+1,".bin");	          
		for (ip=0; ip<lNumProfiles; ip++) {
		  ef->localDiag2davg[id2d][ip]=0.0;      
		}
	  }
    }
   
/* Set the number of 3-d diagnostics */
/*	mobi_diags3d_ini_(&numDiags3d, &minusone, &debugFlag); */
	if (ef->numDiags3d > MAXDIAGS3d) {
	  SETERRQ(PETSC_COMM_WORLD,1,"Number of 3-d diagnostics requested exceeds maximum. You must increase MAXDIAGS3d!");
	}
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Number of 3-d diagnostics requested: %d\n",ef->numDiags3d);CHKERRQ(ierr);	                            	  
	
	if ((ef->numDiags2d == 0) & (ef->numDiags3d == 0)) {
	  SETERRQ(PETSC_COMM_WORLD,1,"You have specified the -calc_diagnostics flag but no diagnostics are being computed!");
    }

	if (ef->numDiags3d>0) {
	  ierr = VecDuplicateVecs(TR1,ef->numDiags3d,&ef->Diag3davg);CHKERRQ(ierr);

	  for (id3d=0; id3d<ef->numDiags3d; id3d++) {
		ierr = VecSet(ef->Diag3davg[id3d],zero);CHKERRQ(ierr);
		ef->Diag3dFile[id3d] = (char *) malloc(PETSC_MAX_PATH_LEN*sizeof(char));
		strcpy(ef->Diag3dFile[id3d],"");
// 		sprintf(Diag3dFile[id3d],"%s%03d%s","diag3d_",id3d+1,".petsc");
	  }
	  ierr = VecGetArrays(ef->Diag3davg,ef->numDiags3d,&ef->localDiag3davg);CHKERRQ(ierr);        
    }
    
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "calcExternalForcing"
PetscErrorCode calcExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *userctx)
{

  PetscInt numTracers;
  const char *prefix;
  void *ctx;

  PetscErrorCode ierr;
  PetscInt itr, ip;
  PetscScalar myTime;
  PetscScalar relyr, day;
  PetscInt itf;
  PetscScalar alpha;
  PetscInt k;
#if defined O_carbon  
  PetscScalar localFocean;
# if defined O_carbon_13_coupled
  PetscScalar localF13ocean;
# endif
#endif /* O_carbon */
#if defined O_sed
PetscInt timeToRunSedModel = 0;
PetscScalar totlocalweathflx;
#endif

  PetscCall(PetscContainerGetPointer(state->extforcctxcontainer, &ctx));
  ExternalForcingContext ef = (ExternalForcingContext)ctx;
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

  myTime = DeltaT*Iter; /* Iter should start at 0 */

#if defined NEEDEMP
  if (useEmP) {
	if (timeDependentBiogeochemForcing) {    
	  ierr = interpTimeDependentProfileSurfaceScalarData(tc,localEmP,timeDependentBiogeochemTimer->numTimes,
														 timeDependentBiogeochemTimer->tdt,localEmPtd,"EmPtd.bin");
	  ierr = dotProdProfileSurfaceScalarData(localEmP,localdA,&EmPglobavg);
	  EmPglobavg=EmPglobavg/totalOceanSurfaceArea;
	} else if (periodicBiogeochemForcing) {   
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
													biogeochemTimer->tdp,localEmPp,"EmP_");
	  ierr = dotProdProfileSurfaceScalarData(localEmP,localdA,&EmPglobavg);
	  EmPglobavg=EmPglobavg/totalOceanSurfaceArea;
	}
  }	
#endif

#if defined O_mobi
# if defined READ_SWRAD
  if (timeDependentBiogeochemForcing) {    
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localswrad,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localswradtd,"swradtd.bin");
  } else if (periodicBiogeochemForcing) {   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localswradp,"swrad_");
  }
# else
   insolation_(&lNumProfiles,&myTime,&locallatitude[0],&localswrad[0],&localtau[0]);
# endif /* READ_SWRAD */

  if (timeDependentBiogeochemForcing) {    
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localdisch,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localdischtd,"dischtd.bin");
  } else if (periodicBiogeochemForcing) {   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localdisch,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localdischp,"disch_");
  }

# if defined O_mobi_iron
  if (timeDependentBiogeochemForcing) {    
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localFe_adep,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localFe_adeptd,"Fe_adeptd.bin");
  } else if (periodicBiogeochemForcing) {   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localFe_adep,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localFe_adepp,"Fe_adep_");
  }                                                
# endif /* O_mobi_iron */

# if defined O_mobi_silicon
  if (periodicBiogeochemForcing) {
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localSi_dep,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localSi_depp,"Si_dep_");
  }
# endif /* O_mobi_silicon */
#endif /* O_mobi */

#if defined O_PaTh
  if (periodicBiogeochemForcing) {
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localdust_adep,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localdust_adepp,"dust_dep_");
  }
#endif /* O_PaTh */

#if defined BGC
  if (timeDependentBiogeochemForcing) {  
    ierr = TimeDependentVecInterp(tc,&Ts,timeDependentBiogeochemTimer->numTimes,timeDependentBiogeochemTimer->tdt,Tstd,"Tstd.petsc");
    ierr = TimeDependentVecInterp(tc,&Ss,timeDependentBiogeochemTimer->numTimes,timeDependentBiogeochemTimer->tdt,Sstd,"Sstd.petsc");
    
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localaice,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localaicetd,"aicetd.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localhice,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localhicetd,"hicetd.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localhsno,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localhsnotd,"hsnotd.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localwind,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localwindtd,"windtd.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localatmosp,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localatmosptd,"atmosptd.bin");
  } else if (periodicBiogeochemForcing) {   
    ierr = PeriodicVecInterp(tc,&Ts,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,Tsp,"Ts_");
    ierr = PeriodicVecInterp(tc,&Ss,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,Ssp,"Ss_");	

    ierr = interpPeriodicProfileSurfaceScalarData(tc,localaice,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localaicep,"aice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhice,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localhicep,"hice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhsno,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localhsnop,"hsno_");                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localwindp,"wind_");   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localatmospp,"atmosp_");	
  }
#endif /* BGC */                                                				                                  

#if defined O_carbon
# if defined O_co2ccn_data
/* Interpolate atmospheric pCO2   */
  if (tc>=TpCO2atm_hist[0]) {
	ierr = calcInterpFactor(numpCO2atm_hist,tc,TpCO2atm_hist,&itf,&alpha); CHKERRQ(ierr);
	pCO2atm = alpha*pCO2atm_hist[itf] + (1.0-alpha)*pCO2atm_hist[itf+1];	  
  } else {
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: time < %10.5f. Assuming pCO2atm=%g\n",TpCO2atm_hist[0],pCO2atm);CHKERRQ(ierr);
  }
# endif /* O_co2ccn_data */

# if defined O_TMM_interactive_atmosphere
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
# endif /* O_TMM_interactive_atmosphere */

# if defined O_c14ccn_data
  if (timeDependentAtmosphericC14) {
	if (tc>=TC14atm_hist[0]) {
	  ierr = calcInterpFactor(numC14atm_hist,tc,TC14atm_hist,&itf,&alpha); CHKERRQ(ierr);
	  dc14ccnnhatm = alpha*C14atmnh_hist[itf] + (1.0-alpha)*C14atmnh_hist[itf+1];
	  dc14ccneqatm = alpha*C14atmeq_hist[itf] + (1.0-alpha)*C14atmeq_hist[itf+1];
	  dc14ccnshatm = alpha*C14atmsh_hist[itf] + (1.0-alpha)*C14atmsh_hist[itf+1];
	} else {
	  ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: time < %10.5f. Assuming DC14atm (NH,EQ,SH)=%g,%g,%g\n",TC14atm_hist[0],dc14ccnnhatm,dc14ccneqatm,dc14ccnshatm);CHKERRQ(ierr);
	}
  }
# endif /* O_c14ccn_data */

# if defined O_c13ccn_data
  if (tc>=TC13atm_hist[0]) {
	ierr = calcInterpFactor(numC13atm_hist,tc,TC13atm_hist,&itf,&alpha); CHKERRQ(ierr);
	dc13atm = alpha*C13atm_hist[itf] + (1.0-alpha)*C13atm_hist[itf+1];
    pC13O2atm = (r13Func(dc13atm)/(1+r13Func(dc13atm)))*pCO2atm; /* mole fraction of atmospheric 13CO2 */
  } else {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Warning: time < %10.5f. Assuming dc13atm=%g\n",TC13atm_hist[0],dc13atm);CHKERRQ(ierr);
  }
# endif /* O_c13ccn_data */
#endif /* O_carbon */

#if defined NEEDEMP
/* compute global means */
  if (useEmP) {
#if !defined O_constant_flux_reference 
	for (itr=0; itr<numTracers; itr++) {    
      ierr = VecDot(surfVolFrac,state->c[itr],&ef->TRglobavg[itr]);CHKERRQ(ierr); /* volume weighted mean surface TR */									              
    }
#endif    
/*    EmPglobavg = 0.0; */ /* this is set to zero above */
  }
#endif

  relyr = myTime/secondsPerYear; /* number of years (and fractional years) of model */
  day = myTime/86400.0 - floor(relyr)*daysPerYear; /* relative day number referenced to the beginning of the current year */

#if defined O_carbon
# if defined O_TMM_interactive_atmosphere
  localFocean = 0.0;  
  Focean = 0.0;
  Fland = 0.0;
#  if defined O_carbon_13_coupled
  localF13ocean = 0.0;
  F13ocean = 0.0; 
#  endif 
# endif /* O_TMM_interactive_atmosphere */
#endif /* O_carbon */

#if defined O_sed
  if ((iLoop % numOceanStepsPerSedStep)==0) { 
   timeToRunSedModel=1;
  }  
#endif

  if (ef->calcDiagnostics) {  
	   if (Iter0+iLoop==(ef->diagTimer->startTimeStep)) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Switching on diagnostics accumulation at step %d\n", Iter0+iLoop);CHKERRQ(ierr);
	     mobi_diags_start_(&debugFlag); /* Switch on diagnostics accumulation at start of cycle */
	   }	
  }                         

  for (itr=0; itr<numTracers; itr++) {    	
	   mobi_copy_data_(&lSize,&lNumProfiles,&itr,ef->localTR[itr],&toModel);
  }  

  mobi_calc_(&lSize,&lNumProfiles,
			 &day,&relyr,
			 localTs,localSs,ef->TRglobavg,
			 localEmP, &EmPglobavg,
#if defined BGC
			 localwind, localaice, localhice, localhsno,
#endif			  			 
#if defined O_carbon
# if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
#  if defined O_carbon_co2_2d
			 pCO2atm,
#  else
			 &pCO2atm,
#  endif                
# endif
# if defined O_c14ccn_data
               &dc14ccnnhatm,&dc14ccneqatm,&dc14ccnshatm,
# endif
# if defined O_c13ccn_data || defined O_carbon_13_coupled
               &pC13O2atm,
# endif               
#endif /* O_carbon */
#if defined O_mobi
			 localswrad,
             &globaldisch, localdisch,
# if defined O_mobi_iron
			 localFe_adep, 
# endif
# if defined O_mobi_silicon
			 localSi_dep, &globalSi_dep, 
# endif
#endif /* O_mobi */
#if defined O_PaTh
			 localdust_adep, 
#endif
#if defined O_carbon
			 ef->localgasexflux, ef->localtotflux, 
# if defined O_carbon_13_coupled
             localc13gasexflux,
# endif
#endif
#if defined O_sed
             &timeToRunSedModel,
             &globalweathflx, &totlocalweathflx,
#endif
			 &debugFlag);

  for (itr=0; itr<numTracers; itr++) {    
	mobi_copy_data_(&lSize,&lNumProfiles,&itr,ef->localJTR[itr],&fromModel);
  }  

  if (ef->calcDiagnostics) {  
	   if (Iter0+iLoop>=(ef->diagTimer->startTimeStep)) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */	
	     mobi_diags_accumulate_(&ef->diagTimer->count, &ef->doAverage, &debugFlag);
	   }	
  }                         

#if defined O_carbon
# if defined O_TMM_interactive_atmosphere
/*  Note-1: following UVic (gasbc.F) we use the gas exchange flux to compute atmospheric CO2 evolution. */
/*  The virtual flux term included in localtotflux should go be zero when integrated over the global ocean surface, */
/*  although not when it includes weathering flux */
/*  Note-2: localdA is in cm^2; we convert it below to m^2 */
  for (ip=0; ip<lNumProfiles; ip++) {
	localFocean = localFocean + ef->localgasexflux[ip]*(localdA[ip]*1.e-4)*(12.0/1.e15)*secondsPerYear; /* PgC/y */
#  if defined O_carbon_13_coupled
/* SPK: we use the same mass as C12 to convert */
	localF13ocean = localF13ocean + localc13gasexflux[ip]*(localdA[ip]*1.e-4)*(12.0/1.e15)*secondsPerYear; /* PgC/y */
#  endif
  }
# endif /* O_TMM_interactive_atmosphere */
#endif /* O_carbon */

  if (ef->calcDiagnostics) {  
	   if (Iter0+iLoop>=(ef->diagTimer->startTimeStep)) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */	
/*    still within same averaging block so increment accumulate count */	
	     ef->diagTimer->count++;
	   }	
  }
    
#if defined O_carbon
# if defined O_TMM_interactive_atmosphere
  MPI_Allreduce(&localFocean, &Focean, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);    
#  if defined O_carbon_13_coupled
  MPI_Allreduce(&localF13ocean, &F13ocean, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
#  endif  
    
  if (useLandModel) {
/*	landsource_(&landState[0],&pCO2atm,&landUseEmission,&deltaTsg,&Fland,&landSource[0]); */ /* returns S and Fl in PgC/y */
/*    time step land */
	for (k=0;k<=2;k++) {
	  landState[k] = landState[k] + atmModelDeltaT*landSource[k];
	}
  }

/* time step atmosphere */
  pCO2atm = pCO2atm + atmModelDeltaT*(fossilFuelEmission + landUseEmission - Focean - Fland)/ppmToPgC;
#  if defined O_carbon_13_coupled
/* SPK: to be consistent with MOBI we use the same conversion factor for C13 */
  pC13O2atm = pC13O2atm + atmModelDeltaT*(- F13ocean)/ppmToPgC;
#  endif  
# endif /* O_TMM_interactive_atmosphere */
#endif /* O_carbon */

#if defined O_sed
# if defined O_sed_weath_diag
  if (timeToRunSedModel == 1) {
	MPI_Allreduce(&totlocalweathflx, &globalweathflx, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Global weathering flux: %d = %g umol C/s\n", Iter, globalweathflx);CHKERRQ(ierr);  
  }	
# endif
#endif /* O_sed */
  
/* Convert to discrete tendency */
  for (itr=0; itr<numTracers; itr++) {
	ierr = VecScale(state->qef[itr],DeltaT);CHKERRQ(ierr);
  }
    
  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "writeExternalForcing"
PetscErrorCode writeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *userctx)
{
  PetscInt numTracers;
  const char *prefix;
  void *ctx;

  PetscErrorCode ierr;

  PetscCall(PetscContainerGetPointer(state->extforcctxcontainer, &ctx));
  ExternalForcingContext ef = (ExternalForcingContext)ctx;
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

/* Note: tc and iLoop are the time and step at the end of the current time step. */

#if defined O_carbon
# if defined O_TMM_interactive_atmosphere
/* accumulate/write instantaneous atmos model state */
  if (Iter0+iLoop>atmWriteTimer->startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
	  atmWriteTimer->count++;
  }

  if (Iter0+iLoop>=atmWriteTimer->startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
   	if ((atmWriteTimer->count==atmWriteTimer->numTimeSteps) || (Iter0+iLoop==atmWriteTimer->startTimeStep)) { /* time to write out */
//    NOTE: Focean and F13ocean were computed at the start of the current time step, whereas pCO2atm and pC13O2atm
//    are the updated values at the end of the time step. Hence the offset of one time step below.
       ierr = PetscPrintf(PETSC_COMM_WORLD,"Focean: %d = %10.5f\n", Iter0+iLoop-1, Focean);CHKERRQ(ierr);
       ierr = PetscPrintf(PETSC_COMM_WORLD,"pCO2atm: %d = %10.5f\n", Iter0+iLoop, pCO2atm);CHKERRQ(ierr);
#  if defined O_carbon_13_coupled
       ierr = PetscPrintf(PETSC_COMM_WORLD,"F13ocean: %d = %10.5f\n", Iter0+iLoop-1, F13ocean);CHKERRQ(ierr);
       ierr = PetscPrintf(PETSC_COMM_WORLD,"pC13O2atm: %d = %10.5f\n", Iter0+iLoop, pC13O2atm);CHKERRQ(ierr); 
       dc13atm=dc13Func(pC13O2atm,pCO2atm);	
       ierr = PetscPrintf(PETSC_COMM_WORLD,"dc13atm: %d = %10.5f\n", Iter0+iLoop, dc13atm);CHKERRQ(ierr);        
#  endif
	  
       ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing atmospheric model output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
       ierr = PetscFPrintf(PETSC_COMM_WORLD,atmfptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);           
       ierr = writeBinaryScalarData("pCO2atm_output.bin",&pCO2atm,1,PETSC_TRUE);
#  if defined O_carbon_13_coupled
       ierr = writeBinaryScalarData("pC13O2atm_output.bin",&pC13O2atm,1,PETSC_TRUE);
       ierr = writeBinaryScalarData("dc13atm_output.bin",&dc13atm,1,PETSC_TRUE);	
#  endif

       if (useEmissions) {
         ierr = PetscPrintf(PETSC_COMM_WORLD,"Cumulative emissions at time %10.5f, step %d = %10.6f PgC\n", tc, Iter0+iLoop, cumulativeEmissions);CHKERRQ(ierr);
       }  

       if (useLandModel) {
/*      write instantaneous land model state */
         ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing land model output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
         ierr = writeBinaryScalarData("land_state_output.bin",landState,3,PETSC_TRUE);
       }

	   }

    if (atmWriteTimer->count==atmWriteTimer->numTimeSteps) {
      ierr = StepTimerUpdate(Iter0+iLoop, ef->atmWriteTimer);CHKERRQ(ierr);
    }
  }

# endif /* O_TMM_interactive_atmosphere */
#endif /* O_carbon */

#if defined O_sed
/*  accumulate/write instantaneous sediment model state */
  if (Iter0+iLoop>sedWriteTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
	  sedWriteTimer.count++;
  }

  if (Iter0+iLoop>=sedWriteTimer.startTimeStep) { /* note: startTimeStep is ABSOLUTE time step */
   	if ((sedWriteTimer.count==sedWriteTimer.numTimeSteps) || (Iter0+iLoop==sedWriteTimer.startTimeStep)) { /* time to write out */
   	  mobi_sed_copy_data_(&lNumProfiles,&localSedMixedTR[0],&localSedBuriedTR[0],&fromModel);
      ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing sediment model output at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
      ierr = VecView(sedMixedTracers,sedmixedfd);CHKERRQ(ierr);
      ierr = VecView(sedBuriedTracers,sedburiedfd);CHKERRQ(ierr);
	   }

	   if (sedWriteTimer.count==sedWriteTimer.numTimeSteps) {
		    ierr = StepTimerUpdate(Iter0+iLoop, ef->sedWriteTimer);CHKERRQ(ierr);
   	}
  }
#endif /* O_sed */

  if (ef->calcDiagnostics) {  
	   if (Iter0+iLoop>=(ef->diagTimer->startTimeStep)) { /* start time averaging (note: startTimeStep is ABSOLUTE time step) */  

   	  if (ef->diagTimer->count==ef->diagTimer->numTimeSteps) { /* time to write averages to file */
        ef->doAverage=1;
        mobi_diags_accumulate_(&ef->diagTimer->count, &ef->doAverage, &debugFlag);
        
        ierr = PetscPrintf(PETSC_COMM_WORLD,"Writing diagnostics time average at time %10.5f, step %d\n", tc, Iter0+iLoop);CHKERRQ(ierr);
        ierr = PetscFPrintf(PETSC_COMM_WORLD,ef->diagfptime,"%d   %10.5f\n",Iter0+iLoop,tc);CHKERRQ(ierr);           

        if (ef->numDiags2d>0) {
          for (id2d=0; id2d<ef->numDiags2d; id2d++) {
            mobi_diags2d_copy_(&id2d,&ef->localDiag2davg[id2d][0],ef->Diag2dFile[id2d],&debugFlag);          		  
            ierr = writeProfileSurfaceScalarData(ef->Diag2dFile[id2d],ef->localDiag2davg[id2d],1,ef->diagAppendOutput);		   
          }  
        }
        
        if (ef->numDiags3d>0) {        
          for (id3d=0; id3d<ef->numDiags3d; id3d++) {
            mobi_diags3d_copy_(&id3d,&ef->localDiag3davg[id3d][0],ef->Diag3dFile[id3d],&debugFlag);
            /* Open file here if first time. The reason for doing it here is that we don't know the file name */
            /* until the first call to mobi_diags3d_copy above */
            if (ef->diagFirstTime) {
                 ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,ef->Diag3dFile[id3d],ef->DIAG_FILE_MODE,&ef->fddiag3dout[id3d]);CHKERRQ(ierr);
            }
            ierr = VecView(ef->Diag3davg[id3d],ef->fddiag3dout[id3d]);CHKERRQ(ierr);	      
		        }  
        }
        
        if (ef->diagFirstTime) {
          ef->diagFirstTime=PETSC_FALSE;
          ef->diagAppendOutput=PETSC_TRUE;
          ef->DIAG_FILE_MODE=FILE_MODE_APPEND;
        }

     		 ierr = StepTimerUpdate(Iter0+iLoop, ef->diagTimer);CHKERRQ(ierr);
		      ef->doAverage = 0;

        if (ef->diagTimer->haveResetStartTimeStep) {
          ierr = PetscPrintf(PETSC_COMM_WORLD,"Switching off diagnostics accumulation at step %d\n", Iter0+iLoop);CHKERRQ(ierr);		
          mobi_diags_stop_(&debugFlag); /* Switch off diagnostics accumulation until next cycle */
        }
	     }
	   }  
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "finalizeExternalForcing"
PetscErrorCode finalizeExternalForcing(PetscScalar tc, PetscInt Iter, TMMState state, void *userctx)
{
  PetscInt numTracers;
  const char *prefix;
  void *ctx;

  PetscErrorCode ierr;

  PetscCall(PetscContainerGetPointer(state->extforcctxcontainer, &ctx));
  ExternalForcingContext ef = (ExternalForcingContext)ctx;
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

/* write final pickup */
#if defined O_carbon
# if defined O_TMM_interactive_atmosphere
/* write instantaneous atmos model state */
  ierr = writeBinaryScalarData("pickup_pCO2atm.bin",&pCO2atm,1,PETSC_FALSE);
#  if defined O_carbon_13_coupled
/*  ierr = writeBinaryScalarData("pickup_pC13O2atm.bin",&pC13O2atm,1,PETSC_FALSE); */
  dc13atm=dc13Func(pC13O2atm,pCO2atm);
  ierr = writeBinaryScalarData("pickup_dc13atm.bin",&dc13atm,1,PETSC_FALSE);
#  endif
  if (useLandModel) {
/*   write instantaneous land model state */
	ierr = writeBinaryScalarData("pickup_land_state.bin",landState,3,PETSC_FALSE);
  }

  ierr = PetscFClose(PETSC_COMM_WORLD,atmfptime);CHKERRQ(ierr);  
# endif /* O_TMM_interactive_atmosphere */
#endif /* O_carbon */

  ierr = VecDestroy(&Ts);CHKERRQ(ierr);
  ierr = VecDestroy(&Ss);CHKERRQ(ierr);

#if defined NEEDEMP
  if (useEmP) {
	if (timeDependentBiogeochemForcing) {
	  ierr = TimeDependentArrayDestroy(&localEmPtd);CHKERRQ(ierr);    
	} else if (periodicBiogeochemForcing) {
	  ierr = PeriodicArrayDestroy(&localEmPp);CHKERRQ(ierr);    
	}
  }
#endif

#if defined O_mobi  
# if defined READ_SWRAD    
  if (timeDependentBiogeochemForcing) {    
	ierr = TimeDependentArrayDestroy(&localswradtd);CHKERRQ(ierr);    
  } else if (periodicBiogeochemForcing) {    
    ierr = PeriodicArrayDestroy(&localswradp);CHKERRQ(ierr);
  }
# endif

  if (timeDependentBiogeochemForcing) {    
	ierr = TimeDependentArrayDestroy(&localdischtd);CHKERRQ(ierr);    
  } else if (periodicBiogeochemForcing) {    
    ierr = PeriodicArrayDestroy(&localdischp);CHKERRQ(ierr);
  }

# if defined O_mobi_iron
  if (timeDependentBiogeochemForcing) {    
	ierr = TimeDependentArrayDestroy(&localFe_adeptd);CHKERRQ(ierr);    
  } else if (periodicBiogeochemForcing) {    
    ierr = PeriodicArrayDestroy(&localFe_adepp);CHKERRQ(ierr);
  }
# endif

# if defined O_mobi_silicon
  if (periodicBiogeochemForcing) {    
    ierr = PeriodicArrayDestroy(&localSi_depp);CHKERRQ(ierr);
  }
# endif
#endif /* O_mobi */

#if defined O_PaTh
  if (periodicBiogeochemForcing) {    
    ierr = PeriodicArrayDestroy(&localdust_adepp);CHKERRQ(ierr);
  }
#endif

#if defined BGC    
  if (timeDependentBiogeochemForcing) {
	ierr = TimeDependentVecDestroy(&Tstd);
	ierr = TimeDependentVecDestroy(&Sstd);
	ierr = TimeDependentArrayDestroy(&localaicetd);CHKERRQ(ierr);    
	ierr = TimeDependentArrayDestroy(&localhicetd);CHKERRQ(ierr);    
	ierr = TimeDependentArrayDestroy(&localhsnotd);CHKERRQ(ierr);    
	ierr = TimeDependentArrayDestroy(&localwindtd);CHKERRQ(ierr);    
	ierr = TimeDependentArrayDestroy(&localatmosptd);CHKERRQ(ierr);    
  } else if (periodicBiogeochemForcing) {
    ierr = PeriodicVecDestroy(&Tsp);CHKERRQ(ierr);
    ierr = PeriodicVecDestroy(&Ssp);CHKERRQ(ierr);
    ierr = PeriodicArrayDestroy(&localaicep);CHKERRQ(ierr);
    ierr = PeriodicArrayDestroy(&localhicep);CHKERRQ(ierr);
    ierr = PeriodicArrayDestroy(&localhsnop);CHKERRQ(ierr);    
    ierr = PeriodicArrayDestroy(&localwindp);CHKERRQ(ierr);    
    ierr = PeriodicArrayDestroy(&localatmospp);CHKERRQ(ierr); 
  }
#endif /* BGC */

#if defined O_sed
  ierr = VecDestroy(&sedMixedTracers);CHKERRQ(ierr); 
  ierr = VecDestroy(&sedBuriedTracers);CHKERRQ(ierr); 
  ierr = PetscViewerDestroy(&sedmixedfd);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&sedburiedfd);CHKERRQ(ierr);    
#endif /* O_sed */
  
  if (ef->calcDiagnostics) {      
    mobi_diags_finalize_(&debugFlag);

	if (ef->numDiags2d>0) {
	  for (id2d=0; id2d<ef->numDiags2d; id2d++) {
		PetscFree(ef->localDiag2davg[id2d]);
	  }	
    }
    
	if (ef->numDiags3d>0) {    
	  ierr = VecDestroyVecs(ef->numDiags3d,&ef->Diag3davg);CHKERRQ(ierr);  
	  for (id3d=0; id3d<ef->numDiags3d; id3d++) {
		ierr = PetscViewerDestroy(&ef->fddiag3dout[id3d]);CHKERRQ(ierr);
	  }
	} 
    ierr = PetscFClose(PETSC_COMM_WORLD,ef->diagfptime);CHKERRQ(ierr);  	   
  }

  return 0;
}

#undef __FUNCT__
#define __FUNCT__ "reInitializeExternalForcing"
PetscErrorCode reInitializeExternalForcing(PetscScalar tc, PetscInt Iter, PetscInt iLoop, TMMState state, void *userctx)
{
  PetscInt numTracers;
  const char *prefix;
  void *ctx;

  PetscErrorCode ierr;
  PetscInt ip;
  PetscScalar myTime;

  PetscCall(PetscContainerGetPointer(state->extforcctxcontainer, &ctx));
  ExternalForcingContext ef = (ExternalForcingContext)ctx;
  
  numTracers=state->numTracers;
  prefix = ((PetscObject)state)->prefix;

  myTime = DeltaT*Iter; /* Iter should start at 0 */

#if defined NEEDEMP
  if (useEmP) {
	if (timeDependentBiogeochemForcing) {    
	  ierr = interpTimeDependentProfileSurfaceScalarData(tc,localEmP,timeDependentBiogeochemTimer->numTimes,
														 timeDependentBiogeochemTimer->tdt,localEmPtd,"EmPtd.bin");
	  ierr = dotProdProfileSurfaceScalarData(localEmP,localdA,&EmPglobavg);
	  EmPglobavg=EmPglobavg/totalOceanSurfaceArea;														 
	} else if (periodicBiogeochemForcing) {   
	  ierr = interpPeriodicProfileSurfaceScalarData(tc,localEmP,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
													biogeochemTimer->tdp,localEmPp,"EmP_");
	  ierr = dotProdProfileSurfaceScalarData(localEmP,localdA,&EmPglobavg);
	  EmPglobavg=EmPglobavg/totalOceanSurfaceArea;													
	}
  }	
#endif
  
#if defined O_mobi
# if defined READ_SWRAD
  if (timeDependentBiogeochemForcing) {    
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localswrad,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localswradtd,"swradtd.bin");
  } else if (periodicBiogeochemForcing) {   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localswrad,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localswradp,"swrad_");
  }
# else
   insolation_(&lNumProfiles,&myTime,&locallatitude[0],&localswrad[0],&localtau[0]);
# endif /* READ_SWRAD */

  if (timeDependentBiogeochemForcing) {    
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localdisch,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localdischtd,"dischtd.bin");
  } else if (periodicBiogeochemForcing) {   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localdisch,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localdischp,"disch_");
  }
  	
# if defined O_mobi_iron
  if (timeDependentBiogeochemForcing) {    
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localFe_adep,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localFe_adeptd,"Fe_adeptd.bin");
  } else if (periodicBiogeochemForcing) {   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localFe_adep,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localFe_adepp,"Fe_adep_");
  }                                                
# endif /* O_mobi_iron */

# if defined O_mobi_silicon
  if (periodicBiogeochemForcing) {
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localSi_dep,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localSi_depp,"Si_dep_");
  }
# endif /* O_mobi_silicon */
#endif /* O_mobi */

#if defined O_PaTh
  if (periodicBiogeochemForcing) {
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localdust_adep,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localdust_adepp,"dust_dep_");
  }
#endif /* O_PaTh */

#if defined BGC
  if (timeDependentBiogeochemForcing) {  
    ierr = TimeDependentVecInterp(tc,&Ts,timeDependentBiogeochemTimer->numTimes,timeDependentBiogeochemTimer->tdt,Tstd,"Tstd.petsc");
    ierr = TimeDependentVecInterp(tc,&Ss,timeDependentBiogeochemTimer->numTimes,timeDependentBiogeochemTimer->tdt,Sstd,"Sstd.petsc");
    
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localaice,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localaicetd,"aicetd.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localhice,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localhicetd,"hicetd.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localhsno,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localhsnotd,"hsnotd.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localwind,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localwindtd,"windtd.bin");
    ierr = interpTimeDependentProfileSurfaceScalarData(tc,localatmosp,timeDependentBiogeochemTimer->numTimes,
                                                         timeDependentBiogeochemTimer->tdt,localatmosptd,"atmosptd.bin");
  } else if (periodicBiogeochemForcing) {   
    ierr = PeriodicVecInterp(tc,&Ts,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,Tsp,"Ts_");
    ierr = PeriodicVecInterp(tc,&Ss,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,biogeochemTimer->tdp,Ssp,"Ss_");	
  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localaice,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localaicep,"aice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhice,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localhicep,"hice_");
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localhsno,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localhsnop,"hsno_");                                                  
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localwind,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localwindp,"wind_");   
    ierr = interpPeriodicProfileSurfaceScalarData(tc,localatmosp,biogeochemTimer->cyclePeriod,biogeochemTimer->numPerPeriod,
                                                  biogeochemTimer->tdp,localatmospp,"atmosp_");	
  }
#endif /* BGC */  
    
  return 0;
}
