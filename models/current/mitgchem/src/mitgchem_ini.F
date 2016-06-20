#include "DIC_OPTIONS.h"

CBOP
C !ROUTINE: INI_MITGCM_BIOGEOCHEM

C !INTERFACE: ==========================================================
      SUBROUTINE MITGCHEM_INI(Nrloc,numTracers,myIter,myTime,
     &                       thetaloc,saltloc,pHloc,SILICAloc,
     &                       hFacCloc,recip_hFacCloc,drFloc,deltaTloc,
     &                       ip)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "DIC_VARS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
C  myIter               :: current timestep
C  myTime               :: current time
C  PTR_DIC              :: dissolced inorganic carbon
C  PTR_ALK              :: alkalinity
C  PTR_PO4              :: phosphate
c  PTR_DOP              :: dissolve organic phosphurous
c  PTR_O2               :: oxygen
C  PTR_FE               :: iron
      INTEGER Nrloc
      INTEGER myIter, numTracers
      real*8 myTime
      real*8 PTR_DIC(Nrloc),PTR_ALK(Nrloc),PTR_PO4(Nrloc),PTR_DOP(Nrloc)
      real*8 PTR_O2(Nrloc)
#ifdef ALLOW_FE
      real*8  PTR_FE(Nrloc)
#endif
      real*8 thetaloc(Nrloc),saltloc(Nrloc),pHloc,SILICAloc
      real*8 hFacCloc(Nrloc),recip_hFacCloc(Nrloc),drFloc(Nr)
      real*8 deltaTloc
      integer ip

C !LOCAL VARIABLES: ====================================================
	  INTEGER myThid
	  PARAMETER (myThid=1)
	  INTEGER k
	  LOGICAL firstTime 
      data firstTime/.TRUE./  ! initialize
      save firstTime
CEOP

C NEW

      IF (firstTime) THEN
C       Set tile indices              
        myBxLo(1)=1
        myBxHi(1)=1
        myByLo(1)=1
        myByHi(1)=1
      
C       The sequence below is similar to that in MITgcm

C       eeboot.F
		standardMessageUnit = 6
		scrUnit1            = 11
		scrUnit2            = 12
		eeDataUnit          = 13
		modelDataUnit       = 14
		errorMessageUnit    = 15

C       ini_parms.F
        nIter0 = 0  ! set this as the default value

        periodicExternalForcing = .FALSE.
        pickupSuff = ' '
        
		do k = 1,Nr
		  dTtracerLev(k) = deltaTloc
		enddo

C       packages_boot.F		
		useGCHEM = .TRUE.

C       packages_readparms.F->
C       ptracers_readparms.F
        PTRACERS_Iter0 = 1 ! this is to avoid pickup files being read
		PTRACERS_numInUse = numTracers
		do k = 1,Nr		
		  PTRACERS_dTLev(k)=dTtracerLev(k)
		enddo  

		CALL GCHEM_READPARMS(myThid)

C       Set vertical grid
C       copy vertical grid spacing to MIT array
		DO k=1,Nr
		 delR(k) = drFloc(k)
        ENDDO		

C       Set some defaults
        Ro_SeaLevel = 0.d0
		rkSign = -1.d0
        
        usingSphericalPolarGrid = .TRUE. ! needed for S/R insol
        usingCurvilinearGrid = .FALSE. ! needed for S/R insol
        rotateGrid = .FALSE. ! needed for S/R insol

        setInterFDr = .TRUE.
        setCenterDr = .FALSE.        
C		CALL READ_R8_FIELD(Nr,1,1,drF,"drF.bin")	   

C       ini_vertical_grid.F
C---  Set Level r-thickness (drF) and Center r-distances (drC)
		
C--     Interface r-distances are defined:		
		DO k=1,Nr
		  drF(k) = delR(k)
		ENDDO
C--     Cell Center r-distances undefined:
C       assume Center at middle between 2 Interfaces
		drC(1)  = 0.5 _d 0 *delR(1)
		DO k=2,Nr
		  drC(k) = 0.5 _d 0 *(delR(k-1)+delR(k))
		ENDDO
		drC(Nr+1) = 0.5 _d 0 *delR(Nr)

C---  Set r-position of  interFace (rF) and cell-Center (rC):
		rF(1)    = Ro_SeaLevel
		DO k=1,Nr
		 rF(k+1) = rF(k)  + rkSign*drF(k)
		ENDDO
		rC(1)   = rF(1)   + rkSign*drC(1)
		DO k=2,Nr
		  rC(k) = rC(k-1) + rkSign*drC(k)
		ENDDO

C-      Calculate reciprol vertical grid spacing :
		DO k=1,Nr+1
		 recip_drC(k)   = 1. _d 0/drC(k)
		ENDDO
		DO k=1,Nr
		 recip_drF(k)   = 1. _d 0/drF(k)
		ENDDO


C       ini_masks_etc.F
        kLowC(1,1,1,1) = Nrloc
        
C       packages_init_fixed.F ->
        CALL GCHEM_INIT_FIXED( myThid )
        
		firstTime = .FALSE.
      ENDIF

C     Copy surface variables
	  SILICA(1,1,1,1)=SILICAloc  ! but this will get overwritten in gchem_init_vari->dic_ini_forcing
C     Note: only need to set surface values, but set entire field to be safe	  
C     Also, we don't really need to set PTR_DOP,PTR_O2, and PTR_FE.
C     Set vertical profile
	  DO k=1,Nrloc
C       Copy masks	  
		maskC(1,1,k,1,1)=1.d0
		hFacC(1,1,k,1,1)=hFacCloc(k)
		recip_hFacC(1,1,k,1,1)=recip_hFacCloc(k)
C       Copy tracers		
        theta(1,1,k,1,1)=thetaloc(k)
        salt(1,1,k,1,1)=saltloc(k)
	  ENDDO

	  DO k=Nrloc+1,Nr
		maskC(1,1,k,1,1)=0.d0	  
		hFacC(1,1,k,1,1)=0.d0
		recip_hFacC(1,1,k,1,1)=0.d0
        theta(1,1,k,1,1)=0.d0
        salt(1,1,k,1,1)=0.d0		
	  ENDDO

C     initialise_varia.F->        
      CALL GCHEM_INIT_VARI( myThid )

C     Send pH back       
      pHloc=pH(1,1,1,1)
      
	  RETURN
	  END
