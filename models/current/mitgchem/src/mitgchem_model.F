#include "DIC_OPTIONS.h"

CBOP
C !ROUTINE: MODEL

C !INTERFACE: ==========================================================
      SUBROUTINE MITGCHEM_MODEL(Nrloc,myIter,myTime,
     &                 thetaloc,saltloc,pHloc,WINDloc,
     &                 ATMOSPloc,pCO2atmloc,
     &                 SILICAloc,FIceloc,
#ifdef READ_PAR     
     &                 PARloc,
#else     
     &                 latloc,
#endif
#ifdef ALLOW_FE     
     &                 InputFeloc,
#endif
#ifdef LIGHT_CHL
     &                 chlloc,
#endif
     &                 alphaloc,rain_ratioloc,
     &                 hFacCloc,recip_hFacCloc,
     &                 pco2loc,co2fluxloc,
     &                 ip         )

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
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
      INTEGER myIter
      real*8 myTime
      real*8 thetaloc(Nrloc),saltloc(Nrloc),pHloc
      real*8 WINDloc,ATMOSPloc,SILICAloc,FIceloc
#ifdef READ_PAR      
      real*8 PARloc
#else
      real*8 latloc
#endif
#ifdef ALLOW_FE
      real*8 InputFeloc
#endif
#ifdef LIGHT_CHL
      real*8 chlloc
#endif
      real*8 pCO2atmloc
      real*8 alphaloc,rain_ratioloc
#ifdef ALLOW_OLD_VIRTUALFLUX      
      real*8 gsm_DICloc,gsm_ALKloc,gsm_Sloc
#endif      
      real*8 hFacCloc(Nrloc),recip_hFacCloc(Nrloc)
C     Output variables
      real*8 pco2loc,co2fluxloc,virtualfluxloc
      integer ip

#ifdef ALLOW_PTRACERS
#ifdef DIC_BIOTIC

C !LOCAL VARIABLES: ====================================================
	  INTEGER myThid,bi,bj
	  PARAMETER (myThid=1,bi=1,bj=1)
	  INTEGER k
CEOP

C       ini_masks_etc.F
      kLowC(1,1,1,1) = Nrloc

C     Copy surface variables
	  pH(1,1,1,1)=pHloc
	  wind(1,1,1,1)=WINDloc
	  AtmosP(1,1,1,1)=ATMOSPloc
      AtmospCO2(1,1,1,1)=pCO2atmloc*1.d-6 ! fix units	  
	  SILICA(1,1,1,1)=SILICAloc
	  FIce(1,1,1,1)=FIceloc
#ifdef READ_PAR
	  PAR(1,1,1,1)=PARloc
#else
      YC(1,1,1,1)=latloc
#endif
#ifdef LIGHT_CHL
      CHL(1,1,1,1)=chlloc
#endif
#ifdef ALLOW_FE
	  InputFe(1,1,1,1)=InputFeloc
#endif	  
#ifdef ALLOW_OLD_VIRTUALFLUX
      gsm_DIC=gsm_DICloc
      gsm_ALK=gsm_ALKloc
      gsm_s=gsm_Sloc
#endif

      if (alphaloc.ge.0.d0) alpha(1,1,1,1)=alphaloc
      if (rain_ratioloc.ge.0.d0) rain_ratio(1,1,1,1)=rain_ratioloc

C     Set vertical profile
	  DO k=1,Nrloc
C       Copy masks	  
		maskC(1,1,k,1,1)=1.d0
		hFacC(1,1,k,1,1)=hFacCloc(k)
		recip_hFacC(1,1,k,1,1)=recip_hFacCloc(k)
C       Copy tracers		
        theta(1,1,k,1,1)=thetaloc(k)
        salt(1,1,k,1,1)=saltloc(k)
C		pTracer(1,1,k,1,1,1)=PTR_DIC(k)
C		pTracer(1,1,k,1,1,2)=PTR_ALK(k)
C		pTracer(1,1,k,1,1,3)=PTR_PO4(k)
C		pTracer(1,1,k,1,1,4)=PTR_DOP(k)
C		pTracer(1,1,k,1,1,5)=PTR_O2(k)
C#ifdef ALLOW_FE         
C		pTracer(1,1,k,1,1,6)=PTR_FE(k)
C#endif                  
	  ENDDO

C      if (ip.eq.0) then
C        write(*,*)'MITGCHEM_MODEL:'
C        k=1
C        write(*,*)theta(1,1,k,1,1),salt(1,1,k,1,1)
C        write(*,*)pTracer(1,1,k,1,1,1),pTracer(1,1,k,1,1,2)
C        write(*,*)pTracer(1,1,k,1,1,3),pTracer(1,1,k,1,1,4)
C        write(*,*)pTracer(1,1,k,1,1,5),pTracer(1,1,k,1,1,6)        
C        write(*,*)InputFe(1,1,1,1)
C	  endif
	  
	  DO k=Nrloc+1,Nr
C        write(*,*)'Copying profile',k
		maskC(1,1,k,1,1)=0.d0	  
		hFacC(1,1,k,1,1)=0.d0
		recip_hFacC(1,1,k,1,1)=0.d0	 
        theta(1,1,k,1,1)=0.d0
        salt(1,1,k,1,1)=0.d0		
C		pTracer(1,1,k,1,1,1)=0.d0
C		pTracer(1,1,k,1,1,2)=0.d0
C		pTracer(1,1,k,1,1,3)=0.d0
C		pTracer(1,1,k,1,1,4)=0.d0
C		pTracer(1,1,k,1,1,5)=0.d0
C#ifdef ALLOW_FE         
C		pTracer(1,1,k,1,1,6)=0.d0
C#endif                  
	  ENDDO

      CALL GCHEM_FORCING_SEP( myTime, myIter, myThid )

C     Send data back
	  pHloc=pH(1,1,1,1)
      pco2loc=pCO2(1,1,1,1)*1.d6 ! convert to ppm
      co2fluxloc=FluxCO2(1,1,1,1)
      
C	  DO k=1,Nrloc
C		GDICloc(k)=(pTracer(1,1,k,1,1,1)-PTR_DIC(k))/PTRACERS_dTLev(k)
C		GALKloc(k)=(pTracer(1,1,k,1,1,2)-PTR_ALK(k))/PTRACERS_dTLev(k)
C		GPO4loc(k)=(pTracer(1,1,k,1,1,3)-PTR_PO4(k))/PTRACERS_dTLev(k)
C		GDOPloc(k)=(pTracer(1,1,k,1,1,4)-PTR_DOP(k))/PTRACERS_dTLev(k)
C		GO2loc(k)=(pTracer(1,1,k,1,1,5)-PTR_O2(k))/PTRACERS_dTLev(k)
C#ifdef ALLOW_FE
C		GFEloc(k)=(pTracer(1,1,k,1,1,6)-PTR_FE(k))/PTRACERS_dTLev(k)
C#endif
C      ENDDO      

C        write(*,*)'MITGCHEM_MODEL: after GCHEM_FORCING_SEP'
C        k=1
C        write(*,*)pTracer(1,1,k,1,1,1),pTracer(1,1,k,1,1,2)
C        write(*,*)pTracer(1,1,k,1,1,3),pTracer(1,1,k,1,1,4)
C        write(*,*)pTracer(1,1,k,1,1,5),pTracer(1,1,k,1,1,6)        
C        write(*,*)GDICloc(k),pTracer(1,1,k,1,1,1),PTR_DIC(k)
      
#endif
#endif

c
	  RETURN
	  END
