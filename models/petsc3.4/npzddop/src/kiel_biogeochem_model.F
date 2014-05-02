CBOP
C !ROUTINE: MODEL

C !INTERFACE: ==========================================================
      SUBROUTINE KIEL_BIOGEOCHEM_MODEL(Nrloc,DeltaT, 
     &                  PTR1,PTR2,PTR3,PTR4,PTR5,PTR6, 
#ifdef CARBON
     &                 PTR7,PTR8,DICglobalavg,localEmP,localpCO2atm,
#endif     
     &                 thetaloc,saltloc,FIceloc,SWRADloc,TAUloc,
     &                 WINDloc,ATMOSPloc,dzloc,
     &                 GPTR1loc,GPTR2loc,GPTR3loc,
     &                 GPTR4loc,GPTR5loc,GPTR6loc,
#ifdef CARBON
     &                  GPTR7loc,GPTR8loc,localph,localco2flux,
#endif     
     &                 useSeparateBiogeochemTS)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE

#include "BGC_PARAMS.h"
#include "BGC_CONTROL.h"
#ifdef CARBON
#include "CAR_PARAMS.h"
#endif

C !INPUT PARAMETERS: ===================================================
      INTEGER Nrloc
      real*8 DeltaT
      real*8 PTR1(Nrloc),PTR2(Nrloc),PTR3(Nrloc)
      real*8 PTR4(Nrloc),PTR5(Nrloc),PTR6(Nrloc)
#ifdef CARBON
      real*8 PTR7(Nrloc),PTR8(Nrloc),DICglobalavg,localEmP,localpCO2atm
#endif      
      real*8 thetaloc(Nrloc),saltloc(Nrloc)
      real*8 FIceloc,SWRADloc,TAUloc,WINDloc,ATMOSPloc,dzloc(Nrloc)
      logical useSeparateBiogeochemTS
C     Output variables
      real*8 GPTR1loc(Nrloc),GPTR2loc(Nrloc),GPTR3loc(Nrloc)
      real*8 GPTR4loc(Nrloc),GPTR5loc(Nrloc),GPTR6loc(Nrloc)
#ifdef CARBON
      real*8 GPTR7loc(Nrloc),GPTR8loc(Nrloc),localph,localco2flux
#endif

C !LOCAL VARIABLES: ====================================================
      INTEGER k,bgc_kloc
      REAL*8 bgc_dz(bgc_ktotal)
      REAL*8 bgc_swr,bgc_tau,bgc_seaice,bgc_wind,bgc_atmosp
      REAL*8 bgc_theta(bgc_ktotal),bgc_salt(bgc_ktotal)
      
C     Copy physics and geometry
      bgc_kloc=Nrloc
      bgc_swr=SWRADloc
      bgc_tau=TAUloc
      bgc_seaice=FIceloc
      bgc_wind=WINDloc
      bgc_atmosp=ATMOSPloc
      DO k=1,Nrloc
        bgc_dz(k)=dzloc(k)
        bgc_theta(k)=thetaloc(k)
        bgc_salt(k)=saltloc(k)
      ENDDO

C Assign the tracers in "wet" grid boxes
      DO k=1,Nrloc

        bgc_tracer(k,1)=PTR1(k)
        bgc_tracer(k,2)=PTR2(k)
        bgc_tracer(k,3)=PTR3(k)
        bgc_tracer(k,4)=PTR4(k)
        bgc_tracer(k,5)=PTR5(k)
        bgc_tracer(k,6)=PTR6(k)
#ifdef CARBON
        bgc_tracer(k,7)=PTR7(k)
        bgc_tracer(k,8)=PTR8(k)
	sph=localph
	emp=localEmP
	pco2atm=localpCO2atm
	dicgave=DICglobalavg
#endif
      ENDDO

      CALL BGC_MODEL(bgc_kloc,bgc_dz,
     &               bgc_swr,bgc_tau,bgc_seaice,bgc_wind,bgc_atmosp,
     &               bgc_theta,bgc_salt)
      
C     Send data back
      if (useSeparateBiogeochemTS) then ! return updated tracer field
        DO k=1,Nrloc
          PTR1(k)=bgc_tracer(k,1)
          PTR2(k)=bgc_tracer(k,2)
          PTR3(k)=bgc_tracer(k,3)
          PTR4(k)=bgc_tracer(k,4)
          PTR5(k)=bgc_tracer(k,5)
          PTR6(k)=bgc_tracer(k,6)
#ifdef CARBON
          PTR7(k)=bgc_tracer(k,7)
          PTR8(k)=bgc_tracer(k,8)
#endif
        ENDDO
      else ! return tracer tendency
        DO k=1,Nrloc
          GPTR1loc(k)=(bgc_tracer(k,1)-PTR1(k))/DeltaT
          GPTR2loc(k)=(bgc_tracer(k,2)-PTR2(k))/DeltaT
          GPTR3loc(k)=(bgc_tracer(k,3)-PTR3(k))/DeltaT
          GPTR4loc(k)=(bgc_tracer(k,4)-PTR4(k))/DeltaT
          GPTR5loc(k)=(bgc_tracer(k,5)-PTR5(k))/DeltaT
          GPTR6loc(k)=(bgc_tracer(k,6)-PTR6(k))/DeltaT
#ifdef CARBON
          GPTR7loc(k)=(bgc_tracer(k,7)-PTR7(k))/DeltaT
          GPTR8loc(k)=(bgc_tracer(k,8)-PTR8(k))/DeltaT
#endif
        ENDDO
      endif
      
#ifdef CARBON
      localph=sph
      localco2flux=co2airseaflux
#endif      

c
      RETURN
      END
