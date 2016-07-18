C$Header: /Users/ikriest/CVS/mops/mops_biogeochem_model.F,v 1.2 2016/06/03 09:28:59 ikriest Exp $
C$Name: mops-2_0 $

CBOP
C !ROUTINE: MODEL

C !INTERFACE: ==========================================================
      SUBROUTINE MOPS_BIOGEOCHEM_MODEL(Nrloc,DeltaT, 
#ifdef CARBON
     &                 DICglobalavg,ALKglobalavg,localEmP,localpCO2atm,
#endif     
     &                 thetaloc,saltloc,FIceloc,SWRADloc,TAUloc,
     &                 WINDloc,ATMOSPloc,dzloc,
#ifdef CARBON
     &                 localph,localco2flux,
#endif     
     &                 localburial, globalrunoff, localrunoffloc,
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
#ifdef CARBON
      real*8 DICglobalavg,ALKglobalavg,localEmP,localpCO2atm
#endif      
      real*8 thetaloc(Nrloc),saltloc(Nrloc)
      real*8 FIceloc,SWRADloc,TAUloc,WINDloc,ATMOSPloc,dzloc(Nrloc)
      logical useSeparateBiogeochemTS
C     Output variables
#ifdef CARBON
      real*8 localph,localco2flux
#endif

      real*8 localburial,globalrunoff,localrunoffloc(Nrloc)

C !LOCAL VARIABLES: ====================================================
      INTEGER k,bgc_kloc
      REAL*8 bgc_dz(bgc_ktotal)
      REAL*8 bgc_swr,bgc_tau,bgc_seaice,bgc_wind,bgc_atmosp
      REAL*8 bgc_theta(bgc_ktotal),bgc_salt(bgc_ktotal)

      REAL*8 bgc_runoffvol(bgc_ktotal),bgc_globalrunoff
      
C     Copy physics and geometry
      bgc_kloc=Nrloc
      bgc_swr=SWRADloc
      bgc_tau=TAUloc
      bgc_seaice=FIceloc
      bgc_wind=WINDloc
      bgc_atmosp=ATMOSPloc
      bgc_globalrunoff = globalrunoff

      DO k=1,Nrloc
        bgc_dz(k)=dzloc(k)
        bgc_theta(k)=thetaloc(k)
        bgc_salt(k)=saltloc(k)
        bgc_runoffvol(k) = localrunoffloc(k)
      ENDDO

C Assign the tracers in "wet" grid boxes
#ifdef CARBON
	  sph=localph
	  emp=localEmP
	  pco2atm=localpCO2atm
	  dicgave=DICglobalavg
C	  alkgave=ALKglobalavg
#endif

      CALL BGC_MODEL(bgc_kloc,bgc_dz,
     &               bgc_runoffvol,bgc_globalrunoff,
     &               bgc_swr,bgc_tau,bgc_seaice,bgc_wind,bgc_atmosp,
     &               bgc_theta,bgc_salt)
      
C     Send data back
#ifdef CARBON
      localph=sph
      localco2flux=co2airseaflux
#endif      

      localburial=flux_bury

      RETURN
      END
