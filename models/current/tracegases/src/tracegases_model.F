#include "TRACEGASES_OPTIONS.h"

C !INTERFACE: ==========================================================
      SUBROUTINE TRACEGASES_MODEL(myIter,myTime,
     &                 PTR,
     &                 thetaloc,saltloc,Vgas660loc,atmosploc,
     &                 gasID,xTRatm,mixingRatioScaleFactor,
     &                 dzsurf,EmPloc,TRemp,
     &                 GTRloc,gasexfluxloc,totfluxloc,TReqloc
     &                          )

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE

C !INPUT PARAMETERS: ===================================================
C  myIter               :: current timestep
C  myTime               :: current time
C  PTR                  :: tracer concentration
      integer myIter
      real*8 myTime
      real*8 PTR
      real*8 thetaloc,saltloc,Vgas660loc,atmosploc
      integer gasID
      real*8 xTRatm, mixingRatioScaleFactor
      real*8 dzsurf,EmPloc,TRemp
C     Output variables
      real*8 GTRloc,gasexfluxloc,totfluxloc,TReqloc

C !LOCAL VARIABLES: ====================================================
      real*8 TR, TReq, pTRatm, Vgas
      real*8 ff, Sc, sol
      real*8 rho0
      parameter (rho0=1024.5d0)
      real*8 N2O_solubility
      real*8 CFC11_solubility, CFC12_solubility
      real*8 SF6_solubility
      real*8 N2O_schmidt_number
      real*8 CFC11_schmidt_number, CFC12_schmidt_number
      real*8 SF6_schmidt_number
CEOP

C     Copy surface variables
      TR=PTR

C     compute solubility and Schmidt number of gas
      if (gasID.eq.1) then
		ff=N2O_solubility(thetaloc,saltloc)
		Sc=N2O_schmidt_number(thetaloc)
      elseif (gasID.eq.2) then
		ff=CFC11_solubility(thetaloc,saltloc)
		Sc=CFC11_schmidt_number(thetaloc)      
      elseif (gasID.eq.3) then
		ff=CFC12_solubility(thetaloc,saltloc)
		Sc=CFC12_schmidt_number(thetaloc)
      elseif (gasID.eq.4) then
		ff=SF6_solubility(thetaloc,saltloc)
		Sc=SF6_schmidt_number(thetaloc)
	  else
        stop 'Error! Unknown gas! 1=N2O, 2=CFC11, 3=CFC12, 4=SF6'
      endif
      
C     ff is solubility of gas in mol/(L*atm) at 1 atm total atmospheric pressure
C     To convert to mol/(m^3*atm), multiply ff by 1000.0, i.e.
	  sol = ff*1000.0d0 ! solubility in mol/(m^3*atm)

C     Equilibrium concentration in mol/m^3 = sol*pTRatm, where
C       pTRatm=local partial pressure of TR [atm]
C        = [atmospheric mixing ratio of TR in dry air] x [total atmospheric pressure in atm]

      pTRatm = xTRatm*atmosploc*mixingRatioScaleFactor
	  TReq = sol*pTRatm ! equilibrium concentration in mol/m^3

      Vgas=Vgas660loc/sqrt(Sc/660.d0)
      
      gasexfluxloc=-Vgas*(TR - TReq)
      totfluxloc=gasexfluxloc + TRemp*EmPloc

C     Send data back      
      GTRloc = totfluxloc/dzsurf ! tendency
      TReqloc = TReq ! equilibrium concentration [mol/m^3]
      
	  RETURN
	  END
