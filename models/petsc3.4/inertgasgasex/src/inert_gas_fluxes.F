C !INTERFACE: ==========================================================
      SUBROUTINE INERT_GAS_FLUXES(myIter,myTime,
     &                 PTRloc,
     &                 thetaloc,saltloc,u10loc,ficeloc,
     &                 atmosploc,gasID,
     &                 Vgas,Finj,Fex,PTReq)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE

C #include "DIAGNOSTICS.h"

C !INPUT PARAMETERS: ===================================================
      INTEGER myIter
      real*8 myTime
      real*8 PTRloc
      real*8 thetaloc,saltloc,u10loc,ficeloc,atmosploc
      integer gasID
C     Output variables
      real*8 Vgas,Finj,Fex,PTReq

C !LOCAL VARIABLES: ====================================================
      real*8 rho0, permil, permeg
      parameter (rho0=1024.5d0, permil=1.d0/rho0, permeg=1.d-6)
      real*8 Ainj, Aex
      parameter (Ainj=2.357d-9, Aex=1.848d-5)      
      real*8 u3,Vgas660,Sc,visc,G_eq,Dgas,T,S,rho
      real*8 gas_mol_fract,solfun,moldiff,sw_dens0
CEOP

      T=thetaloc
      S=saltloc

      rho = sw_dens0(S,T)
      
C       write(*,*)T,S,u10loc,ficeloc,atmosploc
C       write(*,*)gasID
C       write(*,*)PTRloc
            
C     % gas diffusivities
      Dgas = moldiff(S,T,gasID)
      visc=1.d-4*(17.91d0 - 0.5381d0*T + 0.00694d0*(T**2) 
     &     + 0.02305d0*S)/rho
      
      Sc = visc/Dgas

C 	  Vgas660=(1.d0-ficeloc)*(0.222d0*(u10loc**2) + 0.333d0*u10loc)
	  Vgas660=0.222d0*(u10loc**2) + 0.333d0*u10loc
C /* % cm/h to m/s */
	  Vgas=(Vgas660/sqrt(Sc/660.0d0))/(100.d0*60.d0*60.d0)

C       write(*,*)Dgas,Sc,visc,Vgas660,Vgas

      u3=(u10loc-2.27d0)**3
      if (u3.lt.0.d0) u3=0.d0

      G_eq = rho*solfun(S,T,gasID)/1.d6

C       write(*,*)'Geq=',G_eq
C       write(*,*)'u3=',u3
      
      Finj = Ainj*atmosploc*gas_mol_fract(gasID)*u3
C       write(*,*)'Finj=',Finj

      Fex = Aex*sqrt(Dgas)*atmosploc*G_eq*u3
 
      PTReq = G_eq*atmosploc
      
      
C       Fas = -Vgas*(TR - G_eq*atmosploc) + Finj + Fex
C 
C       GTRloc = Fas/dzsurf
      
c
	  RETURN
	  END


