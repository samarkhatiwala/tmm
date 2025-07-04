C !INTERFACE: ==========================================================
      SUBROUTINE INERT_GAS_FLUXES(PTRloc,thetaloc,saltloc,
     &                 u10loc,ficeloc,atmosploc,
     &                 gasID,pistonVelocityCoeff,Vgas,
     &                 Finj,Fex,
     &                 PTReq)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE

C #include "DIAGNOSTICS.h"

C !INPUT PARAMETERS: ===================================================
      real*8 PTRloc
      real*8 thetaloc,saltloc,u10loc,ficeloc,atmosploc
      integer gasID
      real*8 pistonVelocityCoeff ! [cm/hr]*[s^2/m^2]
C     Output variables
      real*8 Vgas,Finj,Fex,PTReq

C !LOCAL VARIABLES: ====================================================
      real*8 rho0, permil, permeg
      parameter (rho0=1024.5d0, permil=1.d0/rho0, permeg=1.d-6)
      real*8 cmph2mps
      parameter (cmph2mps=1.d0/(100.d0*60.d0*60.d0))
      real*8 Ainj, Aex
      parameter (Ainj=1.06d-9, Aex=2.19d-6)      
      real*8 u3,Vgas660,Sc,visc,G_eq,Dgas,T,S,rho
      real*8 gas_mol_fract,solfun,moldiff,sw_dens0
CEOP

      T=thetaloc
      S=saltloc

      rho = sw_dens0(S,T)
      
C     % gas diffusivities
      Dgas = moldiff(S,T,gasID)
      visc=1.d-4*(17.91d0 - 0.5381d0*T + 0.00694d0*(T**2) 
     &     + 0.02305d0*S)/rho
      
      Sc = visc/Dgas

	  Vgas660=(1.d0-ficeloc)*pistonVelocityCoeff*(u10loc**2) 
C /* % cm/h to m/s */
	  Vgas=(Vgas660/sqrt(Sc/660.0d0))*cmph2mps

C     Solubility equilibrium concentration
      G_eq = rho*solfun(S,T,gasID)/1.d6 
      PTReq = G_eq*atmosploc

C   DN: ice term added here too 
      u3= (1.d0-ficeloc)*(u10loc-2.27d0)**3
      if (u3.lt.0.d0) u3=0.d0

      Finj = Ainj*atmosploc*gas_mol_fract(gasID)*u3

      Fex = Aex*sqrt(Dgas)*atmosploc*G_eq*u3
      
C     Total air-sea flux: Fas = -Vgas*(PTR - PTReq) + Finj + Fex
      
	  RETURN
	  END


