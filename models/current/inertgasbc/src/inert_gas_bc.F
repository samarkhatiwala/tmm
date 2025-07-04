      SUBROUTINE INERT_GAS_BC(N, thetaloc, saltloc,
     &                        atmosploc, gasID,
     &                        PTReqloc)

C !USES: ===============================================================
      IMPLICIT NONE

C #include "DIAGNOSTICS.h"

C !INPUT PARAMETERS: ===================================================
      INTEGER N
      real*8 PTRloc
      real*8 thetaloc(N),saltloc(N),atmosploc(N)
      integer gasID
C     Output variables
      real*8 PTReqloc(N)

C !LOCAL VARIABLES: ====================================================
      real*8 rho0, permil, permeg
      parameter (rho0=1024.5d0, permil=1.d0/rho0, permeg=1.d-6)
      real*8 G_eq,T,S,rho
      real*8 gas_mol_fract,solfun,moldiff,sw_dens0
      INTEGER k
CEOP

      DO k=1,N
        T=thetaloc(k)
        S=saltloc(k)        
        rho = sw_dens0(S,T)        
        G_eq = rho*solfun(S,T,gasID)/1.d6
        PTReqloc(k) = G_eq*atmosploc(k)
      ENDDO
      
	  RETURN
	  END


