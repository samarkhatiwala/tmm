
CBOP
C !ROUTINE: MODEL

C !INTERFACE: ==========================================================
      SUBROUTINE INERT_GAS_DIAGNOSTICS(Nrloc,myIter,myTime,
     &                 PTRloc,thetaloc,saltloc,gasID,
     &                 PTReqloc,PTRsatanomloc)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE

C !INPUT PARAMETERS: ===================================================
      INTEGER Nrloc
      INTEGER myIter
      real*8 myTime   
      real*8 PTRloc(Nrloc),thetaloc(Nrloc),saltloc(Nrloc)
      integer gasID
C     Output variables
      real*8 PTReqloc(Nrloc),PTRsatanomloc(Nrloc)
C !LOCAL VARIABLES: ====================================================
      real*8 rho0, permil, permeg
      parameter (rho0=1024.5d0, permil=1.d0/rho0, permeg=1.d-6)
      INTEGER k
      real*8 atmosp,G_eq,T,S,rho
      real*8 solfun,sw_dens0
CEOP
      
      atmosp=1.d0 ! use atmospheric pressure of 1
C     calculate saturation anomaly
      DO k=1,Nrloc
        T=thetaloc(k)
        S=saltloc(k)        
        rho = sw_dens0(S,T)        
        G_eq = rho*solfun(S,T,gasID)/1.d6
        PTReqloc(k) = G_eq*atmosp
        PTRsatanomloc(k)=(PTRloc(k)-PTReqloc(k))/PTReqloc(k)
      ENDDO

	  RETURN
	  END
