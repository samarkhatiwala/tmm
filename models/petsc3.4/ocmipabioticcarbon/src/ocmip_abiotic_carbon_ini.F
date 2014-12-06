CBOP
C !ROUTINE: OCMIP_ABIOTIC_CARBON_INI

C !INTERFACE: ==========================================================
      SUBROUTINE OCMIP_ABIOTIC_CARBON_INI(myIter,myTime,
     &                 PTR_DIC,PTR_ALK,PTR_PO4,PTR_SIO2,
     &                 thetaloc,saltloc,pHloc)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE

C !INPUT PARAMETERS: ===================================================
C  myIter               :: current timestep
C  myTime               :: current time
C  PTR_DIC              :: dissolced inorganic carbon
C  PTR_ALK              :: alkalinity
C  PTR_PO4              :: phosphate
c  PTR_SIO2              :: silicate
      INTEGER myIter
      real*8 myTime
      real*8 PTR_DIC,PTR_ALK,PTR_PO4,PTR_SIO2
      real*8 thetaloc,saltloc,pHloc
C     Output variables

C !LOCAL VARIABLES: ====================================================
      real*8 DIC
      real*8 co2star
      real*8 phlo,phhi
      parameter (phlo=6.d0, phhi=9.d0)      
      real*8 rho0, permil, permeg
      parameter (rho0=1024.5d0, permil=1.d0/rho0, permeg=1.d-6)
      real*8 k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi
	  common/carbconst/k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi      
      real*8 bt,st,ft
      common/species/bt,st,ft      
      real*8 ff,sol,htotal      
      INTEGER it
CEOP

C     Copy surface variables
      DIC=PTR_DIC

	  call carbon_coeffs(thetaloc,saltloc,ff,k0,k1,k2,kb,k1p,k2p,
     &                   k3p,ksi,kw,ks,kf,bt,st,ft,1)

C ---- MIT solver ------------------------------------------
	  pHloc=8.d0
	  do it=1,10 ! iterate to convergence
		call calc_co2_approx(thetaloc,saltloc,DIC,PTR_PO4,PTR_SIO2,
     &             PTR_ALK,k1,k2,k1p,k2p,k3p,ks,kb,kw,ksi,kf,ff,bt,
     &             st,ft,pHloc,co2star)
	  enddo
C ---- MIT solver ------------------------------------------

C ---- Newton solver ------------------------------------------
C      call co2_newton(DIC,PTR_ALK,PTR_PO4,PTR_SIO2,phlo,phhi,
C     &                pHloc,co2star)
C ---- Newton solver ------------------------------------------

	  RETURN
	  END
