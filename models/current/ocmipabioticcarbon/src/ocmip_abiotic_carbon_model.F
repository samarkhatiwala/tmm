C !INTERFACE: ==========================================================
      SUBROUTINE OCMIP_ABIOTIC_CARBON_MODEL(myIter,myTime,
     &                 PTR_DIC,
#ifdef ALLOW_C14
     &                 PTR_DIC14,
#endif
     &                 PTR_ALK,PTR_PO4,PTR_SIO2,
     &                 thetaloc,saltloc,pHloc,Vgasloc,
     &                 atmosploc,pCO2atm_pp,dzsurf,
     &                 EmPloc,DICemp,
     &                 linearChemistryFactor,linearChemistryCO2,
     &                 linearChemistryDIC,
#ifdef ALLOW_C14
     &                 D14Catm,DIC14emp,
#endif     
     &                 GDICloc,gasexfluxloc,totfluxloc,pco2loc
#ifdef ALLOW_C14
     &                 ,GDIC14loc,c14gasexfluxloc,c14totfluxloc
#endif     
     &                          )

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
C  myIter               :: current timestep
C  myTime               :: current time
C  PTR_DIC              :: dissolced inorganic carbon
C  PTR_DIC14            :: DIC-normalized C14 concentration
C  PTR_ALK              :: alkalinity
C  PTR_PO4              :: phosphate
c  PTR_DOP              :: dissolve organic phosphurous
c  PTR_O2               :: oxygen
C  PTR_FE               :: iron
      INTEGER myIter
      real*8 myTime
      real*8 PTR_DIC
#ifdef ALLOW_C14
      real*8 PTR_DIC14
#endif      
      real*8 PTR_ALK,PTR_PO4,PTR_SIO2
      real*8 thetaloc,saltloc,pHloc,Vgasloc,atmosploc,pCO2atm_pp
      real*8 dzsurf,EmPloc,DICemp
      real*8 linearChemistryFactor,linearChemistryCO2
      real*8 linearChemistryDIC
#ifdef ALLOW_C14
      real*8 D14Catm,DIC14emp
#endif                  
C     Output variables
      real*8 GDICloc,gasexfluxloc,totfluxloc,pco2loc
#ifdef ALLOW_C14
      real*8 GDIC14loc,c14gasexfluxloc,c14totfluxloc
#endif            

C !LOCAL VARIABLES: ====================================================
      real*8 DIC,GDIC
      real*8 DIC14,GDIC14
      real*8 Rocn,Ratm
      real*8 pCO2atm
      real*8 co2star,co2starair,Kex
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
#ifdef ALLOW_C14
      DIC14=PTR_DIC14
#endif      

C as per MITgcm pkg/dic
#ifdef DIC_BOUNDS
      DIC = max(0.4d0, min(10.d0, PTR_DIC))
#ifdef ALLOW_C14      
      DIC14 = max(0.4d0, min(10.d0, PTR_DIC14))
#endif      
#endif

C       Kex = Vgasloc/dzsurf
C 
C       co2emp = dicgave*emp
C       
C       co2airseaflux = (co2gasex + co2emp)      

C     pCO2atm=local partial pressure of CO2 (in uatm).
C       Note: pCO2atm = [atmospheric mole fraction of CO2 in dry air (ppmv)] 
C                      x [total atmospheric pressure in atm]	  
      pCO2atm = pCO2atm_pp*atmosploc 
C      call co2calc(thetaloc,saltloc,DIC,PTR_ALK,PTR_PO4,PTR_SIO2,
C     &             pCO2atm,phlo,phhi,co2star,co2starair,pCO2surf,
C     &             dco2ddic,pHloc,1)

	  call carbon_coeffs(thetaloc,saltloc,ff,k0,k1,k2,kb,k1p,k2p,
     &                   k3p,ksi,kw,ks,kf,bt,st,ft,1)

C ---- MIT solver ------------------------------------------
C	  pHloc=8.d0
C	  do it=1,10 ! iterate to convergence
C     make sure the solver is initialized previously
		call calc_co2_approx(thetaloc,saltloc,DIC,PTR_PO4,PTR_SIO2,
     &             PTR_ALK,k1,k2,k1p,k2p,k3p,ks,kb,kw,ksi,kf,ff,bt,
     &             st,ft,pHloc,co2star)
C	  enddo
C ---- MIT solver ------------------------------------------

C ---- Newton solver ------------------------------------------
C      call co2_newton(DIC,PTR_ALK,PTR_PO4,PTR_SIO2,phlo,phhi,
C     &                pHloc,co2star)
C ---- Newton solver ------------------------------------------

C ff is the solubility (computed in S/R CARBON_COEFFS) in mol/(kg*atm)
C To convert to mol/(m^3*uatm), multiply ff by 1e-6*1024.5, i.e.
C sol = ff*permeg*rho0. 
C Equilibrium [CO2]aq in mol/m^3 = sol*pCO2_atm*atmpres, where
C            pCO2_atm = atmospheric mole fraction CO2 in dry air at 1 atm total pres (ppmv)
C            atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)

	  sol = ff*permeg*rho0 ! solubility in mol/(m^3*uatm)
	  co2starair = sol*pCO2atm ! equilibrium CO2aq in mol/m^3
C	  pCO2surf(i) = co2star(i)/sol  ! surface pCO2 in uatm
C	  htotal = 10.0**(-pH(i))
C	  dco2ddic(i) = (htotal/k2)/(dic(i)/co2star(i)) ! approximate d(co2star)/d(dic)

      if (linearChemistryFactor.gt.0.d0) then     
        co2star = linearChemistryCO2 + 
     &            linearChemistryFactor*(DIC-linearChemistryDIC)
      endif
      
      gasexfluxloc=-Vgasloc*(co2star - co2starair)
      totfluxloc=gasexfluxloc+DICemp*EmPloc
      
      GDIC = totfluxloc/dzsurf
      
#ifdef ALLOW_C14      
      Rocn=DIC14/DIC
      Ratm=(1.d0 + D14Catm/1000.d0)

      c14gasexfluxloc=-Vgasloc*(co2star*Rocn - co2starair*Ratm)
      c14totfluxloc=c14gasexfluxloc+DIC14emp*EmPloc
      GDIC14=c14totfluxloc/dzsurf
      
#endif

C     Send data back
      GDICloc=GDIC
      pco2loc=co2star/sol      
#ifdef ALLOW_C14
      GDIC14loc=GDIC14
#endif            

c
	  RETURN
	  END
