C$Header: /Users/ikriest/CVS/mops/BGC_MODEL.F,v 1.2 2015/11/17 14:18:51 ikriest Exp $
C$Name: mops-2_0 $

C The P core of the Kiel BGC model (NPZD-DOP) framework.
C Written by Iris Kriest (modified Aug-Dec 2010).
C Iris Kriest added sediment+river runoff option #ifdef SEDIMENT (Dec. 2010)
C Iris Kriest skipped option #ifdef SEDIMENT and added N-cycle with N-fixation and denitrification (Apr 2011)
C Iris Kriest modified and extended some comments (Jan 2015)

C Basis for air-sea gas exchange of O2: MIT source code, 
C modified by Samar Khatiwala, then by Iris Kriest (Aug-Nov. 2010)

C THIS FILE CONTAINS:
C SR bgc_model		- calculate P core of BGC
C SR o2_surfforcing	- calculate O2 air-sea gas exchange
C FUNCTION phi		- calculate Evans&Garcon light

C CALCULATE THE P CORE (BIOGEOCHEMISTRY) AND CHANGES IN ASS. ELEMENTS (O2,C,N); CALL EXTERNAL ROUTINES FOR ABIOTIC EXCHANGES
C CALCULATE THE P CORE (BIOGEOCHEMISTRY) AND CHANGES IN ASS. ELEMENTS (O2,C,N); CALL EXTERNAL ROUTINES FOR ABIOTIC EXCHANGES
C CALCULATE THE P CORE (BIOGEOCHEMISTRY) AND CHANGES IN ASS. ELEMENTS (O2,C,N); CALL EXTERNAL ROUTINES FOR ABIOTIC EXCHANGES
C CALCULATE THE P CORE (BIOGEOCHEMISTRY) AND CHANGES IN ASS. ELEMENTS (O2,C,N); CALL EXTERNAL ROUTINES FOR ABIOTIC EXCHANGES
C CALCULATE THE P CORE (BIOGEOCHEMISTRY) AND CHANGES IN ASS. ELEMENTS (O2,C,N); CALL EXTERNAL ROUTINES FOR ABIOTIC EXCHANGES

C PREPROCESSOR OPTIONS:
C CARBON	calculate carbon (bio)geochemistry
 
C CALLED BY:	mops_biogeochem_model

C CALLS:	SR o2_surfforcing
C       	SR car_coeffs (#ifdef CARBON)
C        	SR co2_surfforcing (#ifdef CARBON)
C      		function phi

C INPUT/ARGUMENT LIST:	
C bgc_kloc	number of vertical layers for this particular profile
C bgc_dz()	thickness of vertical layers [m]
C bgc_runoffvol()   volume gain for each box due to runoff [m3/day*e-12]
C bgc_globalrunoff  global runoff from burial [Gmol P]
C bgc_swr	daily intgrated surface irradiance [W/m2]
C bgc_tau	daylength [days]
C bgc_seaice	fraction of ice cover
C bgc_wind	wind speed [cm/hr]
C bgc_atmosp	atmospheric pressure
C bgc_theta()	temperature [degs C]
C bgc_salt()	salinity [PSU]

C INPUT/COMMON BLOCK
C biogeochemical parameters. See BGC_PARAMS.h and BGC_INI.F.
C constants for air-sea gas exchange. See BGC_PARAMS.h and BGC_INI.F.
C bgc_dt	time step length [days]. See BGC_CONTROL.h.
C bgc_keuph	number of layers in euphotic zone. See BGC_CONTROL.h.
C sediment (burial) parameters. See BGC_PARAMS.h and BGC_INI.F.
C carbon parameters. See CAR_PARAMS.h and CAR_INI.F
 
C OUTPUT/COMMON BLOCK
C bgc_tracer(,)	updated tracer fields. See BGC_PARAMS.h
C f*_out()	biogeochemical fluxes. See BGC_DIAGNOSTICS.h
C flux_bury	burial in sediment in this call to BGC_MODEL [mmol P/m2]


      SUBROUTINE BGC_MODEL(bgc_kloc,bgc_dz,     
     &            bgc_runoffvol,bgc_globalrunoff,
     &            bgc_swr,bgc_tau,bgc_seaice,bgc_wind,bgc_atmosp,
     &            bgc_theta,bgc_salt)

      implicit none

#include "BGC_PARAMS.h"
#include "BGC_CONTROL.h"
#include "BGC_DIAGNOSTICS.h"
#include "BGC_MISFIT.h"
#include "CAR_PARAMS.h"

! arrays related to physics
      integer bgc_kloc
      real*8 bgc_dz(bgc_ktotal)
      real*8 bgc_swr,bgc_tau,bgc_seaice
      real*8 bgc_theta(bgc_ktotal),bgc_salt(bgc_ktotal)
      real*8 ciz(bgc_keuph)
      real*8 bgc_wind,bgc_atmosp,bgc_pco2atm

! arrays to store the sms and other fluxes
      real*8 topo4(bgc_ktotal)
      real*8 tophy(bgc_ktotal)
      real*8 tozoo(bgc_ktotal)
      real*8 todop(bgc_ktotal)
      real*8 todet(bgc_ktotal)
      real*8 tooxy(bgc_ktotal)
      real*8 todin(bgc_ktotal)
      real*8 flux(bgc_ktotal)
      real*8 fdiv(bgc_ktotal)

! air sea gas exchange in general
      real*8 vgas660
      
! air-sea gas exchange of oxygen
      real*8 surf_oxy,o2gasex
      
! coefficients for biology
      real*8 attlim,atten,flight,flightlim,limnut,fnutlim,parfrac
      real*8 PO4,DOP,OXY,PHY,ZOO,DET,DIN
      real*8 phygrow0,phygrow,phyexu,phyloss,detloss
      real*8 graz0,graz,zooexu,zooloss
      real*8 eflux,flux_u,flux_l
      real*8 tempscale,TACmuphy,TACik
      real*8 remindet,remindop
      real*8 phi,glbygd
      
      integer k,kk,n,it,bgc_keuphloc

      real*8 bgc_runoffvol(bgc_ktotal),bgc_globalrunoff
      real*8 runoff(bgc_ktotal)
      real*8 oxymm,o2req,o2usefrac,fDET

      real*8 dinlim,nfixation,nfixtfac,nfixnfac
      real*8 dinmm,dinreq,dinusefrac
      real*8 denitdet,denitdop
      real*8 ttemp

! carbon chemistry
#ifdef CARBON      
      real*8 co2gasex,co2emp,surf_dic,surf_pho,surf_alk,surf_sil
c      real*8 dicchem(bgc_ktotal)
c      real*8 alkchem(bgc_ktotal)
      real*8 todic(bgc_ktotal)
      real*8 toalk(bgc_ktotal)
#endif

! Reset the diagnostic fluxes per ocean time step.

      bgc_keuphloc = INT(MIN(bgc_kloc,bgc_keuph))
       
      DO K=1,bgc_kloc

        f1_out(k)=0.0d0 
        f2_out(k)=0.0d0 
        f3_out(k)=0.0d0 
        f4_out(k)=0.0d0 
        f5_out(k)=0.0d0 
        f6_out(k)=0.0d0 
        f7_out(k)=0.0d0 

      ENDDO

      DO K=1,bgc_kloc

        m1_out(k)=0.0d0 
        m2_out(k)=0.0d0 
        m3_out(k)=0.0d0 

      ENDDO

! Things for air-sea gas exchange that are common to all tracers 
! (O2,CO2)
! vgas660 = exchange coefficient normalized to a Sc of 660,
! Conversion: m/day (0.24=[.01 m/cm]*[24 hr/day] converts cm/hr to m/d)
      
      vgas660=(0.337d0*bgc_wind**2)*0.24d0*(1.d0-bgc_seaice)
      
#ifdef CARBON

      CALL CAR_COEFFS(bgc_theta(1),bgc_salt(1),1)
      co2airseaflux=0.d0

#endif
      
      flux_bury = 0.0d0

C INTERNAL TIME LOOP FOR BGC

      do it=1,bgc_timesteps
      
! Reset fluxes.

      DO K=1,bgc_kloc

        topo4(k)=0.0d0
        todop(k)=0.0d0
        tooxy(k)=0.0d0
        tophy(k)=0.0d0
        tozoo(k)=0.0d0
        todet(k)=0.0d0
        flux(k) =0.0d0
        fdiv(k) =0.0d0
        todin(k)=0.0d0

#ifdef CARBON
        todic(k)=0.0d0
        toalk(k)=0.0d0
#endif

        runoff(k)=0.0d0 

      ENDDO


! AIR-SEA GAS EXCHANGE OF OXYGEN

      surf_oxy = bgc_tracer(1,ioxy)

      CALL O2_SURFFORCING(vgas660,bgc_atmosp,bgc_theta(1),bgc_salt(1),
     &     surf_oxy,o2gasex)

#ifdef CARBON

! AIR-SEA GAS EXCHANGE OF CO2

      surf_dic = bgc_tracer(1,idic)
      surf_pho = bgc_tracer(1,ipo4)

! Surface total alkalinity from the OCMIP protocol
      surf_alk = ocmip_alkfac*bgc_salt(1)

! Surface silicate from the OCMIP protocol
      surf_sil=ocmip_silfac      
      
      CALL CO2_SURFFORCING(vgas660,bgc_atmosp,
     &     surf_dic,surf_pho,surf_alk,surf_sil,bgc_theta(1),
     &     co2gasex,co2emp)
     
!     co2gasex and co2emp are in mmolC/(m^2 d)
!     co2airseaflux at the end of the internal time stepping loop will be 
!     in mmolC/(m^2 ocean_time_step)
      co2airseaflux = co2airseaflux + (co2gasex + co2emp)*bgc_dt

!      CALL CAR_CHEMISTRY(...,dicchem,alkchem)

#endif

! EUPHOTIC ZONE  AND ITS EXPORT TO DEEPER LAYERS 
! EUPHOTIC ZONE  AND ITS EXPORT TO DEEPER LAYERS 
! EUPHOTIC ZONE  AND ITS EXPORT TO DEEPER LAYERS 
! EUPHOTIC ZONE  AND ITS EXPORT TO DEEPER LAYERS 
! EUPHOTIC ZONE  AND ITS EXPORT TO DEEPER LAYERS 

! Net solar radiation at top of every layer.

      parfrac=0.4d0
      ciz(1)=bgc_swr*(1.0d0-bgc_seaice)*parfrac

      DO K=2,bgc_keuphloc

        attlim=MAX(0.0d0,bgc_tracer(k-1,iphy))
        atten = (ACkw+ACkchl*attlim)*bgc_dz(k-1)
        ciz(k)=ciz(k-1)*exp(-atten)

      ENDDO

! Biogeochemical fluxes in euphotic zone and their export
! Take care of negative tracer concentrations.

      DO K=1,bgc_keuphloc

       PO4=bgc_tracer(k,ipo4)
       DIN=bgc_tracer(k,idin)
       PHY=bgc_tracer(k,iphy)
       ZOO=bgc_tracer(k,izoo)
       attlim=MAX(PHY,0.0d0)

! temperature dependence of phytoplankton growth (Eppley)
! this affects the light-half-saturation constant via acik=acmuphy/alpha
       tempscale = EXP(bgc_theta(k)/TempB)
       TACmuphy = ACmuphy*tempscale
       TACik = ACik*tempscale       
       
! The light limitation function of phytoplankton.
! This function corresponds to Evans and Garcon, 1997.
! Note that the initial slope of the P-I curve, alpha, is ACMuPhy/ACIk
! flightlim thus gives the light limited growth rate, averaged over day 
! and layer, normalised by max. growth rate
       atten = (ACkw+ACkchl*attlim)*bgc_dz(k) !attenuation at bottom of layer
       glbygd = 2.0d0*ciz(k)/(TACik*bgc_tau)   ! 2 * G_L/G_D of EG97
       flightlim = bgc_tau/atten*(phi(glbygd)-phi(glbygd*exp(-atten)))

       if(PHY.gt.0.0d0) then

         limnut = MIN(PO4,DIN/rnp)

         if(limnut.gt.vsafe) then

! The nutrient limitation of phytoplankton
           fnutlim = limnut/(ackpo4+limnut)

! The growth rate of phytoplankton: light*nutrient limitation.
           phygrow0 = TACmuphy*PHY*MIN(flightlim,fnutlim)

! Make sure not to take up more nutrients than available.
           phygrow = MIN(limnut,phygrow0*bgc_dt)/bgc_dt

         else !limnut < vsafe

           phygrow=0.0d0

         endif !limnut

! The exudation of phytoplankton
         phyexu = AClambda * PHY

! Other losses of phytoplankton       
         phyloss = AComni * PHY * PHY

         if(ZOO.gt.0.0d0) then

! Grazing of zooplankton, Holling III
           graz0=ACmuzoo*PHY*PHY/(ACkphy*ACkphy+PHY*PHY)*ZOO

! Make sure not to graze more phytoplankton than available.
           graz = MIN(PHY,graz0*bgc_dt)/bgc_dt

         else !ZOO < 0

           graz=0.0d0

         endif !ZOO

       else !PHY < 0

         phygrow=0.0d0
         phyexu =0.0d0
         phyloss=0.0d0
         graz   =0.0d0

       endif !PHY

       if(ZOO.gt.0.0d0) then

! Zooplankton exudation
          zooexu = AClambdaz * ZOO

! Zooplankton mortality 
          zooloss = AComniz * ZOO * ZOO
       
       else !ZOO < 0

           zooexu = 0.0d0
           zooloss = 0.0d0

       endif !ZOO

! Relaxation of N:P to Redfield values (mimick cyanobacteria)

       if(PO4.gt.vsafe) then

         ttemp = bgc_theta(k)
         nfixtfac = MAX(0.0d0,tf2*ttemp*ttemp + tf1*ttemp + tf0)/tff
         dinlim = MAX(0.0d0,DIN)
         nfixnfac = MAX(0.0d0, 1.0d0-dinlim/(PO4*rnp))
         nfixation = nfixtfac*nfixnfac*nfix

       else

          nfixation = 0.0d0  

       endif  
       
! Photosynthesis stored in this array for diagnostic purposes only.
        f1_out(k) = f1_out(k)+phygrow*bgc_dt
        f2_out(k) = f2_out(k)+graz*bgc_dt
        f6_out(k) = f6_out(k)+nfixation*bgc_dt

! Collect all euphotic zone fluxes in these arrays.
        topo4(k)=-phygrow+zooexu
        todop(k)= graztodop*(1.0d0-ACeff)*graz
     &           +graztodop*(phyexu+zooloss)
     &           +phyloss
        tooxy(k)= tooxy(k)+(phygrow-zooexu)*ro2ut
        tophy(k)= phygrow-graz-phyexu-phyloss
        tozoo(k)= ACeff*graz-zooexu-zooloss
        todet(k) = (1.0d0-graztodop)*(1.0d0-ACeff)*graz
     &        + (1.0d0-graztodop)*(phyexu+zooloss)
        todin(k)=topo4(k)*rnp + nfixation

        ENDDO !loop over euphotic zone


! Explicit sinking of detritus in seperate loop. 
        flux_u = 0.0d0

        do k=1,bgc_kloc-1 !loop over all layers

          DET = MAX(bgc_tracer(k,idet)-alimit*alimit,0.0d0)
          flux_l=wdet(k)*DET
          flux(k)  = flux(k)+flux_u
          fdiv(k) = fdiv(k)+(flux_u-flux_l)/bgc_dz(k)
          flux_u=flux_l          

        enddo

        flux_l = 0.0d0

! account for burial in the sediment	
        DET = MAX(bgc_tracer(bgc_kloc,idet)-alimit*alimit,0.0d0)
        fDET = wdet(bgc_kloc)*DET
        flux_l = MIN(1.0d0,burdige_fac*fDET**burdige_exp)*fDET
        flux_bury = flux_bury + flux_l*bgc_dt

        flux(bgc_kloc) = flux(bgc_kloc)+flux_u
        fdiv(bgc_kloc) = fdiv(bgc_kloc)+(flux_u-flux_l)/bgc_dz(bgc_kloc)

! Store flux for diagnostic purposes.

        f3_out(1)  = flux_bury

        do k=2,bgc_kloc
          f3_out(k)  = f3_out(k)+flux(k)*bgc_dt
        enddo


! PROCESSES AFFECTING THE ENTIRE WATER COLUMN
! PROCESSES AFFECTING THE ENTIRE WATER COLUMN
! PROCESSES AFFECTING THE ENTIRE WATER COLUMN
! PROCESSES AFFECTING THE ENTIRE WATER COLUMN
! PROCESSES AFFECTING THE ENTIRE WATER COLUMN

      DO K=1,bgc_kloc

      DOP = MAX(bgc_tracer(k,idop)-alimit*alimit,0.0d0)
      PHY = MAX(bgc_tracer(k,iphy)-alimit*alimit,0.0d0)
      ZOO = MAX(bgc_tracer(k,izoo)-alimit*alimit,0.0d0)
      DET = MAX(bgc_tracer(k,idet)-alimit*alimit,0.0d0)

c AEROBIC DECAY

c In contrast to the older (Kriest&Oschlies, 2013) version, this option:
c (1) does not degrade OM in the absence of O2, i.e. OM can accumulate 
c (2) uses a Michaelis-Menten Kinetic to slow down bacterial remineralisation under low O2
c (2) takes care not to use more O2 per timestep than available

c Michaelis-Menten limitation for oxic degradation: 

      OXY = MAX(bgc_tracer(k,ioxy)-subox,0.0d0)
      oxymm = OXY*OXY/(OXY*OXY+ACkbaco2*ACkbaco2)

c O2 required for total remineralisation in a time step will then be:
      
      o2req = oxymm*(dlambda*DOP+detlambda*DET)*ro2ut*bgc_dt      
      
c restrict remineralisation to amount of vaialable oxygen
      
      if (o2req.gt.0.0d0) then
         o2usefrac = MIN(OXY,o2req)/o2req
      else
         o2usefrac = 0.0d0
      endif

      remindop = oxymm*dlambda*DOP*o2usefrac
      remindet = oxymm*detlambda*DET*o2usefrac

c ANAEROBIC DECAY INCL. ANAMMOX ETC.

      if(OXY.lt.36.0d0) then
      
      DIN = MAX(bgc_tracer(k,idin)-subdin,0.0d0)
      dinmm = DIN*DIN/(DIN*DIN+ACkbacdin*ACkbacdin)*(1.0d0-oxymm)

c NO3 required for total remineralisation in a time step will then be:
      
      dinreq = dinmm*(dlambda*DOP+detlambda*DET)*rhno3ut*bgc_dt
      
c restrict remineralisation to amount of variable oxygen
      
      if (dinreq.gt.0.0d0) then
         dinusefrac = MIN(DIN,dinreq)/dinreq
      else
         dinusefrac = 0.0d0
      endif

c restrict anaerobic processes to regions with low oxygen concentration

      denitdop = dinmm*dlambda*DOP*dinusefrac
      denitdet = dinmm*detlambda*DET*dinusefrac

      else
      
      denitdop = 0.0d0
      denitdet = 0.0d0

      endif
      
      topo4(k)=topo4(k)+remindop+remindet+denitdop+denitdet
      todop(k)=todop(k)-remindop-denitdop
     &                 +plambda*PHY
     &                 +zlambda*ZOO     
      tooxy(k)=tooxy(k)-(remindop+remindet)*ro2ut
      tophy(k)=tophy(k)-plambda*PHY
      tozoo(k)=tozoo(k)-zlambda*ZOO
      todet(k)=todet(k)-remindet-denitdet
      todin(k)=todin(k)+(remindop+remindet)*rnp
     &                 -(denitdop+denitdet)*rhno3ut 

      f4_out(k) = f4_out(k) + (remindop+remindet)*bgc_dt
      f7_out(k) = f7_out(k) + (denitdop+denitdet)*bgc_dt

      ENDDO


! RESUPPLY OF BURIED MATTER VIA RIVER RUNOFF OR VIA SURFACE
! RESUPPLY OF BURIED MATTER VIA RIVER RUNOFF OR VIA SURFACE
! RESUPPLY OF BURIED MATTER VIA RIVER RUNOFF OR VIA SURFACE
! RESUPPLY OF BURIED MATTER VIA RIVER RUNOFF OR VIA SURFACE
! RESUPPLY OF BURIED MATTER VIA RIVER RUNOFF OR VIA SURFACE

#ifdef RUNOFF
      DO K=1,bgc_kloc      
        runoff(k) = bgc_globalrunoff * bgc_runoffvol(k)
        f5_out(k) = f5_out(k) + runoff(k)*bgc_dt
      ENDDO
#else
      runoff(1) = bgc_globalrunoff/bgc_dz(1)
      f5_out(1) = f5_out(1) + runoff(1)*bgc_dt
#endif
 
! UPDATE MASS CONCENTRATIONS 
! UPDATE MASS CONCENTRATIONS 
! UPDATE MASS CONCENTRATIONS 
! UPDATE MASS CONCENTRATIONS 
! UPDATE MASS CONCENTRATIONS 
! UPDATE MASS CONCENTRATIONS 

      DO K=1,bgc_kloc

! Update tracer concentrations, by adding the fluxes scaled by 
! time step length.

      bgc_tracer(k,ipo4) = bgc_tracer(k,ipo4) +
     &    topo4(k)*bgc_dt + runoff(k)*bgc_dt
      bgc_tracer(k,idop) = bgc_tracer(k,idop) +
     &     todop(k)*bgc_dt
      bgc_tracer(k,ioxy) = bgc_tracer(k,ioxy) +   
     &    tooxy(k)*bgc_dt
      bgc_tracer(k,iphy)= bgc_tracer(k,iphy)  +
     &     tophy(k)*bgc_dt
      bgc_tracer(k,izoo)= bgc_tracer(k,izoo)  +
     &     tozoo(k)*bgc_dt
      bgc_tracer(k,idet)= bgc_tracer(k,idet)  +
     &     (todet(k)+fdiv(k))*bgc_dt
      bgc_tracer(k,idin)= bgc_tracer(k,idin)  +
     &     todin(k)*bgc_dt + rnp*runoff(k)*bgc_dt

#ifdef CARBON
      bgc_tracer(k,idic)= bgc_tracer(k,idic)  +
     &     rcp*topo4(k)*bgc_dt    
     &    + rcp*runoff(k)*bgc_dt
      bgc_tracer(k,ialk)= bgc_tracer(k,ialk)  +
     &       0.0d0
!     &     rcp*topo4(k)*bgc_dt    
#endif

      ENDDO 

      bgc_tracer(1,ioxy)= bgc_tracer(1,ioxy) +   
     &    o2gasex/bgc_dz(1)*bgc_dt

#ifdef CARBON
      bgc_tracer(1,idic)= bgc_tracer(1,idic) +   
     &    (co2emp+co2gasex)/bgc_dz(1)*bgc_dt

!flocalarea is localdx*dy divided by the total atmospheric volume
!      co2_atm=co2_atm+co2gasex*flocalarea*bgc_dt
      
! So far, no carbonate dissolution etc.
!      DO K=1,bgc_kloc
!        bgc_tracer(k,idic)= bgc_tracer(k,idic)   +
!     &       dicchem(k)*bgc_dt    
!        bgc_tracer(k,ialk)= bgc_tracer(k,idic)   +
!             0.0d0
!     &       alkchem(k)*bgc_dt    
!      enddo
#endif

      ENDDO !internal time loop

! For now, I am happy to take tracer concentration at the end of each time loop
! for computation of misfit; I might as well sum it up and divide by the number
! of time steps

      DO K=1,bgc_kloc
      
        m1_out(k) = bgc_tracer(k,ipo4)
        m2_out(k) = bgc_tracer(k,ioxy)
        m3_out(k) = bgc_tracer(k,idin)

      ENDDO

      RETURN

      END

!-----------------------------------------------------------------------
      
      FUNCTION phi(u)
      
      real*8 phi,u
      
!      phi= u*(0.555588d0+0.004926d0*u)/(1.0d0+0.188721d0*u)

      if(u.gt.1.0d-6) then
        phi= LOG(u+SQRT(1.0d0+u*u))-(SQRT(1.0d0+u*u)-1.0d0)/u
      else
        phi=0.0d0
      endif
      
      END

C ----------------------------------------------------------------------

C CALCULATE THE AIR-SEA GAS EXCHANGE OF O2      
C CALCULATE THE AIR-SEA GAS EXCHANGE OF O2      
C CALCULATE THE AIR-SEA GAS EXCHANGE OF O2      
C CALCULATE THE AIR-SEA GAS EXCHANGE OF O2      
C CALCULATE THE AIR-SEA GAS EXCHANGE OF O2      

C PREPROCESSOR OPTIONS:
 
C CALLED BY:	bgc_model

C CALLS:	

C INPUT/ARGUMENT LIST:	
C vgas660	exchange coefficient, depends on wind
C atmosp0       atmopheric pressure
C ttemp		surface temperature
C stemp		surface salinity
C soxy		surface oxygen

C
C INPUT/COMMON BLOCK: 
C Contants for O2 saturation calculation. See BGC_PARAMS.h and BGC_INI.F.

C OUTPUT/ARGUMENT LIST: 
C o2ex		exchange rate [mmol O2/m2/d]

      SUBROUTINE O2_SURFFORCING(vgas660,atmosp,ttemp,stemp,soxy,o2ex)

      implicit none

#include "BGC_PARAMS.h"

      real*8 vgas660,atmosp,ttemp,stemp,soxy,o2ex

! local coefficients

      real*8 SchmidtNoO2,aTT,aTK,aTS,aTS2,aTS3,aTS4,aTS5,
     &       oCnew,o2s,O2sat,Kwexch

      SchmidtNoO2=sox1+sox2*ttemp+sox3*ttemp*ttemp  
     &    + sox4*ttemp*ttemp*ttemp

      KWexch = vgas660/sqrt(SchmidtNoO2/660.0d0)

! Determine saturation O2
! using Garcia and Gordon (1992), L&O (mistake in original???)

      aTT  = 298.15d0 -ttemp
      aTK  = 273.15d0 +ttemp
      aTS  = log(aTT/aTK)
      aTS2 = aTS*aTS
      aTS3 = aTS2*aTS
      aTS4 = aTS3*aTS
      aTS5 = aTS4*aTS
      oCnew  = oA0+oA1*aTS+oA2*aTS2+oA3*aTS3+oA4*aTS4+oA5*aTS5
     &    + stemp*(oB0+oB1*aTS+oB2*aTS2+oB3*aTS3)+oC0*(stemp*stemp)
      o2s = EXP(oCnew)

! Convert from ml/l to mmol/m^3
! Note: o2 in mit is in mol/m3; I use mmol/m3, thus coonvert with 1d3
      O2sat = o2s/22391.6d0 * 1.0d3*1.0d3

! Determine flux, inc. correction for local atmos surface pressure
      o2ex = Kwexch*(atmosp*O2sat-soxy)

      RETURN

      END

