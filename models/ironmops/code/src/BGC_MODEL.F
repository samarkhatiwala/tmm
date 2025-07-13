C$Header: /Users/ikriest/CVS/mops/BGC_MODEL.F,v 1.2 2015/11/17 14:18:51 ikriest Exp $
C$Name: mops-2_0 $

C The P core of the Kiel BGC model (NPZD-DOP) framework.
C Written by Iris Kriest (modified Aug-Dec 2010).
C Iris Kriest added sediment+river runoff option #ifdef SEDIMENT (Dec. 2010)
C Iris Kriest skipped option #ifdef SEDIMENT and added N-cycle with N-fixation and denitrification (Apr 2011)
C Iris Kriest modified and extended some comments (Jan 2015)

C Iris Kriest added the Fe cycle with dissolved (idfe) and particulate (ipfe) iron
C The Fe cycle is largely based on Nickelsen et al. (2015) and Somes et al. (2021)
C Basic features:
C	External source of dFe by dust deposition
C   Internal source of dFe (=sink of pFe) by remineralisation of detritus and parameterised sedimentary release
C   Internal sink of dFe (= source of pFe) by detritus production, scavenging and precipitation
C   External sink of pFe through burial at the sea floor
C   Hence, internal changes in pFe relate / are proportional to changes in detritus; additionally there are the changes related to scavenging and precipitation
C   Hence, detritus has a variable Fe:P ratio while the ratio in all other organic components is constant.
C   Units of pFe and dFe are umol Fe/m3. This is also the unit for ligands.
C Changes with respect to Nickelsen et al. (2015) and Somes et al. (2021) relate to:
C 	- Constant half-saturation constant of dFe uptake by phytoplankton
C   - Fe limitation affects max. growth rate and light affinity (ACik), hence leaving initial slope of the PI-curve untouched
C   - Total phytoplankton limitation is the minimum of light, nutrient and Fe limitation
C   - Otherwise, particulate Fe consists basically of two fraction: one bound to detritus, and an additional part created through scavenging

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
     &            bgc_theta,bgc_salt,
! For Fe
     &            bgc_fedep,
     &            debug)

      implicit none

#include "BGC_PARAMS.h"
#include "BGC_CONTROL.h"
#include "BGC_DIAGNOSTICS.h"
#include "CAR_PARAMS.h"

#ifdef BUDGET
#include "BGC_BUDGET.h"
#endif

! arrays related to physics and boundary conditions
      integer bgc_kloc
      real*8 bgc_dz(bgc_ktotal)
      real*8 bgc_swr,bgc_tau,bgc_seaice
      real*8 bgc_theta(bgc_ktotal),bgc_salt(bgc_ktotal)
      real*8 ciz(bgc_keuph)
      real*8 bgc_wind,bgc_atmosp,bgc_pco2atm
      logical debug

! For Fe: boundary conditions
      real*8 bgc_fedep(bgc_ktotal)

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

! For Fe: arrays to store sms
      real*8 todfe(bgc_ktotal)
      real*8 topfe(bgc_ktotal)
      real*8 feflux(bgc_ktotal)
      real*8 fefdiv(bgc_ktotal)

! air sea gas exchange in general
      real*8 vgas660
      
! air-sea gas exchange of oxygen
      real*8 o2gasex,o2sat
      
! coefficients for biology
      real*8 attlim,atten,flight,flightlim,limnut,fnutlim,parfrac
      real*8 PO4,DOP,OXY,PHY,ZOO,DET,DIN
      real*8 phygrow0,phygrowP,phygrow,phyexu,phyloss,detloss
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
      real*8 co2gasex,co2emp,alkemp,surf_dic,surf_pho,surf_alk,surf_sil
      real*8 caco3_prod,fdiv_caco3
      real*8 todic(bgc_ktotal)
      real*8 toalk(bgc_ktotal)
#endif

! For Fe: auxiliary local variables
      real*8 felim,phygrowFe,ligands,fesed,DFE,PFE,AOU,fPFE,febury
      real*8 fea,feprime,feorgads,feprecip,o2demand

C     if (debug) then
C     DO K=1,bgc_kloc
C       write(*,*)"tr1",k,bgc_tracer(k,ipo4),bgc_tracer(k,iphy),
C    &   bgc_tracer(k,ioxy),bgc_tracer(k,idet)
C       write(*,*)"tr2",k,bgc_tracer(k,idop),bgc_tracer(k,idic),
C    &   bgc_tracer(k,ialk),bgc_tracer(k,idin)
C     ENDDO
C     endif

! Depth of the euphotic zone
      bgc_keuphloc = INT(MIN(bgc_kloc,bgc_keuph))
       
! Reset the diagnostic fluxes per ocean time step.

      DO K=1,bgc_kloc

        f1_out(k)=0.0d0 ! PP
        f2_out(k)=0.0d0 ! Grazing
        f3_out(k)=0.0d0 ! Sedimentation and burial
        f4_out(k)=0.0d0 ! Oxic remineralisation of POP and DOP
        f5_out(k)=0.0d0 ! River runoff or surface deposition of buried material
        f6_out(k)=0.0d0 ! Nitrogen fixation
        f7_out(k)=0.0d0 ! Anoxic remineralisation of POP and DOP
        f8_out(k)=0.0d0 ! Production of calcite (only for option DCARBON)
        f9_out(k)=0.0d0 ! Dissolution of calcite  (only for option DCARBON)
        f10_out(k)=0.0d0 ! Precipitation and organic scavenging of dFe  (only for option DCARBON)
        f11_out(k)=0.0d0 ! Sedimentation of pFe (layer 1: burial)
        f12_out(k)=0.0d0 ! Sedimentary release of dFe

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

#ifdef BUDGET
      budget_o2 = 0.0d0
      budget_n2 = 0.0d0
      budget_h2o = 0.0d0
      budget_co2 = 0.0d0
      budget_fe  = 0.0d0
      budget_c   = 0.0d0
#endif

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

! For Fe       
        todfe(k) = 0.0d0
        topfe(k) = 0.0d0
        feflux(k) =0.0d0
        fefdiv(k) =0.0d0

      ENDDO


! AIR-SEA GAS EXCHANGE OF OXYGEN

      CALL O2SATURATION(bgc_theta(1),bgc_salt(1),o2sat)
      CALL O2_SURFFORCING(vgas660,bgc_atmosp,bgc_theta(1),bgc_salt(1),
     &     bgc_tracer(1,ioxy),o2sat,o2gasex)

#ifdef CARBON

! AIR-SEA GAS EXCHANGE OF CO2

      surf_dic = bgc_tracer(1,idic)
      surf_pho = bgc_tracer(1,ipo4)

! Surface total alkalinity 
      surf_alk = bgc_tracer(1,ialk)

! Surface silicate from the OCMIP protocol
      surf_sil=ocmip_silfac      
      
      CALL CO2_SURFFORCING(vgas660,bgc_atmosp,
     &     surf_dic,surf_pho,surf_alk,surf_sil,bgc_theta(1),
     &     co2gasex,co2emp,alkemp)

C					if (debug) then
C							write(*,*)"surf1 ",o2sat,o2gasex,surf_dic,surf_pho,
C    &      surf_alk,surf_sil
C							write(*,*)"surf2 :",sph,co2gasex,co2emp,alkemp
C     endif
           
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
C       if (debug) then
C         write(*,*)"bio1: ",k,attlim,atten,ciz(k),ACkw,ACkchl,
C    &     bgc_dz(k-1)
C       endif  

      ENDDO

! Biogeochemical fluxes in euphotic zone and their export
! Take care of negative tracer concentrations.

      DO K=1,bgc_keuphloc

       PO4=bgc_tracer(k,ipo4)
       DIN=bgc_tracer(k,idin)
       PHY=bgc_tracer(k,iphy)
       ZOO=bgc_tracer(k,izoo)
       DFE=bgc_tracer(k,idfe)
       attlim=MAX(PHY,0.0d0)

! temperature dependence of phytoplankton growth (Eppley)
! this affects the light-half-saturation constant via acik=acmuphy/alpha
! the extension to dfe means that muphy and ik decline with declining fe, leaving alpha untouched 
       tempscale = EXP(bgc_theta(k)/TempB)
       felim  = DFE/(kfe+DFE)
       TACmuphy = ACmuphy*tempscale*felim
       TACik = ACik*tempscale*felim       
       
! The light limitation function of phytoplankton.
! This function corresponds to Evans and Garcon, 1997.
! Note that the initial slope of the P-I curve, alpha, is ACMuPhy/ACIk
! flightlim thus gives the light limited growth rate, averaged over day 
! and layer, normalised by max. growth rate
       atten = (ACkw+ACkchl*attlim)*bgc_dz(k) !attenuation at bottom of layer
       glbygd = 2.0d0*ciz(k)/(TACik*bgc_tau)   ! 2 * G_L/G_D of EG97
       flightlim = bgc_tau/atten*(phi(glbygd)-phi(glbygd*exp(-atten)))
       
C      if (debug) then
C        write(*,*)"bio2: ",k,tempscale,felim,TACmuphy,TACik,atten,
C    &      glbygd,flightlim
C      endif

       if(PHY.gt.0.0d0) then

         limnut = MIN(PO4,DIN/rnp) ! because I use the same (scaled) half-sat constant for N and P I here determine the limiting nutrient in this simple way

         if(limnut.gt.vsafe) then

! The nutrient limitation of phytoplankton
           fnutlim = limnut/(ackpo4+limnut)

! The growth rate of phytoplankton: minimum of light and nutrients; note that felim is 
! already in TACMuPhy! Thus, here we assume multiplicative limitation

           phygrow0 = TACmuphy*PHY*MIN(flightlim,fnutlim)

! Make sure not to take up more nutrients than available (mmol P/m3/d)
           phygrowFe = MIN(DFE,phygrow0*rfep*bgc_dt)/(bgc_dt*rfep) ! mmol P/m3/d
           phygrowP  = MIN(limnut,phygrow0*bgc_dt)/bgc_dt          ! mmol P/m3/d
           phygrow   = MIN(phygrowFe,phygrowP)                     ! mmol P/m3/d

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

C     if (debug) then
C       write(*,*)"bio3: ",k,phygrow0,phygrowFe,phygrowP,phygrow
C       write(*,*)"bio4: ",k,graz0,graz,zooexu,zooloss
C     endif
      
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
       
C      if (debug) then
C        write(*,*)"bio5: ",k,ttemp,nfixtfac,dinlim,nfixation
C      endif
                     
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

! For Fe:
        todfe(k)=topo4(k)*rfep + todop(k)*rfep    
        topfe(k)=todet(k)*rfep

C       if (debug) then
C         write(*,*)"bio6: ",k,topo4(k),todop(k),tooxy(k),tophy(k)
C         write(*,*)"bio7: ",k,tozoo(k),todet(k),todin(k),todfe(k)
C       endif

#ifdef BUDGET
! Note: in the original MOPS source code I never account for the O2-requirement
!       during N2-Fixation (to form HNO3). Thus, here it is an external source.
        budget_o2 = budget_o2 + 1.25d0*nfixation*bgc_dt*bgc_dz(k)
        budget_n2 = budget_n2 + nfixation*bgc_dt*bgc_dz(k)
        budget_h2o = budget_h2o 
     &             - 0.5d0*nfixation*bgc_dt*bgc_dz(k)
     &             + orh2o*(zooexu-phygrow)*bgc_dt*bgc_dz(k)
        budget_co2 = budget_co2 
     &             + rcp*(zooexu-phygrow)*bgc_dt*bgc_dz(k)
#endif


        ENDDO !loop over euphotic zone

! Explicit sinking of detritus in seperate loop. 
        flux_u = 0.0d0

        do k=1,bgc_kloc-1 !loop over all layers except last one, which accounts for burial (see below)

          DET = MAX(bgc_tracer(k,idet)-alimit*alimit,0.0d0)
          flux_l=wdet(k)*DET
          flux(k)  = flux(k)+flux_u
          fdiv(k) = fdiv(k)+(flux_u-flux_l)/bgc_dz(k)
          flux_u=flux_l          

C									if (debug) then
C											write(*,*)"flux1: ",k,DET,flux_l,flux(k),fdiv(k),flux_u
C									endif

        enddo

        flux_l = 0.0d0
        
! account for burial in the sediment	
        DET = MAX(bgc_tracer(bgc_kloc,idet)-alimit*alimit,0.0d0)
        fDET = wdet(bgc_kloc)*DET
        flux_l = MIN(1.0d0,burdige_fac*fDET**burdige_exp)*fDET
        flux_bury = flux_bury + flux_l*bgc_dt

        flux(bgc_kloc) = flux(bgc_kloc)+flux_u
        fdiv(bgc_kloc) = fdiv(bgc_kloc)+(flux_u-flux_l)/bgc_dz(bgc_kloc)

! Store flux for diagnostic purposes. Flux is at the top of a box; for box # 1 it is the burial

        f3_out(1)  = flux_bury

        do k=2,bgc_kloc
          f3_out(k)  = f3_out(k)+flux(k)*bgc_dt
        enddo

! For Fe: Explicit sinking of pFe in seperate loop 
! Sedimentary release of dFe depends on respiration divided by available oxygen
! In contrast to Dale et al. (2015) and Somes et al. (2021) I don't calculate the dependence from the ratio of Cox and O2, but rather from O2-consumption and O2
! In principle, o2demand/OXY is essentially the height of the water column [m] that would be stripped of oxygen within this time step because of the non-buried organic matter
! This lowers the value for dFe release a bit.
! 2024-04-16: Changed to omitting multiplication by bgc_dt. This is now the oxidative demand over a day. 

        o2demand = (fDET-flux_l)*ro2ut  ! the oxidative demand of all organic matter that is NOT buried within a time step [mmol O2/m2]
        OXY = MAX(bgc_tracer(bgc_kloc,ioxy)-alimit*alimit,0.0d0) ! bottom water oxygen

        flux_u = 0.0d0

        do k=1,bgc_kloc !loop over all layers - assumes no boundary parameterisation for burial (in contrast to detritus)

          PFE = MAX(bgc_tracer(k,ipfe)-alimit*alimit,0.0d0)
          flux_l=wdet(k)*PFE
          feflux(k)  = feflux(k)+flux_u
          fefdiv(k)  = fefdiv(k)+(flux_u-flux_l)/bgc_dz(k)
          flux_u=flux_l         
          
C         if (debug) then
C           write(*,*),"flux2: ",k,PFE,wdet(k),flux_l,feflux(k),
C    &        fefdiv(k)
C         endif
        enddo

! sedimentary release of dFe: in contrast to Somes et al. (2021) do NOT release more dFE than would be buried
        fesed = MIN(flux_l,fesedmax*tanh(o2demand/OXY))

#ifdef BUDGET
        budget_fe = budget_fe - flux_l*bgc_dt
#endif
    
! Store flux for diagnostic purposes. Flux is at the top of a box; for box # 1 it is the burial

        f11_out(1)  = f11_out(1) + flux_l*bgc_dt

        do k=2,bgc_kloc
          f11_out(k)  = f11_out(k)+feflux(k)*bgc_dt
        enddo

        flux_l = 0.0d0


#ifdef CARBON

! effect of CaCO3 production on alkalinity and DIC
       caco3_prod = 0.0d0
       do k=1,bgc_kloc ! IK 2023-03-16: This should be bgc_kloc!
         caco3_prod = caco3_prod + todet(k)*rcp*frac_caco3*bgc_dz(k)
         toalk(k) = toalk(k) - 2.0d0 * todet(k)*rcp*frac_caco3 ! mmol Alk/m3/d
         todic(k) = todic(k) - 1.0d0 * todet(k)*rcp*frac_caco3 ! mmol Alk/m3/d
C        if (debug) then
C          write(*,*)"toalk1: ",k,caco3_prod,toalk(k),todic(k),todet(k)
C        endif
         f9_out(k) = f9_out(k) + todet(k)*rcp*frac_caco3*bgc_dt
       enddo

! effect of CaCO3 flux divergence in aphotic layers on alkalinity and CaCO3
       do k=1,bgc_kloc-1
         fdiv_caco3 = caco3_prod*(fcaco3(k)-fcaco3(k+1))/bgc_dz(k)
         toalk(k) = toalk(k) + 2.0d0 * fdiv_caco3
         todic(k) = todic(k) + 1.0d0 * fdiv_caco3
C        if (debug) then
C          write(*,*)"toalk2: ",k,caco3_prod,toalk(k),todic(k),todet(k)
C        endif
         f8_out(k) = f8_out(k) + fdiv_caco3*bgc_dt
       enddo
       
! The bottom box is closed = no exchange with the sediment
         fdiv_caco3 = caco3_prod*fcaco3(bgc_kloc)/bgc_dz(bgc_kloc)
         toalk(bgc_kloc) = toalk(bgc_kloc) + 2.0d0 * fdiv_caco3
         todic(bgc_kloc) = todic(bgc_kloc) + 1.0d0 * fdiv_caco3
C        if (debug) then
C          write(*,*)"toalk3: ",caco3_prod,toalk(k),todic(k),todet(k)
C        endif
         f8_out(bgc_kloc) = f8_out(bgc_kloc) + fdiv_caco3*bgc_dt

#endif

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
      DFE = MAX(bgc_tracer(k,idfe)-alimit*alimit,0.0d0)
      PFE = MAX(bgc_tracer(k,ipfe)-alimit*alimit,0.0d0)

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

! For Fe: release of Fe from within detritus (O2 and NO3 dependent) and additional release
      todfe(k)=todfe(k)+(remindet+denitdet)*rfep    
     &                 +(plambda*PHY+zlambda*ZOO )*rfep
!     &                 +(remindop+denitdop)*rfep   
     &                 +detlambda*MAX(0.0d0,PFE-DET*rfep)    
      topfe(k)=topfe(k)-(remindet+denitdet)*rfep
     &                 -detlambda*MAX(0.0d0,PFE-DET*rfep)    

! For Fe: scavenging by organic matter and precipipation

      CALL O2SATURATION(bgc_theta(k),bgc_salt(k),o2sat)
      AOU = max(40.0d0,o2sat-OXY)
      ligands = max(fealpha*AOU**0.8d0+febeta*DOP**0.8d0,0.5d0)
      fea = 1.0d0+kfeleq*(ligands-DFE)
      feprime = ((fea**2.0d0+4.0d0*kfeleq*DFE)**0.5d0-fea)
     &            /(2.0d0*kfeleq)
      feorgads = kfeorg*feprime*DET**0.58d0
      feprecip = kfepre*feprime**2.0d0*tanh(max(OXY,0.0d0)) ! not in the paper, but in the code
      
      todfe(k) = todfe(k) - feorgads - feprecip
      topfe(k) = topfe(k) + feorgads + feprecip

#ifdef CARBON
      f10_out(k) = f10_out(k) + (feprecip+feorgads)*bgc_dt
#else
      f8_out(k)  = ligands
      f9_out(k)  = f9_out(k) + feorgads*bgc_dt
      f10_out(k) = f10_out(k) + feprecip*bgc_dt
#endif

#ifdef BUDGET
      budget_h2o = budget_h2o 
     &           + orh2o*(remindop+remindet)*bgc_dt*bgc_dz(k)
      budget_h2o = budget_h2o 
     &           + srh2o*(denitdop+denitdet)*bgc_dt*bgc_dz(k)
      budget_co2 = budget_co2 
     &           + rcp*(remindop+remindet)*bgc_dt*bgc_dz(k)
      budget_co2 = budget_co2 
     &           + rcp*(denitdop+denitdet)*bgc_dt*bgc_dz(k)
      budget_n2 = budget_n2 
     &          - (denitdop+denitdet)*(rhno3ut+rnp)*bgc_dt*bgc_dz(k)
#endif

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
     &     rcp*topo4(k)*bgc_dt ! this is from organic production and decay   
     &    + todic(k)*bgc_dt ! this is from CaCO3 production and dissolution
     &    + rcp*runoff(k)*bgc_dt ! this is DIC runoff
      bgc_tracer(k,ialk)= bgc_tracer(k,ialk)  -
     &     (todin(k) + topo4(k))*bgc_dt    ! this is from organic production and decay    
     &     + toalk(k)*bgc_dt    ! this is from CaCO3 production and dissolution
C        if (debug) then
C          write(*,*)"bgc1: ",k,bgc_tracer(k,idic),bgc_tracer(k,ialk),
C    &      todic(k),toalk(k)
C        endif
#endif

! For Fe - including sources from deposition and hydrothermal fluxes
      bgc_tracer(k,idfe)= bgc_tracer(k,idfe)  +
     &     (todfe(k)+bgc_fedep(k))*bgc_dt
      bgc_tracer(k,ipfe)= bgc_tracer(k,ipfe)  +
     &     (topfe(k)+fefdiv(k))*bgc_dt 

      ENDDO 

      bgc_tracer(1,ioxy)= bgc_tracer(1,ioxy) +   
     &    o2gasex/bgc_dz(1)*bgc_dt

#ifdef BUDGET
      budget_o2 = budget_o2 + o2gasex*bgc_dt
#endif

! For Fe: sedimentary release
      bgc_tracer(bgc_kloc,idfe)= bgc_tracer(bgc_kloc,idfe) +   
     &    fesed/bgc_dz(bgc_kloc)*bgc_dt

      f12_out(bgc_kloc) = f12_out(bgc_kloc) + fesed*bgc_dt 
     
#ifdef BUDGET
      do k=1,bgc_kloc
        budget_fe = budget_fe + bgc_fedep(k)*bgc_dz(k)*bgc_dt
      enddo
      budget_fe = budget_fe + fesed*bgc_dt
#endif

! Air-sea gas exchange of CO2 and E-P
#ifdef CARBON
      bgc_tracer(1,idic)= bgc_tracer(1,idic) +   
     &    (co2emp+co2gasex)/bgc_dz(1)*bgc_dt

      bgc_tracer(1,ialk)= bgc_tracer(1,ialk) +   
     &    (alkemp)/bgc_dz(1)*bgc_dt

C					if (debug) then
C							write(*,*)"bgc2: ",bgc_tracer(1,idic),bgc_tracer(1,ialk),
C    &    co2emp,co2gasex,alkemp
C					endif

#ifdef BUDGET
      budget_c = budget_c + (co2emp+co2gasex)*bgc_dt
#endif

#endif

      ENDDO !internal time loop

C     if (debug) then
C     DO K=1,bgc_kloc
C       write(*,*)"tr3",k,bgc_tracer(k,ipo4),bgc_tracer(k,iphy),
C    &   bgc_tracer(k,ioxy),bgc_tracer(k,idet)
C       write(*,*)"tr4",k,bgc_tracer(k,idop),bgc_tracer(k,idic),
C    &   bgc_tracer(k,ialk),bgc_tracer(k,idin)
C     ENDDO
C     endif

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
C ssat     oxygen saturation

C
C INPUT/COMMON BLOCK: 
C Contants for O2 saturation calculation. See BGC_PARAMS.h and BGC_INI.F.

C OUTPUT/ARGUMENT LIST: 
C o2ex		exchange rate [mmol O2/m2/d]

      SUBROUTINE O2_SURFFORCING(vgas660,atmosp,ttemp,stemp,
     &                          soxy,ssat,o2ex)

      implicit none

#include "BGC_PARAMS.h"

      real*8 vgas660,atmosp,ttemp,stemp,soxy,ssat,o2ex

! local coefficients

      real*8 SchmidtNoO2,Kwexch

      SchmidtNoO2=sox1+sox2*ttemp+sox3*ttemp*ttemp  
     &    + sox4*ttemp*ttemp*ttemp

      KWexch = vgas660/sqrt(SchmidtNoO2/660.0d0)

! Determine flux, inc. correction for local atmos surface pressure
      o2ex = Kwexch*(atmosp*ssat-soxy)

      RETURN

      END

C ----------------------------------------------------------------------

C CALCULATE O2 SATURATION      
C CALCULATE O2 SATURATION      
C CALCULATE O2 SATURATION      
C CALCULATE O2 SATURATION      
C CALCULATE O2 SATURATION      

C PREPROCESSOR OPTIONS:
 
C CALLED BY:	bgc_model

C CALLS:	

C INPUT/ARGUMENT LIST:	
C ttemp		temperature
C stemp		salinity

C
C INPUT/COMMON BLOCK: 
C Contants for O2 saturation calculation. See BGC_PARAMS.h and BGC_INI.F.

      SUBROUTINE o2saturation(ttemp,stemp,o2sat)

      implicit none

#include "BGC_PARAMS.h"

      real*8 ttemp,stemp,o2sat

! local coefficients

      real*8 aTT,aTK,aTS,aTS2,aTS3,aTS4,aTS5,oCnew,o2s

! Determine saturation O2 using Garcia and Gordon (1992), L&O (mistake in original???)

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
      o2sat = o2s/22391.6d0 * 1.0d3*1.0d3

      RETURN

      END

