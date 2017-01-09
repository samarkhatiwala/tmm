C$Header: /Users/ikriest/CVS/mops/BGC_INI.F,v 1.3 2016/06/03 09:28:59 ikriest Exp $
C$Name: mops-2_0 $

C Initialization for P core of the Kiel BGC model (NPZD-DOP) framework.
C Written by Iris Kriest (modified Aug-Nov 2010).
C Iris Kriest added sediment+river runoff option #ifdef SEDIMENT (Dec 2010)
C Iris Kriest skipped option #ifdef SEDIMENT and added N-cycle with N-fixation and denitrification (Apr 2011)
C Iris Kriest modified and extended some comments (Jan 2015)
C Iris Kriest added -DIMPRO option for use of implicit profiles, as in Kriest and Oschlies, 2011 (April 2015)
C Iris Kriest modified slightly for use within optimization (May 2015)

C Basis for air-sea gas exchange of O2: MIT source code, 
C modified by Samar Khatiwala, then by Iris Kriest (Aug-Nov 2010)

C THIS FILE CONTAINS:
C SR bgc_ini	- calculate initilization of P core of BGC

C INITIALIZE CALCULATE THE P CORE (BIOGEOCHEMISTRY) 
C INITIALIZE CALCULATE THE P CORE (BIOGEOCHEMISTRY) 
C INITIALIZE CALCULATE THE P CORE (BIOGEOCHEMISTRY) 
C INITIALIZE CALCULATE THE P CORE (BIOGEOCHEMISTRY) 
C INITIALIZE CALCULATE THE P CORE (BIOGEOCHEMISTRY) 

C PREPROCESSOR OPTIONS:

C CALLED BY:	mops_biogeochem_ini

C CALLS:	
		
C INPUT/ARGUMENT LIST:	
C bgc_zu()	depth of upper boundary of vertical layers
C setDefaults	switch for parameter reinitialization

C INPUT/COMMON BLOCK:
C bgc_kmax	nominal number of vertical layers (by mops_biogeochem_ini)

C OUTPUT/COMMON BLOCK:
C biogeochemical parameters. See BGC_PARAMS.h and BGC_INI.F.
C contants for air-sea gas exchange. See BGC_PARAMS.h and BGC_INI.F.

      SUBROUTINE BGC_INI(bgc_zu,setDefaults)

      implicit none

#include "BGC_PARAMS.h"
#include "BGC_CONTROL.h"

      real*8 bgc_zu(bgc_ktotal+1)

#ifdef IMPRO
      real*8 anafacz, anafac0
#endif
     
      integer k,kk
      logical setDefaults
      
      if (setDefaults) then
      
C AIR-SEA GAS EXCHANGE OF O2, taken from MIT model

C define Schmidt no. coefficients for O2
C based on Keeling et al [GBC, 12, 141, (1998)]
        sox1 = 1638.0d0
        sox2 = -81.83d0
        sox3 = 1.483d0
        sox4 = -0.008004d0

C coefficients for determining saturation O2
        oA0=  2.00907d0
        oA1=  3.22014d0
        oA2=  4.05010d0
        oA3=  4.94457d0
        oA4= -2.56847d-1
        oA5=  3.88767d0
        oB0= -6.24523d-3
        oB1= -7.37614d-3
        oB2= -1.03410d-2
        oB3= -8.17083d-3
        oC0= -4.88682d-7

C NUMERICAL ISSUES RELATED TO BIOGEOCHEMISTRY

        vsafe = 1.0d-6
  
C BGC PARAMETERS FOR SURFACE OCEAN BIOLOGY      

! Stoichiometry
        rcp=117.0d0 !redfield ratio C:P
        rnp=16.0d0 !redfield ratio N:P
        ro2ut = 170.0d0  !redfield -O2:P ratio 
        subox = 1.0d0 !no oxic degradation below this level
        subdin = 4.0d0 !no denitrification below this level

! N2-Fixatioon
! Factors tf2, tf1 and tf0 are a polynomial (2nd order) 
! approximation to the functional relationship by Breitbarth et al. (2007),
! for temperature dependence of Trichodesmium growth, 
! their eq. (2), assuming that their powers relate to "e" (10^), and
! that the last but one sign has to be reversed.
! The relation will be scaled to their max. growth rate.
! Note that the second order approx. is basically similar to their
! function 2 for T-dependent nitrogen fixation multiplied by 4 
! (2 [N atoms per mole] * 12 [light hrs per day]/6 [C-atoms per N-atoms])
	    tf2 = -0.0042d0
	    tf1 = 0.2253d0
	    tf0 = -2.7819d0
	    tff = 0.2395d0  
	    nfix = 2.0d-3 ![mmol N/m3/d]; if PO4>0 and DIN = 0 and T=26.82 this corresponds 
	                  !to an integral over 100 m of ~ 200 umol N/m2/d, a 
		              !rather high value acc. to Mahaffey et al., 2005
		              !alternatively, 2 umol/m3/d = 0.08 nmol/L/hr,
		              !a rather high value acc. to Mulholland (2007)

! Phytoplankton parameters
        TempB = 15.65d0 !ref. temperature for T-dependent growth
        ACmuphy = 0.6d0 !max. growth rate [1/day]
        ACik = 9.653d0 !Light half-saturation constant [W/m2]
        ACkpo4 = 0.4995d0 !half-saturation constant for PO4 uptake [mmol P/m3]
        ACkchl  = 0.03d0*rnp !att. of Phy [1/(m*mmol P/m3)]
        ACkw = 0.04d0 !attenuation of water [1/m]
        AClambda = 0.03d0      !exudation rate [1/day]
        AComni = 0.0d0      !density dep. loss rate [m3/(mmol P * day)]
        plambda = 0.01d0 !phytoplankton mortality

! Zooplankton parameters
        ACmuzoo=1.893d0 !max. grazing rate [1/d]
        ACkphy=SQRT(ACmuzoo/1.0d0)/rnp !zoo half-saturation constant [mmol P]
        AClambdaz=0.03d0 !zooplankton excretion [1/d]
        AComniz=4.548d0 !zooplankton mortality [m3/(mmol P * day)]
        ACeff=0.75d0 !assimilation efficiency
        zlambda=0.01d0 !zooplankton mortality

! DOPparameters
        graztodop = 0.15d0 ! fraction grazing that goes into DOP
        dlambda  = 0.17d0/360.0d0 !DOP remineralization rate [1/day]  (SLOW recycling)

c A minimum value for the degradation of PHY and DOP
        alimit = 1.0d-3

! Detritus parameters
        detlambda = 0.05d0 !detritus remineralisation rate [1/d]    
        detmartin = 0.8580d0 ! Exponent for Martin curve

! Parameters specific to MOPS: (Kriest and Oschlies 2013, 2015)
        burdige_fac = 1.6828d0 ! factor for sediment burial (see Kriest and Oschlies, 2013)
        burdige_exp = 0.799d0  ! exponent for sediment burial (see Kriest and Oschlies, 2013)
        ACkbaco2 = 2.0d0  !Half sat.-constant for oxic degradation (see Kriest and Oschlies, 2015)
        ACkbacdin = 8.0d0 !Half sat.-constant for suboxic degradation (see Kriest and Oschlies, 2015)
        
      endif !setDefaults; anything beyond this line will be affected during optimization 
      
C Derived parameters

        rhno3ut = 0.8d0*ro2ut - rnp ! -HNO3:P ratio for denitrification

        detwa = detlambda/detmartin !w = a*z+b for aphotic layers
        detwb = 0.0d0 !offset for detritus sinking

C REMINERALISATION LENGTH SCALES

      do k=1,bgc_kmax
        wdet(k) = detwb + (bgc_zu(k)+bgc_zu(k+1))*0.5d0*detwa
      enddo

#ifdef IMPRO

C Use implicit profiles for particle sinking to overcome numerical diffusion,
C as explained in Kriest and Oschlies, 2011, Ocean Model., 39, 275-283.
C Only do this for deeper layers, as many other processes (physical and biological)
C beside sinking and remineralization might play a role in the euphotic zone.

      anafac0 = (1.0d0+detwa*bgc_dt)**(detlambda/detwa)-1.0d0
      do k=bgc_keuph+1,bgc_kmax
        anafacz = (bgc_zu(k)/bgc_zu(k+1))**(detlambda/detwa)
        wdet(k) = ((bgc_zu(k+1)-bgc_zu(k))*anafac0*
     &            anafacz/(1.0d0-anafacz))/bgc_dt
      enddo

#endif

      return
      end
