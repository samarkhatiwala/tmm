C Initialization for P core of the Kiel BGC model (NPZD-DOP) framework.
C Written by Iris Kriest (modified Aug-Nov 2010).

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
 
C CALLED BY:	kiel_biogeochem_ini

C CALLS:	
		
C INPUT/ARGUMENT LIST:	
C bgc_zu()	depth of upper boundary of vertical layers
C setDefaults	switch for parameter reinitialization

C INPUT/COMMON BLOCK:
C bgc_kmax	nominal number of vertical layers (by kiel_biogeochem_ini)

C OUTPUT/COMMON BLOCK:
C biogeochemical parameters. See BGC_PARAMS.h and BGC_INI.F.
C contants for air-sea gas exchange. See BGC_PARAMS.h and BGC_INI.F.

      SUBROUTINE BGC_INI(bgc_zu,setDefaults)

      implicit none

#include "BGC_PARAMS.h"
#include "BGC_CONTROL.h"

      real*8 bgc_zu(bgc_ktotal+1)
     
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
        rcp=106.0d0 !redfield ratio C:P
        rnp=16.0d0 !redfield ratio N:P
        ro2ut = 170.0d0  !redfield -02:P ratio
        subox = 4.0d0 !no oxic degradation below this level

! Phytoplankton parameters
        TempB = 15.65d0 !ref. temperature for T-dependent growth
        ACmuphy = 0.6d0 !max. growth rate [((mmol P/m3)/day)]
        ACik = 24.0d0 !Light half-saturation constant [W/m2]
        ACkpo4 = 0.5d0/rnp !half-saturation constant for PO4 uptake [mmol P/m3]
        ACkchl  = 0.03d0*rnp !att. of Phy [1/(m*mmol P/m3)]
        ACkw = 0.04d0 !attenuation of water [1/m]
        AClambda = 0.03d0      !exudation rate [1/day]
        AComni = 0.0d0      !density dep. loss rate [m3/(mmol P * day)]
! Note: 0.01 mmol P = 1.22 mmol C = 14.64 g C = ~ 0.3-0.4 g Chl
! If mortality is to be around 0.1 at 2 g Chl ~0.05 mmol P we have to choose 
! AComni = 0.1/0.05 = 2
        plambda = 0.01d0 !phytoplankton mortality

! Zooplankton parameters
        ACmuzoo=2.0d0 !max. grazing rate [1/d]
        AClambdaz=0.03d0 !zooplankton excretion [1/d]
        AComniz=0.2d0*rnp !zooplankton mortality [m3/(mmol P * day)]
        ACeff=0.75d0 !assimilation efficiency
        zlambda=0.01d0 !zooplankton mortality

! DOPparameters
        graztodop = 0.15d0 ! fraction grazing that goes into DOP
        dlambda  = 0.17d0/360.0d0 !DOP remineralization rate [1/day]  (SLOW recycling)

c A minimum value for the degradation of PHY and DOP
        alimit = 1.0d-3

! Detritus parameters
        detlambda = 0.05d0 !detritus remineralisation rate [1/d]    

      endif !setDefaults; anything beyond this line will be affected during optimization 
      
C Derived parameters

        ACkphy=SQRT(ACmuzoo/1.0d0)/rnp !zoo half-saturation constant [mmol P]

        detwa = detlambda/1.0725d0 !w = a*z+b for aphotic layers
        detwb = 0.0d0 !offset for detritus sinking

C REMINERALISATION LENGTH SCALES

      do k=1,bgc_kmax
        wdet(k) = detwb + (bgc_zu(k)+bgc_zu(k+1))*0.5d0*detwa
      enddo

      return
      end
