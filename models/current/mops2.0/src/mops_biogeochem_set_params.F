C$Header: /Users/ikriest/CVS/mops/mops_biogeochem_set_params.F,v 1.3 2016/06/03 09:28:59 ikriest Exp $
C$Name: mops-2_0 $

      SUBROUTINE MOPS_BIOGEOCHEM_SET_PARAMS(npar,parvec)

      implicit none

#include "BGC_PARAMS.h"

      integer npar
      real*8 parvec(npar)
      
C BGC PARAMETERS FOR SURFACE OCEAN BIOLOGY      
C THESE MAY CHANGE DURING OPTIMISATION
C Theoretically, if I wanted to optimize ALL (well, many) parameters, I'd have to do the 
C following for MOPS:
C ! Stoichiometry
C       rnc         = parvec(1)      !redfield ratio N:P
C       rnp         = parvec(2)      !redfield ratio N:P
C       ro2ut       = parvec(3)      !redfield -02:P ratio
C       subox       = parvec(4)      !no oxic degradation below this level
C       subdin      = parvec(5)      !no oxic degradation below this level
C 
C ! Phytoplankton parameters
C       TempB       = parvec(6)
C       ACmuphy     = parvec(7)      !max. growth rate [1/day]
C       ACik        = parvec(8)      !Light half-saturation constant [W/m2]
C       ACkpo4      = parvec(9)      !half-saturation constant for PO4 uptake [mmol P/m3]
C       AClambda    = parvec(10)     !exudation rate [1/day]
C       AComni      = parvec(11)     !density dep. loss rate [m3/(mmol P * day)]
C       plambda     = parvec(12)     !phytoplankton mortality [1/day]
C  
C        ACmuzoo    =parvec(13) !max. grazing rate [1/d]
C        AClambdaz  =parvec(14) !zooplankton excretion [1/d]
C        AComniz    =parvec(15)*rnp !zooplankton mortality [m3/(mmol P * day)]
C        ACeff      =parvec(16) !assimilation efficiency
C        zlambda    =parvec(17) !zooplankton mortality
C
! DOPparameters
C        graztodop  = parvec(18) ! fraction grazing that goes into DOP
C        dlambda    = parvec(19) !DOP remineralization rate [1/day]  (SLOW recycling)
C
! Detritus parameters
C        detlambda  = parvec(20) !detritus remineralisation rate [1/d]    
C        detmartin  = parvec(21) !w = a*z+b
C
! Parameters specific to MOPS: (Kriest and Oschlies 2013, 2015)
C        burdige_fac = parvec(22) ! factor for sediment burial
C        burdige_exp = parvec(23) ! exponent for sediment burial
C        ACkbaco2 = parvec(25) !Half sat.-constant for oxic degradation.
C        ACkbacdin = parvec(26)  ! Half sat.-constant for suboxic degradation
 
C In reality, I'll probably never want/try to optimize all of these, but rather stick
C to 1-5 parameters. The following is an example of a parameter vector parvec, 
C that consists of npar=three parameters ACMUphy, graztodop and detmartin:

! Third version for Volkmar's optimization tests, 16 October 2015
         ro2ut = parvec(1)
         ACik = parvec(2)
         ACkpo4 = parvec(3)
         ACmuzoo =parvec(4)
         AComniz = parvec(5)
         detmartin = parvec(6)
        
      RETURN
      END
