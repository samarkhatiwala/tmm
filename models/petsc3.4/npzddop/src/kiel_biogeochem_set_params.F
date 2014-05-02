      SUBROUTINE kiel_biogeochem_set_params(npar,parvec)

      implicit none

#include "BGC_PARAMS.h"

      integer npar
      real*8 parvec(npar)
      
C BGC PARAMETERS FOR SURFACE OCEAN BIOLOGY      
C THESE MAY CHANGE DURING OPTIMISATION
C 
C ! Stoichiometry
C       rnp         = parvec(1)      !redfield ratio N:P
C       ro2ut       = parvec(2)      !redfield -02:P ratio
C       subox       = parvec(3)      !no oxic degradation below this level
C 
C ! Phytoplankton parameters
C       TempB = parvec(4)
C       ACmuphy     = parvec(5)      !max. growth rate [1/day]
C       ACik        = parvec(6)      !Light half-saturation constant [W/m2]
C       ACkpo4      = parvec(7)      !half-saturation constant for PO4 uptake [mmol P/m3]
C       ACkchl      = parvec(8)*rnp  !att. of Phy [1/(m*mmol P/m3)]
C       ACkw        = parvec(9)      !attenuation of water [1/m]
C       AClambda    = parvec(10)     !exudation rate [1/day]
C       AComni      = parvec(11)     !density dep. loss rate [m3/(mmol P * day)]
C       plambda     = parvec(12)     !phytoplankton mortality [1/day]
C  
C        ACmuzoo=parvec(13) !max. grazing rate [1/d]
C        AClambdaz=parvec(14) !zooplankton excretion [1/d]
C        AComniz=parvec(15)*rnp !zooplankton mortality [m3/(mmol P * day)]
C        ACeff=parvec(16) !assimilation efficiency
C        zlambda=parvec(17) !zooplankton mortality

! DOPparameters
C        graztodop = parvec(18) ! fraction grazing that goes into DOP
C        dlambda  = parvec(19) !DOP remineralization rate [1/day]  (SLOW recycling)

! Detritus parameters
c        detlambda = parvec(20) !detritus remineralisation rate [1/d]    
c        detwa = parvec(21) !w = a*z+b
c        detwmin = parvec(22) !minimum sinking speed of detritus [m/d]

      RETURN
      END
