C$Header: /Users/ikriest/CVS/mops/mops_biogeochem_set_params.F,v 1.3 2016/06/03 09:28:59 ikriest Exp $
C$Name: mops-2_0 $

      SUBROUTINE MOPS_BIOGEOCHEM_SET_PARAMS(nlen,parvec,parname)

      implicit none

#include "BGC_PARAMS.h"
      integer nlen
      real*8 parvec
      character(len=nlen) :: parname
      
C EDITED 15-AUG-2018: Sophy Oliver: Sets each parameter value one at a time according to parameter name (parname) and parameter value (parvec), so this script needs to be called for each parameter to be set
C BGC PARAMETERS FOR SURFACE OCEAN BIOLOGY      
C THESE MAY CHANGE DURING OPTIMISATION
C Theoretically, if I wanted to optimize ALL (well, many) parameters, I'd have to do the 
C following for MOPS:

! Stoichiometry
        if (parname == 'rcp') then
            rcp = parvec                !redfield ratio C:P
            print*, "The parameter rcp has been set to ", parvec
        else if (parname == 'rnp') then
            rnp = parvec                !redfield ratio N:P
            print*, "The parameter rnp has been set to ", parvec
        else if (parname == 'ro2ut') then
            ro2ut = parvec                !redfield -02:P ratio
            print*, "The parameter ro2ut has been set to ", parvec
        else if (parname == 'subox') then
            subox = parvec                !no oxic degradation below this level
            print*, "The parameter subox has been set to ", parvec
        else if (parname == 'subdin') then
            subdin = parvec                !no oxic degradation below this level
            print*, "The parameter subdin has been set to ", parvec
 
! Phytoplankton parameters
        else if (parname == 'TempB') then
            TempB = parvec
            print*, "The parameter TempB has been set to ", parvec
        else if (parname == 'ACmuphy') then
            ACmuphy = parvec            !max. growth rate [1/day]
            print*, "The parameter ACmuphy has been set to ", parvec
        else if (parname == 'ACik') then
            ACik = parvec                !Light half-saturation constant [W/m2]
            print*, "The parameter ACik has been set to ", parvec
        else if (parname == 'ACkpo4') then
            ACkpo4 = parvec                !half-saturation constant for PO4 uptake [mmol P/m3]
            print*, "The parameter ACkpo4 has been set to ", parvec
        else if (parname == 'AClambda') then
            AClambda = parvec            !exudation rate [1/day]
            print*, "The parameter AClambda has been set to ", parvec
        else if (parname == 'AComni') then
            AComni = parvec             !density dep. loss rate [m3/(mmol P * day)]
            print*, "The parameter AComni' has been set to ", parvec
        else if (parname == 'plambda') then
            plambda = parvec            !phytoplankton mortality [1/day]
            print*, "The parameter plambda has been set to ", parvec

        else if (parname == 'ACmuzoo') then
            ACmuzoo = parvec            !max. grazing rate [1/d]
            print*, "The parameter ACmuzoo has been set to ", parvec
        else if (parname == 'AClambdaz') then
            AClambdaz = parvec            !zooplankton excretion [1/d]
            print*, "The parameter AClambdaz has been set to ", parvec
        else if (parname == 'AComniz') then
            AComniz = parvec	        !zooplankton mortality [m3/(mmol P * day)]
            print*, "The parameter AComniz has been set to ", parvec
        else if (parname == 'ACeff') then
            ACeff = parvec                !assimilation efficiency
            print*, "The parameter ACeff has been set to ", parvec
        else if (parname == 'zlambda') then
            zlambda = parvec            !zooplankton mortality
            print*, "The parameter zlambda has been set to ", parvec

! DOPparameters
        else if (parname == 'graztodop') then
            graztodop = parvec            !fraction grazing that goes into DOP
            print*, "The parameter graztodop has been set to ", parvec
        else if (parname == 'dlambda') then
            dlambda = parvec            !DOP remineralization rate [1/day]  (SLOW recycling)
            print*, "The parameter dlambda has been set to ", parvec

! Detritus parameters
        else if (parname == 'detlambda') then
            detlambda = parvec            !detritus remineralisation rate [1/d]    
            print*, "The parameter detlambda has been set to ", parvec
        else if (parname == 'detmartin') then
            detmartin = parvec            !w = a*z+b
            print*, "The parameter detmartin has been set to ", parvec

! Parameters specific to MOPS: (Kriest and Oschlies 2013, 2015)
        else if (parname == 'burdige_fac') then
            burdige_fac = parvec        !factor for sediment burial
            print*, "The parameter burdige_fac has been set to ", parvec
        else if (parname == 'burdige_exp') then
            burdige_exp = parvec        !exponent for sediment burial
            print*, "The parameter burdige_exp has been set to ", parvec
        else if (parname == 'ACkbaco2') then
            ACkbaco2 = parvec            !Half sat.-constant for oxic degradation.
            print*, "The parameter ACkbaco2 has been set to ", parvec
        else if (parname == 'ACkbacdin') then
            ACkbacdin = parvec            !Half sat.-constant for suboxic degradation
            print*, "The parameter ACkbacdin has been set to ", parvec
        end if
 
C In reality, I'll probably never want/try to optimize all of these, but rather stick
C to 1-5 parameters. The following is an example of a parameter vector parvec, 
C that consists of npar=three parameters ACMUphy, graztodop and detmartin:

! Third version for Volkmar's optimization tests, 16 October 2015
 
C        ACkbaco2 = parvec(2)
C        ACkbacdin =parvec(3)
C        subdin = parvec(4)
C        nfix = parvec(5)
C        detmartin = parvec(6)
        
      RETURN
      END
