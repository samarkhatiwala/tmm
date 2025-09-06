C$Header: /Users/ikriest/CVS/mops/mops_biogeochem_set_params.F,v 1.3 2016/06/03 09:28:59 ikriest Exp $
C$Name: mops-2_0 $

      SUBROUTINE MOPS_BIOGEOCHEM_SET_PARAMS(npar,parvec)

      implicit none

#include "BGC_PARAMS.h"

      integer npar
      real*8 parvec(npar)

! Third version for Volkmar's optimization tests, 16 October 2015
         ACik      = parvec(1)
         ACmuzoo   = parvec(2)
         ro2ut     = parvec(3)
         detmartin = parvec(4)
        
      RETURN
      END
