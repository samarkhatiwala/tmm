C$Header: /Users/ikriest/CVS/mops/mops_biogeochem_misfit.F,v 1.1 2015/11/17 14:18:51 ikriest Exp $
C$Name: mops-2_0 $

C !INTERFACE: ==========================================================
      SUBROUTINE MOPS_BIOGEOCHEM_MISFIT(Nrloc,
     &           mbgc1loc,mbgc2loc,mbgc3loc)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE

#include "BGC_PARAMS.h"
#include "BGC_MISFIT.h"

C !INPUT PARAMETERS: ===================================================
      INTEGER Nrloc
C     Output variables
      real*8 mbgc1loc(Nrloc),mbgc2loc(Nrloc),mbgc3loc(Nrloc)
C !LOCAL VARIABLES: ====================================================
          INTEGER k

C     Copy diagnostic variables
      DO k=1,Nrloc
        mbgc1loc(k)=m1_out(k)
        mbgc2loc(k)=m2_out(k)
        mbgc3loc(k)=m3_out(k)
      ENDDO

      RETURN
      END
