C$Header: /Users/ikriest/CVS/mops/mops_biogeochem_diagnostics.F,v 1.2 2016/06/03 09:28:59 ikriest Exp $
C$Name: mops-2_0 $


CBOP
C !ROUTINE: MODEL

C !INTERFACE: ==========================================================
      SUBROUTINE MOPS_BIOGEOCHEM_DIAGNOSTICS(Nrloc,
     &           fbgc1loc,fbgc2loc,fbgc3loc,fbgc4loc,fbgc5loc,
     &           fbgc6loc,fbgc7loc)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE

#include "BGC_PARAMS.h"
#include "BGC_DIAGNOSTICS.h"

C !INPUT PARAMETERS: ===================================================
      INTEGER Nrloc
C     Output variables
      real*8 fbgc1loc(Nrloc),fbgc2loc(Nrloc),fbgc3loc(Nrloc),
     &       fbgc4loc(Nrloc),fbgc5loc(Nrloc),
     &       fbgc6loc(Nrloc),fbgc7loc(Nrloc)
C !LOCAL VARIABLES: ====================================================
          INTEGER k
CEOP

C     Copy diagnostic variables
      DO k=1,Nrloc
        fbgc1loc(k)=f1_out(k)
        fbgc2loc(k)=f2_out(k)
        fbgc3loc(k)=f3_out(k)
        fbgc4loc(k)=f4_out(k)
        fbgc5loc(k)=f5_out(k)
        fbgc6loc(k)=f6_out(k)
        fbgc7loc(k)=f7_out(k)
      ENDDO

      RETURN
      END
