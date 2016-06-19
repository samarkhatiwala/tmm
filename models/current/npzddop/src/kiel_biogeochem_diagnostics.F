
CBOP
C !ROUTINE: MODEL

C !INTERFACE: ==========================================================
      SUBROUTINE KIEL_BIOGEOCHEM_DIAGNOSTICS(Nrloc,
     &           fbgc1loc,fbgc2loc,fbgc3loc,fbgc4loc,fbgc5loc)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE

#include "BGC_PARAMS.h"
#include "BGC_DIAGNOSTICS.h"

C !INPUT PARAMETERS: ===================================================
      INTEGER Nrloc
C     Output variables
      real*8 fbgc1loc(Nrloc),fbgc2loc(Nrloc),fbgc3loc(Nrloc),
     &       fbgc4loc(Nrloc),fbgc5loc(Nrloc)
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
      ENDDO

	  RETURN
	  END
