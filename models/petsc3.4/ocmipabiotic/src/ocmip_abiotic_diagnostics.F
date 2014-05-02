#define DIAGROUTINE

CBOP
C !ROUTINE: MODEL

C !INTERFACE: ==========================================================
      SUBROUTINE OCMIP_ABIOTIC_DIAGNOSTICS(Nrloc,myIter,myTime,pco2loc)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE

#include "DIAGNOSTICS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
C  myIter               :: current timestep
C  myTime               :: current time
      INTEGER Nrloc
      INTEGER myIter
      real*8 myTime      
C     Output variables
      real*8 pco2loc
   
C !LOCAL VARIABLES: ====================================================
CEOP

C     Copy diagnostic variables
      pco2loc=co2star/sol
      
c
	  RETURN
	  END
