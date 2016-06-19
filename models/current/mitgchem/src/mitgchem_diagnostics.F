#define ISDIAGROUTINE
#include "DIC_OPTIONS.h"

CBOP
C !ROUTINE: MODEL

C !INTERFACE: ==========================================================
      SUBROUTINE MITGCHEM_DIAGNOSTICS(Nrloc,myIter,myTime,
     &                 bioacloc)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "DYNVARS.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DIC_VARS.h"
#include "TMM_MITGCHEM_DIAGS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid               :: thread number
C  myIter               :: current timestep
C  myTime               :: current time
      INTEGER Nrloc
      INTEGER myIter
      real*8 myTime      
C     Output variables
      real*8 bioacloc(Nrloc)
   
#ifdef ALLOW_PTRACERS
#ifdef DIC_BIOTIC

C !LOCAL VARIABLES: ====================================================
	  INTEGER myThid,bi,bj
	  PARAMETER (myThid=1,bi=1,bj=1)
	  INTEGER k
CEOP

C     Copy diagnostic variables
      DO k=1,Nrloc
        bioacloc(k)=bioac(1,1,k)
      ENDDO      
#endif
#endif

c
	  RETURN
	  END
