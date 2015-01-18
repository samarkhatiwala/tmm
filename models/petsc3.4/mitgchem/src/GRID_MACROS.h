C $Header: /u/gcmpack/MITgcm/model/inc/GRID_MACROS.h,v 1.13 2013/08/11 14:27:37 jmc Exp $
C $Name:  $
C
CBOP
C    !ROUTINE: GRID_MACROS.h
C    !INTERFACE:
C    include GRID_MACROS.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | GRID_MACROS.h
C     *==========================================================*
C     | These macros are used to substitute definitions for
C     | GRID.h variables for particular configurations.
C     | In setting these variables the following convention
C     | applies.
C     | undef  phi_CONST   - Indicates the variable phi is fixed
C     |                      in X, Y and Z.
C     | undef  phi_FX      - Indicates the variable phi only
C     |                      varies in X (i.e.not in X or Z).
C     | undef  phi_FY      - Indicates the variable phi only
C     |                      varies in Y (i.e.not in X or Z).
C     | undef  phi_FXY     - Indicates the variable phi only
C     |                      varies in X and Y ( i.e. not Z).
C     *==========================================================*
C     \ev
CEOP

#undef    HFACC_CONST
#undef    HFACC_FX
#undef    HFACC_FY
#undef    HFACC_FXY
#include "HFACC_MACROS.h"

#undef    RECIP_HFACC_CONST
#undef    RECIP_HFACC_FX
#undef    RECIP_HFACC_FY
#undef    RECIP_HFACC_FXY
#include "RECIP_HFACC_MACROS.h"
