
      SUBROUTINE READ_REC_XY_RS( fName,field,iRec,myIter,myThid )

      IMPLICIT NONE
C Global
C     !INPUT/OUTPUT PARAMETERS:

      CHARACTER*(*) fName
      real*8 field
      INTEGER iRec
      INTEGER myIter
      INTEGER myThid

      END
      
      SUBROUTINE EXCH_XY_RS(
     U                       phi,
     I                       myThid )
      
      IMPLICIT NONE

      real phi
      INTEGER myThid
      
      END
      
      SUBROUTINE GLOBAL_SUM_TILE_RL(
     I                       phiTile,
     O                       sumPhi,
     I                       myThid )

      IMPLICIT NONE

      real*8     phiTile
      real*8     sumPhi
      INTEGER myThid
      
      END

      SUBROUTINE BARRIER(myThid)     

      IMPLICIT NONE

      INTEGER myThid
      
      END
       