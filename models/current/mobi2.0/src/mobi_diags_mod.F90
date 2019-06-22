      MODULE mobi_diags_mod

      implicit none
      
      INTEGER :: MAXDIAGS
      PARAMETER (MAXDIAGS=100)
      
      INTEGER, SAVE :: writeFlag = 0
      INTEGER, SAVE :: numProfiles = 0
      INTEGER, SAVE :: totNumPoints = 0
      INTEGER, SAVE :: numTracerFluxDiags = 0
      INTEGER, SAVE :: num2dDiags = 0
      INTEGER, SAVE :: num3dDiags = 0

      INTEGER, SAVE :: diagsLogFileUnit=-1
      
	  CHARACTER(len=30), save :: diag2dFileNames(MAXDIAGS)
	  CHARACTER(len=30), save :: diag3dFileNames(MAXDIAGS)
	  CHARACTER(len=30), save :: diag2dNames(MAXDIAGS)
	  CHARACTER(len=30), save :: diag3dNames(MAXDIAGS)
      
      REAL, DIMENSION (:,:), ALLOCATABLE, SAVE :: diags2d
      REAL, DIMENSION (:,:), ALLOCATABLE, SAVE :: diags3d
      
      INTEGER, DIMENSION (:,:), ALLOCATABLE, SAVE :: kmtdiags

      END MODULE mobi_diags_mod
