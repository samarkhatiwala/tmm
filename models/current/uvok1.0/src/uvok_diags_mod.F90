#ifdef O_TMM
      MODULE uvok_diags_mod

      implicit none
      
      INTEGER :: MAXDIAGS
      PARAMETER (MAXDIAGS=100)
      
      INTEGER, SAVE :: writeFlag = 0
      INTEGER, SAVE :: numProfiles = 0
      INTEGER, SAVE :: totNumPoints = 0
      INTEGER, SAVE :: numTracerFluxDiags = 0
      INTEGER, SAVE :: num2dDiags = 0
      INTEGER, SAVE :: num3dDiags = 0
      INTEGER, SAVE :: id_F_dic=-1,id_F_dic13=-1,id_F_alk=-1
      INTEGER, SAVE :: id_F_o2=-1,id_F_po4=-1,id_F_dop=-1,id_F_no3=-1
      INTEGER, SAVE :: id_F_dfe=-1,id_F_don=-1,id_F_din15=-1
      INTEGER, SAVE :: id_F_don15=-1,id_F_doc13=-1,id_F_c14=-1
      INTEGER, SAVE :: id_O_phsur=-1,id_O_co3sur=-1,id_O_ocalcsur=-1
      INTEGER, SAVE :: id_O_oaragsur=-1,id_O_pco2sur=-1
      INTEGER, SAVE :: id_O_caco3pro=-1,id_O_sedrr=-1
      INTEGER, SAVE :: id_O_dc14
      INTEGER, SAVE :: id_O_phytnpp=-1,id_O_phytnpp_dop=-1
      INTEGER, SAVE :: id_O_phytgraz=-1,id_O_zoograz=-1
      INTEGER, SAVE :: id_O_detgraz=-1,id_O_phytmort=-1
      INTEGER, SAVE :: id_O_phytrecy=-1,id_O_zoopmort=-1
      INTEGER, SAVE :: id_O_excret=-1,id_O_avej=-1,id_O_avej_D=-1
      INTEGER, SAVE :: id_O_gmax=-1,id_O_no3P=-1,id_O_po4P=-1
      INTEGER, SAVE :: id_O_po4_D=-1,id_O_diaznpp=-1
      INTEGER, SAVE :: id_O_diaznpp_dop=-1,id_O_diazgraz=-1
      INTEGER, SAVE :: id_O_diazmort=-1,id_O_diazrecy=-1,id_O_nfix=-1
      INTEGER, SAVE :: id_O_wcdeni=-1,id_O_bdeni=-1,id_O_detrfeexpo=-1
      INTEGER, SAVE :: id_O_detrferemi=-1,id_O_feorgads=-1
      INTEGER, SAVE :: id_O_deffe=-1,id_O_feprime=-1,id_O_fesed=-1
      INTEGER, SAVE :: id_O_bfe=-1,id_O_fecol=-1,id_O_detrremi=-1
      INTEGER, SAVE :: id_O_detrexp=-1,id_O_caco3exp=-1

      INTEGER, SAVE :: diagsLogFileUnit=-1
      
	  CHARACTER(len=30), save :: diag2dFileNames(MAXDIAGS)
	  CHARACTER(len=30), save :: diag3dFileNames(MAXDIAGS)
      
      REAL, DIMENSION (:,:), ALLOCATABLE, SAVE :: diags2d
      REAL, DIMENSION (:,:), ALLOCATABLE, SAVE :: diags3d

      END MODULE uvok_diags_mod
#endif