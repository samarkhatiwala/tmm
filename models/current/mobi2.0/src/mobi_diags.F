      SUBROUTINE MOBI_DIAGS_INI(lNumProfiles, lTotNumPoints, 
     &  lNum2dDiags, lNum3dDiags, debugFlag)

      use mobi_diags_mod 
       
      implicit none

#include "size.h"

	  integer :: lNumProfiles, lTotNumPoints, debugFlag
	  integer :: lNum2dDiags, lNum3dDiags

#include "pconst.h"
#include "mw.h"
#include "levind.h"
#include "npzd.h"
#include "timeavgs.h"
#include "switch.h"
#include "mobi_diags.h"

	  integer :: ierr, ip, idd, idn, n
	  integer :: it(10), iou
      integer id_xt, id_xu, id_yt, id_yu, id_zt, id_zw, id_zl, id_xt_e

      real c100, c500, c1e3, c1e4, c1e6, c1e20

      c100 = 100.
      c500 = 500.
      c1e3 = 1.e3
      c1e4 = 1.e4
      c1e6 = 1.e6
      c1e20 = 1.e20
      
!   First store some TMM variables we'll need
      writeFlag=debugFlag
	  numProfiles=lNumProfiles
	  totNumPoints=lTotNumPoints

	  ALLOCATE( kmtdiags(imt,jmt) )
      kmtdiags(:,:) = kmt(:,:)
      
!     Set some flags so UVic doesn't accumulate diagnostics
      euler2 = .false.
      timavgperts = .false.
      tsiperts = .false.

	  CALL FINDUNIT(diagsLogFileUnit)
	  OPEN(diagsLogFileUnit, file='available_diagnostics.txt',
     &     status='unknown', form='formatted')

C     Code based on S/R mom_tavg_def in mom_tavg.F
C     Define 2-d diagnostics (what is called 3d in mom_tavg.F)
C     Surface fluxes of tracers 
      iou=1 ! flag for 2-d diagnostic
      it(1)=0  ! index for diagnostic (incremented within defvar)
      do n=3,nt
        if (trim(mapt(n)) .eq. 'temp') then
          call defvar ('F_heat', iou, 3, it, -c1e6, c1e6, ' ', 'F'
     &,     'surface downward heat flux', ' ', 'W m-2')
        elseif (trim(mapt(n)) .eq. 'salt') then
          call defvar ('F_salt', iou,3, it, -c100, c100, ' ', 'F'
     &,     'surface downward salt flux', ' ', 'kg m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'dic') then
          call defvar ('F_dic', iou,3, it, -c100, c100, ' ', 'F'
     &,   'surface downward carbon flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'dic13') then
          call defvar ('F_dic13', iou,3, it, -c100, c100, ' ', 'F'
     &,   'surface downward carbon 13 flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'alk') then
          call defvar ('F_alk', iou,3, it, -c100, c100, ' ', 'F'
     &,   'surface downward alkalinity flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'o2') then
          call defvar ('F_o2', iou,3, it, -c100, c100, ' ', 'F'
     &,   'surface downward oxygen flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'po4') then
          call defvar ('F_po4', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward phosphate flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'dop') then
          call defvar ('F_dop', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward DOP flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'no3') then
          call defvar ('F_no3', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward nitrate flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'din15') then
          call defvar ('F_din15', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward nitrate 15 flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'dfe') then
          call defvar ('F_dfe', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward iron flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'don') then
          call defvar ('F_don', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward DON flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'don15') then
          call defvar ('F_don15', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward DON 15 flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'doc13') then
          call defvar ('F_doc13', iou,3, it, -c100, c100, ' ', 'F'
     &,   'surface downward DOC 13 flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'c14') then
          call defvar ('F_c14', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward carbon 14 flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'cfc11') then
          call defvar ('F_cfc11', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward CFC11 flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'cfc12') then
          call defvar ('F_cfc12', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward CFC12 flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'phyt') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'phytn15') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'phytc13') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'zoop') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'zoopn15') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'zoopc13') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'diaz') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'diazn15') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'diazc13') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detr') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detrn15') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detrc13') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detrfe') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'cocc') then
!         skip since no surface flux
        elseif (trim(mapt(n)) .eq. 'coccn15') then
!         skip since no surface flux
        elseif (trim(mapt(n)) .eq. 'coccc13') then
!         skip since no surface flux
        elseif (trim(mapt(n)) .eq. 'caco3') then
!         skip since no surface flux
        elseif (trim(mapt(n)) .eq. 'caco3c13') then
!         skip since no surface flux
        elseif (trim(mapt(n)) .eq. 'detr_B') then
!         skip since no surface flux
        else
C         if (n .lt. 1000) write(a3,'(i3)') n
C         if (n .lt. 100) write(a3,'(i2)') n
C         if (n .lt. 10) write(a3,'(i1)') n
C         call defvar ('F_'//trim(a3), iou ,3, it, -c1e6, c1e6, ' '
C    &,     'F', 'tracer flux '//trim(a3)
C    &,     'tracer_flux_'//trim(a3), 'unknown')
        endif
      enddo

	  numTracerFluxDiags = it(1)

C     Other 2-d variables:
# if defined O_save_carbon_carbonate_chem
      call defvar ('O_phsur', iou, 3, it, c0, c100, ' ', 'F'
     &, 'sea surface pH', ' ', '1')
      call defvar ('O_co3sur', iou, 3, it, c0, c100, ' ', 'F'
     &, 'sea surface CO3', ' ', 'mol m-3')
      call defvar ('O_ocalcsur', iou, 3, it, c0, c100, ' ', 'F'
     &, 'sea surface Omega calcite', ' ', '1')
      call defvar ('O_oaragsur', iou, 3, it, c0, c100, ' ', 'F'
     &, 'sea surface Omega aragonite', ' ', '1')
      call defvar ('O_pco2sur', iou, 3, it, c0, c100, ' ', 'F'
     &, 'sea surface pCO2', ' ', '1')
# endif
# if defined O_save_npzd
#  if !defined O_npzd_caco3
      call defvar ('O_caco3pro', iou, 3, it, c0, c100, ' ', 'F'
     &, 'calcite production', ' ', 'mol C m-2 s-1')
#  endif
      call defvar ('O_sedrr', iou, 3, it, -c1, c100, ' ', 'F'
     &, 'detrital sediment rain rate', ' ', 'mol N m-2 s-1')
# endif

	  num2dDiags = it(1) ! total number of 2-d diagnostics
	  lNum2dDiags = num2dDiags

	  ALLOCATE( diags2d(numProfiles,num2dDiags) )	
	  diags2d(:,:) = 0.0

C     Define 3-d diagnostics (what is called 4d in mom_tavg.F)
      iou=2 ! flag for 3-d diagnostic
      it(1)=0  ! index for diagnostic (incremented within defvar)
#if defined O_carbon_14
          call defvar ('O_dc14', iou, 4, it, -c1e3, c1e3, ' ', 'F'
     &,     'delta carbon 14', ' ', 'permil')
#endif 

# if defined O_save_npzd
      it(3) = id_zl
# else
      it(3) = id_zt
# endif     
# if defined O_save_npzd
      call defvar ('O_phytnpp', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'ocean net primary production rate', ' ', 'mol N m-3 s-1')
#  if defined O_npzd_nitrogen
      call defvar ('O_phytnpp_dop', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'ocean net primary production rate from DOP', ' '
     &, 'mol N m-3 s-1')
#   if defined O_npzd_caco3
      call defvar ('O_calnpp_dop', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'ocean net primary production rate from DOP', ' '
     &, 'mol N m-3 s-1')
#   endif      
#  endif
      call defvar ('O_phytgraz', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'grazing rate on P', ' ', 'mol N m-3 s-1')
      call defvar ('O_zoograz', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'grazing rate on Z', ' ', 'mol N m-3 s-1')
      call defvar ('O_detgraz', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'grazing rate on detritus', ' ', 'mol N m-3 s-1')
      call defvar ('O_phytmort', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'phytoplankton mortality rate', ' ', 'mol N m-3 s-1')
      call defvar ('O_phytrecy', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'microbial fast recycling rate', ' ', 'mol N m-3 s-1')
      call defvar ('O_zoopmort', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'zooplankton mortality rate', ' ', 'mol N m-3 s-1')
      call defvar ('O_excret', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'excretion rate', ' ', 'mol N m-3 s-1')
#  if defined O_npzd_caco3
      call defvar ('O_coccnpp', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'ocean cocc npp rate', ' ', 'mol N m-3 s-1')
      call defvar ('O_coccgraz', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'cocc grazing rate', ' ', 'mol N m-3 s-1')
      call defvar ('O_coccmort', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'cocc mortality rate', ' ', 'mol N m-3 s-1')
      call defvar ('O_coccrecy', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'cocc microbial fast recycling rate', ' ', 'mol N m-3 s-1')
      call defvar ('O_caco3pro', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'calcite production',' ', 'mol C m-3 s-1')
      call defvar ('O_caco3att', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'attached calcite',' ', 'mol C m-3 s-1')
#  endif
#  if defined O_kk_ballast
      call defvar ('O_detgraz_B', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'grazing rate on CaCO3 ballast detritus', ' ', 'mol N m-3 s-1')
#  endif
#  if defined O_npzd_extra_diagnostics
      call defvar ('O_avej', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'light-dependant phyt. growth rate', ' ', 'd-1')
      call defvar ('O_gmax', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'temp-dependant zoo growth rate', ' ', 'd-1')
#   if defined O_npzd_nitrogen      
      call defvar ('O_no3P', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'no3-dependant phyt. growth rate', ' ', 'd-1')
      call defvar ('O_avej_D', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'light-dependant Diaz growth rate', ' ', 'd-1')
      call defvar ('O_po4_D', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'po4-dependant Diaz growth rate', ' ', 'd-1')
#   endif      
      call defvar ('O_po4P', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'po4-dependant phyt. growth rate', ' ', 'd-1')
#  endif
#  if defined O_npzd_nitrogen
      call defvar ('O_diaznpp', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'net primary production diazotrophs', ' ', 'mol N m-3 s-1')
      call defvar ('O_diaznpp_dop', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'net primary production diazotrophs from DOP', ' '
     &, 'mol N m-3 s-1')
      call defvar ('O_diazgraz', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'grazing rate on diazotrophs', ' ', 'mol N m-3 s-1')
      call defvar ('O_diazmort', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'mortality diazotrophs', ' ', 'mol N m-3 s-1')
      call defvar ('O_diazrecy', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'fast-recycling diazotrophs', ' ', 'mol N m-3 s-1')
      call defvar ('O_nfix', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'nitrogen fixation', ' ', 'mol N m-3 s-1')
      it(3) = id_zt
#   if defined O_npzd_o2      
      call defvar ('O_wcdeni', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'water column denitrification', ' ', 'mol N m-3 s-1')
      call defvar ('O_bdeni', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'benthic denitrification', ' ', 'mol N m-3 s-1')
#   endif      
#  endif
#  if defined O_npzd_iron
      call defvar ('O_detrfeexpo', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'particulate Fe export rate', ' ', 'mol Fe m-3 s-1')
      call defvar ('O_detrferemi', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'Fe remineralization rate', ' ', '-')
#   if defined O_npzd_iron_diagnostics
      call defvar ('O_feorgads', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'iron organic scavenging rate', ' ', 'mol Fe m-3 s-1')
      call defvar ('O_deffe', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'iron deficiency', ' ', '-')
      call defvar ('O_feprime', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'free iron concentration', ' ', 'mol Fe m-3')
      call defvar ('O_fesed', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'iron sediment release (without iron from OMZs)', ' '
     &, 'mol Fe m-3 s-1')
      call defvar ('O_bfe', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'iron sediment release rate (with iron from OMZs)', ' '
     &, 'mol Fe m-3 s-1')
      call defvar ('O_fecol', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'iron inorganic scavenging rate', ' ', 'mol Fe m-3 s-1')
#    if defined O_npzd_caco3
      call defvar ('O_deffe_C', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'iron deficiency coccs', ' ', '-')
#    endif      
#   endif
#  endif
      it(3) = id_zt
      call defvar ('O_detrremi', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'remineralisation rate', ' ', 'mol N m-3 s-1')
      it(3) = id_zt
      call defvar ('O_detrexp', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'detrital export rate', ' ', 'mol N m-3 s-1')
      it(3) = id_zt
      call defvar ('O_caco3exp', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'calcite export rate', ' ', 'mol C m-3 s-1')
#  if defined O_kk_ballast
      it(3) = id_zt
      call defvar ('O_detrremi_B', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'CaCO3 ballasted remin rate', ' ', 'mol N m-3 s-1')
      it(3) = id_zt
      call defvar ('O_detrexp_B', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'CaCO3 ballasted detrital export rate', ' ', 'mol N m-3 s-1')
#  endif
#  if defined O_npzd_caco3
      it(3) = id_zt
      call defvar ('O_caco3dissl', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'calcite dissolution rate', ' ', 'mol C m-3 s-1')
      it(3) = id_zt
#  endif
# endif

	  num3dDiags = it(1) ! total number of 3-d diagnostics
      lNum3dDiags = num3dDiags
      
	  ALLOCATE( diags3d(totNumPoints, num3dDiags) )
	  diags3d(:,:) = 0.0

      CLOSE(diagsLogFileUnit)

      CALL MOBI_DIAGS_RESET(debugFlag)

	  END

      SUBROUTINE MOBI_DIAGS_START(debugFlag)

      use mobi_diags_mod 
       
      implicit none

#include "size.h"

	  INTEGER :: debugFlag

#include "pconst.h"
#include "mw.h"
#include "npzd.h"
#include "timeavgs.h"
#include "switch.h"
#include "mobi_diags.h"

!     Set some flags so UVic accumulates diagnostics
      euler2 = .false.
      timavgperts = .true.
      tsiperts = .true.

      CALL MOBI_DIAGS_RESET(debugFlag)
    
	  END

      SUBROUTINE MOBI_DIAGS_RESET(debugFlag)

      use mobi_diags_mod 
       
      implicit none

#include "size.h"

	  INTEGER :: debugFlag

#include "pconst.h"
#include "mw.h"
#include "npzd.h"
#include "timeavgs.h"
#include "switch.h"
#include "mobi_diags.h"

!     Reset arrays to zero
CSPK  Based on timeavgs.F

#  if defined O_carbon && defined O_carbon_14
      ta_dc14(:,:,:) = c0
#  endif
#  if defined O_save_npzd
      ta_rnpp(:,:,:) = c0
      ta_rgraz(:,:,:) = c0
      ta_rgraz_Z(:,:,:) = c0
      ta_rgraz_Det(:,:,:) = c0
      ta_rmorp(:,:,:) = c0
      ta_rmorpt(:,:,:) = c0
      ta_rmorz(:,:,:) = c0
      ta_rexcr(:,:,:) = c0
      ta_rremi(:,:,:) = c0
      ta_rexpo(:,:,:) = c0
      ta_rexpocal(:,:,:) = c0
      ta_rsedrr(:,:) = c0
#   if defined O_npzd_nitrogen
      ta_rnpp_dop(:,:,:) = c0
      ta_rnpp_D(:,:,:) = c0
      ta_rnpp_D_dop(:,:,:) = c0
      ta_rgraz_D(:,:,:) = c0
      ta_rmorp_D(:,:,:) = c0
      ta_rmorpt_D(:,:,:) = c0
      ta_rwcdeni(:,:,:) = c0
      ta_rbdeni(:,:,:) = c0
      ta_rnfix(:,:,:) = c0
#    if defined O_npzd_caco3
      ta_rnpp_C_dop(:,:,:) = c0
#    endif      
#   endif
#   if defined O_npzd_iron
      ta_rexpofe(:,:,:) = c0
      ta_rremife(:,:,:) = c0
#    if defined O_npzd_iron_diagnostics
      ta_rfeorgads(:,:,:) = c0
      ta_rdeffe(:,:,:) = c0
      ta_rfeprime(:,:,:) = c0
      ta_rfesed(:,:,:) = c0
      ta_rbfe(:,:,:) = c0 
      ta_rfecol(:,:,:) = c0 
#     if defined O_npzd_caco3
      ta_rdeffe_C(:,:,:) = c0
#     endif      
#    endif
#   endif
#   if defined O_kk_ballast
      ta_rgraz_Det_B(:,:,:) = c0
      ta_rexpo_B(:,:,:) = c0
      ta_rremi_B(:,:,:) = c0
#   endif
#   if defined O_npzd_caco3
      ta_rcalatt(:,:,:) = c0
      ta_rdissl(:,:,:) = c0
      ta_rexpocaco3(:,:,:) = c0
      ta_rnpp_C(:,:,:) = c0
      ta_rgraz_C(:,:,:) = c0
      ta_rmorp_C(:,:,:) = c0
      ta_rmorpt_C(:,:,:) = c0
      ta_rprocal(:,:,:) = c0
#   else
      ta_rprocal(:,:) = c0
#   endif
#   if defined O_npzd_extra_diagnostics
      ta_ravej(:,:,:) = c0
      ta_ravej_D(:,:,:) = c0
      ta_rgmax(:,:,:) = c0
      ta_rno3P(:,:,:) = c0
      ta_rpo4P(:,:,:) = c0
      ta_rpo4_D(:,:,:) = c0
#   endif
#  endif
#  if defined O_save_carbon_carbonate_chem
      ta_sspH(:,:) = c0
      ta_ssCO3(:,:) = c0
      ta_ssOc(:,:) = c0
      ta_ssOa(:,:) = c0
      ta_sspCO2(:,:) = c0
#  endif
          
	  END

      SUBROUTINE MOBI_DIAGS_STOP(debugFlag)

      use mobi_diags_mod 
       
      implicit none

#include "size.h"

	  INTEGER :: debugFlag

#include "pconst.h"
#include "mw.h"
#include "npzd.h"
#include "timeavgs.h"
#include "switch.h"
#include "mobi_diags.h"

!     Set some flags so UVic doesn't accumulate diagnostics
      euler2 = .false.
      timavgperts = .false.
      tsiperts = .false.

	  END
        
      SUBROUTINE MOBI_DIAGS_ACCUMULATE(numAvg, doAverage, debugFlag)

      use mobi_diags_mod
      
      implicit none

#include "size.h"

      integer :: numAvg, doAverage, debugFlag

#include "pconst.h"
#include "mw.h"
#include "levind.h"
#include "npzd.h"
#include "diaga.h"
#include "timeavgs.h"
#include "mobi_diags.h"

      integer is, ie, jj

CSPK  These ta_XX variables are set in diago.F. But as this routine
C     is not called by TMM-MOBI we set them here.
	  ta_sspH(1:numProfiles,1) = sspH(1:numProfiles,1)
	  ta_ssCO3(1:numProfiles,1) = ssCO3(1:numProfiles,1)
	  ta_ssOc(1:numProfiles,1) = ssOc(1:numProfiles,1)
	  ta_ssOa(1:numProfiles,1) = ssOa(1:numProfiles,1)
	  ta_sspCO2(1:numProfiles,1) = sspCO2(1:numProfiles,1)

CSPK  The array ta_rsedrr does not seem to be set anywhere in the 
C     current version of MOBI. I think this is a bug and it should 
C     be set to rsedrr for each horizontal location in tracer.F. 
C     Until this is fixed it will be zero in the output file.

CSPK  This is based on timeavgs.F 
        is=1
        ie=imt
        jj=0

        call mom_tavg_out (numAvg, doAverage, debugFlag
     &,   is, ie, jj+1, jj+1, imt, jmt, km, nt
     &,   mapt
     &,   t, stf
#  if defined O_carbon && defined O_carbon_14
     &,   ta_dc14(is,1,jj+1)
#  endif
#  if defined O_save_npzd
#   if defined O_save_carbon_totals
     &,   redctn
#   endif
     &,   kpzd, ta_rnpp(is,1,jj+1)
     &,   ta_rgraz(is,1,jj+1)
     &,   ta_rmorp(is,1,jj+1), ta_rmorpt(is,1,jj+1)
     &,   ta_rmorz(is,1,jj+1), ta_rremi(is,1,jj+1)
     &,   ta_rexcr(is,1,jj+1), ta_rexpo(is,1,jj+1)
     &,   ta_rexpocal(is,1,jj+1)
#    if !defined O_npzd_caco3 
     &,   ta_rprocal(is,jj+1)
#    else
     &,   ta_rprocal(is,1,jj+1)        
#    endif        
     &,   ta_rsedrr(is,jj+1)
     &,   ta_rgraz_Z(is,1,jj+1), ta_rgraz_Det(is,1,jj+1)
#   if defined O_kk_ballast
     &,   ta_rexpo_B(is,1,jj+1), ta_rremi_B(is,1,jj+1)
     &,   ta_rgraz_Det_B(is,1,jj+1)
#   endif
#   if defined O_npzd_nitrogen
     &,   ta_rnpp_D(is,1,jj+1), ta_rnpp_D_dop(is,1,jj+1)
     &,   ta_rgraz_D(is,1,jj+1), ta_rnpp_dop(is,1,jj+1)
     &,   ta_rmorp_D(is,1,jj+1), ta_rmorpt_D(is,1,jj+1)
     &,   ta_rnfix(is,1,jj+1), ta_rwcdeni(is,1,jj+1)
     &,   ta_rbdeni(is,1,jj+1) 
#    if defined O_npzd_caco3
     &,   ta_rnpp_C_dop(is,1,jj+1)
#    endif        
#   endif
#   if defined O_npzd_caco3
     &,   ta_rnpp_C(is,1,jj+1), ta_rgraz_C(is,1,jj+1)
     &,   ta_rmorp_C(is,1,jj+1), ta_rmorpt_C(is,1,jj+1)
     &,   ta_rdissl(is,1,jj+1),   ta_rexpocaco3(is,1,jj+1)
     &,   ta_rcalatt(is,1,jj+1)
#   endif
#   if defined O_npzd_iron
     &,   ta_rexpofe(is,1,jj+1), ta_rremife(is,1,jj+1)
#    if defined O_npzd_iron_diagnostics
     &,   ta_rfeorgads(is,1,jj+1)
     &,   ta_rdeffe(is,1,jj+1), ta_rfeprime(is,1,jj+1)
     &,   ta_rfesed(is,1,jj+1), ta_rbfe(is,1,jj+1)
     &,   ta_rfecol(is,1,jj+1)
#     if defined O_npzd_caco3
     &,   ta_rdeffe_C(is,1,jj+1)
#     endif
#    endif
#   endif
#   if defined O_npzd_extra_diagnostics
     &,    ta_ravej(is,1,jj+1), ta_ravej_D(is,1,jj+1)
     &,    ta_rgmax(is,1,jj+1), ta_rno3P(is,1,jj+1)
     &,    ta_rpo4P(is,1,jj+1), ta_rpo4_D(is,1,jj+1)
#   endif 
#  endif
#  if defined O_save_carbon_carbonate_chem
     &,   ta_sspH(is,jj+1), ta_ssCO3(is,jj+1), ta_ssOc(is,jj+1)
     &,   ta_ssOa(is,jj+1), ta_sspCO2(is,jj+1)
#  endif
     &,   kmt(is,jj+1))

      CALL MOBI_DIAGS_RESET(debugFlag)

      END
      
      subroutine mom_tavg_out (numAvg, doAverage, debugFlag
     &,                        ids, ide, jds, jde, imt, jmt, km, nt
     &,                        mapt
     &,                        t, stf
# if defined O_carbon && defined O_carbon_14
     &,                        dc14
# endif
# if defined O_save_npzd
#  if defined O_save_carbon_totals
     &,                        redctn
#  endif
     &,                        kpzd, npp, graz, morp, morpt
     &,                        morz, remi
     &,                        excr, expo, expocal, procal, sedrr 
     &,                        graz_Z, graz_Det
#  if defined O_kk_ballast
     &,                        expo_B,remi_B,graz_Det_B
#  endif
#  if defined O_npzd_nitrogen
     &,                        npp_D, npp_D_dop, graz_D, npp_dop, morp_D
     &,                        morpt_D
     &,                        nfix, wcdeni, bdeni
#   if defined O_npzd_caco3
     &,                        npp_C_dop
#   endif      
#  endif
#  if defined O_npzd_caco3
     &,                        npp_C, graz_C, morp_C, morpt_C
     &,                        dissl, expocaco3, calatt
#  endif
#  if defined O_npzd_iron
     &,                        expofe, remife
#     if defined O_npzd_iron_diagnostics
     &,                        feorgads
     &,                        deffe, feprime
     &,                        fesed, bfe, fecol
#      if defined O_npzd_caco3
     &,                        deffe_C
#      endif
#     endif
#  endif
#  if defined O_npzd_extra_diagnostics
     &,                        avej, avej_D, gmax, no3P, po4P, po4_D
#  endif
# endif
# if defined O_save_carbon_carbonate_chem
     &,                        sspH, ssCO3, ssOc, ssOa, sspCO2
# endif
     &,                        kmt)

      use mobi_diags_mod
      
      implicit none

      integer numAvg, doAverage, debugFlag
      integer iou(5), j, ln, n, ntrec, imt, jmt, km, nt, kpzd, ids, ide
      integer jds, jde, igs, ige, ig, jgs, jge, jg, kgs, kge, kg, lgs
      integer lge, lg, ils, ile, jls, jle, kls, kle, lls, lle, ib(10)
      integer ic(10), kmt(ids:ide,jds:jde)

      character(10) :: mapt(nt)

      real t(ids:ide,jds:jde,km,nt)
      real stf(ids:ide,jds:jde,nt)
# if defined O_carbon && defined O_carbon_14
      real dc14(ids:ide,jds:jde,km)
# endif
# if defined O_save_npzd
      real npp(ids:ide,jds:jde,kpzd)
      real graz_Z(ids:ide,jds:jde,kpzd) , graz(ids:ide,jds:jde,kpzd)
      real graz_Det(ids:ide,jds:jde,kpzd)
      real morp(ids:ide,jds:jde,kpzd), morz(ids:ide,jds:jde,kpzd)
      real excr(ids:ide,jds:jde,kpzd), morpt(ids:ide,jds:jde,kpzd)
      real remi(ids:ide,jds:jde,km), expo(ids:ide,jds:jde,km)
      real expocal(ids:ide,jds:jde,km)
      real sedrr(ids:ide,jds:jde)
#  if defined O_kk_ballast
      real remi_B(ids:ide,jds:jde,km), expo_B(ids:ide,jds:jde,km)
      real graz_Det_B(ids:ide,jds:jde,kpzd)
#  endif
#  if defined O_npzd_caco3
      real npp_C(ids:ide,jds:jde,kpzd), graz_C(ids:ide,jds:jde,kpzd)
      real morp_C(ids:ide,jds:jde,kpzd), morpt_C(ids:ide,jds:jde,kpzd)
      real dissl(ids:ide,jds:jde,km), expocaco3(ids:ide,jds:jde,km)
      real calatt(ids:ide,jds:jde,kpzd)
      real procal(ids:ide,jds:jde,kpzd)
#  else      
      real procal(ids:ide,jds:jde)
#  endif
#  if defined O_npzd_nitrogen
      real npp_D(ids:ide,jds:jde,kpzd), npp_D_dop(ids:ide,jds:jde,kpzd)
      real npp_dop(ids:ide,jds:jde,kpzd)
      real morpt_D(ids:ide,jds:jde,kpzd), nfix(ids:ide,jds:jde,kpzd)
      real wcdeni(ids:ide,jds:jde,km), bdeni(ids:ide,jds:jde,km)
      real morp_D(ids:ide,jds:jde,kpzd), graz_D(ids:ide,jds:jde,kpzd)
#   if defined O_npzd_caco3
      real npp_C_dop(ids:ide,jds:jde,kpzd)
#   endif      
#  endif
#  if defined O_npzd_iron
      real expofe(ids:ide,jds:jde,kpzd), remife(ids:ide,jds:jde,kpzd)
#    if defined O_npzd_iron_diagnostics
      real feorgads(ids:ide,jds:jde,kpzd)
      real deffe(ids:ide,jds:jde,kpzd), feprime(ids:ide,jds:jde,kpzd)
      real fesed(ids:ide,jds:jde,kpzd), bfe(ids:ide,jds:jde,kpzd)
      real fecol(ids:ide,jds:jde,kpzd)
#     if defined O_npzd_caco3
      real deffe_C(ids:ide,jds:jde,kpzd)
#     endif      
#    endif
#  endif
#  if defined O_npzd_extra_diagnostics
      real avej(ids:ide,jds:jde,kpzd), avej_D(ids:ide,jds:jde,kpzd)
      real gmax(ids:ide,jds:jde,kpzd), no3P(ids:ide,jds:jde,kpzd)
      real po4P(ids:ide,jds:jde,kpzd), po4_D(ids:ide,jds:jde,kpzd)
#  endif
# endif
# if defined O_save_carbon_carbonate_chem
      real sspH(ids:ide,jds:jde), ssCO3(ids:ide,jds:jde)
      real ssOc(ids:ide,jds:jde), ssOa(ids:ide,jds:jde)
      real sspCO2(ids:ide,jds:jde)
# endif

      real :: tm(imt,jmt,km)
      real :: tmpij(imt,jmt)
      real :: tmpijm(imt,jmt)
      real :: tmpijk(imt,jmt,km)
      real :: tmpijkm(imt,jmt,km)
# if defined O_save_npzd
      real :: tmpijl(imt,jmt,km)
      real :: tmpijlm(imt,jmt,km)
# endif

      real :: avgFac
      real c0, c1, c100, c500, c1e3, c1e4, c1e20
      real c1e5, c1e6, p1, p001, p035, C2K, cal2J

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c500 = 500.
      c1e3 = 1.e3
      c1e4 = 1.e4
      c1e5 = 1.e5
      c1e6 = 1.e6
      p1 = 0.1
      p001 = 0.001
      p035 = 0.035
      C2K = 273.15
      cal2J = 2.389e-05
      
      ils=1
      ile=numProfiles
      jls=1
      jle=jmt
      lls=1
      kls=1
      lle=km
      kle=km
      
	  if (doAverage.eq.0) then ! keep accumulating
		tm(:,:,:)=1.0
        iou(1)=1 ! flag for accumulating
	  else  ! time average
		avgFac = 1.0/numAvg
		tm(:,:,:)=avgFac
        iou(1)=2 ! flag for averaging
	  endif

C 2-d diagnostics (this is called 3d data in mom_tavg.F)
      iou(2)=1 ! flag for 2-d diagnostic
      iou(3)=imt
      iou(4)=jmt
      iou(5)=1

      ib(1) = ils-igs+1
      ic(1) = ile-ils+1
      ib(2) = jls-jgs+1
      ic(2) = jle-jls+1
      ib(3) = ntrec
      ic(3) = 1
      ln = ic(1)*ic(2)*ic(3)      
      tmpijm(ils:ile,jls:jle) = tm(ils:ile,jls:jle,1)
      do n=3,nt
        tmpij(ils:ile,jls:jle) = stf(ils:ile,jls:jle,n)
        if (trim(mapt(n)) .eq. 'temp') then
          call putvaramsk('F_heat', iou, ln, ib, ic, tmpij, tmpijm
     &,     cal2J, c0)
        elseif (trim(mapt(n)) .eq. 'salt') then
          call putvaramsk('F_salt', iou, ln, ib, ic, tmpij, tmpijm
     &,     p1, c0)
        elseif (trim(mapt(n)) .eq. 'dic') then
          call putvaramsk('F_dic', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'dic13') then
          call putvaramsk('F_dic13', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'alk') then
          call putvaramsk('F_alk', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'o2') then
          call putvaramsk('F_o2', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'po4') then
          call putvaramsk('F_po4', iou, ln, ib, ic, tmpij, tmpijm
     &,     c1e5, c0)
       elseif (trim(mapt(n)) .eq. 'dop') then
          call putvaramsk('F_dop', iou, ln, ib, ic, tmpij, tmpijm
     &,     c1e5, c0)
        elseif (trim(mapt(n)) .eq. 'no3') then
          call putvaramsk('F_no3', iou, ln, ib, ic, tmpij, tmpijm
     &,     c1e5, c0)
       elseif (trim(mapt(n)) .eq. 'dfe') then
          call putvaramsk('F_dfe', iou, ln, ib, ic, tmpij, tmpijm
     &,     c1e5, c0)
        elseif (trim(mapt(n)) .eq. 'don') then
          call putvaramsk('F_don', iou, ln, ib, ic, tmpij, tmpijm
     &,     c1e5, c0)
        elseif (trim(mapt(n)) .eq. 'din15') then
          call putvaramsk('F_din15', iou, ln, ib, ic, tmpij, tmpijm
     &,     c1e5, c0)
        elseif (trim(mapt(n)) .eq. 'don15') then
          call putvaramsk('F_don15', iou, ln, ib, ic, tmpij, tmpijm
     &,     c1e5, c0)
        elseif (trim(mapt(n)) .eq. 'doc13') then
          call putvaramsk('F_doc13', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'c14') then
          call putvaramsk('F_c14', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'cfc11') then
          call putvaramsk('F_cfc11', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'cfc12') then
          call putvaramsk('F_cfc12', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'phyt') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'phytn15') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'phytc13') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'zoop') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'zoopn15') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'zoopc13') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'diaz') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'diazn15') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'diazc13') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detr') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detrn15') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detrc13') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detrfe') then
!         skip since so surface flux  
        elseif (trim(mapt(n)) .eq. 'cocc') then
!         skip since no surface flux
        elseif (trim(mapt(n)) .eq. 'coccn15') then
!         skip since no surface flux
        elseif (trim(mapt(n)) .eq. 'coccc13') then
!         skip since no surface flux
        elseif (trim(mapt(n)) .eq. 'caco3') then
!         skip since no surface flux
        elseif (trim(mapt(n)) .eq. 'caco3c13') then
!         skip since no surface flux
        elseif (trim(mapt(n)) .eq. 'detr_B') then
!         skip since no surface flux
        else
C         if (n .lt. 1000) write(a3, '(i3)') n
C         if (n .lt. 100) write(a3, '(i2)') n
C         if (n .lt. 10) write(a3, '(i1)') n
C         call putvaramsk('F_'//trim(a3), iou, ln, ib, ic, tmpij
C    &,     tmpijm, c1, c0)
        endif
      enddo

# if defined O_save_carbon_carbonate_chem
      tmpijm(ils:ile,jls:jle) = tm(ils:ile,jls:jle,1)
      where (sspH(ils:ile,jls:jle) .le. 0.) tmpijm(ils:ile,jls:jle) = 0.
      tmpij(ils:ile,jls:jle) = sspH(ils:ile,jls:jle)
      call putvaramsk ('O_phsur', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      tmpij(ils:ile,jls:jle) = ssCO3(ils:ile,jls:jle)
      call putvaramsk ('O_co3sur', iou, ln, ib, ic, tmpij, tmpijm, c1
     &, c0)
      tmpij(ils:ile,jls:jle) = ssOc(ils:ile,jls:jle)
      call putvaramsk ('O_ocalcsur', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      tmpij(ils:ile,jls:jle) = ssOa(ils:ile,jls:jle)
      call putvaramsk ('O_oaragsur', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      tmpij(ils:ile,jls:jle) = sspCO2(ils:ile,jls:jle)
      call putvaramsk ('O_pco2sur', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
# endif
# if defined O_save_npzd
      tmpijm(ils:ile,jls:jle) = tm(ils:ile,jls:jle,1)
#  if !defined O_npzd_caco3
      tmpij(ils:ile,jls:jle) = procal(ils:ile,jls:jle)
      call putvaramsk ('O_caco3pro', iou, ln, ib, ic, tmpij, tmpijm
     &, c100, c0)
#  endif
      tmpij(ils:ile,jls:jle) = sedrr(ils:ile,jls:jle)
      call putvaramsk ('O_sedrr', iou, ln, ib, ic, tmpij, tmpijm
     &, c1e5, c0)
# endif

C 3-d diagnostics (this is called 4d data in mom_tavg.F)
      iou(2)=2 ! flag for 3-d diagnostic
      iou(3)=imt
      iou(4)=jmt
      iou(5)=km

      tmpijkm(ils:ile,jls:jle,kls:kle) = tm(ils:ile,jls:jle,kls:kle)      
# if defined O_carbon && defined O_carbon_14
      tmpijk(ils:ile,jls:jle,kls:kle) = dc14(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_dc14', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c1, c0)
# endif

# if defined O_save_npzd
      tmpijlm(ils:ile,jls:jle,lls:lle) = tm(ils:ile,jls:jle,lls:lle)
      ib(3) = lls-lgs+1
      ic(3) = lle-lls+1
      ln = ic(1)*ic(2)*ic(3)*ic(4)
      tmpijl(ils:ile,jls:jle,lls:lle) = npp(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_phytnpp', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
#  if defined O_npzd_caco3
      tmpijl(ils:ile,jls:jle,lls:lle) = npp_C(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_coccnpp', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = graz_C(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_coccgraz', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = morp_C(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_coccmort', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = morpt_C(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_coccrecy', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = procal(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_caco3pro', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = calatt(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_caco3att', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
#  endif
#  if defined O_kk_ballast
      tmpijl(ils:ile,jls:jle,lls:lle) =
     & graz_Det_B(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_detgraz_B', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
#  endif
#  if defined O_npzd_nitrogen
      tmpijl(ils:ile,jls:jle,lls:lle) = npp_dop(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_phytnpp_dop', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
#   if defined O_npzd_caco3
      tmpijl(ils:ile,jls:jle,lls:lle) =
     &                       npp_C_dop(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_calnpp_dop', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
#   endif      
#  endif
      tmpijl(ils:ile,jls:jle,lls:lle) = graz(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_phytgraz', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = graz_Z(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_zoograz', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = 
     & graz_Det(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_detgraz', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = morp(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_phytmort', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = morpt(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_phytrecy', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = morz(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_zoopmort', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = excr(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_excret', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
#  if defined O_npzd_extra_diagnostics
      tmpijl(ils:ile,jls:jle,lls:lle) = avej(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_avej', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = gmax(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_gmax', iou, ln, ib, ic, tmpijl, tmpijlm 
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = po4P(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_po4P', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
#   if defined O_npzd_nitrogen      
      tmpijl(ils:ile,jls:jle,lls:lle) = po4_D(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_po4_D', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = avej_D(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_avej_D', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = no3P(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_no3P', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
#   endif      
#  endif
#  if defined O_npzd_nitrogen
      tmpijl(ils:ile,jls:jle,lls:lle) = npp_D(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_diaznpp', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle)=npp_D_dop(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_diaznpp_dop', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = graz_D(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_diazgraz', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = morp_D(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_diazmort', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = morpt_D(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_diazrecy', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = nfix(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_nfix', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(3) = kls-kgs+1
      ic(3) = kle-kls+1
      ln = ic(1)*ic(2)*ic(3)*ic(4)
      tmpijkm(ils:ile,jls:jle,kls:kle) = tm(ils:ile,jls:jle,kls:kle)
#   if defined O_npzd_o2      
      tmpijk(ils:ile,jls:jle,kls:kle) = wcdeni(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_wcdeni', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c1e3, c0)
      tmpijk(ils:ile,jls:jle,kls:kle) = bdeni(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_bdeni', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c1e3, c0)
#   endif      
#  endif
#  if defined O_kk_ballast
      tmpijk(ils:ile,jls:jle,kls:kle) = remi_B(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_detrremi_B', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c1e3, c0)
      tmpijk(ils:ile,jls:jle,kls:kle) = expo_B(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_detrexp_B', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c1e3, c0)
#  endif
#  if defined O_npzd_caco3
      tmpijk(ils:ile,jls:jle,kls:kle) = dissl(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_caco3dissl', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c1e3, c0)
      tmpijk(ils:ile,jls:jle,kls:kle) =
     &          expocaco3(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_caco3exp', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c1e3, c0)
#  endif
#  if defined O_npzd_iron
      tmpijl(ils:ile,jls:jle,lls:lle) = expofe(ils:ile,jls:jle
     &, lls:lle)
      call putvaramsk ('O_detrfeexpo', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = remife(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_detrferemi', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0) 
#   if defined O_npzd_iron_diagnostics 
      tmpijl(ils:ile,jls:jle,lls:lle) = feorgads(ils:ile,jls:jle
     &, lls:lle)
      call putvaramsk ('O_feorgads', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = deffe(ils:ile,jls:jle
     &, lls:lle)
      call putvaramsk ('O_deffe', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = feprime(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_feprime', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = fesed(ils:ile,jls:jle
     &, lls:lle)
      call putvaramsk ('O_fesed', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = bfe(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_bfe', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      tmpijl(ils:ile,jls:jle,lls:lle) = fecol(ils:ile,jls:jle
     &, lls:lle)
      call putvaramsk ('O_fecol', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
#    if defined O_npzd_caco3
      tmpijl(ils:ile,jls:jle,lls:lle) = deffe_C(ils:ile,jls:jle
     &, lls:lle)
      call putvaramsk ('O_deffe_C', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1, c0)
#    endif      
#   endif
#  endif
      ib(3) = kls-kgs+1
      ic(3) = kle-kls+1
      ln = ic(1)*ic(2)*ic(3)*ic(4)
      tmpijkm(ils:ile,jls:jle,kls:kle) = tm(ils:ile,jls:jle,kls:kle)
      tmpijk(ils:ile,jls:jle,kls:kle) = remi(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_detrremi', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c1e3, c0)
      tmpijk(ils:ile,jls:jle,kls:kle) = expo(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_detrexp', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c1e3, c0)
#  if !defined O_npzd_caco3      
      tmpijk(ils:ile,jls:jle,kls:kle) = expocal(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_caco3exp', iou, ln, ib, ic, tmpijk, tmpijkm
     &     , c1, c0)
#  endif      
# endif

      END
      
      SUBROUTINE MOBI_DIAGS2d_COPY(id, arr2d, cfilename, debugFlag)

	  use iso_c_binding, only: C_CHAR, c_null_char
      use mobi_diags_mod
      
      implicit none

#include "size.h"

      integer :: id
      real :: arr2d(numProfiles)
	  character (kind=c_char, len=1), intent (inout) :: cfilename(*)
      integer :: debugFlag

#include "pconst.h"
#include "mw.h"
#include "npzd.h"

      logical, save :: firstTime = .true.
      integer :: n

	  n=id+1  ! id starts at 0

      if (firstTime) then
        call F_to_C_string(trim(diag2dFileNames(n)),cfilename)
        if (n.eq.num2dDiags) firstTime=.false.
      endif
        
      arr2d(:) = diags2d(:,n)
      diags2d(:,n)=0.0 ! reset for next averaging block
            
      END

      SUBROUTINE MOBI_DIAGS3d_COPY(id, arr3d, cfilename, debugFlag)

	  use iso_c_binding, only: C_CHAR, c_null_char
      use mobi_diags_mod
      
      implicit none

#include "size.h"

      integer :: id
      real :: arr3d(totNumPoints)
	  character (kind=c_char, len=1), intent (inout) :: cfilename(*)      
      integer :: debugFlag

#include "pconst.h"
#include "mw.h"
#include "npzd.h"

      logical, save :: firstTime = .true.
      integer n

	  n=id+1  ! id starts at 0

      if (firstTime) then
        call F_to_C_string(trim(diag3dFileNames(n)),cfilename)
        if (n.eq.num3dDiags) firstTime=.false.
      endif
      
      arr3d(:) = diags3d(:,n)
      diags3d(:,n)=0.0 ! reset for next averaging block
            
      END

      SUBROUTINE MOBI_DIAGS_FINALIZE(debugFlag)

      use mobi_diags_mod
      
      implicit none

#include "size.h"

      integer :: debugFlag

#include "pconst.h"
#include "mw.h"
#include "npzd.h"
#include "mobi_diags.h"

      deallocate(diags2d)
      deallocate(diags3d)
      
      END

	  subroutine F_to_C_string(fstr,cstr)

!     Copy a fortran string fstr to a C string cstr

	  use iso_c_binding, only: C_CHAR, c_null_char
	  implicit none

	  character (kind=c_char, len=1), intent (inout) :: cstr(*)
	  character(len=*), intent(in) :: fstr

	  integer :: i, nchar

	  nchar=len(fstr)
	  do i=1,nchar
		cstr(i)=fstr(i:i)
	  enddo
	  cstr(nchar+1)=c_null_char

	  return

	  end subroutine F_to_C_string

C=============================================================
C Find a free unit for I/O
C=============================================================
      SUBROUTINE FINDUNIT( iounit )
C OUT:
C     iounit   integer - unit number
C
C MDSFINDUNIT returns a valid, unused unit number for f77 I/O
C The routine stops the program is an error occurs in the process
C of searching the I/O channels.
C
C Created: 03/20/99 adcroft@mit.edu

      implicit none

C Arguments
      integer iounit
C Local
      integer ii
      logical op
      integer ios
C     ------------------------------------------------------------------

C Sweep through a valid range of unit numbers
      iounit=-1
      do ii=9,99
        if (iounit.eq.-1) then
          inquire(unit=ii,iostat=ios,opened=op)
          if (ios.ne.0) then
            stop 'ABNORMAL END: S/R MDSFINDUNIT'
          endif
          if (.NOT. op) then
            iounit=ii
          endif
        endif
      enddo

C Was there an available unit number
      if (iounit.eq.-1) then
        stop ' MDSFINDUNIT: could not find an available unit number!'
      endif

C     ------------------------------------------------------------------
      return
      END
      