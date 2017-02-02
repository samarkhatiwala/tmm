#ifdef O_TMM
      SUBROUTINE UVOK_DIAGS_INI(lNumProfiles, lTotNumPoints, 
     &  lNum2dDiags, lNum3dDiags, debugFlag)

      use uvok_diags_mod 
       
      implicit none

#include "size.h"

	  integer :: lNumProfiles, lTotNumPoints, debugFlag
	  integer :: lNum2dDiags, lNum3dDiags

#include "pconst.h"
#include "mw.h"
#include "npzd.h"
#include "timeavgs.h"
#include "switch.h"
#include "uvok_diags.h"

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

!     Set some flags so UVOK doesn't accumulate diagnostics
      euler2 = .false.
      timavgperts = .false.
      tsiperts = .false.

	  CALL FINDUNIT(diagsLogFileUnit)
	  OPEN(diagsLogFileUnit, file='available_diagnostics.txt',
     &     status='unknown', form='formatted')

C     Code based on S/R mom_tavg_def in mom_tavg.F
C     Define 2-d diagnostics    
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
          id_F_dic=it(1)
        elseif (trim(mapt(n)) .eq. 'dic13') then
          call defvar ('F_dic13', iou,3, it, -c100, c100, ' ', 'F'
     &,   'surface downward carbon 13 flux', ' ', 'mol m-2 s-1')
          id_F_dic13=it(1)
        elseif (trim(mapt(n)) .eq. 'alk') then
          call defvar ('F_alk', iou,3, it, -c100, c100, ' ', 'F'
     &,   'surface downward alkalinity flux', ' ', 'mol m-2 s-1')
          id_F_alk=it(1)     
        elseif (trim(mapt(n)) .eq. 'o2') then
          call defvar ('F_o2', iou,3, it, -c100, c100, ' ', 'F'
     &,   'surface downward oxygen flux', ' ', 'mol m-2 s-1')
          id_F_o2=it(1)     
        elseif (trim(mapt(n)) .eq. 'po4') then
          call defvar ('F_po4', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward phosphate flux', ' ', 'mol m-2 s-1')
          id_F_po4=it(1)     
        elseif (trim(mapt(n)) .eq. 'dop') then
          call defvar ('F_dop', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward DOP flux', ' ', 'mol m-2 s-1')
          id_F_dop=it(1)     
        elseif (trim(mapt(n)) .eq. 'no3') then
          call defvar ('F_no3', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward nitrate flux', ' ', 'mol m-2 s-1')
          id_F_no3=it(1)     
        elseif (trim(mapt(n)) .eq. 'dfe') then
          call defvar ('F_dfe', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward iron flux', ' ', 'mol m-2 s-1')
          id_F_dfe=it(1)     
        elseif (trim(mapt(n)) .eq. 'don') then
          call defvar ('F_don', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward DON flux', ' ', 'mol m-2 s-1')
          id_F_don=it(1)     
        elseif (trim(mapt(n)) .eq. 'din15') then
          call defvar ('F_din15', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward nitrate 15 flux', ' ', 'mol m-2 s-1')
          id_F_din15=it(1)     
        elseif (trim(mapt(n)) .eq. 'don15') then
          call defvar ('F_don15', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward DON 15 flux', ' ', 'mol m-2 s-1')
          id_F_don15=it(1)     
       elseif (trim(mapt(n)) .eq. 'doc13') then
         call defvar ('F_doc13', iou,3, it, -c100, c100, ' ', 'F'
     &,   'surface downward DOC 13 flux', ' ', 'mol m-2 s-1')
         id_F_doc13=it(1)     
        elseif (trim(mapt(n)) .eq. 'c14') then
          call defvar ('F_c14', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward carbon 14 flux', ' ', 'mol m-2 s-1')
          id_F_c14=it(1)     
        elseif (trim(mapt(n)) .eq. 'cfc11') then
          call defvar ('F_cfc11', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward CFC11 flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'cfc12') then
          call defvar ('F_cfc12', iou,3, it, -c1e4, c1e4, ' ', 'F'
     &,   'surface downward CFC12 flux', ' ', 'mol m-2 s-1')
        elseif (trim(mapt(n)) .eq. 'phyt') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'zoop') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'diaz') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detr') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'phytn15') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'zoopn15') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detrn15') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'diazn15') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'phytc13') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'zoopc13') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detrc13') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'diazc13') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detrfe') then
!         skip since so surface flux
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
      id_O_phsur=it(1)
      call defvar ('O_co3sur', iou, 3, it, c0, c100, ' ', 'F'
     &, 'sea surface CO3', ' ', 'mol m-3')
      id_O_co3sur=it(1)     
      call defvar ('O_ocalcsur', iou, 3, it, c0, c100, ' ', 'F'
     &, 'sea surface Omega calcite', ' ', '1')
      id_O_ocalcsur=it(1)     
      call defvar ('O_oaragsur', iou, 3, it, c0, c100, ' ', 'F'
     &, 'sea surface Omega aragonite', ' ', '1')
      id_O_oaragsur=it(1)     
      call defvar ('O_pco2sur', iou, 3, it, c0, c100, ' ', 'F'
     &, 'sea surface pCO2', ' ', '1')
      id_O_pco2sur=it(1)     
# endif
# if defined O_save_npzd
      call defvar ('O_caco3pro', iou, 3, it, c0, c100, ' ', 'F'
     &, 'calcite production', ' ', 'mol C m-2 s-1')
      id_O_caco3pro=it(1)     
      call defvar ('O_sedrr', iou, 3, it, -c1, c100, ' ', 'F'
     &, 'detrital sediment rain rate', ' ', 'mol N m-2 s-1')
      id_O_sedrr=it(1)     
# endif

	  num2dDiags = it(1) ! total number of 2-d diagnostics
	  lNum2dDiags = num2dDiags

	  ALLOCATE( diags2d(numProfiles,num2dDiags) )	
	  diags2d(:,:) = 0.0

C     Define 3-d diagnostics    
      iou=2 ! flag for 3-d diagnostic
      it(1)=0  ! index for diagnostic (incremented within defvar)
#if defined O_carbon_14
	  call defvar ('O_dc14', iou, 4, it, -c1e3, c1e3, ' ', 'F'
     &,     'delta carbon 14', ' ', 'permil')
      id_O_dc14=it(1)
#endif      
# if defined O_save_npzd
      call defvar ('O_phytnpp', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'ocean net primary production rate', ' ', 'mol N m-3 s-1')
      id_O_phytnpp=it(1)
      call defvar ('O_phytnpp_dop', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'ocean net primary production rate from DOP', ' '
     &, 'mol N m-3 s-1')
      id_O_phytnpp_dop=it(1)     
      call defvar ('O_phytgraz', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'grazing rate on P', ' ', 'mol N m-3 s-1')
      id_O_phytgraz=it(1)     
      call defvar ('O_zoograz', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'grazing rate on Z', ' ', 'mol N m-3 s-1')
      id_O_zoograz=it(1)     
      call defvar ('O_detgraz', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'grazing rate on detritus', ' ', 'mol N m-3 s-1')
      id_O_detgraz=it(1)     
      call defvar ('O_phytmort', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'phytoplankton mortality rate', ' ', 'mol N m-3 s-1')
      id_O_phytmort=it(1)     
      call defvar ('O_phytrecy', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'microbial fast recycling rate', ' ', 'mol N m-3 s-1')
      id_O_phytrecy=it(1)     
      call defvar ('O_zoopmort', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'zooplankton mortality rate', ' ', 'mol N m-3 s-1')
      id_O_zoopmort=it(1)     
      call defvar ('O_excret', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'excretion rate', ' ', 'mol N m-3 s-1')
      id_O_excret=it(1)     
# if defined O_npzd_extra_diagnostics
      call defvar ('O_avej', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'light-dependant phyt. growth rate', ' ', 'd-1')
      id_O_avej=it(1)     
      call defvar ('O_avej_D', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'light-dependant Diaz growth rate', ' ', 'd-1')
      id_O_avej_D=it(1)     
      call defvar ('O_gmax', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'temp-dependant zoo growth rate', ' ', 'd-1')
      id_O_gmax=it(1)     
      call defvar ('O_no3P', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'no3-dependant phyt. growth rate', ' ', 'd-1')
      id_O_no3P=it(1)     
      call defvar ('O_po4P', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'po4-dependant phyt. growth rate', ' ', 'd-1')
      id_O_po4P=it(1)     
      call defvar ('O_po4_D', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'po4-dependant Diaz growth rate', ' ', 'd-1')
      id_O_po4_D=it(1)     
# endif
#  if defined O_npzd_nitrogen
      call defvar ('O_diaznpp', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'net primary production diazotrophs', ' ', 'mol N m-3 s-1')
      id_O_diaznpp=it(1)     
      call defvar ('O_diaznpp_dop', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'net primary production diazotrophs from DOP', ' '
     &, 'mol N m-3 s-1')
      id_O_diaznpp_dop=it(1)     
      call defvar ('O_diazgraz', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'grazing rate on diazotrophs', ' ', 'mol N m-3 s-1')
      id_O_diazgraz=it(1)     
      call defvar ('O_diazmort', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'mortality diazotrophs', ' ', 'mol N m-3 s-1')
      id_O_diazmort=it(1)     
      call defvar ('O_diazrecy', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'fast-recycling diazotrophs', ' ', 'mol N m-3 s-1')
      id_O_diazrecy=it(1)     
      call defvar ('O_nfix', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'nitrogen fixation', ' ', 'mol N m-3 s-1')
      id_O_nfix=it(1)     
      call defvar ('O_wcdeni', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'water column denitrification', ' ', 'mol N m-3 s-1')
      id_O_wcdeni=it(1)     
      call defvar ('O_bdeni', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'benthic denitrification', ' ', 'mol N m-3 s-1')
      id_O_bdeni=it(1)     
#  endif
#  if defined O_npzd_iron
      call defvar ('O_detrfeexpo', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'particulate Fe export rate', ' ', 'mol Fe m-3 s-1')
      id_O_detrfeexpo=it(1)     
      call defvar ('O_detrferemi', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'Fe remineralization rate', ' ', '-')
      id_O_detrferemi=it(1)     
#  if defined O_npzd_iron_diagnostics
      call defvar ('O_feorgads', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'iron organic scavenging rate', ' ', 'mol Fe m-3 s-1')
      id_O_feorgads=it(1)     
      call defvar ('O_deffe', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'iron deficiency', ' ', '-')
      id_O_deffe=it(1)     
      call defvar ('O_feprime', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'free iron concentration', ' ', 'mol Fe m-3')
      id_O_feprime=it(1)     
      call defvar ('O_fesed', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'iron sediment release (without iron from OMZs)', ' '
     &, 'mol Fe m-3 s-1')
      id_O_fesed=it(1)     
      call defvar ('O_bfe', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'iron sediment release rate (with iron from OMZs)', ' '
     &, 'mol Fe m-3 s-1')
      id_O_bfe=it(1)     
      call defvar ('O_fecol', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'iron inorganic scavenging rate', ' ', 'mol Fe m-3 s-1')
      id_O_fecol=it(1)     
# endif
# endif
# if defined O_npzd_chl
      call defvar ('O_chl', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'diagnostic chlorophyll concentration', ' '
     &, 'g Chl m-3')
      id_O_chl=it(1)
#   if defined O_npzd_nitrogen
      call defvar ('O_chl_diaz', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'chlorophyll concentration diazotrophs', ' '
     &, 'g Chl m-3')
      id_O_chl_diaz=it(1)
#   endif
#  endif
      call defvar ('O_detrremi', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'remineralisation rate', ' ', 'mol N m-3 s-1')
      id_O_detrremi=it(1)     
      call defvar ('O_detrexp', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'detrital export rate', ' ', 'mol N m-3 s-1')
      id_O_detrexp=it(1)     
      call defvar ('O_caco3exp', iou, 4, it, -c1, c100, ' ', 'F'
     &, 'calcite export rate', ' ', 'mol C m-3 s-1')
      id_O_caco3exp=it(1)     
# endif

	  num3dDiags = it(1) ! total number of 3-d diagnostics
      lNum3dDiags = num3dDiags
      
	  ALLOCATE( diags3d(totNumPoints, num3dDiags) )
	  diags3d(:,:) = 0.0

      CLOSE(diagsLogFileUnit)

      CALL UVOK_DIAGS_RESET(debugFlag)

	  END

      SUBROUTINE UVOK_DIAGS_START(debugFlag)

      use uvok_diags_mod 
       
      implicit none

#include "size.h"

	  INTEGER :: debugFlag

#include "pconst.h"
#include "mw.h"
#include "npzd.h"
#include "timeavgs.h"
#include "switch.h"
#include "uvok_diags.h"

!     Set some flags so UVOK accumulates diagnostics
      euler2 = .false.
      timavgperts = .true.
      tsiperts = .true.

      CALL UVOK_DIAGS_RESET(debugFlag)
    
	  END

      SUBROUTINE UVOK_DIAGS_RESET(debugFlag)

      use uvok_diags_mod 
       
      implicit none

#include "size.h"

	  INTEGER :: debugFlag

#include "pconst.h"
#include "mw.h"
#include "npzd.h"
#include "timeavgs.h"
#include "switch.h"
#include "uvok_diags.h"

!     Reset arrays to zero
      
#if defined O_carbon_14
      ta_dc14(:,:,:)=0.0
#endif
#if defined O_save_npzd
      ta_rnpp(:,:,:)=0.0
      ta_rgraz(:,:,:)=0.0
      ta_rmorp(:,:,:)=0.0
      ta_rmorpt(:,:,:)=0.0
      ta_rmorz(:,:,:)=0.0
      ta_rexcr(:,:,:)=0.0
      ta_rremi(:,:,:)=0.0
      ta_rexpo(:,:,:)=0.0
      ta_rnpp_D(:,:,:)=0.0
      ta_rgraz_D(:,:,:)=0.0
!      ta_rmorpt_D(:,:,:)=0.0
      ta_rnfix(:,:,:)=0.0
      ta_rdeni(:,:,:)=0.0
      ta_rexpocal(:,:,:)=0.0
      ta_rprocal(:,:)=0.0
      ta_rgraz_Z(:,:,:)=0.0
      ta_rgraz_Det(:,:,:)=0.0
      ta_ravej(:,:,:)=0.0
      ta_ravej_D(:,:,:)=0.0
      ta_rgmax(:,:,:)=0.0
      ta_rno3P(:,:,:)=0.0
      ta_rpo4P(:,:,:)=0.0
      ta_rpo4_D(:,:,:)=0.0
!      ta_rsedrr(:,:)=0.0
!      ta_rbdeni(:,:,:)=0.0
      ta_rmorp_D(:,:,:)=0.0
!      ta_rnpp_dop(:,:,:)=0.0
!      ta_rnpp_D_dop(:,:,:)=0.0
#endif
#   if defined O_chl 
       ta_rchl(:,:,:) = 0.0 
#    if defined O_npzd_nitrogen
       ta_rchl_D(:,:,:) = 0.0
#    endif
#   endif
#if defined O_save_carbon_carbonate_chem
      ta_sspH(:,:)=0.0
      ta_ssCO3(:,:)=0.0
      ta_ssOc(:,:)=0.0
      ta_ssOa(:,:)=0.0
      ta_sspCO2(:,:)=0.0
#endif
#if defined O_npzd_iron
      ta_rfeorgads(:,:,:)=0.0
      ta_rfecol(:,:,:)=0.0
      ta_rdeffe(:,:,:)=0.0
      ta_rremife(:,:,:)=0.0
      ta_rexpofe(:,:,:)=0.0
      ta_rfeprime(:,:,:)=0.0
      ta_rfesed(:,:,:)=0.0
      ta_rbfe(:,:,:)=0.0
#endif          
	  END

      SUBROUTINE UVOK_DIAGS_STOP(debugFlag)

      use uvok_diags_mod 
       
      implicit none

#include "size.h"

	  INTEGER :: debugFlag

#include "pconst.h"
#include "mw.h"
#include "npzd.h"
#include "timeavgs.h"
#include "switch.h"
#include "uvok_diags.h"

!     Set some flags so UVOK doesn't accumulate diagnostics
      euler2 = .false.
      timavgperts = .false.
      tsiperts = .false.

	  END
        
      SUBROUTINE UVOK_DIAGS_ACCUMULATE(ipro, kl, nzloc, numAvg, 
     &                                 doAverage, debugFlag)

      use uvok_diags_mod
      
      implicit none

#include "size.h"

      integer :: ipro, kl, nzloc, numAvg, doAverage, debugFlag

#include "pconst.h"
#include "mw.h"
#include "npzd.h"
#include "diaga.h"
#include "timeavgs.h"
#include "uvok_diags.h"

      integer :: idd, idn, n
      integer :: ip, kstart, kend
      integer :: ils, ile, jls, jle, lls, kls
      parameter (ils=1, ile=1, jls=1, jle=1, lls=1, kls=1)
      real :: tmpij(ils:ile,jls:jle)
      real :: tmpijm(ils:ile,jls:jle)=1.0
      real :: tmpijl(1,1:km,1)
      real :: tmpijlm(1,1:km,1)=1.0
      real :: tmpijk(1,1:km,1)
      real :: tmpijkm(1,1:km,1)=1.0
      integer :: lle, kle
      integer ic(10), ib(10), iou, ln
      real :: avgFac
      real c100, c500, c1e3, c1e4, c1e20
      real c1e5, c1e6, p1, p001, p035, C2K, cal2J

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
      
C     tmpijm(:,:)=1.0
C     tmpijlm(:,:,:)=1.0
      
      ip = ipro + 1 ! ipro starts at 0
      kstart = kl +  1 ! kl starts at 0
      kend = kstart + nzloc - 1
      lle=nzloc
      kle=nzloc
      
	  if (doAverage.eq.0) then ! keep accumulating
C	tmpijm(:,:)=1.0
C	tmpijlm(:,:,:)=1.0
C	tmpijkm(:,:,:)=1.0		
        ic(1)=1 ! flag for accumulating
	  else  ! time average
		avgFac = 1.0/numAvg
		tmpijm(:,:)=avgFac
		tmpijlm(:,:,:)=avgFac
		tmpijkm(:,:,:)=avgFac
        ic(1)=2 ! flag for averaging
	  endif
     	   
      iou=1 ! flag for 2-d diagnostic
      ib(2)=ip ! profile number
      ib(3)=ip
	  ln=1 ! size of 2-d diagnostic
      do n=3,nt
        tmpij(ils:ile,jls:jle) = stf(ils:ile,jls:jle,n)
        if (trim(mapt(n)) .eq. 'temp') then
          call putvaramsk('F_heat', iou, ln, ib, ic, tmpij, tmpijm
     &,     cal2J, c0)
        elseif (trim(mapt(n)) .eq. 'salt') then
          call putvaramsk('F_salt', iou, ln, ib, ic, tmpij, tmpijm
     &,     p1, c0)
        elseif (trim(mapt(n)) .eq. 'dic') then
          ib(1)=id_F_dic
          call putvaramsk('F_dic', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'dic13') then
          ib(1)=id_F_dic13
          call putvaramsk('F_dic13', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'alk') then
          ib(1)=id_F_alk        
          call putvaramsk('F_alk', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'o2') then
          ib(1)=id_F_o2        
          call putvaramsk('F_o2', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'po4') then
          ib(1)=id_F_po4        
          call putvaramsk('F_po4', iou, ln, ib, ic, tmpij, tmpijm
     &,     c1e5, c0)
        elseif (trim(mapt(n)) .eq. 'dop') then
          ib(1)=id_F_dop        
          call putvaramsk('F_dop', iou, ln, ib, ic, tmpij, tmpijm
     &,     c1e5, c0)
        elseif (trim(mapt(n)) .eq. 'no3') then
          ib(1)=id_F_no3        
          call putvaramsk('F_no3', iou, ln, ib, ic, tmpij, tmpijm
     &,     c1e5, c0)
        elseif (trim(mapt(n)) .eq. 'dfe') then
          ib(1)=id_F_dfe        
          call putvaramsk('F_dfe', iou, ln, ib, ic, tmpij, tmpijm
     &,     c1e5, c0)
        elseif (trim(mapt(n)) .eq. 'don') then
          ib(1)=id_F_don        
          call putvaramsk('F_don', iou, ln, ib, ic, tmpij, tmpijm
     &,     c1e5, c0)
        elseif (trim(mapt(n)) .eq. 'din15') then
          ib(1)=id_F_din15        
          call putvaramsk('F_din15', iou, ln, ib, ic, tmpij, tmpijm
     &,     c1e5, c0)
        elseif (trim(mapt(n)) .eq. 'don15') then
          ib(1)=id_F_don15        
          call putvaramsk('F_don15', iou, ln, ib, ic, tmpij, tmpijm
     &,     c1e5, c0)
        elseif (trim(mapt(n)) .eq. 'doc13') then
          ib(1)=id_F_doc13        
          call putvaramsk('F_doc13', iou, ln, ib, ic, tmpij, tmpijm
     &,     c100, c0)
        elseif (trim(mapt(n)) .eq. 'c14') then
          ib(1)=id_F_c14        
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
        elseif (trim(mapt(n)) .eq. 'zoop') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'diaz') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detr') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'diaz') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'phytn15') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'zoopn15') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detrn15') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'diazn15') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'phytc13') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'zoopc13') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detrc13') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'diazc13') then
!         skip since so surface flux
        elseif (trim(mapt(n)) .eq. 'detrfe') then
!         skip since so surface flux  
!
        else
C         if (n .lt. 1000) write(a3, '(i3)') n
C         if (n .lt. 100) write(a3, '(i2)') n
C         if (n .lt. 10) write(a3, '(i1)') n
C         call putvaramsk('F_'//trim(a3), iou, ln, ib, ic, tmpij
C    &,     tmpijm, c1, c0)
        endif
      enddo

# if defined O_save_carbon_carbonate_chem
C     tmpijm(ils:ile,jls:jle) = tm(ils:ile,jls:jle,1)
C     where (sspH(ils:ile,jls:jle) .le. 0.) tmpijm(ils:ile,jls:jle) = 0.
      ib(1)=id_O_phsur
      tmpij(ils:ile,jls:jle) = sspH(ils:ile,jls:jle)
      call putvaramsk ('O_phsur', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      ib(1)=id_O_co3sur     
      tmpij(ils:ile,jls:jle) = ssCO3(ils:ile,jls:jle)
      call putvaramsk ('O_co3sur', iou, ln, ib, ic, tmpij, tmpijm, c1
     &, c0)
      ib(1)=id_O_ocalcsur     
      tmpij(ils:ile,jls:jle) = ssOc(ils:ile,jls:jle)
      call putvaramsk ('O_ocalcsur', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      ib(1)=id_O_oaragsur     
      tmpij(ils:ile,jls:jle) = ssOa(ils:ile,jls:jle)
      call putvaramsk ('O_oaragsur', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      ib(1)=id_O_pco2sur     
      tmpij(ils:ile,jls:jle) = sspCO2(ils:ile,jls:jle)
      call putvaramsk ('O_pco2sur', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
# endif
# if defined O_save_npzd
C     tmpijm(ils:ile,jls:jle) = tm(ils:ile,jls:jle,1)
      ib(1)=id_O_caco3pro
      tmpij(ils:ile,jls:jle) = ta_rprocal(ils:ile,jls:jle)
      call putvaramsk ('O_caco3pro', iou, ln, ib, ic, tmpij, tmpijm
     &, c100, c0)
!      ib(1)=id_O_sedrr  !   
!      tmpij(ils:ile,jls:jle) = 0.0 ! sedrr(ils:ile,jls:jle)  this variable is not set anywhere
!      call putvaramsk ('O_sedrr', iou, ln, ib, ic, tmpij, tmpijm
!     &, c1e5, c0)
# endif

C 3-d diagnostics
      iou=2 ! flag for 3-d diagnostic
      ib(2)=kstart
      ib(3)=kend
      ln=nzloc ! size of 3-d diagnostic
# if defined O_carbon && defined O_carbon_14
      ib(1)=id_O_dc14
      tmpijk(1,kls:kle,1) = ta_dc14(1,kls:kle,1) ! dc14(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_dc14', iou, ln, ib, ic, tmpijk, tmpijkm
     &, c1, c0)
# endif      
# if defined O_save_npzd
      ib(1)=id_O_phytnpp
      tmpijl(1,lls:lle,1) = ta_rnpp(1,lls:lle,1) ! npp(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_phytnpp', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_phytnpp_dop
!      tmpijl(1,lls:lle,1) = ta_rnpp_dop(1,lls:lle,1) ! npp_dop(ils:ile,jls:jle,lls:lle)
!      call putvaramsk ('O_phytnpp_dop', iou, ln, ib, ic, tmpijl, tmpijlm
!     &, c1e3, c0)
      ib(1)=id_O_phytgraz
      tmpijl(1,lls:lle,1) = ta_rgraz(1,lls:lle,1) ! graz(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_phytgraz', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_zoograz
      tmpijl(1,lls:lle,1) = ta_rgraz_Z(1,lls:lle,1) ! graz_Z(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_zoograz', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_detgraz
      tmpijl(1,lls:lle,1) = 
     & ta_rgraz_Det(1,lls:lle,1) ! graz_Det(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_detgraz', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_phytmort
      tmpijl(1,lls:lle,1) = ta_rmorp(1,lls:lle,1) ! morp(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_phytmort', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_phytrecy
      tmpijl(1,lls:lle,1) = ta_rmorpt(1,lls:lle,1) ! morpt(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_phytrecy', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_zoopmort
      tmpijl(1,lls:lle,1) = ta_rmorz(1,lls:lle,1) ! morz(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_zoopmort', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_excret
      tmpijl(1,lls:lle,1) = ta_rexcr(1,lls:lle,1) ! excr(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_excret', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)      
# if defined O_npzd_extra_diagnostics
      ib(1)=id_O_avej
      tmpijl(1,lls:lle,1) = ta_ravej(1,lls:lle,1) ! avej(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_avej', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_avej_D
      tmpijl(1,lls:lle,1) = ta_ravej_D(1,lls:lle,1) ! avej_D(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_avej_D', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_gmax
      tmpijl(1,lls:lle,1) = ta_rgmax(1,lls:lle,1) ! gmax(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_gmax', iou, ln, ib, ic, tmpijl, tmpijlm 
     &, c1e3, c0)
      ib(1)=id_O_no3P
      tmpijl(1,lls:lle,1) = ta_rno3P(1,lls:lle,1) ! no3P(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_no3P', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_po4P
      tmpijl(1,lls:lle,1) = ta_rpo4P(1,lls:lle,1) ! po4P(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_po4P', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_po4_D
      tmpijl(1,lls:lle,1) = ta_rpo4_D(1,lls:lle,1) ! po4_D(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_po4_D', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
# endif
#  if defined O_npzd_nitrogen
      ib(1)=id_O_diaznpp
      tmpijl(1,lls:lle,1) = ta_rnpp_D(1,lls:lle,1) ! npp_D(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_diaznpp', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
!      ib(1)=id_O_diaznpp_dop
!      tmpijl(1,lls:lle,1)=ta_rnpp_D_dop(1,lls:lle,1) ! npp_D_dop(ils:ile,jls:jle,lls:lle)
!      call putvaramsk ('O_diaznpp_dop', iou, ln, ib, ic, tmpijl, tmpijlm
!     &, c1e3, c0)
      ib(1)=id_O_diazgraz
      tmpijl(1,lls:lle,1) = ta_rgraz_D(1,lls:lle,1) ! graz_D(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_diazgraz', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_diazmort
      tmpijl(1,lls:lle,1) = ta_rmorp_D(1,lls:lle,1) ! morp_D(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_diazmort', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_diazrecy
!      tmpijl(1,lls:lle,1) = ta_rmorpt_D(1,lls:lle,1) ! morpt_D(ils:ile,jls:jle,lls:lle)
!      call putvaramsk ('O_diazrecy', iou, ln, ib, ic, tmpijl, tmpijlm
!     &, c1e3, c0)
      ib(1)=id_O_nfix
      tmpijl(1,lls:lle,1) = ta_rnfix(1,lls:lle,1) ! nfix(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_nfix', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
C     ib(3) = kls-kgs+1
C     ic(3) = kle-kls+1
C     ln = ic(1)*ic(2)*ic(3)*ic(4)
C     tmpijkm(ils:ile,jls:jle,kls:kle) = tm(ils:ile,jls:jle,kls:kle)
      ib(1)=id_O_wcdeni
      tmpijl(1,lls:lle,1) = ta_rdeni(1,lls:lle,1) ! wcdeni(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_wcdeni', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
!      ib(1)=id_O_bdeni
!      tmpijk(1,kls:kle,1) = ta_rbdeni(1,kls:kle,1) ! bdeni(ils:ile,jls:jle,kls:kle)
!      call putvaramsk ('O_bdeni', iou, ln, ib, ic, tmpijk, tmpijkm
!     &, c1e3, c0)
#  endif
#  if defined O_npzd_iron
      ib(1)=id_O_detrfeexpo
      tmpijl(1,lls:lle,1) = ta_rexpofe(1,lls:lle,1) ! expofe(ils:ile,jls:jle
!     &, lls:lle)
      call putvaramsk ('O_detrfeexpo', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_detrferemi
      tmpijl(1,lls:lle,1) = ta_rremife(1,lls:lle,1) ! remife(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_detrferemi', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0) 
#  if defined O_npzd_iron_diagnostics 
      ib(1)=id_O_feorgads
      tmpijl(1,lls:lle,1) = ta_rfeorgads(1,lls:lle,1) ! feorgads(ils:ile,jls:jle
!     &, lls:lle)
      call putvaramsk ('O_feorgads', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_deffe
      tmpijl(1,lls:lle,1) = ta_rdeffe(1,lls:lle,1) ! deffe(ils:ile,jls:jle
!     &, lls:lle)
      call putvaramsk ('O_deffe', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1, c0)
      ib(1)=id_O_feprime
      tmpijl(1,lls:lle,1) = ta_rfeprime(1,lls:lle,1) ! feprime(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_feprime', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_fesed
      tmpijl(1,lls:lle,1) = ta_rfesed(1,lls:lle,1) ! fesed(ils:ile,jls:jle
!     &, lls:lle)
      call putvaramsk ('O_fesed', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_bfe
      tmpijl(1,lls:lle,1) = ta_rbfe(1,lls:lle,1) ! bfe(ils:ile,jls:jle,lls:lle)
      call putvaramsk ('O_bfe', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_fecol
      tmpijl(1,lls:lle,1) = ta_rfecol(1,lls:lle,1) ! fecol(ils:ile,jls:jle
!     &, lls:lle)
      call putvaramsk ('O_fecol', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
#  endif
#  endif
#  if defined O_npzd_chl
      ib(1)=id_O_chl
      tmpijl(1,lls:lle,1) = ta_rchl(1,lls:lle,1)
      call putvaramsk ('O_chl', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1, c0)
#   if defined O_npzd_nitrogen
      ib(1)=id_O_chl_D
      tmpijl(1,lls:lle,1) = ta_rchl_D(1,lls:lle,1)
      call putvaramsk ('O_chl_diaz', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1, c0)
#   endif
#  endif
C     ib(3) = kls-kgs+1
C     ic(3) = kle-kls+1
C     ln = ic(1)*ic(2)*ic(3)*ic(4)
C     tmpijkm(ils:ile,jls:jle,kls:kle) = tm(ils:ile,jls:jle,kls:kle)
      ib(1)=id_O_detrremi
      tmpijl(1,lls:lle,1) = ta_rremi(1,lls:lle,1) ! remi(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_detrremi', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_detrexp
      tmpijl(1,lls:lle,1) = ta_rexpo(1,lls:lle,1) ! expo(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_detrexp', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1e3, c0)
      ib(1)=id_O_caco3exp
      tmpijl(1,lls:lle,1) = ta_rexpocal(1,lls:lle,1) ! expocal(ils:ile,jls:jle,kls:kle)
      call putvaramsk ('O_caco3exp', iou, ln, ib, ic, tmpijl, tmpijlm
     &, c1, c0)
# endif

      CALL UVOK_DIAGS_RESET(debugFlag)

      END
      
      SUBROUTINE UVOK_DIAGS2d_COPY(id, arr2d, cfilename, debugFlag)

	  use iso_c_binding, only: C_CHAR, c_null_char
      use uvok_diags_mod
      
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

      SUBROUTINE UVOK_DIAGS3d_COPY(id, arr3d, cfilename, debugFlag)

	  use iso_c_binding, only: C_CHAR, c_null_char
      use uvok_diags_mod
      
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

      SUBROUTINE UVOK_DIAGS_FINALIZE(debugFlag)

      use uvok_diags_mod
      
      implicit none

#include "size.h"

      integer :: debugFlag

#include "pconst.h"
#include "mw.h"
#include "npzd.h"
#include "uvok_diags.h"

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
#endif
