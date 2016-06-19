#ifdef O_TMM
#include "MOBI_TMM_OPTIONS.h"
#endif

      subroutine gosbc (is, ie, js, je)

#if defined O_mom || defined O_embm
!=======================================================================
!     calculate the average fluxes for next ocean time step
!=======================================================================

      implicit none

      integer ie, is, je, js, i, j, nc
      integer iem1, isp1, jem1, jsp1, k

      real f1, f1a, f1l, fh, fs, fwcflx, fwaflx, time
      real area, tarea, tsflx, rsocn, tmp, fg, fgs, sulphfac

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "calendar.h"
#include "csbc.h"
#include "coord.h"
#include "grdvar.h"
#include "tmngr.h"
#include "switch.h"
# if defined O_embm
#include "cembm.h"
#include "atm.h"
# endif
# if defined O_mom
#include "mw.h"
# endif
# if defined O_ice
#  if defined O_ice_cpts
#include "cpts.h"
#  endif
#include "ice.h"
# endif
# if defined O_mtlm
#include "mtlm.h"
# endif
#include "levind.h"
# if defined O_fwa
#include "fwa.h"
#include "cregin.h"
# endif
# if defined O_sed && !defined O_sed_uncoupled
#include "sed.h"
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
#include "insolation.h"
      real tot, twt, cz(1), zero(1)
# endif

# if defined O_embm_read_sflx || defined O_embm_write_sflx
      integer ntrec
      save ntrec
      data ntrec /0/
# endif

#ifndef O_TMM
      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1
#else
      isp1 = is
      iem1 = ie
      jsp1 = js
      jem1 = je
#endif
#ifndef O_TMM
# if defined O_mom && defined O_embm
      f1 = 1./atatm
      fh = 2.389e-8/atatm
      fs = -socn/atatm

!-----------------------------------------------------------------------
!     calculate average net fluxes. convert heat flux to cal/cm**2/s
!     from kW/m**2 and fresh water flux (cm/s) to an apparent salt
!     flux (g/cm**2/s) using global ocean average salinity, socn
!-----------------------------------------------------------------------

      do j=2,jmtm1
        do i=2,imtm1
          if (tmsk(i,j) .ge. 0.5) then
#  if defined O_plume
#   if !defined O_plume_brine
            subflux(i,j,1) = flux(i,j,isat)
            if (subflux(i,j,1) .gt. 0.) subflux(i,j,1) = 0.
            flux(i,j,isat) = flux(i,j,isat) - subflux(i,j,1)
            subflux(i,j,1) = fh*subflux(i,j,1)

            subflux(i,j,2) = flux(i,j,ishum)
#   endif
            if (subflux(i,j,2) .gt. 0.) subflux(i,j,2) = 0.
            flux(i,j,ishum) = flux(i,j,ishum) - subflux(i,j,2)
            subflux(i,j,2) = fs*subflux(i,j,2)
#  endif
#  if defined O_convect_brine
            cba0(i,j) = 0.
            do nc=0,ncat
              cba(i,j,nc) = f1*cba(i,j,nc)
              if (cba(i,j,nc) .gt. 0.) then
                cbf(i,j,nc) = fs*cbf(i,j,nc)/cba(i,j,nc)
                cba0(i,j) = cba0(i,j) + cba(i,j,nc)
              else
                cbf(i,j,nc) = 0.
                cba(i,j,nc) = 0.
              endif
            enddo
            if (cba0(i,j) .gt. 1.) then
              if (cba0(i,j) .gt. 1.000001) then
                print*, "==> Warning: ice area > 1: ", cba0(i,j)
                endif
              cba0(i,j) = 1.
            endif
            cba0(i,j) = 1. - cba0(i,j)
#  endif
            sbc(i,j,ihflx) = sbc(i,j,ihflx) + fh*flux(i,j,isat)
!           add virtual fluxes of salinity
            sbc(i,j,isflx) = sbc(i,j,isflx) + fs*flux(i,j,ishum)

          else
#  if defined O_plume
            subflux(i,j,1) = 0.
            subflux(i,j,2) = 0.
#  endif
#  if defined O_convect_brine
            cbf(i,j,:) = 0.
            cba(i,j,:) = 0.
            cba0(i,j) = 1.
#  endif
            sbc(i,j,ihflx) = 0.
            sbc(i,j,isflx) = 0.
          endif
          if (umsk(i,j) .ge. 0.5) then
#  if defined O_ice_evp || defined O_embm_awind
            sbc(i,j,itaux) = f1*flux(i,j,nat+1)
            sbc(i,j,itauy) = f1*flux(i,j,nat+2)
#  endif
          else
            sbc(i,j,itaux) = 0.
            sbc(i,j,itauy) = 0.
          endif
        enddo
      enddo

      call setbcx (sbc(1,1,ihflx), imt, jmt)
      call setbcx (sbc(1,1,isflx), imt, jmt)
      call setbcx (sbc(1,1,itaux), imt, jmt)
      call setbcx (sbc(1,1,itauy), imt, jmt)
#  if defined O_convect_brine
      do nc=0,ncat
        call setbcx (cbf(1,1,nc), imt, jmt)
        call setbcx (cba(1,1,nc), imt, jmt)
      enddo
      call setbcx (cba0, imt, jmt)
#  endif
# endif

!-----------------------------------------------------------------------
!     update boundary conditions from the land model
!     do this now instead of in gasbc so fields can be written out
!-----------------------------------------------------------------------

      f1l = 0.
      f1a = 0.
# if defined O_embm
      if (atatm .ne. 0.) f1a = 1.0/atatm
# endif
# if defined O_mtlm
      if (atlnd .ne. 0.) f1l = 1.0/atlnd
      do j=2,jmtm1
        do i=2,imtm1
          if (land_map(i,j) .ne. 0) then
            sbc(i,j,iro) = sbc(i,j,iro)*f1l
            sbc(i,j,isca) = sbc(i,j,isca)*f1l
            sbc(i,j,ievap) = sbc(i,j,ievap)*f1l
            sbc(i,j,ilwr) = sbc(i,j,ilwr)*f1l
            sbc(i,j,isens) = sbc(i,j,isens)*f1l
#  if defined O_carbon
            sbc(i,j,inpp) =  sbc(i,j,inpp)*f1l
            sbc(i,j,isr) =  sbc(i,j,isr)*f1l
#  endif
#  if defined O_mtlm_carbon_13
            sbc(i,j,inpp13) =  sbc(i,j,inpp13)*f1l
            sbc(i,j,isr13) =  sbc(i,j,isr13)*f1l
            sbc(i,j,iburn13) =  sbc(i,j,iburn13)*f1l
#  endif
          else
            sbc(i,j,iro) = sbc(i,j,iro)*f1a
            sbc(i,j,ievap) = 0.
            sbc(i,j,ilwr) = 0.
            sbc(i,j,isens) = 0.
#  if defined O_mtlm && defined O_carbon
            sbc(i,j,inpp) = 0.
            sbc(i,j,isr) = 0.
#  endif
#  if defined O_mtlm_carbon_13
            sbc(i,j,inpp13) = 0.
            sbc(i,j,isr13) = 0.
            sbc(i,j,iburn13) = 0.
#  endif
#  if defined O_mtlm_carbon_14
            sbc(i,j,inpp14) = 0.
            sbc(i,j,isr14) = 0.
            sbc(i,j,iburn14) = 0.
#  endif
          endif
        enddo
      enddo
      call setbcx (sbc(1,1,isca), imt, jmt)
      call setbcx (sbc(1,1,ievap), imt, jmt)
      call setbcx (sbc(1,1,ilwr), imt, jmt)
      call setbcx (sbc(1,1,isens), imt, jmt)
#  if defined O_carbon
      call setbcx (sbc(1,1,inpp), imt, jmt)
      call setbcx (sbc(1,1,isr), imt, jmt)
#  endif
#  if defined O_mtlm_carbon_13
      call setbcx (sbc(1,1,inpp13), imt, jmt)
      call setbcx (sbc(1,1,isr13), imt, jmt)
      call setbcx (sbc(1,1,iburn13), imt, jmt)
#  endif
#  if defined O_mtlm_carbon_14
      call setbcx (sbc(1,1,inpp14), imt, jmt)
      call setbcx (sbc(1,1,isr14), imt, jmt)
      call setbcx (sbc(1,1,iburn14), imt, jmt)
#  endif
# else
      sbc(:,:,iro) = sbc(:,:,iro)*f1a
# endif
      call setbcx (sbc(1,1,iro), imt, jmt)

# if defined O_embm
!-----------------------------------------------------------------------
!     zero diagnostic for river discharge and call river model
!-----------------------------------------------------------------------
      disch(:,:) = 0.
      call rivmodel

#  if defined O_sed && !defined O_sed_uncoupled
!-----------------------------------------------------------------------
!     apply CaCO3 weathering flux proportional to discharge
!-----------------------------------------------------------------------
#   if defined O_sed_constrain_weath
      if (weathflx .lt. 0.) weathflx = 0.
#   endif
#   if defined O_save_carbon_totals
      dicwflx =  weathflx
#   endif
      tmp = 0.
      do j=2,jmtm1
        do i=2,imtm1
          if (kmt(i,j) .ne. 0.) then
            tmp = tmp + disch(i,j)*dxt(i)*dyt(j)*cst(j)
          endif
        enddo
      enddo
      tmp = weathflx/tmp
      do j=2,jmtm1
#   if defined O_global_sums || defined O_save_carbon_totals
        fgs = dyt(j)*cst(j)*segtim/secday
#   endif
        do i=2,imtm1
          if (kmt(i,j) .ne. 0.) then
#   if defined O_carbon
            sbc(i,j,idicflx) = sbc(i,j,idicflx) + disch(i,j)*tmp
#    if defined O_carbon_13
              sbc(i,j,idic13flx) = sbc(i,j,idic13flx) + dic13ocn*tmp
#    endif
#   if defined O_npzd_alk
            sbc(i,j,ialkflx) = sbc(i,j,ialkflx) + disch(i,j)*2.*tmp
#   endif
#    if defined O_save_carbon_totals
            carblith = carblith - disch(i,j)*dxt(i)*tmp*fgs
     &              - sbc(i,j,ibdicfx)*dxt(i)*fgs*dzt(kmt(i,j))
#    endif
#    if defined O_global_sums
            dtoic = dtoic - disch(i,j)*dxt(i)*tmp*fgs
#    endif
#   endif
          endif
        enddo
      enddo

#  endif
# endif
# if defined O_sealev_data_transient && defined O_sealev_salinity
!-----------------------------------------------------------------------
!     add in flux from sea level change
!-----------------------------------------------------------------------
      do j=2,jmtm1
        do i=2,imtm1
          if (kmt(i,j) .gt. 0)
     &      sbc(i,j,isflx) = sbc(i,j,isflx) + fs*dsealev
        enddo
      enddo

# endif
#endif ! O_TMM
# if defined O_mom && defined O_embm
#  if defined O_fwa
!-----------------------------------------------------------------------
!     add additional fresh water flux anomaly
!-----------------------------------------------------------------------

      time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
      if (time .ge. fwayri .and. time .le. fwayrf) then

        if (areafwa .gt. 0) then
!         fwaflxi is in Sv (1e6 m3 s-1) and fwarate is in Sv/yr
          fwaflx = fwaflxi + (time - fwayri)*fwarate
!         convert to flux in g salt cm-2 s-1
          fwaflx = -socn*fwaflx*1.e12/areafwa
#   if defined O_fwa_precip
          call areaavg (precip, fwawt, tmp)
          if (tmp .gt. 0) fwaflx = fwaflx/tmp
#   endif
          do j=2,jmtm1
            do i=2,imtm1
              sbc(i,j,isflx) = sbc(i,j,isflx) + fwaflx*fwawt(i,j)
#   if defined O_fwa_precip
     &                         *precip(i,j)
#   endif
            enddo
          enddo
        endif

        if (compensate .and. areafwc .gt. 0) then
!         fwaflxi is in Sv (1e6 m3 s-1) and fwarate is in Sv/yr
          fwcflx = fwaflxi + (time - fwayri)*fwarate
!         convert to opposite flux in g salt cm-2 s-1
          fwcflx = socn*fwcflx*1.e12/areafwc
#   if defined O_fwa_compevap
          call areaavg (evap, fwcwt, tmp)
          if (tmp .gt. 0) fwcflx = fwcflx/tmp
#   endif
          do j=2,jmtm1
            do i=2,imtm1
              sbc(i,j,isflx) = sbc(i,j,isflx) + fwcflx*fwcwt(i,j)
#   if defined O_fwa_compevap
     &                         *evap(i,j)
#   endif
            enddo
          enddo
        endif

      endif

#  endif
#  if defined O_carbon || defined O_npzd_alk || defined O_npzd_o2 || defined O_npzd || defined O_cfcs_data || defined O_cfcs_data_transient
!-----------------------------------------------------------------------
!     add normalized virtual fluxes to other tracers
!-----------------------------------------------------------------------
      tarea = 0.
      tsflx = 0.
      rsocn = 1./socn
      do j=2,jmtm1
        do i=2,imtm1
          if (tmsk(i,j) .ge. 0.5) then
            area = dxt(i)*dyt(j)*cst(j)
            tarea = tarea + area
            tsflx = tsflx + sbc(i,j,isflx)*area
          endif
        enddo
      enddo
      tsflx = tsflx/tarea
      do j=2,jmtm1
        do i=2,imtm1
          if (tmsk(i,j) .ge. 0.5) then
            tmp = (sbc(i,j,isflx) - tsflx)*rsocn
            vflux(i,j) = tmp
#   if defined O_carbon
            sbc(i,j,idicflx) = sbc(i,j,idicflx) + gaost(idic)*tmp
#    if defined O_carbon_13
            sbc(i,j,idic13flx) = sbc(i,j,idic13flx) + gaost(idic13)*tmp
#    endif
#    if defined O_carbon_14
            sbc(i,j,ic14flx) = sbc(i,j,ic14flx) + gaost(ic14)*tmp
#    endif
#   endif
#   if defined O_npzd_alk
            sbc(i,j,ialkflx) = sbc(i,j,ialkflx) + gaost(ialk)*tmp
#   endif
#   if defined O_npzd_o2
            sbc(i,j,io2flx) = sbc(i,j,io2flx) + gaost(io2)*tmp
#   endif
#   if defined O_npzd
            sbc(i,j,ipo4flx) = sbc(i,j,ipo4flx) + gaost(ipo4)*tmp
            sbc(i,j,idopflx) = sbc(i,j,idopflx) + gaost(idop)*tmp
#    if !defined O_npzd_no_vflux
            sbc(i,j,iphytflx) = sbc(i,j,iphytflx) + gaost(iphyt)*tmp
            sbc(i,j,izoopflx) = sbc(i,j,izoopflx) + gaost(izoop)*tmp
            sbc(i,j,idetrflx) = sbc(i,j,idetrflx) + gaost(idetr)*tmp
#    endif
#    if defined O_npzd_nitrogen
            sbc(i,j,ino3flx) = sbc(i,j,ino3flx) + gaost(ino3)*tmp
            sbc(i,j,idonflx) = sbc(i,j,idonflx) + gaost(idon)*tmp
#     if !defined O_npzd_no_vflux
            sbc(i,j,idiazflx) = sbc(i,j,idiazflx) + gaost(idiaz)*tmp
#     endif
#     if defined O_npzd_nitrogen_15
            sbc(i,j,idin15flx) = sbc(i,j,idin15flx) + gaost(idin15)*tmp
            sbc(i,j,idon15flx) = sbc(i,j,idon15flx) + gaost(idon15)*tmp
#      if !defined O_npzd_no_vflux 
            sbc(i,j,iphytn15flx) = sbc(i,j,iphytn15flx) 
     &                          + gaost(iphytn15)*tmp
            sbc(i,j,izoopn15flx) = sbc(i,j,izoopn15flx) 
     &                          + gaost(izoopn15)*tmp
            sbc(i,j,idetrn15flx) = sbc(i,j,idetrn15flx) 
     &                          + gaost(idetrn15)*tmp
            sbc(i,j,idiazn15flx) = sbc(i,j,idiazn15flx) 
     &                          + gaost(idiazn15)*tmp
#      endif
#     endif
#    endif
#   endif
#   if defined O_carbon_13
            sbc(i,j,idoc13flx) = sbc(i,j,idoc13flx) + gaost(idoc13)*tmp
#   endif
#   if defined O_cfcs_data || defined O_cfcs_data_transient
            sbc(i,j,icfc11flx) = sbc(i,j,icfc11flx) + gaost(icfc11)*tmp
            sbc(i,j,icfc12flx) = sbc(i,j,icfc12flx) + gaost(icfc12)*tmp
#   endif
          endif
        enddo
      enddo
#   if defined O_carbon
      call setbcx (sbc(1,1,idicflx), imt, jmt)
#    if defined O_carbon_13
      call setbcx (sbc(1,1,idic13flx), imt, jmt)
#    endif
#    if defined O_carbon_14
      call setbcx (sbc(1,1,ic14flx), imt, jmt)
#    endif
#   endif
#   if defined O_npzd_alk
      call setbcx (sbc(1,1,ialkflx), imt, jmt)
#   endif
#   if defined O_npzd_o2
      call setbcx (sbc(1,1,io2flx), imt, jmt)
#   endif
#   if defined O_npzd
      call setbcx (sbc(1,1,ipo4flx), imt, jmt)
      call setbcx (sbc(1,1,idopflx), imt, jmt)
#    if !defined O_npzd_no_vflux
      call setbcx (sbc(1,1,iphytflx), imt, jmt)
      call setbcx (sbc(1,1,izoopflx), imt, jmt)
      call setbcx (sbc(1,1,idetrflx), imt, jmt)
#    endif
#    if defined O_npzd_nitrogen
      call setbcx (sbc(1,1,ino3flx), imt, jmt)
      call setbcx (sbc(1,1,idonflx), imt, jmt)
#     if !defined O_npzd_no_vflux
      call setbcx (sbc(1,1,idiazflx), imt, jmt)
#     endif
#     if defined O_npzd_nitrogen_15
      call setbcx (sbc(1,1,idin15flx), imt, jmt)
      call setbcx (sbc(1,1,idon15flx), imt, jmt)
#      if !defined O_npzd_no_vflux
      call setbcx (sbc(1,1,iphytn15flx), imt, jmt)
      call setbcx (sbc(1,1,izoopn15flx), imt, jmt)
      call setbcx (sbc(1,1,idetrn15flx), imt, jmt)
      call setbcx (sbc(1,1,idiazn15flx), imt, jmt)
#      endif
#     endif
#     if defined O_carbon_13
      call setbcx (sbc(1,1,idoc13flx), imt, jmt)
#      if !defined O_npzd_no_vflux
      call setbcx (sbc(1,1,iphytc13flx), imt, jmt)
      call setbcx (sbc(1,1,izoopc13flx), imt, jmt)
      call setbcx (sbc(1,1,idetrc13flx), imt, jmt)
#       if defined O_npzd_nitrogen
      call setbcx (sbc(1,1,idiazc13flx), imt, jmt)
#       endif
#      endif
#     endif
#    endif
#   endif
#   if defined O_cfcs_data || defined O_cfcs_data_transient
      call setbcx (sbc(1,1,icfc11flx), imt, jmt)
      call setbcx (sbc(1,1,icfc12flx), imt, jmt)
#   endif
#  endif
# endif
# if defined O_embm_read_sflx
      call read_var ("F_salt.nc", "F_salt", sbc(1,1,isflx), ntrec)
# elif defined O_embm_write_sflx
      call write_var ("F_salt.nc", "F_salt", sbc(1,1,isflx), ntrec)
# endif
#endif
# if defined O_embm && defined O_sulphate_data || defined O_sulphate_data_transient

!-----------------------------------------------------------------------
!     set anthropogenic sulphate forcing (increase in surface albedo)
!     for the next atmospheric step from updated surface albedo
!-----------------------------------------------------------------------
      call sulphdata

      zero(1) = 0.
      cz(1) = 0.
!     full direct and indirect sulphate forcing
      sulphfac = 0.29
#  if defined O_sulphate_data_direct && !defined O_sulphate_data_indirect
!     reduce sulphate forcing to just the direct effect (0.4/1.1)
      sulphfac = sulphfac*.4/1.1
#  endif
#  if defined O_sulphate_data_indirect && !defined O_sulphate_data_direct
!     reduce sulphate forcing to just the indirect effect (0.7/1.1)
      sulphfac = sulphfac*.7/1.1
#  endif
      do j=jsp1,jem1
        do i=isp1,iem1
          if (sulph(i,j,2) .gt. 0. .and. solins(i,j) .gt. 0.) then
            tot = 0.
            twt = 0.
!           calculate sunlight weighted daily average zenith angle
            do k=1,24
                call zenith (1, (k-1)*3600., 3600., daylen, tlat(i,j)
     &,                      zero(1), sindec, cz(1))
              tot = tot + cz(1)*cz(1)
              twt = twt + cz(1)
            enddo
            cz(1) = 0.
            if (twt .gt. 0) cz(1) = tot/twt
!           convert from optical depth to albedo
            if (cz(1) .gt. 0.05) then
!             reonstruct surface coalbedo
              tmp = solins(i,j) - solins(i,j)*(1. - sbc(i,j,iaca))
     &            - solins(i,j)*sbc(i,j,iaca)*scatter
              if (tmp .gt. 0.) then
                tmp = dnswr(i,j)/tmp
                sulph(i,j,2) = sulphfac*sulph(i,j,2)*tmp**2/cz(1)
              else
                sulph(i,j,2) = 0.
              endif
            else
              sulph(i,j,2) = 0.
            endif
          endif
        enddo
      enddo
# endif

      return
      end

#if defined O_embm_read_sflx
      subroutine read_var (fname, vname, var, ntrec)

      implicit none

      character (*) :: fname, vname

      integer ib(10), ic(10), iou, ntrec

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "tmngr.h"

      real var(imt,jmt), tmpij(imtm2,jmtm2)

      ntrec = ntrec + 1
      call openfile (trim(fname), iou)
      ib(:) = 1
      ib(3) = ntrec
      ic(:) = 1
      ic(1) = imtm2
      ic(2) = jmtm2
      call getvara (trim(vname), iou, imtm2*jmtm2, ib, ic, tmpij
     &, 1., 0.)
      var(2:imtm1,2:jmtm1) = tmpij(1:imtm2:1:jmtm2)
      embmbc (var)

      return
      end
#endif

#if defined O_embm_write_sflx
      subroutine write_var (fname, vname, var, ntrec)

      implicit none

      character (*) :: fname, vname

      integer ib(10), ic(10), it(10), id_time, id_xt, id_yt, iou, ntrec
      integer nyear, nmonth, nday, nhour, nmin, nsec

      real tmp

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "coord.h"
#include "tmngr.h"

      real var(imt,jmt), tmpij(imtm2,jmtm2), tmpi(imtm2), tmpj(jmtm2)

      ntrec = ntrec + 1
      call openfile (trim(fname), iou)
      if (ntrec .eq. 1) then
        call redef (iou)

        call defdim ('time', iou, 0, id_time)
        call defdim ('longitude', iou, imtm2, id_xt)
        call defdim ('latitude', iou, jmtm2, id_yt)

        it(:) = id_time
        call defvar ('time', iou, 1, it, c0, c0, 'T', 'D'
     &,   'time', 'time', 'years since 0-1-1')
        call putatttext (iou, 'time', 'calendar', calendar)

        it(1) = id_xt
        call defvar ('longitude', iou, 1, it, 0., 0., 'X', 'D'
     &,   'longitude', 'longitude', 'degrees_east')
        it(1) = id_yt
        call defvar ('latitude', iou, 1, it, 0., 0., 'Y', 'D'
     &,   'latitude', 'latitude', 'degrees_north')

        it(:) = id_xt
        it(2) = id_yt
        it(3) = id_time
        call defvar (trim(vname), iou, 3, it, -1.e20, 1.e20, ' ', 'D'
     &,   ' ',' ', ' ')

        call enddef (iou)

!-----------------------------------------------------------------------
!       write 1d data (t)
!-----------------------------------------------------------------------
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)

        ib(:) = 1
        ic(:) = imtm2
        tmpi(1:imtm2) = xt(2:imtm1)
        call putvara ('longitude', iou, imtm2, ib, ic, tmpi, c1, c0)
        ic(:) = jmtm2
        tmpj(1:jmtm2) = xt(2:jmtm1)
        call putvara ('latitude', iou, jmtm2, ib, ic, tmpj, c1, c0)

      endif
      call putvars ('time', iou, ntrec, relyr, c1, c0)
      ib(:) = 1
      ib(3) = ntrec
      ic(:) = 1
      ic(1) = imtm2
      ic(2) = jmtm2
      tmpij(1:imtm2:1:jmtm2) = var(2:imtm1,2:jmtm1)
      call putvara (trim(vname), iou, imtm2*jmtm2, ib, ic, tmpij
     &, 1., 0.)

      return
      end
#endif
