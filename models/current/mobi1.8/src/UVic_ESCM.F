#ifdef O_TMM
#include "MOBI_TMM_OPTIONS.h"
#else
      program UVic_ESCM

!=======================================================================

!     UNIVERSITY OF VICTORIA EARTH SYSTEM CLIMATE MODEL (UVic ESCM)

!     The UVic ESCM is a climate model developed by researchers in the
!     Climate Modelling Group, within the School of Earth and Ocean
!     Sciences, located at the University of Victoria, Victoria,
!     British Columbia, Canada.

!     Many people have contributed to the development of this code.
!     This is a collective effort and although individual contributions
!     are appreciated, individual authorship is not indicated.

!     Some of the people who have contributed to the development of the
!     model are:
!       D. Archer, C. Avis, A. Berger, R. Betts, A. Biastoch, C. Bitz,
!       J. Brauch, C. Brennan, B. Bryan, J. Burton, S. Carto,
!       M. Cottet-Puinel, M. Cox, P. Cox, G. Danabasoglu, K. Dixon,
!       J. Dukowicz, M. Eby, R. Essery, T. Ewen, A. Fanning, J. Fyke,
!       R. Gerdes, A. Gnanadesikan, C. Goldberg, D. Goldberg,
!       D. Gregory, J. Gregory, S. Griffies, R. Hanson, R. Hetherington,
!       H. Hickey, M. Holland, G. Holloway, T. Huck, T. Hughes,
!       E. Hunke, W. Hurlin, W. Ingram, R. Key, E. Kluzek, C. Koeberle,
!       K. Kvale, C. Lawson, J. Lewis, M. Loutre, R. Malone,
!       D. Matthews, K. Meissner, A. Montenegro, A. Mouchet, T. Murdock,
!       P. Myers, A. Oschlies, R. Pacanowski, M. Pahlow, P. Poussart
!       S. Rahmstorf, R. Redler,, D. Robitaille, A. Rosati, M. Roth
!       C. Sabine, O. Saenko, A. Schmittner, B. Semtner, W. Sijp
!       T. Silva, H. Simmons, A. Skvortsov, R. Smith, P. Spence
!       D. Stone, R. Tonkonojenkov, C. Tricot, S. Valcke, A. Weaver
!       E. Wiebe, J. Willebrandt, M. Yoshimori, K. Zickfeld

!     Please direct problems or questions to the code contact
!     person at: http://climate.uvic.ca/model

!     Requirements:

!     Standard fortran 90 is used

!     Disclaimer:

!     The UVic Earth System Climate Model is a climate modelling
!     research tool developed at the University of Victoria. Others may
!     use it freely but we assume no responsibility for problems or
!     incorrect use. It is left to the user to ensure that a particular
!     configuration is working correctly.
!=======================================================================

      implicit none

      integer i, j, numots, numats, numseg, n, loop

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "coord.h"
#include "csbc.h"
#include "iounit.h"
#include "emode.h"
#include "levind.h"
#include "scalar.h"
#include "switch.h"
#include "tmngr.h"
#include "cembm.h"
#if defined O_embm
#include "atm.h"
#endif
#if defined O_mom
#include "mw.h"
#endif
      print*, '== UNIVERSITY OF VICTORIA EARTH SYSTEM CLIMATE MODEL =='
      print*, '                                            '

!-----------------------------------------------------------------------
!     initialize i/o units
!-----------------------------------------------------------------------

      call ioinit

!-----------------------------------------------------------------------
!     setup file renaming
!-----------------------------------------------------------------------

      call file_names

!-----------------------------------------------------------------------
!     Initialize S.B.C. indices
!-----------------------------------------------------------------------

      call sbc_init

!-----------------------------------------------------------------------
!     Initialize tracers
!-----------------------------------------------------------------------

      call tracer_init

!-----------------------------------------------------------------------
!     read namelist variables
!-----------------------------------------------------------------------

      call read_namelist

!-----------------------------------------------------------------------
!     read grid
!-----------------------------------------------------------------------

      call grids

!-----------------------------------------------------------------------
!     read topography
!-----------------------------------------------------------------------

      call topog (kmt, kmu, map, xt, yt, zt, xu, yu, zw, imt, jmt, km
     &,           sg_bathy)

      call isleperim (kmt, map, iperm, jperm, iofs, nippts, nisle, imt
     &,               jmt, km, mnisle, maxipp, xu, yu, zw)

!-----------------------------------------------------------------------
!     common setup
!-----------------------------------------------------------------------

      call setcom (1, imt, 1, jmt)

#if defined O_mom
!-----------------------------------------------------------------------
!     ocean setup
!-----------------------------------------------------------------------

      call setmom (1, imt, 1, jmt)

# if defined O_sed
!-----------------------------------------------------------------------
!     sediment setup (requires mom)
!-----------------------------------------------------------------------

      call setsed (1, imt, 1, jmt)

# endif
#endif
#if defined O_embm
!-----------------------------------------------------------------------
!     atmosphere setup
!-----------------------------------------------------------------------
      call setembm (1, imt, 1, jmt)

# if defined O_ism
!-----------------------------------------------------------------------
!     ice sheet setup (requires embm)
!-----------------------------------------------------------------------

      call setism (1, imt, 1, jmt)

# endif
# if defined O_mtlm
!-----------------------------------------------------------------------
!     land setup (requires embm)
!-----------------------------------------------------------------------
      call setmtlm (1, imt, 1, jmt)

# endif
#endif
!-----------------------------------------------------------------------
!     initialize data reading routine
!-----------------------------------------------------------------------

      call setdata (1, imt, 1, jmt)

!-----------------------------------------------------------------------
!     compute the number of ocean time steps "numots" for this run and
!     the number of ocean time steps per ocean segment "ntspos".
!     compute the number of atmos time steps "numats" for this run and
!     the number of atmos time steps per atmos segment "ntspas".
!     divide the integration time "days" into "numseg" segments.
!     each will be length "segtim" days. Surface boundary conditions
!     are supplied every "segtim" days.
!-----------------------------------------------------------------------

      numots = nint(rundays/(dtocn*secday))
      ntspos = nint(segtim/(dtocn*secday))
      numats = nint(rundays/(dtatm*secday))
      ntspas = nint(segtim/(dtatm*secday))
      numseg = numots/ntspos
#if defined O_mtlm
# if defined O_mtlm_segday
      if (segtim .gt. 1.) then
        ntspls = nint(c1/(dtlnd*secday))
      else
        ntspls = nint(segtim/(dtlnd*secday))
      endif
# else
      ntspls = nint(segtim/(dtlnd*secday))
# endif
#endif

#if !defined O_embm && defined O_mom
!     load the tracers (SST & SSS) for each row "j". (zero on land)
!     load from the MW if fully opened otherwise load from disk

      if (wide_open_mw) then
        do j=1,jmt
          do i=1,imt
# if !defined O_replacst
            if (isst .ne. 0) sbc(i,j,isst) = t(i,1,j,itemp,taup1)
            if (isss .ne. 0) sbc(i,j,isss) = t(i,1,j,isalt,taup1)
# endif
# if defined O_carbon
            if (issdic .ne. 0) sbc(i,j,issdic) = t(i,1,j,idic,taup1)
#  if defined O_carbon_13
            if (issdic13 .ne. 0) sbc(i,j,issdic13) = 
     &           t(i,1,j,idic13,taup1)
#  endif
#  if defined O_carbon_14
            if (issc14 .ne. 0) sbc(i,j,issc14) = t(i,1,j,ic14,taup1)
#  endif
# endif
# if defined O_npzd_alk
            if (issalk .ne. 0) sbc(i,j,issalk) = t(i,1,j,ialk,taup1)
# endif
# if defined O_npzd_o2
            if (isso2 .ne. 0) sbc(i,j,isso2) = t(i,1,j,io2,taup1)
# endif
# if defined O_npzd
            if (isspo4 .ne. 0) sbc(i,j,isspo4) = t(i,1,j,ipo4,taup1)
            if (issdop .ne. 0) sbc(i,j,issdon) = t(i,1,j,idop,taup1)
            if (issphyt .ne. 0) sbc(i,j,issphyt) = t(i,1,j,iphyt,taup1)
            if (isszoop .ne. 0) sbc(i,j,isszoop) = t(i,1,j,izoop,taup1)
            if (issdetr .ne. 0) sbc(i,j,issdetr) = t(i,1,j,idetr,taup1)
#  if defined O_npzd_nitrogen
            if (issno3 .ne. 0) sbc(i,j,issno3) = t(i,1,j,ino3,taup1)
            if (issdon .ne. 0) sbc(i,j,issdon) = t(i,1,j,idon,taup1)
            if (issdiaz .ne. 0) sbc(i,j,issdiaz) = t(i,1,j,idiaz,taup1)
#   if defined O_npzd_nitrogen_15
            if (issdin15 .ne. 0) sbc(i,j,issdin15) 
     &           = t(i,1,j,idin15,taup1)
            if (issdon15 .ne. 0) sbc(i,j,issdon15) 
     &           = t(i,1,j,idon15,taup1)
            if (issphytn15 .ne. 0) sbc(i,j,issphytn15) 
     &           = t(i,1,j,iphytn15,taup1)
            if (isszoopn15 .ne. 0) sbc(i,j,isszoopn15) 
     &           = t(i,1,j,izoopn15,taup1)
            if (issdetrn15 .ne. 0) sbc(i,j,issdetrn15) 
     &           = t(i,1,j,idetrn15,taup1)
            if (issdiazn15 .ne. 0) sbc(i,j,issdiazn15) 
     &           = t(i,1,j,idiazn15,taup1)
#   endif
#  endif
#  if defined O_carbon_13
            if (issdoc13 .ne. 0) sbc(i,j,issdoc13) = 
     &           t(i,1,j,idoc13,taup1)
            if (issphytc13 .ne. 0) sbc(i,j,issphytc13) = 
     &           t(i,1,j,iphytc13,taup1)
            if (isszoopc13 .ne. 0) sbc(i,j,isszoopc13) = 
     &           t(i,1,j,izoopc13,taup1)
            if (issdetrc13 .ne. 0) sbc(i,j,issdetrc13) = 
     &           t(i,1,j,idetrc13,taup1)
#   if defined O_npzd_nitrogen
            if (issdiazc13 .ne. 0) sbc(i,j,issdiazc13) = 
     &           t(i,1,j,idiazc13,taup1)
#   endif
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
            if (isscfc11 .ne. 0) sbc(i,j,isscfc11) =
     &        t(i,1,j,icfc11,taup1)
            if (isscfc12 .ne. 0) sbc(i,j,isscfc12) =
     &        t(i,1,j,icfc12,taup1)
# endif
          enddo
        enddo
      else
        do j=1,jmt
# if !defined O_replacst
          if (isst .ne. 0) call getst (j, sbc(1,1,isst), itemp)
          if (isss .ne. 0) call getst (j, sbc(1,1,isss), isalt)
# endif
# if defined O_carbon
          if (issdic .ne. 0) call getst (j, sbc(1,1,issdic), idic)
#  if defined O_carbon_13
          if (issdic13 .ne. 0) call getst (j, sbc(1,1,issdic13), idic13)
#  endif
#  if defined O_carbon_14
          if (issc14 .ne. 0) call getst (j, sbc(1,1,issc14), ic14)
#  endif
# endif
# if defined O_npzd_alk
          if (issalk .ne. 0) call getst (j, sbc(1,1,issalk), ialk)
# endif
# if defined O_npzd_o2
          if (isso2 .ne. 0) call getst (j, sbc(1,1,isso2), io2)
# endif
# if defined O_npzd
            if (isspo4 .ne. 0) call getst (j, sbc(1,1,isspo4), ipo4)
            if (issdop .ne. 0) call getst (j, sbc(1,1,issdop), idop)
            if (issphyt .ne. 0) call getst (j, sbc(1,1,issphyt), iphyt)
            if (isszoop .ne. 0) call getst (j, sbc(1,1,isszoop), izoop)
            if (issdetr .ne. 0) call getst (j, sbc(1,1,issdetr), idetr)
#  if defined O_npzd_nitrogen
            if (issno3 .ne. 0) call getst (j, sbc(1,1,issno3), ino3)
            if (issdon .ne. 0) call getst (j, sbc(1,1,issdon), idon)
            if (issdiaz .ne. 0) call getst (j, sbc(1,1,issdiaz), idiaz)
#   if defined O_npzd_nitrogen_15
            if (issdin15 .ne. 0) 
     &           call getst (j, sbc(1,1,issdin15), idin15)
            if (issdon15 .ne. 0) 
     &           call getst (j, sbc(1,1,issdon15), idon15)
            if (issphytn15 .ne. 0) 
     &           call getst (j, sbc(1,1,issphytn15), iphytn15)
            if (isszoopn15 .ne. 0) 
     &           call getst (j, sbc(1,1,isszoopn15), izoopn15)
            if (issdetrn15 .ne. 0) 
     &           call getst (j, sbc(1,1,issdetrn15), idetrn15)
            if (issdiazn15 .ne. 0) 
     &           call getst (j, sbc(1,1,issdiazn15), idiazn15)
#   endif
#  endif
#  if defined O_carbon_13
            if (issdoc13 .ne. 0) call getst (j, sbc(1,1,issdoc13)
     &,       idoc13)
            if (issphytc13 .ne. 0) call getst (j, sbc(1,1,issphytc13)
     &,       iphytc13)
            if (isszoopc13 .ne. 0) call getst (j, sbc(1,1,isszoopc13)
     &,       izoopc13)
            if (issdetrc13 .ne. 0) call getst (j, sbc(1,1,issdetrc13)
     &,       idetrc13)
#   if defined O_npzd_nitrogen
            if (issdiazc13 .ne. 0) call getst (j, sbc(1,1,issdiazc13)
     &,       idiazc13)
#   endif
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
            if (isscfc11 .ne. 0) call getst (j, sbc(1,1,isscfc11)
     &,       icfc11)
            if (isscfc12 .ne. 0) call getst (j, sbc(1,1,isscfc12)
     &,       icfc12)
# endif
        enddo
      endif

#endif
!-----------------------------------------------------------------------
!     check for consistency in the S.B.C. setup
!-----------------------------------------------------------------------

      call chkcpl

#if defined O_global_sums || defined O_co2emit_diag
!-----------------------------------------------------------------------
!     get global sums at the start of the run
!-----------------------------------------------------------------------

      call globalsum (1)

#endif
!-----------------------------------------------------------------------
!     S T A R T    S E G M E N T    L O O P
!-----------------------------------------------------------------------

      do n=1,numseg

!-----------------------------------------------------------------------
!       get the atmospheric S.B.C.
!-----------------------------------------------------------------------

        call gasbc (1, imt, 1, jmt)

!-----------------------------------------------------------------------
!       call the atmospheric model once for each time step until one
!       segment of "segtim" days is complete. hold atmos S.B.C. fixed
!       during each segment and predict average S.B.C. for ocean
!-----------------------------------------------------------------------

        do loop=1,ntspas
#if defined O_embm
          call embm (1, imt, 1, jmt)

# if !defined O_mom
#  if defined O_mtlm
          call mtlmout (1, imt, 1, jmt)
#  endif
          if (tsits .and. iotsi .ne. stdout .and. iotsi .gt. 0) then
            write (*,'(1x, a3, i7, 1x, a32)') 'ts=',itt, stamp
          endif
# endif
#else
          call data (1,imt,1,jmt)
#endif
        enddo

#if defined O_mtlm
!-----------------------------------------------------------------------
!       get land S.B.C.s
!-----------------------------------------------------------------------

        call glsbc (1, imt, 1, jmt)

!----------------------------------------------------------------------
!       call the land-surface and vegetation  model once for each time
!       step until one segment of "segtim" days is complete.
!-----------------------------------------------------------------------

        do loop=1,ntspls
          call mtlm (1, imt, 1, jmt)
        enddo

#endif
!-----------------------------------------------------------------------
!       get ocean S.B.C.s
!-----------------------------------------------------------------------

        call gosbc (1, imt, 1, jmt)

#if defined O_mom
!-----------------------------------------------------------------------
!       call the ocean model once for each time step until one
!       segment of "segtim" days is complete. hold ocean S.B.C. fixed
!       during each segment and predict average S.B.C. for atmos
!-----------------------------------------------------------------------

        do loop=1,ntspos
          call mom
# if defined O_embm
          call embmout (1, imt, 1, jmt)
# endif
# if defined O_ism
          if (abs(relyr/dtismyr - nint(relyr/dtismyr)) .lt. 1.e-6)
     &      call ism (1, imt, 1, jmt)
          call ismout (1, imt, 1, jmt)
# endif
# if defined O_mtlm
          call mtlmout (1, imt, 1, jmt)
# endif
# if defined O_sed
          if (abs(relyr/dtsedyr - nint(relyr/dtsedyr)) .lt. 1.e-6)
     &      call sed (1, imt, 1, jmt)
          call sedout (1, imt, 1, jmt)
# endif
          if (tsits .and. iotsi .ne. stdout .and. iotsi .gt. 0) then
            write (*,'(1x, a3, i7, 1x, a32)') 'ts=',itt, stamp
          endif
        enddo

#endif
#if defined O_global_sums || defined O_co2emit_diag
!-----------------------------------------------------------------------
!       write change in global sums for heat and fresh water
!-----------------------------------------------------------------------

        if (tsits) call globalsum (2)

#endif
!-----------------------------------------------------------------------
!       close any open netcdf files to flush buffers
!-----------------------------------------------------------------------

        call closeall

      enddo

!-----------------------------------------------------------------------
!     E N D    S E G M E N T    L O O P
!-----------------------------------------------------------------------

#if defined O_global_sums || defined O_co2emit_diag
!-----------------------------------------------------------------------
!     get global sums at the end of the run
!-----------------------------------------------------------------------

      call globalsum (3)

#endif
      print*, ' ==>  UVIC_ESCM integration is complete.'

      call closeall

      call release_all

      stop
      end

#endif ! O_TMM

#ifndef O_TMM
      subroutine chkcpl

      implicit none

      integer jrow, j, i, k

      logical errorc

      real critv, tmp, t1, r1, r2, r3, r4, r5, r6

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "csbc.h"
#include "cembm.h"
#include "calendar.h"
#include "switch.h"
#include "tmngr.h"

#if defined O_sealev_salinity && O_mom
#include "mw.h"
#include "iounit.h"
#include "levind.h"
#include "grdvar.h"

!-----------------------------------------------------------------------
!     adjust salinity for changes in sea level
!-----------------------------------------------------------------------

      if (dsealev .ne. 0.) then

        tmp = -socn*tcella(1)*dsealev/tcellv
        print*, ' '
        print*, '==> Warning: salinity changed due to sea level change'
        print*, '             sea level change (m):', dsealev*0.01
        print*, '             salinity change (psu):', tmp*1000.
        print*, ' '

        do jrow=1,jmt
          if (wide_open_mw) then
            j = jrow
          else
            j = jmw
            call getrow (latdisk(taudisk), nslab, jrow
     &,                  u(1,1,j,1,tau), t(1,1,j,1,tau))
            call getrow (latdisk(taup1disk), nslab, jrow
     &,                  u(1,1,j,1,taup1), t(1,1,j,1,taup1))
          endif
          do i=2,imtm1
            do k=1,kmt(i,jrow)
              t(i,k,j,isalt,tau) = t(i,k,j,isalt,tau) + tmp
              t(i,k,j,isalt,taup1) = t(i,k,j,isalt,taup1) + tmp
            enddo
          enddo
!         initialize every latitude
          if (wide_open_mw) then
!           do nothing since "tau" data is in place in the MW
          else
            call putrow (latdisk(taudisk), nslab, jrow
     &,                  u(1,1,j,1,tau), t(1,1,j,1,tau))
            call putrow (latdisk(taup1disk), nslab, jrow
     &,                  u(1,1,j,1,taup1), t(1,1,j,1,taup1))
          endif
        enddo
      endif
      dsealev = 0.
#endif

!-----------------------------------------------------------------------
!     adjust output intervals for forcing acceleration
!-----------------------------------------------------------------------

      if (accel .le. 1) then
        accel = 1.
      else
        print*, ' '
        print*, '==> Warning: any forcing that changes on time scales'
     &        ,             ' longer than 1 year has been accelerated'
        print*, '             acceleration factor: ', accel
        print*, '             starts at relyr: ', accel_yr0
        print*, '             current year: ', year0 + relyr
        tmp = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        print*, '             current accelerated year: ', tmp
        print*, ' '
        print*, '==> Warning: adjusting output intervals for'
     &,  ' long term forcing acceleration'
        print*, ' '
        if (tsiint .ge. yrlen*accel) tsiint = tsiint/accel
        if (tavgint .ge. yrlen*accel) tavgint = tavgint/accel
        if (tmbint .ge. yrlen*accel) tmbint = tmbint/accel
        if (stabint .ge. yrlen*accel) stabint = stabint/accel
        if (zmbcint .ge. yrlen*accel) zmbcint = zmbcint/accel
        if (glenint .ge. yrlen*accel) glenint = glenint/accel
        if (trmbint .ge. yrlen*accel) trmbint = trmbint/accel
        if (vmsfint .ge. yrlen*accel) vmsfint = vmsfint/accel
        if (gyreint .ge. yrlen*accel) gyreint = gyreint/accel
        if (extint .ge. yrlen*accel) extint = extint/accel
        if (prxzint .ge. yrlen*accel) prxzint = prxzint/accel
        if (exconvint .ge. yrlen*accel) exconvint = exconvint/accel
        if (dspint .ge. yrlen*accel) dspint = dspint/accel
        if (timavgint .ge. yrlen*accel) timavgint = timavgint/accel
        if (cmixint .ge. yrlen*accel) cmixint = cmixint/accel
        if (xbtint .ge. yrlen*accel) xbtint = xbtint/accel
        if (crossint .ge. yrlen*accel) crossint = crossint/accel
        if (densityint .ge. yrlen*accel) densityint = densityint/accel
        if (fctint .ge. yrlen*accel) fctint = fctint/accel
        if (tyzint .ge. yrlen*accel) tyzint = tyzint/accel
        if (restint .ge. yrlen*accel) restint = restint/accel
        if (tbtint .ge. yrlen*accel) tbtint = tbtint/accel
      endif

      if (mod(timavgint, segtim) .gt. 1.e-6 .and. timavgint .gt. 0.)
     &  then
        t1 = nint(timavgint/segtim)*segtim
        if (t1 .lt. segtim) t1 = segtim
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: "timavgint" does not contain an integral number'
     &,   '              of coupling time steps "segtim".              '
        write (stdout,*) '              (changed "timavgint" from '
     &, timavgint,' days to ', t1,' days to insure this condition)'
        timavgint = t1
      endif
      if (timavgint .eq. 0.) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: averaging interval "timavgint" = 0. implies no '
     &,   '             averaging when "time_averages" is enabled      '
      endif
      if (timavgint .gt. timavgper) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: the interval "timavgint" exceeds the averaging '
     &,   '             period "timavgper" for option "time_averages"  '
      endif
      if (timavgint .lt. timavgper) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: averaging period "timavgper" exceeds interval  '
     &,   '             "timavgint". Setting timavgper = timavgint     '
        timavgper = timavgint
      endif
      write (stdout,'(/,1x,a,f10.2,a,/,1x,a,f10.2,a)')
     &  '==> Time averages written every ', timavgint, ' days, '
     &, '    with an averaging period of ', timavgper, ' days. '

      if (mod(tsiint, segtim) .gt. 1.e-6 .and. tsiint .gt. 0.)
     &  then
        t1 = nint(tsiint/segtim)*segtim
        if (t1 .lt. segtim) t1 = segtim
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: "tsiint" does not contain an integral number'
     &,   '              of coupling time steps "segtim".              '
        write (stdout,*) '              (changed "tsiint" from '
     &, tsiint,' days to ', t1,' days to insure this condition)'
        tsiint = t1
      endif
      if (tsiint .eq. 0.) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: averaging interval "tsiint" = 0. implies no '
     &,   '             averaging when "time_averages" is enabled      '
      endif
      if (tsiint .gt. tsiper) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: the interval "tsiint" exceeds the averaging '
     &,   '             period "tsiper"'
      endif
      if (tsiint .lt. tsiper) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: averaging period "tsiper" exceeds interval  '
     &,   '             "tsiint". Setting tsiper = tsiint     '
        tsiper = tsiint
      endif
      write (stdout,'(/,1x,a,f10.2,a,/,1x,a,f10.2,a)')
     &  '==> Time step integrals written every ', tsiint, ' days, '
     &, '    with an averaging period of       ', tsiper, ' days. '

# if defined O_even_fluxes
      if (mod(namix,2) .ne. mod(nint(segtim*daylen/dtatm),2)) then
        write(*,*) '==> Error: time steps between mixing and coupling'
        write(*,*) '           in the atmosphere must both be even to'
        write(*,*) '           use even_fluxes.'
        stop '=>UVic_ESCM'
      endif

      if (mod(nmix,2) .ne. mod(nint(segtim*daylen/dtocn),2)) then
        write(*,*) '==> Error: time steps between mixing and coupling'
        write(*,*) '           in the ocean must both be even to use'
        write(*,*) '           even_fluxes.'
        stop '=>UVic_ESCM'
      endif

#endif
      if (mod(namix,2) .ne. mod(nats,2)) then
        write(*,*) '==> Warning: restart was not saved with even flux'
        write(*,*) '             averaging. Starting with a mixing time'
        write(*,*) '             step in atmosphere.'
        nats = namix
      endif

      if (mod(nmix,2) .ne. mod(nots,2)) then
        write(*,*) '==> Warning: restart was not saved with even flux'
        write(*,*) '             averaging. Starting with a mixing time'
        write(*,*) '             step in ocean.'
        nots = nmix
      endif

!-----------------------------------------------------------------------
!     do consistency checks before allowing model to continue
!-----------------------------------------------------------------------

      errorc = .false.
      write (stdout,*) ' '
      write (stdout,*) '    (checking setup)'

#if defined O_embm && !defined O_mom && !defined O_replacst
      write (stdout,9000)
     &  '==> Error: use option replacst with embm and not mom'
      errorc = .true.
#endif
#if defined O_redi_diffusion || defined O_gent_mcwilliams
# if !defined O_isopycmix
      write (stdout,9000)
     &  '==> Error: use option isopycmix with redi_diffusion '
     &, 'or gent_mcwilliams'
      errorc = .true.
# endif
#endif
      if (dtatm .eq. c0) then
          write (stdout,9000)
     & '==> Error: the atmospheric time step must be set in "setatm"'
          errorc = .true.
          dtatm = 1.e-6
      endif
!      critv = 1.e-6
      critv = 1.e-4
      if (segtim .ne. c0) then
        r1 = rundays/segtim
      else
        r1 = 0.5
      endif
      if (segtim .eq. c0) then
          write (stdout,9000)
     & '==> Error: coupling period "segtim" must be specified'
          errorc = .true.
      endif
      if (abs(r1-nint(r1)) .gt. critv) then
          write (stdout,9000)
     & '==> Error: there must be an integral number of segments '
     &,'"segtim"  within "rundays" (the length of the run)'
          errorc = .true.
      endif
#if defined O_mom
      r2 = segtim/(dtocn*secday)
      if (abs(r2-nint(r2)) .gt. critv) then
        write (stdout,9000)
     & '==> Error: there must be an integral number of density time '
     &,'steps "dtocn"  within "segtim" (the segment time)'
        errorc = .true.
      endif
#endif
#if defined O_embm
      r3 = segtim/(dtatm*secday)
      if (abs(r3-nint(r3)) .gt. critv) then
        write (stdout,9000)
     & '==> Error: there must be an integral number of atmos time '
     &,'steps "dtatm"  within "segtim" (the segment time)'
        errorc = .true.
      endif
#endif
#if defined O_ism
      r4 = (dtism*secday)/segtim
      if (abs(r4-nint(r4)) .gt. critv) then
        write (stdout,9000)
     & '==> Error: there must be an integral number of "segtim" steps '
     &,'    within "dtism"'
        errorc = .true.
      endif
#endif
#if defined O_mtlm
      r5 = segtim/(dtlnd*secday)
      if (abs(r5-nint(r5)) .gt. critv) then
        write (stdout,9000)
     & '==> Error: there must be an integral number of land time '
     &,'steps "dtlnd"  within "segtim" (the segment time)'
        errorc = .true.
      endif
#endif
#if defined O_sed
      r6 = (dtsed*secday)/segtim
      if (abs(r6-nint(r6)) .gt. critv) then
        write (stdout,9000)
     & '==> Error: there must be an integral number of "segtim" steps '
     &,'within "dtsed"'
        errorc = .true.
      endif
#endif
#if defined O_mom && defined O_embm && defined O_restorst
      write (stdout,9000)
     & '==> Warning: restoring to surface tracers (restorst) '
     &,'    with coupled model'
#endif
#if defined O_mom && defined O_embm && defined O_replacst
      write (stdout,9000)
     & '==> Warning: replacing surface tracers (replacst) '
     &,'with coupled model'
#endif
#if defined O_restorst && defined O_replacst
      write (stdout,9000)
     & '==> Warning: both restoring and replacing surface tracers '
     &,'(restorst, replacst ). Replacing takes precidence'
#endif
#if defined O_co2ccn_data_transient && !defined O_co2ccn_data
      write (stdout,9000)
     & '==> Error: transient option needs O_co2ccn_data'
      errorc = .true.
#endif
#if defined O_co2emit_data_transient && !defined O_co2emit_data
      write (stdout,9000)
     & '==> Error: transient option needs O_co2emit_data'
      errorc = .true.
#endif
#if defined O_co2emit_track_sat_transient && !defined O_co2emit_track_sat
      write (stdout,9000)
     & '==> Error: transient option needs O_co2emit_track_sat'
      errorc = .true.
#endif
#if defined O_co2emit_track_co2_transient && !defined O_co2emit_track_co2
      write (stdout,9000)
     & '==> Error: transient option needs O_co2emit_track_co2'
      errorc = .true.
#endif
#if defined O_agric_data_transient && !defined O_agric_data
      write (stdout,9000)
     & '==> Error: transient option needs O_agric_data'
      errorc = .true.
#endif
#if defined O_crop_data_transient && !defined O_crop_data
      write (stdout,9000)
     & '==> Error: transient option needs O_crop_data'
      errorc = .true.
#endif
#if defined O_landice_data_transient && !defined O_landice_data
      write (stdout,9000)
     & '==> Error: transient option needs O_landice_data'
      errorc = .true.
#endif
#if defined O_solar_data_transient && !defined O_solar_data
      write (stdout,9000)
     & '==> Error: transient option needs O_solar_data'
      errorc = .true.
#endif
#if defined O_volcano_data_transient && !defined O_volcano_data
      write (stdout,9000)
     & '==> Error: transient option needs O_volcano_data'
      errorc = .true.
#endif
#if defined O_sulphate_data_transient && !defined O_sulphate_data
      write (stdout,9000)
     & '==> Error: transient option needs O_sulphate_data'
      errorc = .true.
#endif
#if defined O_aggfor_data_transient && !defined O_aggfor_data
      write (stdout,9000)
     & '==> Error: transient option needs O_aggfor_data'
      errorc = .true.
#endif
#if defined O_cfcs_data_transient && !defined O_cfcs_data
      write (stdout,9000)
     & '==> Error: transient option needs O_cfcs_data'
      errorc = .true.
#endif
#if defined O_c14ccn_data_transient && !defined O_c14ccn_data
      write (stdout,9000)
     & '==> Error: transient option needs O_c14ccn_data'
      errorc = .true.
#endif
#if defined O_co2ccn_data|| defined O_co2emit_data || defined O_co2emit_track_co2 || defined O_co2emit_track_sat

          tmp = 0
# if defined O_co2ccn_data
     &        + 1
# endif
# if defined O_co2emit_data
     &        + 1
# endif
# if defined O_co2emit_track_co2
     &        + 1
# endif
# if defined O_co2emit_track_sat
     &        + 1
# endif
# if defined O_co2ccn_user
     &        + 1
# endif
      if (tmp .gt. 1) then
        write (stdout,9000)
     & '==> Error: use only one of O_co2ccn_data, O_co2emit_data, '
     &,'O_co2emit_track_co2, O_co2emit_track_sat or O_co2ccn_user'
        errorc = .true.
      endif
#endif

      write (stdout,*) ' '
      write (stdout,*) '    (End of checks) '
      write (stdout,*) ' '
      if (errorc) stop '=> ERRORS found in chkcpl'

9000  format (/,(1x,a))
      return
      end

#endif ! O_TMM

      subroutine set (index, num, name, text, inc)

!-----------------------------------------------------------------------
!     increment counter, set index and text
!-----------------------------------------------------------------------

      character(*) :: name, text

      name = text
      index = num
      inc = index + 1

      return
      end

#ifndef O_TMM
#if defined O_mom
      subroutine getst (jrow, ocnout, ntabc)

!-----------------------------------------------------------------------
!     read surface tracers from disk row "jrow"
!-----------------------------------------------------------------------

      implicit none

      integer i, jrow, ntabc

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "iounit.h"
#include "mw.h"
#include "tmngr.h"

      real ocnout(imt,jmt)

      call getrow (latdisk(taup1disk), nslab, jrow
     &,          u(1,1,jmw,1,taup1), t(1,1,jmw,1,taup1))
      do i=1,imt
        ocnout(i,jrow) = t(i,1,jmw,ntabc,taup1)
      enddo

      return
      end
#endif

#endif ! O_TMM

      subroutine sbc_init

!-----------------------------------------------------------------------
!     Initialize S.B.C. indices
!-----------------------------------------------------------------------

      implicit none

      integer m

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "csbc.h"

      sbc(:,:,:) = 0.0
      mapsbc(:) = " "
      itaux = 0
      itauy = 0
      iws = 0
      iaca = 0
      isca = 0
      ihflx = 0
      isflx = 0
      isst = 0
      isss = 0
      iro = 0
      iwa = 0
      iwxq = 0
      iwyq = 0
      iwxt = 0
      iwyt = 0
      iwxc = 0
      iwyc = 0
      ipsw = 0
      isu = 0
      isv = 0
      igu = 0
      igv = 0
      issdic = 0
      idicflx = 0
      issdic13 = 0
      idic13flx = 0
      issalk = 0
      ialkflx = 0
      isso2 = 0
      io2flx = 0
      isspo4 = 0
      ipo4flx = 0
      issdop = 0
      idopflx = 0
      issphyt = 0
      iphytflx = 0
      isszoop = 0
      izoopflx = 0
      issdetr = 0
      idetrflx = 0
      issno3 = 0
      ino3flx = 0
      issdon = 0
      idonflx = 0
      issdiaz = 0
      idiazflx = 0
      issdin15 = 0
      idin15flx = 0
      issdon15 = 0
      idon15flx = 0
      issphytn15 = 0
      iphytn15flx = 0
      isszoopn15 = 0
      izoopn15flx = 0
      issdetrn15 = 0
      idetrn15flx = 0
      issdiazn15 = 0
      idiazn15flx = 0
      issdoc13 = 0
      idoc13flx = 0
      issphytc13 = 0
      iphytc13flx = 0
      isszoopc13 = 0
      izoopc13flx = 0
      issdetrc13 = 0
      idetrc13flx = 0
      issdiazc13 = 0
      idiazc13flx = 0
      issc14 = 0
      ic14flx = 0
      isscfc11 = 0
      icfc11flx = 0
      isscfc12 = 0
      icfc12flx = 0
      iat = 0
      irh = 0
      ipr = 0
      ips = 0
      iaws = 0
      iswr = 0
      ilwr = 0
      isens = 0
      ievap = 0
      idtr = 0
      isr = 0
      inpp = 0
      iburn = 0
      isr13 = 0
      inpp13 = 0
      iburn13 = 0
      isr14 = 0
      inpp14 = 0
      iburn14 = 0
      ibtemp = 0
      ibsalt = 0
      ibdic = 0
      ibdicfx = 0
      ibalk = 0
      ibalkfx = 0
      ibo2 = 0
      ircal = 0
      irorg = 0
      ibtemp = 0
      ibsalt = 0
      ibo2 = 0
      ibalk = 0
      ibdic = 0
      ibdicfx = 0
      ibalkfx  = 0

      m = 1
      call set (itaux, m, mapsbc(m), 'taux', m)
      call set (itauy, m, mapsbc(m), 'tauy', m)
      call set (iws, m, mapsbc(m), 'ws', m)
      call set (iaca, m, mapsbc(m), 'a_calb', m)
      call set (isca, m, mapsbc(m), 's_calb', m)
      call set (ihflx, m, mapsbc(m), 'hflx', m)
      call set (isflx, m, mapsbc(m), 'sflx', m)
      call set (isst, m, mapsbc(m), 'sst', m)
      call set (isss, m, mapsbc(m), 'sss', m)
      call set (iro, m, mapsbc(m), 'ro', m)
#if defined O_embm_awind && defined O_embm
      call set (iwa, m, mapsbc(m), 'wa', m)
#endif
#if defined O_embm
      call set (iwxq, m, mapsbc(m), 'wx_q', m)
      call set (iwyq, m, mapsbc(m), 'wy_q', m)
      call set (iwxt, m, mapsbc(m), 'wx_t', m)
      call set (iwyt, m, mapsbc(m), 'wy_t', m)
# if defined O_carbon_co2_2d
      call set (iwxc, m, mapsbc(m), 'wx_c', m)
      call set (iwyc, m, mapsbc(m), 'wy_c', m)
# endif
#endif
#if defined O_shortwave
      call set (ipsw, m, mapsbc(m), 'psw', m)
#endif
#if defined O_ice_evp
      call set (isu, m, mapsbc(m), 'su', m)
      call set (isv, m, mapsbc(m), 'sv', m)
      call set (igu, m, mapsbc(m), 'gu', m)
      call set (igv, m, mapsbc(m), 'gv', m)
#endif
#if defined O_carbon
      call set (issdic, m, mapsbc(m), 'ssdic', m)
      call set (idicflx, m, mapsbc(m), 'dicflx', m)
# if defined O_carbon_13
      call set (issdic13, m, mapsbc(m), 'ssdic13', m)
      call set (idic13flx, m, mapsbc(m), 'dic13flx', m)
# endif
# if defined O_carbon_14
      call set (issc14, m, mapsbc(m), 'ssc14', m)
      call set (ic14flx, m, mapsbc(m), 'c14flx', m)
# endif
#endif
#if defined O_npzd_alk
      call set (issalk, m, mapsbc(m), 'ssalk', m)
      call set (ialkflx, m, mapsbc(m), 'alkflx', m)
#endif
#if defined O_npzd_o2
      call set (isso2, m, mapsbc(m), 'sso2', m)
      call set (io2flx, m, mapsbc(m), 'o2flx', m)
#endif
#if defined O_npzd
      call set (isspo4, m, mapsbc(m), 'sspo4', m)
      call set (ipo4flx, m, mapsbc(m), 'po4flx', m)
      call set (issdop, m, mapsbc(m), 'ssdop', m)
      call set (idopflx, m, mapsbc(m), 'dopflx', m)
# if !defined O_npzd_no_vflux
      call set (issphyt, m, mapsbc(m), 'ssphyt', m)
      call set (iphytflx, m, mapsbc(m), 'phytflx', m)
      call set (isszoop, m, mapsbc(m), 'sszoop', m)
      call set (izoopflx, m, mapsbc(m), 'zoopflx', m)
      call set (issdetr, m, mapsbc(m), 'ssdetr', m)
      call set (idetrflx, m, mapsbc(m), 'detrflx', m)
# endif
# if defined O_npzd_nitrogen
      call set (issno3, m, mapsbc(m), 'ssno3', m)
      call set (ino3flx, m, mapsbc(m), 'no3flx', m)
      call set (issdon, m, mapsbc(m), 'ssdon', m)
      call set (idonflx, m, mapsbc(m), 'donflx', m)
#  if !defined O_npzd_no_vflux
      call set (issdiaz, m, mapsbc(m), 'ssdiaz', m)
      call set (idiazflx, m, mapsbc(m), 'diazflx', m)
#  endif
#  if defined O_npzd_nitrogen_15
      call set (issdin15, m, mapsbc(m), 'ssdin15', m)
      call set (idin15flx, m, mapsbc(m), 'din15flx', m)
      call set (issdon15, m, mapsbc(m), 'ssdon15', m)
      call set (idon15flx, m, mapsbc(m), 'don15flx', m)
#   if !defined O_npzd_no_vflux
      call set (issphytn15, m, mapsbc(m), 'ssphytn15', m)
      call set (iphytn15flx, m, mapsbc(m), 'phytn15flx', m)
      call set (isszoopn15, m, mapsbc(m), 'sszoopn15', m)
      call set (izoopn15flx, m, mapsbc(m), 'zoopn15flx', m)
      call set (issdetrn15, m, mapsbc(m), 'ssdetrn15', m)
      call set (idetrn15flx, m, mapsbc(m), 'detrn15flx', m)
      call set (issdiazn15, m, mapsbc(m), 'ssdiazn15', m)
      call set (idiazn15flx, m, mapsbc(m), 'diazn15flx', m)
#   endif
#  endif
#  if defined O_carbon_13
      call set (issdoc13, m, mapsbc(m), 'ssdoc13', m)
      call set (idoc13flx, m, mapsbc(m), 'doc13flx', m)
#   if !defined O_npzd_no_vflux
      call set (issphytc13, m, mapsbc(m), 'ssphytc13', m)
      call set (iphytc13flx, m, mapsbc(m), 'phytc13flx', m)
      call set (isszoopc13, m, mapsbc(m), 'sszoopc13', m)
      call set (izoopc13flx, m, mapsbc(m), 'zoopc13flx', m)
      call set (issdetrc13, m, mapsbc(m), 'ssdetrc13', m)
      call set (idetrc13flx, m, mapsbc(m), 'detrc13flx', m)
#    if defined O_npzd_nitrogen
      call set (issdiaz, m, mapsbc(m), 'ssdiazc13', m)
      call set (idiazflx, m, mapsbc(m), 'diazc13flx', m)
#    endif
#   endif
#  endif
# endif
#endif
#if defined O_cfcs_data || defined O_cfcs_data_transient
      call set (isscfc11, m, mapsbc(m), 'sscfc11', m)
      call set (icfc11flx, m, mapsbc(m), 'cfc11flx', m)
      call set (isscfc12, m, mapsbc(m), 'sscfc12', m)
      call set (icfc12flx, m, mapsbc(m), 'cfc12flx', m)
#endif
#if defined O_mtlm
      call set (iat, m, mapsbc(m), 'at', m)
      call set (irh, m, mapsbc(m), 'rh', m)
      call set (ipr, m, mapsbc(m), 'pr', m)
      call set (ips, m, mapsbc(m), 'ps', m)
      call set (iaws, m, mapsbc(m), 'aws', m)
      call set (iswr, m, mapsbc(m), 'swr', m)
      call set (ilwr, m, mapsbc(m), 'lwr', m)
      call set (isens, m, mapsbc(m), 'sens', m)
      call set (ievap, m, mapsbc(m), 'evap', m)
      call set (idtr, m, mapsbc(m), 'dtr', m)
#endif
#if defined O_mtlm && defined O_carbon
      call set (isr, m, mapsbc(m), 'sr', m)
      call set (inpp, m, mapsbc(m), 'npp', m)
      call set (iburn, m, mapsbc(m), 'burn', m)
#endif
#if defined O_mtlm_carbon_13
      call set (isr13, m, mapsbc(m), 'sr13', m)
      call set (inpp13, m, mapsbc(m), 'npp13', m)
      call set (iburn13, m, mapsbc(m), 'burn13', m)
#endif
#if defined O_mtlm_carbon_14
      call set (isr14, m, mapsbc(m), 'sr14', m)
      call set (inpp14, m, mapsbc(m), 'npp14', m)
      call set (iburn14, m, mapsbc(m), 'burn14', m)
#endif
#if defined O_sed
      call set (ircal, m, mapsbc(m), 'rcal', m)
      call set (irorg, m, mapsbc(m), 'rorg', m)
      call set (ibtemp, m, mapsbc(m), 'btemp', m)
      call set (ibsalt, m, mapsbc(m), 'bsalt', m)
      call set (ibo2, m, mapsbc(m), 'bo2', m)
      call set (ibalk, m, mapsbc(m), 'balk', m)
      call set (ibdic, m, mapsbc(m), 'bdic', m)
# if defined O_carbon
      call set (ibdicfx, m, mapsbc(m), 'bdicfx', m)
# endif
# if defined O_npzd_alk
      call set (ibalkfx, m, mapsbc(m), 'balkfx', m)
# endif
#endif

      if ( m-1 .gt. numsbc) then
        print*, '==> Error: increase numsbc in csbc.h to ', m-1
        stop '=>UVic_ESCM'
      endif

      return
      end

      subroutine tracer_init

!-----------------------------------------------------------------------
!     Initialize ocean tracer names
!-----------------------------------------------------------------------

      implicit none

      integer m

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "atm.h"
#if defined O_mom
#include "mw.h"

!-----------------------------------------------------------------------
!     Initialize ocean tracer names
!-----------------------------------------------------------------------

      mapt(:) = " "
      itemp = 0
      isalt = 0
      idic = 0
      idic13 = 0
      ic14 = 0
      icfc11 = 0
      icfc12 = 0
      io2 = 0
      ialk = 0
      ipo4 = 0
      idop = 0
      iphyt = 0
      izoop = 0
      idetr = 0
      ino3 = 0
      idon = 0
      idiaz = 0
      idin15 = 0
      idon15 = 0
      iphytn15 = 0
      izoopn15 = 0
      idetrn15 = 0
      idiazn15 = 0
      idoc13 = 0
      iphytc13 = 0
      izoopc13 = 0
      idetrc13 = 0
      idiazc13 = 0

      m = 1
      call set (itemp, m, mapt(m), 'temp', m)
      call set (isalt, m, mapt(m), 'salt', m)
# if defined O_carbon
      call set (idic, m, mapt(m), 'dic', m)
#  if defined O_carbon_13
      call set (idic13, m, mapt(m), 'dic13', m)
#  endif
#  if defined O_carbon_14
      call set (ic14, m, mapt(m), 'c14', m)
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
      call set (icfc11, m, mapt(m), 'cfc11', m)
      call set (icfc12, m, mapt(m), 'cfc12', m)
# endif
# if defined O_npzd_o2
      call set (io2, m, mapt(m), 'o2', m)
# endif
# if defined O_npzd_alk
      call set (ialk, m, mapt(m), 'alk', m)
# endif
# if defined O_npzd
      call set (ipo4, m, mapt(m), 'po4', m)
      call set (idop, m, mapt(m), 'dop', m)
      call set (iphyt, m, mapt(m), 'phyt', m)
      call set (izoop, m, mapt(m), 'zoop', m)
      call set (idetr, m, mapt(m), 'detr', m)
#  if defined O_npzd_nitrogen
      call set (ino3, m, mapt(m), 'no3', m)
      call set (idon, m, mapt(m), 'don', m)
      call set (idiaz, m, mapt(m), 'diaz', m)
#   if defined O_npzd_nitrogen_15
      call set (idin15, m, mapt(m), 'din15', m)
      call set (idon15, m, mapt(m), 'don15', m)
      call set (iphytn15, m, mapt(m), 'phytn15', m)
      call set (izoopn15, m, mapt(m), 'zoopn15', m)
      call set (idetrn15, m, mapt(m), 'detrn15', m)
      call set (idiazn15, m, mapt(m), 'diazn15', m)
#   endif
#  endif
#  if defined O_carbon_13
      call set (idoc13, m, mapt(m), 'doc13', m)
      call set (iphytc13, m, mapt(m), 'phytc13', m)
      call set (izoopc13, m, mapt(m), 'zoopc13', m)
      call set (idetrc13, m, mapt(m), 'detrc13', m)
#   if defined O_npzd_nitrogen
      call set (idiazc13, m, mapt(m), 'diazc13', m)
#   endif
#  endif
# endif
      if ( m-1 .gt. nt) then
        print*, '==> Error: increase nt for tracers in size.h'
        stop '=>UVic_ESCM'
      endif

!-----------------------------------------------------------------------
!     Initialize ocean tracer sourcr names, must have equivalent tracer
!-----------------------------------------------------------------------

      mapst(:) = " "
      itrc(:) = 0

      m = 1
# if defined O_carbon && defined O_carbon_14
      call set (isc14, m, mapst(m), 'sc14', m)
      itrc(ic14) = m-1
# endif
# if defined O_npzd
#  if defined O_carbon
      call set (isdic, m, mapst(m), 'sdic', m)
      itrc(idic) = m-1
#   if defined O_carbon_13
      call set (isdic13, m, mapst(m), 'sdic13', m)
      itrc(idic13) = m-1
#   endif
#  endif
#  if defined O_npzd_o2
      call set (iso2, m, mapst(m), 'so2', m)
      itrc(io2) = m-1
#  endif
#  if defined O_npzd_alk
      call set (isalk, m, mapst(m), 'salk', m)
      itrc(ialk) = m-1
#  endif
      call set (ispo4, m, mapst(m), 'spo4', m)
      itrc(ipo4) = m-1
      call set (isdop, m, mapst(m), 'sdop', m)
      itrc(idop) = m-1
      call set (isphyt, m, mapst(m), 'sphyt', m)
      itrc(iphyt) = m-1
      call set (iszoop, m, mapst(m), 'szoop', m)
      itrc(izoop) = m-1
      call set (isdetr, m, mapst(m), 'sdetr', m)
      itrc(idetr) = m-1
#  if defined O_npzd_nitrogen
      call set (isno3, m, mapst(m), 'sno3', m)
      itrc(ino3) = m-1
      call set (isdon, m, mapst(m), 'sdon', m)
      itrc(idon) = m-1
      call set (isdiaz, m, mapst(m), 'sdiaz', m)
      itrc(idiaz) = m-1
#   if defined O_npzd_nitrogen_15
      call set (isdin15, m, mapst(m), 'sdin15', m)
      itrc(idin15) = m-1
      call set (isdon15, m, mapst(m), 'sdon15', m)
      itrc(idon15) = m-1
      call set (isphytn15, m, mapst(m), 'sphytn15', m)
      itrc(iphytn15) = m-1
      call set (iszoopn15, m, mapst(m), 'szoopn15', m)
      itrc(izoopn15) = m-1
      call set (isdetrn15, m, mapst(m), 'sdetrn15', m)
      itrc(idetrn15) = m-1
      call set (isdiazn15, m, mapst(m), 'sdiazn15', m)
      itrc(idiazn15) = m-1
#   endif
#  endif
#  if defined O_carbon_13
      call set (isdoc13, m, mapst(m), 'sdoc13', m) 
      itrc(idoc13) = m-1
      call set (isphytc13, m, mapst(m), 'sphytc13', m) 
      itrc(iphytc13) = m-1
      call set (iszoopc13, m, mapst(m), 'szoopc13', m) 
      itrc(izoopc13) = m-1
      call set (isdetrc13, m, mapst(m), 'sdetrc13', m) 
      itrc(idetrc13) = m-1
#   if defined O_npzd_nitrogen
      call set (isdiazc13, m, mapst(m), 'sdiazc13', m) 
      itrc(idiazc13) = m-1
#   endif
#  endif
# endif
      if ( m-1 .gt. nt) then
        print*, '==> Error: increase nsrc for tracer sources in size.h'
        stop '=>UVic_ESCM'
      endif

#endif
#if defined O_embm
!-----------------------------------------------------------------------
!     Initialize atmosphere tracer names
!-----------------------------------------------------------------------

      mapat(:) = " "
      isat = 0
      ishum = 0
      ico2 = 0

      m = 1
      call set (isat, m, mapat(m), 'sat', m)
      call set (ishum, m, mapat(m), 'shum', m)
# if defined O_carbon && defined O_carbon_co2_2d
      call set (ico2, m, mapat(m), 'co2', m)
# endif

      if ( m-1 .gt. nat) then
        print*, '==> Error: increase nat in size.h'
        stop '=>UVic_ESCM'
      endif
#endif

      return
      end

#ifndef O_TMM
      subroutine read_namelist

!-----------------------------------------------------------------------
!     read all model namelist variables
!-----------------------------------------------------------------------

      implicit none

      character (120) :: fname, new_file_name, logfile

      integer iotraj, i, j, k, ip, kr, jq, n, ioun, num_processors

      logical exists, annlevobc, annlev, initpt

      real dtatms, snapint, snapls, snaple, snapde, ahbkg, runstep
      real slmx, s_dm, afkph, dfkph, sfkph, zfkph, ahs, ahb, trajint
      real crops_yr

#if !defined O_mom || !defined O_time_step_monitor
      integer isot1, ieot1, isot2, ieot2, jsot, jeot, ksot, keot, mrot
#endif
#if !defined O_mom || !defined O_isopycmix
      real ahisop, del_dm, slmxr, s_dmr
#endif
#if !defined O_mom || !defined O_isopycmix || !defined O_gent_mcwilliams
      real athkdf
#endif

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "accel.h"
#include "calendar.h"
#include "cembm.h"
#include "coord.h"
#include "cnep.h"
#include "csbc.h"
#include "cprnts.h"
#include "diag.h"
#include "emode.h"
#include "fwa.h"
#include "hmixc.h"
#include "insolation.h"
#include "iounit.h"
#if defined O_mom && defined O_isopycmix
#include "isopyc.h"
#endif
#include "mtlm.h"
#include "npzd.h"
#include "scalar.h"
#include "sed.h"
#include "stab.h"
#include "switch.h"
#include "tmngr.h"
#include "veg.h"
#include "vmixc.h"

!-----------------------------------------------------------------------
!     read all model namelist variables
!-----------------------------------------------------------------------

      namelist /contrl/ init, runlen, rununits, restrt, initpt
     &,                 num_processors, runstep
      namelist /tsteps/ dtts, dtuv, dtsf, dtatm, dtatms, namix, segtim
     &,                 daylen
      namelist /riglid/ mxscan, sor, tolrsf, tolrsp, tolrfs
      namelist /mixing/ am, ah, ahbkg, ambi, ahbi, kappa_m, kappa_h
     &,                 cdbot, spnep, senep, aidif, ncon, nmix, eb
     &,                 acor, dampts, dampdz, annlev, annlevobc
      namelist /diagn/  tsiint, tsiper, tavgint, itavg, tmbint, tmbper
     &,                 itmb, stabint, zmbcint, glenint, trmbint, itrmb
     &,                 vmsfint, gyreint,igyre, extint, prxzint, trajint
     &,                 exconvint, dspint, dspper, snapint, snapls
     &,                 snaple, snapde, timavgint, timavgper, cmixint
     &,                 prlat, prslon, prelon, prsdpt, predpt
     &,                 slatxy, elatxy, slonxy, elonxy
     &,                 cflons, cflone, cflats, cflate, cfldps, cfldpe
     &,                 maxcfl, xbtint, xbtper, crossint, densityint
     &,                 fctint, tyzint, restint, tbtint, tbtper
     &,                 accel, accel_yr0
      namelist /io/     expnam, logfile, runstamp, iotavg, iotmb, iotrmb
     &,                 ioglen, iovmsf, iogyre, ioprxz, ioext, iodsp
     &,                 iotsi, iozmbc, iotraj, ioxbt, isot1, ieot1
     &,                 isot2, ieot2, jsot, jeot, ksot, keot, mrot
      namelist /ictime/ year0, month0, day0, hour0, min0, sec0, ryear
     &,                 rmonth, rday, rhour, rmin, rsec, refrun, refinit
     &,                 refuser, eqyear, eqmon, monlen, init_time
     &,                 init_time_in, init_time_out, omega
      namelist /fwa/    isfwa1, iefwa1, isfwa2, iefwa2, jsfwa, jefwa
     &,                 mrfwa, fwaflxi, fwayri, fwayrf, fwarate
     &,                 compensate
      namelist /blmix/  Ahv, Ahh
      namelist /hlmix/  hl_depth, hl_back, hl_max
      namelist /isopyc/ slmx, ahisop, athkdf, del_dm, s_dm
      namelist /ppmix/  wndmix, fricmx, diff_cbt_back, visc_cbu_back
     &,                 visc_cbu_limit, diff_cbt_limit
      namelist /smagnl/ diff_c_back
      namelist /sed/    dtsedyr, nsedacc, weath, kmin
      namelist /embm/   rlapse, rf1, rf2, scatter, rhmax, vcsref, vcsfac
     &,                 vcsyri, aggfor_os, adiff
      namelist /carbon/ co2ccn, co2emit, co2for, co2_yr, c14_yr, dc14ccn
     &,                 c14prod, dc13ccn
      namelist /paleo/  orbit_yr, pyear, eccen, obliq, mvelp, sealev
     &,                 sealev_yr, volcano_yr, sulph_yr, aggfor_yr
     &,                 cfcs_yr
      namelist /ice/    niats, nivts, dampice, tsno, hsno_max, ice_yr
     &,                 landice_yr, ice_calb, sno_calb
      namelist /veg/    veg_alb, veg_rl, veg_rs, veg_smd, agric_yr
     &,                 crops_yr, iagric, icrops, idesert, iice
      namelist /solar/  solarconst, solar_yr
      namelist /mtlm/   TIMESTEP, INT_VEG, VEG_EQUIL, DAY_TRIF, DAY_PHEN
     &,                 BF

!     physical constants

      rho0 = 1.035
      rho0r = c1/rho0
      grav = 980.6
      radius = 6370.0e5
      pi = atan(1.0) * 4.0

!     set defaults for namelist contrl
      init = .true.
      runlen = 365.0
      rununits = 'days'
      restrt = .true.
      initpt = .false.
      num_processors = 1
      runstep = -1.0

!     set defaults for namelist tsteps
      dtts = 43200.0
      dtuv = 600.0
      dtsf = 600.0
      dtatm = 43200.0
      dtatms = 43200.0
      namix = 16
      segtim = 1.0
      daylen = 86400.0

!     set defaults for namelist riglid
      mxscan = 300
      sor = 1.60
      tolrsf = 5.0e8
      tolrsp = 1.0e4
      tolrfs = 1.0e4

!     set defaults for namelist mixing
      am = 2.0e9
      ah = 2.0e7
      ahbkg = 0.
      ambi = 1.0e23
      ahbi = 5.0e22
      kappa_m = 10.0
      kappa_h = 1.0
      cdbot = 1.3e-3
      spnep = 3.0e5
      senep = 12.0e5
#if defined O_implicitvmix || defined O_isopycmix || defined O_redi_diffusion
      aidif = 1.0
#else
      aidif = 0.0
#endif
      ncon = 1
      nmix = 16
      eb = .false.
      acor = 0.0
      do n=1,nt
        dampts(n) = 50.0
        dampdz(n) = 26.575e2
      enddo
      annlev = .false.
      annlevobc = .false.

!     set defaults for namelist diagn
      tsiint = 1.0
      tsiper = 1.0
      tavgint = -36500.0
      itavg = .true.
      tmbint = -36500.0
      tmbper = 365.0
      itmb = .true.
      stabint = -36500.0
      zmbcint = -36500.0
      glenint = -36500.0
      trmbint = -36500.0
      itrmb = .true.
      vmsfint = -36500.0
      gyreint = -36500.0
      igyre = .true.
      extint = -36500.0
      prxzint = -36500.0
      trajint = 0.0
      exconvint = -36500.0
      dspint = -36500.0
      dspper = -365.0
      snapint = 0.0
      snapls = 0.0
      snaple = 0.0
      snapde = 0.0
      timavgint = 36500.0
      timavgper = 365.0
      cmixint = -36500.0
      do n=1, nlatpr
        prlat(n) = 100.0
        prslon(n) = 0.0
        prelon(n) = 0.0
        prsdpt(n) = 0.0
        predpt(n) = 6000.0e2
        if (n. le. 4) then
          prslon(n) = 180.0
          prelon(n) = 250.0
        endif
      enddo
      prlat(1) = -60.0
      prlat(2) = 0.0
      prlat(3) = 27.0
      prlat(4) = 55.0
      slatxy = -90.0
      elatxy = 90.0
      slonxy = 3.0
      elonxy = 357.0
      cflons = 0.0
      cflone = 360.0
      cflats = -90.0
      cflate = 90.0
      cfldps = 0.0
      cfldpe = 6000.0e2
      maxcfl = 3
      xbtint = -36500.0
      xbtper = -365.0
      crossint = 365000.0
      densityint = -36500.0
      fctint = -36500.0
      tyzint = -36500.0
      restint = 36500.0
      tbtint = -36500.0
      tbtper = -365.0
      accel = 1.
      accel_yr0 = 0.

!     set defaults for namelist io
      expnam = ' '
      logfile = 'machine.log'
      runstamp = ' '
      restrt = .false.
      iotavg = -1
      iotmb = -1
      iotrmb = -1
      ioglen = -1
      iovmsf = -1
      iogyre = -1
      ioprxz = -1
      ioext = -1
      iodsp = -1
      iotsi = -1
      iozmbc = -1
      iotraj = -1
      ioxbt = -1
      isot1 = 2
      ieot1 = imtm1
      isot2 = 2
      ieot2 = 1
      jsot = 2
      jeot = jmtm1
      ksot = 1
      keot = km
      mrot = 0

!     set defaults for namelist ictime
      year0 = 1
      month0 = 1
      day0 = 1
      hour0 = 0
      min0 = 0
      sec0 = 0
      ryear = 1
      rmonth = 1
      rday = 1
      rhour = 0
      rmin = 0
      rsec = 0
      refrun = .false.
      refinit = .true.
      refuser = .false.
      eqyear = .true.
      eqmon = .false.
      monlen = 30
      init_time = .false.
      init_time_in = .false.
      init_time_out = .false.
      omega = pi/43082.0

!     set defaults for namelist fwa
      isfwa1 = 1
      iefwa1 = imt
      isfwa2 = 0
      iefwa2 = 0
      jsfwa = 1
      jefwa = jmt
      mrfwa = 1
      fwaflxi = 0.
      fwayri = -1.e20
      fwayrf = 1.e20
      fwarate = 0.
      compensate = .false.

!     set defaults for namelist blmix
!     Reference:
!     A Water Mass Model of the World Ocean  K. Bryan, L.J. Lewis
!     JGR, vol 84, No. C5, May 20, 1979
      afkph = 0.8
      dfkph = 1.05
      sfkph = 4.5e-5
      zfkph = 2500.0e2
!     Use Bryan & Lewis values for vertical tracer diffusion
!     Ahv range of 0.3 to 1.3, crossover at 2500m.
!     compute depth dependent vertical diffusion coefficients for
!     tracers using the relationship of Bryan and Lewis
      pi = 4.0 * atan(1.0)
# if defined O_bryan_lewis_vertical
      do k=1,km
        Ahv(k) = (afkph + (dfkph/pi)*(atan(sfkph*(zw(k) - zfkph))))
      enddo
# endif
!     compute depth dependent horizontal diffusion coefficients for
!     tracers using the relationship of Bryan and Lewis
      ahs = 5.0e+3
      ahb = 1.0e+3
# if defined O_bryan_lewis_vertical
      do k=1,km
        Ahh(k) = (ahb + (ahs - ahb)*exp(-zt(k)/50000.0))*1.0e4
      enddo
# endif

!     set defaults for namelist hlmix
      hl_depth = 500.0e2
      hl_back  = 1.e4
      hl_max   = 1.e9

!     set defaults for namelist isopyc
      slmx  = 1.0/100.0   ! maximum isopycnal slope
      ahisop = 1.e7       ! isopycnal diffusion coefficient
      athkdf = 1.0e7      ! isopycnal thickness diffusion coefficient
      del_dm = 4.0/1000.0 ! transition for scaling diffusion coefficient
      s_dm = 1.0/1000.0   ! half width scaling for diffusion coefficient

!     set defaults for namelist ppmix
      wndmix    = 10.0
      fricmx    = 50.0
      diff_cbt_back =  0.1
      visc_cbu_back =  1.0
# if defined O_implicitvmix
!     simulate convective adjustment with large mixing coefficient limits
      visc_cbu_limit = fricmx
      diff_cbt_limit = 1.0e6
# else
!     in regions of gravitational instability set mixing limits to the
!     maximum consistant with the "cfl" criterion. convective adjustment
!     will also act on the instability.
      visc_cbu_limit = fricmx
      diff_cbt_limit = fricmx
# endif

!     set defaults for namelist smagnl
      diff_c_back = c0

!     set defaults for namelist sed
      dtsedyr = 1.   ! time step of sediment model in years
      nsedacc = 1    ! number of steps for accelerating sediments
      weath = 2.e20  ! weathering flux in Pg year-1
                     ! set high to check for constant namelist input
      kmin = 8       ! minimum ocean level for sediments

!     set defaults namelist embm
      rlapse     = 5.e-5
      rf1        = 0.3
      rf2        = 3.e5
      scatter    = 0.23
      rhmax      = 0.85
      vcsref     = 0.
      vcsfac     = 0.
      vcsyri     = 0.
      aggfor_os  = 0.
      adiff      = 0.

!     set defaults namelist carbon
      co2ccn     = 280.
      co2emit    = 0.
      co2for     = 5.35e03
      co2_yr     = 2.e20
      c14_yr     = 2.e20
      dc14ccn    = 0.
      c14prod    = 0.
      dc13ccn    = -6.5

!     set defaults namelist paleo
      pyear      = 2.e20
      orbit_yr   = 2.e20
!     default user orbit is for 1950
      eccen      = 1.672393E-02
      obliq      = 23.44627
      mvelp      = 102.0391
      sealev     = 0.
      sealev_yr  = 2.e20
      volcano_yr = 2.e20
      sulph_yr   = 2.e20
      aggfor_yr  = 2.e20
      cfcs_yr    = 2.e20

!     set defaults namelist ice
      niats      = 1
      nivts      = 1
      dampice    = 5.
      tsno       = 0.
      hsno_max   = 1000.
# if defined O_embm_snow_transient
!     no limit
      hsno_max   = 1.e30
# endif
      ice_yr     = 2.e20
      landice_yr = 2.e20
      ice_calb   = 0.25
      sno_calb   = 0.2

!     set defaults namelist veg
      veg_alb(1:7) = (/0.17, 0.17, 0.22, 0.22, 0.22, 0.30, 0.60/)
      veg_rl(1:7) = (/1.5, 1.0, 0.1, 0.3, 0.01, 0.01, 0.005/)
      veg_rs(1:7) = (/0.75, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0/)
      veg_smd(1:7) = (/2.0, 3.0, 0.1, 0.1, 0.5, 0.01, 0.01/)
      agric_yr = 2.e20
      crops_yr = 2.e20
      iagric = 3
      icrops = 3
      idesert = 6
      iice = 7

!     set defaults namelist solar
      solarconst = 1.368e6  ! mw m-2
      solar_yr = 2.e20

!     set defaults namelist mtlm
      TIMESTEP = 3600.
      INT_VEG = .true.
      VEG_EQUIL = .false.
      DAY_TRIF = 30
      DAY_PHEN = 1
      BF = 1.
    

!-----------------------------------------------------------------------
!     provide for change in above presets using "namelist"
!-----------------------------------------------------------------------

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, contrl, end=101)
101   continue
      runstep = float(int(runstep))
      write (stdout, contrl)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, tsteps, end=102)
102   continue
      write (stdout, tsteps)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, riglid, end=103)
103   continue
      write (stdout, riglid)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, mixing, end=104)
104   continue
      write (stdout, mixing)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, diagn, end=105)
105   continue
      write (stdout, diagn)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, io, end=106)
106   continue
      write (stdout, io)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, ictime, end=107)
107   continue
      write (stdout, ictime)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, fwa, end=109)
109   continue
      write (stdout, fwa)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, blmix, end=110)
110   continue
      write (stdout, blmix)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read (ioun, hlmix, end=111)
      write (stdout, hlmix)
      call relunit (ioun)
111   continue

      call getunit (ioun, 'control.in', 'f s r')
      read (ioun, isopyc, end=112)
112   continue
      write (stdout, isopyc)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read (ioun, ppmix, end=113)
113   continue
      write (stdout, ppmix)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read (ioun, smagnl, end=114)
114   continue
      write (stdout, smagnl)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, sed, end=115)
115   continue
      write (stdout, sed)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, embm, end=116)
116   continue
      write (stdout, embm)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, carbon, end=117)
117   continue
      write (stdout, carbon)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, paleo, end=118)
118   continue
      write (stdout, paleo)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, ice, end=119)
119   continue
      write (stdout, ice)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, veg, end=120)
120   continue
      write (stdout, veg)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, solar, end=121)
121   continue
      write (stdout, solar)
      call relunit (ioun)

      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, mtlm, end=122)
122   continue
      write (stdout, mtlm)
      call relunit (ioun)

!     limit or convert some variables

# if defined O_embm
!     ensure pass is between zero and one.
      pass =  min(max((1. - scatter), 0.), 1.)
# else
      dtatm = dtts
      dtatms = dtts
      namix = nmix
# endif
# if defined O_mom
      dtocn = dtts
      nots = nmix
      isot1 = max(isot1, 1)
      ieot1 = min(ieot1, imt)
      isot2 = max(isot1, 1)
      ieot2 = min(ieot1, imt)
      jsot = max(jsot, 1)
      jeot = min(jeot, jmt)
      ksot = max(ksot, 1)
      keot = min(keot, km)
      do k=1,km
        dtxcel(k) = 1.0
      enddo
#  if defined O_isopycmix
      ah = ahbkg
      slmxr = 0.
      if (slmx .gt. 0.) slmxr = c1/slmx
      s_dm = 0.
      if (s_dm .gt. 0.) s_dmr = c1/s_dm
#  endif
# else
!     if not coupled set ocean and coupling time to atmosphere
      dtts = dtatm
      dtocn = dtatm
      nots = namix
# endif
# if defined O_carbon_14

!     calculate c14ccn from dc14ccn and co2ccn
      c14ccn = (1 + dc14ccn*0.001)*rc14std*co2ccn
     &      /(1+(1 + dc14ccn*0.001)*rc14std)
# endif
# if defined O_carbon_13

!     calculate c13ccn from dc13ccn and co2ccn
      c13ccn = (1 + dc13ccn*0.001)*rc13std*co2ccn
     &      /(1+(1 + dc13ccn*0.001)*rc13std)
# endif

!     sort out forcing years
      if (pyear  .gt. 1.e20) pyear = 1800.
      if (co2_yr .gt. 1.e20) co2_yr = pyear
      if (c14_yr .gt. 1.e20) c14_yr = pyear
      if (orbit_yr .gt. 1.e20) orbit_yr = pyear
      if (volcano_yr .gt. 1.e20) volcano_yr = pyear
      if (sulph_yr .gt. 1.e20) sulph_yr = pyear
      if (aggfor_yr .gt. 1.e20) aggfor_yr = pyear
      if (cfcs_yr .gt. 1.e20) cfcs_yr = pyear
      if (sealev_yr .gt. 1.e20) sealev_yr = pyear
      if (solar_yr .gt. 1.e20) solar_yr = pyear
      if (agric_yr .gt. 1.e20 .and. crops_yr .lt. 1.e20)
     &  agric_yr = crops_yr
      if (agric_yr .gt. 1.e20) agric_yr = pyear
      crops_yr = agric_yr
      if (landice_yr .gt. 1.e20 .and. ice_yr .lt. 1.e20)
     &  landice_yr = ice_yr
      if (landice_yr .gt. 1.e20) landice_yr = pyear
      ice_yr = landice_yr

      tsiper = min(tsiper, tsiint)
      tbtper = min(tbtper, tbtint)
      if (accel .le. 1) accel = 1.
      runlen = runlen/accel
      runstep = runstep/accel
      if (runstep .gt. 0.0) runlen = min(runstep, runlen)
      if (init_time) then
        init_time_in = .true.
        init_time_out = .true.
      endif

      monlen = 30
# if defined O_calendar_360_day
      eqyear = .true.
      eqmon = .true.
      calendar = '360_day'
# elif defined O_calendar_gregorian
      eqyear = .false.
      eqmon = .false.
      calendar = 'gregorian'
# else
      eqyear = .true.
      eqmon = .false.
      calendar = 'noleap'
# endif

!-----------------------------------------------------------------------
!     set runstamp
!-----------------------------------------------------------------------
      if (runstamp .eq. ' ') then
!       read runstamp from last line of logfile if it exists
        fname = new_file_name (logfile)
        inquire (file=trim(fname), exist=exists)
        if (exists) then
          call getunit (ioun, fname, 'f s r')
          do while (exists)
            read  (ioun, '(a130)', end=130) runstamp
          enddo
130       continue
          call relunit (ioun)
          print*, " "
          print*, "runstamp: ", trim(runstamp)
        endif
      endif

      return
      end

#endif ! O_TMM
