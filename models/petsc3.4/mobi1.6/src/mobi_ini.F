#include "MOBI_TMM_OPTIONS.h"
      SUBROUTINE MOBI_INI(z, drF, DeltaT,
     &                    S_surf_glob, TR_surf_glob, debugFlag)

      implicit none
#include "size.h"
#include "npzd.h"
#include "calendar.h"
#include "coord.h"
#include "grdvar.h"
#include "pconst.h"
#include "scalar.h"
#include "csbc.h"
#include "mw.h"
#include "atm.h"
#include "cembm.h"

      real z(km), drF(km), DeltaT
      real S_surf_glob, TR_surf_glob(nsrc)
      integer debugFlag
      
      integer m, n, k, itr

! set some constants
      pi = 4.0 * atan(1.0)
      radian = 360./(2.*pi)
      rho0 = 1.035
      rho0r = c1/rho0
      socn = 0.03475
Csetmom.F:        socn = gaost(isalt)

      
!-----------------------------------------------------------------------
!     Initialize ocean tracer names
!-----------------------------------------------------------------------
! SPK: taken from S/R tracer_init

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

      if (debugFlag.eq.1) then
		do n=1,nt
		  print*,'n=',n,mapt(n)
		enddo
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

!-----------------------------------------------------------------------
!     Initialize S.B.C. indices
!-----------------------------------------------------------------------
! SPK: taken from S/R sbc_init

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

CSPK taken from setmom.F
C set initial surface values
      gaost(:) = 0.
      if (isalt .ne. 0 .and. isss .ne. 0) then
        gaost(isalt) = S_surf_glob
C        call areaavg (sbc(1,1,isss), dmsk, gaost(isalt))
        gaost(isalt) = gaost(isalt) + 0.035
        socn = gaost(isalt)
      endif
      itr = 0
      do n=3,nt
        if (itrc(n) .ne. 0) then
         itr = itr + 1
		 gaost(n) = TR_surf_glob(itr)		 
        endif
      enddo
CSPK check
      if (itr.ne.nsrc) then
        print*,'ERROR in mobi_ini: number of tracers incorrect!'
        stop
      endif  
      
      if (debugFlag.eq.1) then
      if (idic .ne. 0 .and. issdic .ne. 0) then
C        gaost(idic) = 
C        call areaavg (sbc(1,1,issdic), dmsk, gaost(idic))
        print*, "global average sea surface dic (mol m-3) = "
     &,   gaost(idic)
      endif
      if (idic13 .ne. 0 .and. issdic13 .ne. 0) then
C        gaost(idic13) = 
C        call areaavg (sbc(1,1,issdic13), dmsk, gaost(idic13))
        print*, "global average sea surface dic 13 (mol m-3) = "
     &,   gaost(idic13)
      endif
      if (io2 .ne. 0 .and. isso2 .ne. 0) then
C        gaost(io2) = 
C        call areaavg (sbc(1,1,isso2), dmsk, gaost(io2))
        print*, "global average sea surface oxygen (mol m-3) = "
     &,   gaost(io2)
      endif
      if (ialk .ne. 0 .and. issalk .ne. 0) then
C        gaost(ialk) = 
C        call areaavg (sbc(1,1,issalk), dmsk, gaost(ialk))
        print*, "global average sea surface alkalinity (mol m-3) = "
     &,   gaost(ialk)
      endif
      if (ipo4 .ne. 0 .and. isspo4 .ne. 0) then
C        gaost(ipo4) = 
C        call areaavg (sbc(1,1,isspo4), dmsk, gaost(ipo4))
        print*, "global average sea surface phosphate (mol m-3) = "
     &,   gaost(ipo4)*0.001
      endif
      if (idop .ne. 0 .and. issdop .ne. 0) then
C        gaost(idop) = 
C        call areaavg (sbc(1,1,issdop), dmsk, gaost(idop))
        print*, "global average sea surface DOP (mol m-3) = "
     &,   gaost(idop)*0.001
      endif
# if !defined O_npzd_no_vflux
      if (iphyt .ne. 0 .and. issphyt .ne. 0) then
C        gaost(iphyt) = 
C        call areaavg (sbc(1,1,issphyt), dmsk, gaost(iphyt))
        print*, "global average sea surface phytoplankton (mol m-3) = "
     &,   gaost(iphyt)*0.001
      endif
      if (izoop .ne. 0 .and. isszoop .ne. 0) then
C        gaost(izoop) = 
C        call areaavg (sbc(1,1,isszoop), dmsk, gaost(izoop))
        print*, "global average sea surface zooplankton (mol m-3) = "
     &,   gaost(izoop)*0.001
      endif
      if (idetr .ne. 0 .and. issdetr .ne. 0) then
C        gaost(idetr) = 
C        call areaavg (sbc(1,1,issdetr), dmsk, gaost(idetr))
        print*, "global average sea surface detritus (mol m-3) = "
     &,   gaost(idetr)*0.001
      endif
# endif
      if (ino3 .ne. 0 .and. issno3 .ne. 0) then
C        gaost(ino3) = 
C        call areaavg (sbc(1,1,issno3), dmsk, gaost(ino3))
        print*, "global average sea surface nitrate (mol m-3) = "
     &,   gaost(ino3)*0.001
      endif
      if (idon .ne. 0 .and. issdon .ne. 0) then
C        gaost(idon) = 
C        call areaavg (sbc(1,1,issdon), dmsk, gaost(idon))
        print*, "global average sea surface DON (mol m-3) = "
     &,   gaost(idon)*0.001
      endif
# if !defined O_npzd_no_vflux
      if (idiaz .ne. 0 .and. issdiaz .ne. 0) then
C        gaost(idiaz) = 
C        call areaavg (sbc(1,1,issdiaz), dmsk, gaost(idiaz))
        print*, "global average sea surface diazotrophs (mol m-3) = "
     &,   gaost(idiaz)*0.001
      endif
# endif
      if (idin15 .ne. 0 .and. issdin15 .ne. 0) then
C        gaost(idin15) = 
C        call areaavg (sbc(1,1,issdin15), dmsk, gaost(idin15))
        print*, "global average sea surface nitrate 15 (mol m-3) = "
     &,   gaost(idin15)*0.001
      endif
      if (idon15 .ne. 0 .and. issdon15 .ne. 0) then
C        gaost(idon15) = 
C        call areaavg (sbc(1,1,issdon15), dmsk, gaost(idon15))
        print*, "global average sea surface DON15 (mol m-3) = "
     &,   gaost(idon15)*0.001
      endif
# if !defined O_npzd_no_vflux
      if (iphytn15 .ne. 0 .and. issphytn15 .ne. 0) then
C        gaost(iphytn15) = 
C        call areaavg (sbc(1,1,issphytn15), dmsk, gaost(iphytn15))
        print*, "global average sea surface phytoplankton n15 " 
     &,         "(mol m-3) = ", gaost(iphytn15)*0.001
      endif
      if (izoopn15 .ne. 0 .and. isszoopn15 .ne. 0) then
C        gaost(izoopn15) = 
C        call areaavg (sbc(1,1,isszoopn15), dmsk, gaost(izoopn15))
        print*, "global average sea surface zooplankton n15 (mol m-3) 
     &           = ", gaost(izoopn15)*0.001
      endif
      if (idetrn15 .ne. 0 .and. issdetrn15 .ne. 0) then
C        gaost(idetrn15) = 
C        call areaavg (sbc(1,1,issdetrn15), dmsk, gaost(idetrn15))
        print*, "global average sea surface detritus n15 (mol m-3) 
     &           = ", gaost(idetrn15)*0.001
      endif
      if (idiazn15 .ne. 0 .and. issdiazn15 .ne. 0) then
C        gaost(idiazn15) = 
C        call areaavg (sbc(1,1,issdiazn15), dmsk, gaost(idiazn15))
        print*, "global average sea surface diazotrophs n15 (mol m-3) 
     &           = ", gaost(idiazn15)*0.001
      endif
# endif
      if (idoc13 .ne. 0 .and. issdoc13 .ne. 0) then
C        gaost(idoc13) = 
C        call areaavg (sbc(1,1,issdoc13), dmsk, gaost(idoc13))
        print*, "global average sea surface DOC13"
     &,         " (mol m-3) = ", gaost(idoc13)*0.001
      endif
# if !defined o_npzd_no_vflux
      if (iphytc13 .ne. 0 .and. issphytc13 .ne. 0) then
C        gaost(iphytc13) = 
C        call areaavg (sbc(1,1,issphytc13), dmsk, gaost(iphytc13))
        print*, "global average sea surface phytoplankton C13"
     &,         " (mol m-3) = ", gaost(iphytc13)*0.001
      endif
      if (izoopc13 .ne. 0 .and. isszoopc13 .ne. 0) then
C        gaost(izoopc13) = 
C        call areaavg (sbc(1,1,isszoopc13), dmsk, gaost(izoopc13))
        print*, "global average sea surface zooplankton C13"
     &,         " (mol m-3) = ", gaost(izoopc13)*0.001
      endif
      if (idetrc13 .ne. 0 .and. issdetrc13 .ne. 0) then
C        gaost(idetrc13) = 
C        call areaavg (sbc(1,1,issdetrc13), dmsk, gaost(idetrc13))
        print*, "global average sea surface detritus c13"
     &,         " (mol m-3) = ", gaost(idetrc13)*0.001
      endif
      if (idiazc13 .ne. 0 .and. issdiazc13 .ne. 0) then
C        gaost(idiazc13) =
C        call areaavg (sbc(1,1,issdiazc13), dmsk, gaost(idiazc13))
        print*, "global average sea surface diazotrophs c13"
     &,         " (mol m-3) = ", gaost(idiazc13)*0.001
      endif
# endif
      if (ic14 .ne. 0 .and. issc14 .ne. 0) then
C        gaost(ic14) = 
C        call areaavg (sbc(1,1,issc14), dmsk, gaost(ic14))
        print*, "global average sea surface carbon 14 (mol m-3) = "
     &,   gaost(ic14)
      endif
      endif ! debugFlag

CSPK copy over data/set constants
      dzt(:) = c0
      dzw(:) = c0
      zt(:) = c0
      zw(:) = c0

      do k=1,km
        dzt(k) = drF(k)
        zt(k) = z(k)      
      enddo

      zw(1) = drF(1)      
      do k=2,km
        zw(k) = zw(k-1) + drF(k)
      enddo
      
      do k=1,km-1
        dzw(k) = zt(k+1) - zt(k)
      enddo
      dzw(0) = zt(1)
      dzw(km) = zw(km) - zt(km)

      do k=1,km
        c2dzt(k) = c2*dzt(k)
        dzt2r(k) = c1/c2dzt(k)
      enddo

      dzwr(km)  = c1/dzw(km)
      dzw2r(km) = p5/dzw(km)

      do k=1,km
        dzwr(k-1)    = c1/dzw(k-1)
        dzw2r(k-1)   = p5/dzw(k-1)
      enddo

      do k=1,km
        dztr(k)  = c1/dzt(k)
      enddo

      dtts=DeltaT  ! make sure this matches the value in control.in
      daylen=86400.0
      c2dtts  = dtts
      tmsk(:,:) = 1.0 ! land/ocean mask
C      
      yrlen = 360    ! the length of a typical (non-leap) year in days
      daylen = 86400.0 ! the length of a day in seconds      
      dc14ccn = 0.0  ! atmospheric dC14 (permil)
      dc13ccn    = -6.5 ! atmospheric dC13
CCCCC

CSPK Now initialize the main MOBI code     
      call mobii

      if (dtts.ne.DeltaT) then
         print*,'ERROR: MOBI_INI: dtts does not match DeltaT'
         stop
      endif
      
      if (debugFlag.eq.1) then
        print*,'MOBI time steps: dtts=',dtts,' c2dtts=',c2dtts,
     &         ' dtnpzd=',dtnpzd
#ifdef O_npzd_fe_limitation
        print*,'Iron limitation depth index kmfe=',kmfe
#endif
      endif     

      END
