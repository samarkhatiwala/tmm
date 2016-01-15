#include "MOBI_TMM_OPTIONS.h"
      SUBROUTINE MOBI_INI(nzmax, z, drF, DeltaT,
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

      integer nzmax
      real z(km), drF(km), DeltaT
      real S_surf_glob, TR_surf_glob(nsrc)
      integer debugFlag
      
      integer m, n, k, itr
      real dmsk(imt,jmt)  ! SPK: setmom.F
      
! set some constants
      pi = 4.0 * atan(1.0)
      radian = 360./(2.*pi)
      rho0 = 1.035
      rho0r = c1/rho0
      socn = 0.03475
Csetmom.F:        socn = gaost(isalt)


	  if (nzmax.ne.km) then
		print*,'ERROR in MOBI_INI: number of depth levels incorrect!'
		stop
	  endif  

      call tracer_init

      if (debugFlag.eq.1) then
		do n=1,nt
		  print*,'n=',n,mapt(n)
		enddo
      endif

      call sbc_init

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
        print*,'ERROR in MOBI_INI: number of tracers incorrect!'
        stop
      endif  
      
      if (debugFlag.eq.1) then
      if (idic .ne. 0 .and. issdic .ne. 0) then
        call areaavg (sbc(1,1,issdic), dmsk, gaost(idic))
        print*, "global average sea surface dic (mol m-3) = "
     &,   gaost(idic)
      endif
      if (idic13 .ne. 0 .and. issdic13 .ne. 0) then
        call areaavg (sbc(1,1,issdic13), dmsk, gaost(idic13))
        print*, "global average sea surface dic 13 (mol m-3) = "
     &,   gaost(idic13)
      endif
      if (io2 .ne. 0 .and. isso2 .ne. 0) then
        call areaavg (sbc(1,1,isso2), dmsk, gaost(io2))
        print*, "global average sea surface oxygen (mol m-3) = "
     &,   gaost(io2)
      endif
      if (ialk .ne. 0 .and. issalk .ne. 0) then
        call areaavg (sbc(1,1,issalk), dmsk, gaost(ialk))
        print*, "global average sea surface alkalinity (mol m-3) = "
     &,   gaost(ialk)
      endif
      if (ipo4 .ne. 0 .and. isspo4 .ne. 0) then
        call areaavg (sbc(1,1,isspo4), dmsk, gaost(ipo4))
        print*, "global average sea surface phosphate (mol m-3) = "
     &,   gaost(ipo4)*0.001
      endif
      if (idop .ne. 0 .and. issdop .ne. 0) then
        call areaavg (sbc(1,1,issdop), dmsk, gaost(idop))
        print*, "global average sea surface DOP (mol m-3) = "
     &,   gaost(idop)*0.001
      endif
# if !defined O_npzd_no_vflux
      if (iphyt .ne. 0 .and. issphyt .ne. 0) then
        call areaavg (sbc(1,1,issphyt), dmsk, gaost(iphyt))
        print*, "global average sea surface phytoplankton (mol m-3) = "
     &,   gaost(iphyt)*0.001
      endif
      if (izoop .ne. 0 .and. isszoop .ne. 0) then
        call areaavg (sbc(1,1,isszoop), dmsk, gaost(izoop))
        print*, "global average sea surface zooplankton (mol m-3) = "
     &,   gaost(izoop)*0.001
      endif
      if (idetr .ne. 0 .and. issdetr .ne. 0) then
        call areaavg (sbc(1,1,issdetr), dmsk, gaost(idetr))
        print*, "global average sea surface detritus (mol m-3) = "
     &,   gaost(idetr)*0.001
      endif
# endif
      if (ino3 .ne. 0 .and. issno3 .ne. 0) then
        call areaavg (sbc(1,1,issno3), dmsk, gaost(ino3))
        print*, "global average sea surface nitrate (mol m-3) = "
     &,   gaost(ino3)*0.001
      endif
      if (idon .ne. 0 .and. issdon .ne. 0) then
        call areaavg (sbc(1,1,issdon), dmsk, gaost(idon))
        print*, "global average sea surface DON (mol m-3) = "
     &,   gaost(idon)*0.001
      endif
# if !defined O_npzd_no_vflux
      if (idiaz .ne. 0 .and. issdiaz .ne. 0) then
        call areaavg (sbc(1,1,issdiaz), dmsk, gaost(idiaz))
        print*, "global average sea surface diazotrophs (mol m-3) = "
     &,   gaost(idiaz)*0.001
      endif
# endif
      if (idin15 .ne. 0 .and. issdin15 .ne. 0) then
        call areaavg (sbc(1,1,issdin15), dmsk, gaost(idin15))
        print*, "global average sea surface nitrate 15 (mol m-3) = "
     &,   gaost(idin15)*0.001
      endif
      if (idon15 .ne. 0 .and. issdon15 .ne. 0) then
        call areaavg (sbc(1,1,issdon15), dmsk, gaost(idon15))
        print*, "global average sea surface DON15 (mol m-3) = "
     &,   gaost(idon15)*0.001
      endif
# if !defined O_npzd_no_vflux
      if (iphytn15 .ne. 0 .and. issphytn15 .ne. 0) then
        call areaavg (sbc(1,1,issphytn15), dmsk, gaost(iphytn15))
        print*, "global average sea surface phytoplankton n15 " 
     &,         "(mol m-3) = ", gaost(iphytn15)*0.001
      endif
      if (izoopn15 .ne. 0 .and. isszoopn15 .ne. 0) then
        call areaavg (sbc(1,1,isszoopn15), dmsk, gaost(izoopn15))
        print*, "global average sea surface zooplankton n15 (mol m-3) 
     &           = ", gaost(izoopn15)*0.001
      endif
      if (idetrn15 .ne. 0 .and. issdetrn15 .ne. 0) then
        call areaavg (sbc(1,1,issdetrn15), dmsk, gaost(idetrn15))
        print*, "global average sea surface detritus n15 (mol m-3) 
     &           = ", gaost(idetrn15)*0.001
      endif
      if (idiazn15 .ne. 0 .and. issdiazn15 .ne. 0) then
        call areaavg (sbc(1,1,issdiazn15), dmsk, gaost(idiazn15))
        print*, "global average sea surface diazotrophs n15 (mol m-3) 
     &           = ", gaost(idiazn15)*0.001
      endif
# endif
      if (idoc13 .ne. 0 .and. issdoc13 .ne. 0) then
        call areaavg (sbc(1,1,issdoc13), dmsk, gaost(idoc13))
        print*, "global average sea surface DOC13"
     &,         " (mol m-3) = ", gaost(idoc13)*0.001
      endif
# if !defined o_npzd_no_vflux
      if (iphytc13 .ne. 0 .and. issphytc13 .ne. 0) then
        call areaavg (sbc(1,1,issphytc13), dmsk, gaost(iphytc13))
        print*, "global average sea surface phytoplankton C13"
     &,         " (mol m-3) = ", gaost(iphytc13)*0.001
      endif
      if (izoopc13 .ne. 0 .and. isszoopc13 .ne. 0) then
        call areaavg (sbc(1,1,isszoopc13), dmsk, gaost(izoopc13))
        print*, "global average sea surface zooplankton C13"
     &,         " (mol m-3) = ", gaost(izoopc13)*0.001
      endif
      if (idetrc13 .ne. 0 .and. issdetrc13 .ne. 0) then
        call areaavg (sbc(1,1,issdetrc13), dmsk, gaost(idetrc13))
        print*, "global average sea surface detritus c13"
     &,         " (mol m-3) = ", gaost(idetrc13)*0.001
      endif
      if (idiazc13 .ne. 0 .and. issdiazc13 .ne. 0) then
        call areaavg (sbc(1,1,issdiazc13), dmsk, gaost(idiazc13))
        print*, "global average sea surface diazotrophs c13"
     &,         " (mol m-3) = ", gaost(idiazc13)*0.001
      endif
# endif
      if (ic14 .ne. 0 .and. issc14 .ne. 0) then
        call areaavg (sbc(1,1,issc14), dmsk, gaost(ic14))
        print*, "global average sea surface carbon 14 (mol m-3) = "
     &,   gaost(ic14)
      endif
      if (icfc11 .ne. 0 .and. isscfc11 .ne. 0) then
        call areaavg (sbc(1,1,isscfc11), dmsk, gaost(icfc11))
        print*, "global average sea surface cfc 11 (mol m-3) = "
     &,   gaost(icfc11)
      endif
      if (icfc12 .ne. 0 .and. isscfc12 .ne. 0) then
        call areaavg (sbc(1,1,isscfc12), dmsk, gaost(icfc12))
        print*, "global average sea surface cfc 12 (mol m-3) = "
     &,   gaost(icfc12)
      endif
      print*, " "
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
