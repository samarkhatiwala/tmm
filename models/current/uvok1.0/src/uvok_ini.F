#ifdef O_TMM
      SUBROUTINE UVOK_INI(nzmax, z, drF, DeltaT,
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
#include "levind.h"
#include "switch.h"

      integer nzmax
      real z(km), drF(km), DeltaT
      real S_surf_glob, TR_surf_glob(nsrc)
      integer debugFlag
      
      integer m, n, k, itr
      real dmsk(imt,jmt)  ! SPK: setmom.F

#ifndef O_co2ccn_user
      THIS LINE IS HERE TO GIVE AN ERROR WHEN COMPILING.
      And IN CASE THAT FAILS
	  print*,'ERROR: you must define O_co2ccn_user in UVOK_TMM_OPTIONS.h'
	  stop
#endif
      
! set some constants
      pi = 4.0 * atan(1.0)
      radian = 360./(2.*pi)
      rho0 = 1.035
      rho0r = c1/rho0
      socn = 0.03475
Csetmom.F:        socn = gaost(isalt)

	  if (nzmax.ne.km) then
		print*,'ERROR in UVOK_INI: number of depth levels incorrect!'
		stop
	  endif  

      kmt(1,1) = nzmax ! this will be set to the local value in uvok_calc

!SPK     Initialization sequence based on UVic_ESCM.F

!-----------------------------------------------------------------------
!     Initialize S.B.C. indices
!-----------------------------------------------------------------------

      call sbc_init

!-----------------------------------------------------------------------
!     Initialize tracers
!-----------------------------------------------------------------------

      call tracer_init

CSPK  check      
      if (debugFlag.eq.1) then
		do n=1,nt
		  print*,'n=',n,mapt(n)
		enddo
      endif
C

!-----------------------------------------------------------------------
!     read namelist variables
!-----------------------------------------------------------------------

      call read_namelist

CSPK  overwrite some values
      dtts=DeltaT  ! make sure this matches the value in control.in
      daylen=86400.0
      c2dtts  = dtts
C

!-----------------------------------------------------------------------
!     read grid
!-----------------------------------------------------------------------

C     call grids

!-----------------------------------------------------------------------
!     read topography
!-----------------------------------------------------------------------

C      call topog (kmt, kmu, map, xt, yt, zt, xu, yu, zw, imt, jmt, km
C     &,           sg_bathy)

C      call isleperim (kmt, map, iperm, jperm, iofs, nippts, nisle, imt
C     &,               jmt, km, mnisle, maxipp, xu, yu, zw)


CSPK copy over/set grid variables (based on grids.F)
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

CKK from S/R setmom.F
#if defined O_npzd
!     taken from setmom.F
!     convert units of NPZD parameters to MOM units
      redctn = redctn*1.e-3
      redotn = redotn*1.e-3
      redotp = redotn/redptn
      redctp = redctn/redptn
      redntp = 1./redptn
      k1p   = k1n*redptn
      kw = kw*1.e-2
      kc = kc*1.e-2
      ki = ki*1.e-2
      wd0 = wd0*1.e2
      alpha = alpha/daylen
      abio = abio/daylen
      nup = nup/daylen
      nupt0 = nupt0/daylen
      gbio = gbio/daylen
      nuz = nuz/daylen
      nud0 = nud0/daylen
#endif

      do k=1,km
!     linear increase wd0-200m with depth
         wd(k) = (wd0+6.0e-2*zt(k))/daylen/dzt(k) ! [s-1]
         rkwz(k) = 1./(kw*dzt(k))
      enddo
      ztt(1)=0.0
      do k=1,km-1
         ztt(k+1)=(-1)*zw(k)
      enddo
      ztt(km)=(-1)*zw(km-1)
# if defined O_carbon || defined O_npzd_alk

!---------------------------------------------------------------------
!     calculate variables used in calcite remineralization
!---------------------------------------------------------------------

      rcak(1) = -(exp(-zw(1)/dcaco3)-1.0)/dzt(1)
      rcab(1) = 1./dzt(1)
      do k=2,km
        rcak(k) = -(exp(-zw(k)/dcaco3))/dzt(k)
     &          + (exp(-zw(k-1)/dcaco3))/dzt(k)
        rcab(k) = (exp(-zw(k-1)/dcaco3))/dzt(k)
      enddo
# endif
!-----------------------------------------------------------------------
!     common setup
!-----------------------------------------------------------------------

C      call setcom (1, imt, 1, jmt)

!-----------------------------------------------------------------------
!     ocean setup
!-----------------------------------------------------------------------

C      call setmom (1, imt, 1, jmt)

CSPK we set these to all point to the same place
      tau   = 0
      taum1 = 0
      taup1 = 0
C
      
CSPK copy initial surface values
      gaost(:) = 0.
      if (isalt .ne. 0 .and. isss .ne. 0) then
        gaost(isalt) = S_surf_glob
C        write(*,*)'Salt: ',gaost(isalt)
      endif
      itr = 0
      do n=3,nt
        if (itrc(n) .ne. 0) then
         itr = itr + 1
		 gaost(n) = TR_surf_glob(itr)	
C		 write(*,*)mapt(n),gaost(n)
        endif
      enddo
C     check
      if (itr.ne.nsrc) then
        print*,'ERROR in UVOK_INI: number of tracers incorrect!'
        stop
      endif  
C
      
      if (debugFlag.eq.1) then
CSPK taken from setmom.F, S/R rowi
!-----------------------------------------------------------------------
!     find inital average surface references
!-----------------------------------------------------------------------
      print*, " "
      print*, "inital average surface references: "
      dmsk(:,:) = 1.
      where (kmt(:,:) .eq. 0) dmsk(:,:) = 0.
C     gaost(:) = 0.
      if (isalt .ne. 0 .and. isss .ne. 0) then
        call areaavg (sbc(1,1,isss), dmsk, gaost(isalt))
        gaost(isalt) = gaost(isalt) + 0.035
        socn = gaost(isalt)
        print*, "global average sea surface salinity (psu) = "
     &,   gaost(isalt)*1000.
      endif
      if (idic .ne. 0 .and. issdic .ne. 0) then
        call areaavg (sbc(1,1,issdic), dmsk, gaost(idic))
        print*, "global average sea surface dic (mol m-3) = "
     &,   gaost(idic)
      endif
C      if (idic13 .ne. 0 .and. issdic13 .ne. 0) then
C        call areaavg (sbc(1,1,issdic13), dmsk, gaost(idic13))
C        print*, "global average sea surface dic 13 (mol m-3) = "
C     &,   gaost(idic13)
C      endif
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
C      if (idop .ne. 0 .and. issdop .ne. 0) then
C        call areaavg (sbc(1,1,issdop), dmsk, gaost(idop))
C        print*, "global average sea surface DOP (mol m-3) = "
C     &,   gaost(idop)*0.001
C      endif
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
C      if (idon .ne. 0 .and. issdon .ne. 0) then
C        call areaavg (sbc(1,1,issdon), dmsk, gaost(idon))
C        print*, "global average sea surface DON (mol m-3) = "
C     &,   gaost(idon)*0.001
C      endif
# if !defined O_npzd_no_vflux
      if (idiaz .ne. 0 .and. issdiaz .ne. 0) then
        call areaavg (sbc(1,1,issdiaz), dmsk, gaost(idiaz))
        print*, "global average sea surface diazotrophs (mol m-3) = "
     &,   gaost(idiaz)*0.001
      endif
# endif
C      if (idin15 .ne. 0 .and. issdin15 .ne. 0) then
C        call areaavg (sbc(1,1,issdin15), dmsk, gaost(idin15))
C        print*, "global average sea surface nitrate 15 (mol m-3) = "
C     &,   gaost(idin15)*0.001
C      endif
C      if (idon15 .ne. 0 .and. issdon15 .ne. 0) then
C        call areaavg (sbc(1,1,issdon15), dmsk, gaost(idon15))
C        print*, "global average sea surface DON15 (mol m-3) = "
C     &,   gaost(idon15)*0.001
C      endif
# if !defined O_npzd_no_vflux
C      if (iphytn15 .ne. 0 .and. issphytn15 .ne. 0) then
C        call areaavg (sbc(1,1,issphytn15), dmsk, gaost(iphytn15))
C        print*, "global average sea surface phytoplankton n15 " 
C     &,         "(mol m-3) = ", gaost(iphytn15)*0.001
C      endif
C      if (izoopn15 .ne. 0 .and. isszoopn15 .ne. 0) then
C        call areaavg (sbc(1,1,isszoopn15), dmsk, gaost(izoopn15))
C        print*, "global average sea surface zooplankton n15 (mol m-3) 
C     &           = ", gaost(izoopn15)*0.001
C      endif
C      if (idetrn15 .ne. 0 .and. issdetrn15 .ne. 0) then
C        call areaavg (sbc(1,1,issdetrn15), dmsk, gaost(idetrn15))
C        print*, "global average sea surface detritus n15 (mol m-3) 
C     &           = ", gaost(idetrn15)*0.001
C      endif
C      if (idiazn15 .ne. 0 .and. issdiazn15 .ne. 0) then
C        call areaavg (sbc(1,1,issdiazn15), dmsk, gaost(idiazn15))
C        print*, "global average sea surface diazotrophs n15 (mol m-3) 
C     &           = ", gaost(idiazn15)*0.001
C      endif
# endif
C      if (idoc13 .ne. 0 .and. issdoc13 .ne. 0) then
C        call areaavg (sbc(1,1,issdoc13), dmsk, gaost(idoc13))
C        print*, "global average sea surface DOC13"
C     &,         " (mol m-3) = ", gaost(idoc13)*0.001
C      endif
# if !defined O_npzd_no_vflux
C      if (iphytc13 .ne. 0 .and. issphytc13 .ne. 0) then
C        call areaavg (sbc(1,1,issphytc13), dmsk, gaost(iphytc13))
C        print*, "global average sea surface phytoplankton C13"
C     &,         " (mol m-3) = ", gaost(iphytc13)*0.001
C      endif
C      if (izoopc13 .ne. 0 .and. isszoopc13 .ne. 0) then
C        call areaavg (sbc(1,1,isszoopc13), dmsk, gaost(izoopc13))
C        print*, "global average sea surface zooplankton C13"
C     &,         " (mol m-3) = ", gaost(izoopc13)*0.001
C      endif
C      if (idetrc13 .ne. 0 .and. issdetrc13 .ne. 0) then
C        call areaavg (sbc(1,1,issdetrc13), dmsk, gaost(idetrc13))
C        print*, "global average sea surface detritus c13"
C     &,         " (mol m-3) = ", gaost(idetrc13)*0.001
C      endif
C      if (idiazc13 .ne. 0 .and. issdiazc13 .ne. 0) then
C        call areaavg (sbc(1,1,issdiazc13), dmsk, gaost(idiazc13))
C        print*, "global average sea surface diazotrophs c13"
C     &,         " (mol m-3) = ", gaost(idiazc13)*0.001
C      endif
# endif
C      if (idfe .ne. 0 .and. issdfe .ne. 0) then
C        call areaavg (sbc(1,1,issdfe), dmsk, gaost(idfe))
C        print*, "global average sea surface iron (mol m-3) = "
C     &,   gaost(idfe)*0.001
C      endif
# if !defined O_npzd_no_vflux
C      if (idetrfe .ne. 0 .and. issdetrfe .ne. 0) then
C        call areaavg (sbc(1,1,issdetrfe), dmsk, gaost(idetrfe))
C        print*, "global average sea surface particulate iron (mol m-3)="
C     &,   gaost(idetrfe)*0.001
C      endif
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
CSPK
      endif ! debugFlag

CSPK from setmom.F, S/R setmom
C

!-----------------------------------------------------------------------
!     atmosphere setup
!-----------------------------------------------------------------------
C      call setembm (1, imt, 1, jmt)

CSPK based on setembm.F
      tmsk(:,:) = 1.0 ! land/ocean mask
C

CSPK this is generally set in S/R loadmw called from S/R mom
      tmask(:,:,:) = c1
C

CSPK set some flags so that UVic doesn't start accumulating time averages
      euler2 = .false.
      timavgperts = .false.
      tsiperts = .false.

# if defined O_npzd_fe_limitation
      do k=km,1,-1
         if (zw(k) .gt. 24000) kmfe = k-1
      enddo
# endif
      
CSPK  Check
      if (dtts.ne.DeltaT) then
         print*,'ERROR: UVOK_INI: dtts does not match DeltaT'
         stop
      endif
      
      if (debugFlag.eq.1) then
        print*,'UVOK time steps: dtts=',dtts,' c2dtts=',c2dtts,
     &         ' dtnpzd=',dtnpzd
#ifdef O_npzd_fe_limitation
        print*,'Iron limitation depth index kmfe=',kmfe
#endif
      endif     

      END
#endif