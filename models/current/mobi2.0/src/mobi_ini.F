
      SUBROUTINE MOBI_INI(numTracers, lSize, nx_loc, nzmax, nz_loc_arr,
     &                    z, drF, DeltaT, tlat_loc,
     &                    sg_bathy_loc,
     &                    S_surf_glob, TR_surf_glob, 
     &                    debugFlag)

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
#include "accel.h"

      integer numTracers, lSize, nx_loc, nzmax, nz_loc_arr(nx_loc)
      real z(km), drF(km), DeltaT, tlat_loc(nx_loc)
      real sg_bathy_loc(lSize)
      real S_surf_glob, TR_surf_glob(nsrc)
      integer debugFlag
      
      integer m, n, i, j, k, ii, itr
      real dmsk(imt,jmt)  ! SPK: setmom.F

#ifndef O_co2ccn_user
      THIS LINE IS HERE TO GIVE AN ERROR WHEN COMPILING.
      And IN CASE THAT FAILS
	  print*,'ERROR: you must define O_co2ccn_user in MOBI_TMM_OPTIONS.h'
	  stop
#endif
      
! set some constants
      pi = 4.0 * atan(1.0)
      radian = 360./(2.*pi)
      rho0 = 1.035
      rho0r = c1/rho0
      socn = 0.03475
Csetmom.F:        socn = gaost(isalt)

	  if (nx_loc.gt.imt) then
		print*,'ERROR in MOBI_INI: number of local profiles exceeds imt!'
		stop
	  endif  

	  if (nzmax.ne.km) then
		print*,'ERROR in MOBI_INI: number of depth levels incorrect!'
		stop
	  endif  

      j=1
      kmt(:,:) = 0
	  kmt(1:nx_loc,j) = nz_loc_arr(1:nx_loc)
	  tlat(1:nx_loc,j) = tlat_loc(1:nx_loc)

      ii=1
      do i=1,nx_loc
        do k=1,kmt(i,j)
	      sg_bathy(i,j,k) = sg_bathy_loc(ii)
	      ii=ii+1
	    enddo  
	  enddo

!SPK     Initialization sequence based on UVic_ESCM.F

!-----------------------------------------------------------------------
!     Initialize tracers
!-----------------------------------------------------------------------

      call tracer_init

!-----------------------------------------------------------------------
!     Initialize S.B.C. indices
!-----------------------------------------------------------------------

      call sbc_init

      if (numTracers.ne.(nt-2)) then
        write(*,*),'ERROR in MOBI_INI: number of tracers incorrect!'
        write(*,*),'                   numTracers not equal to nt-2'
        stop
      endif  
      
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

!-----------------------------------------------------------------------
!     common setup
!-----------------------------------------------------------------------

C      call setcom (1, imt, 1, jmt)

!-----------------------------------------------------------------------
!     ocean setup
!-----------------------------------------------------------------------
      do k=1,km
        dtxcel(k)=1.0
      enddo

CSPK we set these to all point to the same place
      tau   = 0
      taum1 = 0
      taup1 = 0
C
      
CSPK  Copy initial surface global mean values
CSPK  Note that gaost(isalt) is not subsequently changed
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
        print*,'ERROR in MOBI_INI: number of tracers incorrect!'
        stop
      endif  
C

      call setmom (1, nx_loc, 1, jmt)
C
CC     if (debugFlag.eq.1) then
CCSPK taken from setmom.F, S/R rowi
C!-----------------------------------------------------------------------
C!     find inital average surface references
C!-----------------------------------------------------------------------
C     print*, " "
C     print*, "inital average surface references: "
C     dmsk(:,:) = 1.
C     where (kmt(:,:) .eq. 0) dmsk(:,:) = 0.
CC     gaost(:) = 0.
C     if (isalt .ne. 0 .and. isss .ne. 0) then
C       call areaavg (sbc(1,1,isss), dmsk, gaost(isalt))
C       gaost(isalt) = gaost(isalt) + 0.035
C       socn = gaost(isalt)
C       print*, "global average sea surface salinity (psu) = "
C    &,   gaost(isalt)*1000.
C     endif
C     if (idic .ne. 0 .and. issdic .ne. 0) then
C       call areaavg (sbc(1,1,issdic), dmsk, gaost(idic))
C       print*, "global average sea surface dic (mol m-3) = "
C    &,   gaost(idic)
C     endif
C     if (idic13 .ne. 0 .and. issdic13 .ne. 0) then
C       call areaavg (sbc(1,1,issdic13), dmsk, gaost(idic13))
C       print*, "global average sea surface dic 13 (mol m-3) = "
C    &,   gaost(idic13)
C     endif
C     if (io2 .ne. 0 .and. isso2 .ne. 0) then
C       call areaavg (sbc(1,1,isso2), dmsk, gaost(io2))
C       print*, "global average sea surface oxygen (mol m-3) = "
C    &,   gaost(io2)
C     endif
C     if (ialk .ne. 0 .and. issalk .ne. 0) then
C       call areaavg (sbc(1,1,issalk), dmsk, gaost(ialk))
C       print*, "global average sea surface alkalinity (mol m-3) = "
C    &,   gaost(ialk)
C     endif
C     if (ipo4 .ne. 0 .and. isspo4 .ne. 0) then
C       call areaavg (sbc(1,1,isspo4), dmsk, gaost(ipo4))
C       print*, "global average sea surface phosphate (mol m-3) = "
C    &,   gaost(ipo4)*0.001
C     endif
C     if (idop .ne. 0 .and. issdop .ne. 0) then
C       call areaavg (sbc(1,1,issdop), dmsk, gaost(idop))
C       print*, "global average sea surface DOP (mol m-3) = "
C    &,   gaost(idop)*0.001
C     endif
C# if !defined O_npzd_no_vflux
C     if (iphyt .ne. 0 .and. issphyt .ne. 0) then
C       call areaavg (sbc(1,1,issphyt), dmsk, gaost(iphyt))
C       print*, "global average sea surface phytoplankton (mol m-3) = "
C    &,   gaost(iphyt)*0.001
C     endif
C#  if defined O_npzd_caco3
C     if (icocc .ne. 0 .and. isscocc .ne. 0) then
C       call areaavg (sbc(1,1,isscocc), dmsk, gaost(icocc))
C       print*, "global average sea surface coccs (mol m-3) = "
C    &,   gaost(icocc)*0.001
C     endif
C     if (icaco3 .ne. 0 .and. isscaco3 .ne. 0) then
C       call areaavg (sbc(1,1,isscaco3), dmsk, gaost(icaco3))
C       print*, "global average sea surface calcite (mol C m-3) = "
C    &,   gaost(icaco3)*0.001
C     endif
C#  endif
C#  if defined O_kk_ballast
C     if (idetr_B .ne. 0 .and. issdetr_B .ne. 0) then
C       call areaavg (sbc(1,1,issdetr_B), dmsk, gaost(idetr_B))
C       print*, "globl avg sea surf ballasted detritus (mol m-3) = "
C    &,   gaost(idetr_B)*0.001
C     endif
C#  endif
C     if (izoop .ne. 0 .and. isszoop .ne. 0) then
C       call areaavg (sbc(1,1,isszoop), dmsk, gaost(izoop))
C       print*, "global average sea surface zooplankton (mol m-3) = "
C    &,   gaost(izoop)*0.001
C     endif
C     if (idetr .ne. 0 .and. issdetr .ne. 0) then
C       call areaavg (sbc(1,1,issdetr), dmsk, gaost(idetr))
C       print*, "global average sea surface detritus (mol m-3) = "
C    &,   gaost(idetr)*0.001
C     endif
C# endif
C     if (ino3 .ne. 0 .and. issno3 .ne. 0) then
C       call areaavg (sbc(1,1,issno3), dmsk, gaost(ino3))
C       print*, "global average sea surface nitrate (mol m-3) = "
C    &,   gaost(ino3)*0.001
C     endif
C     if (idon .ne. 0 .and. issdon .ne. 0) then
C       call areaavg (sbc(1,1,issdon), dmsk, gaost(idon))
C       print*, "global average sea surface DON (mol m-3) = "
C    &,   gaost(idon)*0.001
C     endif
C# if !defined O_npzd_no_vflux
C     if (idiaz .ne. 0 .and. issdiaz .ne. 0) then
C       call areaavg (sbc(1,1,issdiaz), dmsk, gaost(idiaz))
C       print*, "global average sea surface diazotrophs (mol m-3) = "
C    &,   gaost(idiaz)*0.001
C     endif
C# endif
C     if (idin15 .ne. 0 .and. issdin15 .ne. 0) then
C       call areaavg (sbc(1,1,issdin15), dmsk, gaost(idin15))
C       print*, "global average sea surface nitrate 15 (mol m-3) = "
C    &,   gaost(idin15)*0.001
C     endif
C     if (idon15 .ne. 0 .and. issdon15 .ne. 0) then
C       call areaavg (sbc(1,1,issdon15), dmsk, gaost(idon15))
C       print*, "global average sea surface DON15 (mol m-3) = "
C    &,   gaost(idon15)*0.001
C     endif
C# if !defined O_npzd_no_vflux
C     if (iphytn15 .ne. 0 .and. issphytn15 .ne. 0) then
C       call areaavg (sbc(1,1,issphytn15), dmsk, gaost(iphytn15))
C       print*, "global average sea surface phytoplankton n15 " 
C    &,         "(mol m-3) = ", gaost(iphytn15)*0.001
C     endif
C#  if defined O_npzd_caco3
C     if (icoccn15 .ne. 0 .and. isscoccn15 .ne. 0) then
C       call areaavg (sbc(1,1,isscoccn15), dmsk, gaost(icoccn15))
C       print*, "global average sea surface calcifiers n15 " 
C    &,         "(mol m-3) = ", gaost(icoccn15)*0.001
C     endif
C#  endif      
C     if (izoopn15 .ne. 0 .and. isszoopn15 .ne. 0) then
C       call areaavg (sbc(1,1,isszoopn15), dmsk, gaost(izoopn15))
C       print*, "global average sea surface zooplankton n15 (mol m-3) 
C    &           = ", gaost(izoopn15)*0.001
C     endif
C     if (idetrn15 .ne. 0 .and. issdetrn15 .ne. 0) then
C       call areaavg (sbc(1,1,issdetrn15), dmsk, gaost(idetrn15))
C       print*, "global average sea surface detritus n15 (mol m-3) 
C    &           = ", gaost(idetrn15)*0.001
C     endif
C     if (idiazn15 .ne. 0 .and. issdiazn15 .ne. 0) then
C       call areaavg (sbc(1,1,issdiazn15), dmsk, gaost(idiazn15))
C       print*, "global average sea surface diazotrophs n15 (mol m-3) 
C    &           = ", gaost(idiazn15)*0.001
C     endif
C# endif
C     if (idoc13 .ne. 0 .and. issdoc13 .ne. 0) then
C       call areaavg (sbc(1,1,issdoc13), dmsk, gaost(idoc13))
C       print*, "global average sea surface DOC13"
C    &,         " (mol m-3) = ", gaost(idoc13)*0.001
C     endif
C# if !defined O_npzd_no_vflux
C     if (iphytc13 .ne. 0 .and. issphytc13 .ne. 0) then
C       call areaavg (sbc(1,1,issphytc13), dmsk, gaost(iphytc13))
C       print*, "global average sea surface phytoplankton C13"
C    &,         " (mol m-3) = ", gaost(iphytc13)*0.001
C     endif
C#  if defined O_npzd_caco3
C     if (icoccc13 .ne. 0 .and. isscoccc13 .ne. 0) then
C       call areaavg (sbc(1,1,isscoccc13), dmsk, gaost(icoccc13))
C       print*, "global average sea surface calcifiers C13"
C    &,         " (mol m-3) = ", gaost(icoccc13)*0.001
C     endif
C     if (icaco3c13 .ne. 0 .and. isscaco3c13 .ne. 0) then
C       call areaavg (sbc(1,1,isscaco3c13), dmsk, gaost(icaco3c13))
C       print*, "global average sea surface CaCO3 C13"
C    &,         " (mol m-3) = ", gaost(icaco3c13)*0.001
C     endif
C#  endif      
C     if (izoopc13 .ne. 0 .and. isszoopc13 .ne. 0) then
C       call areaavg (sbc(1,1,isszoopc13), dmsk, gaost(izoopc13))
C       print*, "global average sea surface zooplankton C13"
C    &,         " (mol m-3) = ", gaost(izoopc13)*0.001
C     endif
C     if (idetrc13 .ne. 0 .and. issdetrc13 .ne. 0) then
C       call areaavg (sbc(1,1,issdetrc13), dmsk, gaost(idetrc13))
C       print*, "global average sea surface detritus c13"
C    &,         " (mol m-3) = ", gaost(idetrc13)*0.001
C     endif
C     if (idiazc13 .ne. 0 .and. issdiazc13 .ne. 0) then
C       call areaavg (sbc(1,1,issdiazc13), dmsk, gaost(idiazc13))
C       print*, "global average sea surface diazotrophs c13"
C    &,         " (mol m-3) = ", gaost(idiazc13)*0.001
C     endif
C# endif
C     if (idfe .ne. 0 .and. issdfe .ne. 0) then
C       call areaavg (sbc(1,1,issdfe), dmsk, gaost(idfe))
C       print*, "global average sea surface iron (mol m-3) = "
C    &,   gaost(idfe)*0.001
C     endif
C# if !defined O_npzd_no_vflux
C     if (idetrfe .ne. 0 .and. issdetrfe .ne. 0) then
C       call areaavg (sbc(1,1,issdetrfe), dmsk, gaost(idetrfe))
C       print*, "global average sea surface particulate iron (mol m-3)="
C    &,   gaost(idetrfe)*0.001
C     endif
C# endif
C     if (ic14 .ne. 0 .and. issc14 .ne. 0) then
C       call areaavg (sbc(1,1,issc14), dmsk, gaost(ic14))
C       print*, "global average sea surface carbon 14 (mol m-3) = "
C    &,   gaost(ic14)
C     endif
C     if (icfc11 .ne. 0 .and. isscfc11 .ne. 0) then
C       call areaavg (sbc(1,1,isscfc11), dmsk, gaost(icfc11))
C       print*, "global average sea surface cfc 11 (mol m-3) = "
C    &,   gaost(icfc11)
C     endif
C     if (icfc12 .ne. 0 .and. isscfc12 .ne. 0) then
C       call areaavg (sbc(1,1,isscfc12), dmsk, gaost(icfc12))
C       print*, "global average sea surface cfc 12 (mol m-3) = "
C    &,   gaost(icfc12)
C     endif
C     print*, " "
CCSPK
CC     endif ! debugFlag
C
CCSPK from setmom.F, S/R setmom
C     call mobi_init
CC

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
      
CSPK  Check
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
