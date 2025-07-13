#include "CPP_OPTIONS_CHECK.h"

      SUBROUTINE MOBI_INI(numTracers, lSize, nx_loc, kmmax, nz_loc_arr,
     &                    z, drF, DeltaT, tlat_loc, area_loc,
     &                    S_surf_glob, TR_surf_glob,
#if defined O_carbon
     &                    co2atm, 
# if defined O_carbon_13
     &                    dc13atm, 
# endif     
# if defined O_carbon_14     
     &                    dc14atm,      
# endif
#endif ! O_carbon
#if defined O_mobi
     &                    sg_bathy_loc,
# if defined O_mobi_iron
     &                    fe_hydr_loc, 
# endif
# if defined O_mobi_silicon
     &                    si_hydr_loc, 
# endif     
#endif
#if defined O_PaTh 
     &                    lith_loc, wLith_loc, 
# if !defined O_mobi
     &                    pom_loc, caco3_loc, opal_loc, 
     &                    wPOM_loc, wCaCO3_loc, wOpal_loc, 
# endif     
#endif
#if defined O_sed
     &                    numOceanStepsPerSedStep,
     &                    nzmaxSed, ibmaxSed, numSedMixedTracers, 
     &                    numSedBuriedTracers, globalweathflx,
     &                    sedsa_loc, sedmask_loc,
#endif ! O_sed   
     &                    debugFlag)

      implicit none
#include "size.h"
#if defined O_mobi
#include "mobi.h"
#endif
#if defined O_PaTh
#include "protac_thor.h"
#endif
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
#if defined O_sed
#include "sed.h"
#endif

      integer numTracers, lSize, nx_loc, kmmax, nz_loc_arr(nx_loc)
      real z(km), drF(km), DeltaT, tlat_loc(nx_loc), area_loc(nx_loc)
      real S_surf_glob, TR_surf_glob(nsrc)
#if defined O_carbon
      real co2atm   
# if defined O_carbon_13
      real dc13atm
# endif     
# if defined O_carbon_14     
      real dc14atm
# endif
#endif ! O_carbon
#if defined O_mobi
      real sg_bathy_loc(lSize)
# if defined O_mobi_iron
      real fe_hydr_loc(lSize)
# endif
# if defined O_mobi_silicon
      real si_hydr_loc(lSize)
# endif      
#endif
#if defined O_PaTh
      real lith_loc(lSize), wLith_loc(km)
# if !defined O_mobi
      real pom_loc(lSize), caco3_loc(lSize), opal_loc(lSize)
      real wPOM_loc(km), wCaCO3_loc(km), wOpal_loc(km)
# endif      
#endif
#if defined O_sed
      integer numOceanStepsPerSedStep
      integer nzmaxSed, ibmaxSed, numSedMixedTracers
      integer numSedBuriedTracers
      real globalweathflx
      real sedsa_loc, sedmask_loc(nx_loc)
#endif  
      integer debugFlag
      
      integer m, n, i, j, k, ii, itr
      real dmsk(imt,jmt)  ! SPK: setmom.F

! set some constants
      pi = 4.0 * atan(1.0)
      radian = 360./(2.*pi)
      rho0 = 1.035
      rho0r = c1/rho0
      socn = 0.03475

	  if (nx_loc.gt.imt) then
		print*,'ERROR in MOBI_INI: number of local profiles exceeds imt!'
		stop
	  endif  

	  if (kmmax.ne.km) then
		print*,'ERROR in MOBI_INI: number of depth levels incorrect!'
		stop
	  endif  

      j=1
      kmt(:,:) = 0
	  kmt(1:nx_loc,j) = nz_loc_arr(1:nx_loc)
	  tlat(1:nx_loc,j) = tlat_loc(1:nx_loc)

#if defined O_mobi
      ii=1
      do i=1,nx_loc
        do k=1,kmt(i,j)
	      sg_bathy(i,j,k) = sg_bathy_loc(ii)
# if defined O_mobi_iron
          fe_hydr(i,j,k) = fe_hydr_loc(ii)*1000.0 ! input is molFe/m^3/s; we convert it to mmolFe/m^3/s as per code in mobi.F where it is otherwise done when read from file
# endif
# if defined O_mobi_silicon
          si_hydr(i,j,k) = si_hydr_loc(ii) ! input is  mol Si per m^3 per sec (otherwise read from file in mobi.F)
# endif
	      ii=ii+1
	    enddo  
	  enddo
#endif


#if defined O_PaTh
      PaTh_lith(:,:,:) = 0.0
# if !defined O_mobi
      PaTh_pom(:,:,:) = 0.0
      PaTh_caco3(:,:,:) = 0.0
      PaTh_opal(:,:,:) = 0.0
#endif      
      ii=1
      do i=1,nx_loc
		do k=1,kmt(i,j)
          PaTh_lith(i,j,k) = lith_loc(ii) ! kg/m^3
# if !defined O_mobi
          PaTh_pom(i,j,k) = pom_loc(ii) ! kg/m^3
          PaTh_caco3(i,j,k) = caco3_loc(ii) ! kg/m^3
          PaTh_opal(i,j,k) = opal_loc(ii) ! kg/m^3
# endif      
          ii=ii+1
        enddo
      enddo

	  wLith(:) = wLith_loc(:) ! m/s
# if !defined O_mobi
	  wPOM(:) = wPOM_loc(:) ! m/s
	  wCaCO3(:) = wCaCO3_loc(:) ! m/s
	  wOpal(:) = wOpal_loc(:) ! m/s
# endif      
#endif


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
        write(*,*)'ERROR in MOBI_INI: number of tracers incorrect!'
        write(*,*)'                   numTracers not equal to nt-2'
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

!     Set the grid spacings. Only the product of these three matters so 
!     we set dxt to the surface area of that cell.
      dxt(1:nx_loc)=area_loc(1:nx_loc) ! cm^2
      dyt(:)=1.0 ! cm
      cst(:)=1.0

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
CSPK  Note that gaost(isalt) and socn are subsequently changed in S/R setmom
      gaost(:) = 0.
      if (isalt .ne. 0 .and. isss .ne. 0) then
        gaost(isalt) = S_surf_glob
C        write(*,*)'Salt: ',gaost(isalt)
        sbc(1:nx_loc,j,isflx) = 0.0 ! we initialize this here to be safe
      endif
      itr = 0
      do n=3,nt
        if (itrc(n) .ne. 0) then
         itr = itr + 1
		 gaost(n) = TR_surf_glob(itr)	
		 write(*,*)n,mapt(n),gaost(n)	 
        endif
      enddo
C     check
      if (itr.ne.nsrc) then
        print*,'ERROR in MOBI_INI: number of tracers incorrect!'
        stop
      endif  
C

      call setmom (1, nx_loc, 1, jmt)

#if defined O_sed
!-----------------------------------------------------------------------
!     sediment setup (requires mom)
!-----------------------------------------------------------------------

      call setsed (1, nx_loc, 1, jmt)

#endif

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

CSPK return some namelist variables
#if defined O_carbon
	  co2atm = co2ccn
# if defined O_carbon_13	  
	  dc13atm = dc13ccn
# endif
# if defined O_carbon_14
	  dc14atm = dc14ccn
# endif
#endif
	  	
CSPK set some flags so that UVic doesn't start accumulating time averages
      euler2 = .false.
      timavgperts = .false.
      tsiperts = .false.

#if defined O_sed
CSPK  Set some flags needed by the sediment model
      addflxo = .true.
      eots = .true.
!     overwrite values set by S/R setsed
      dtsed = numOceanStepsPerSedStep*dtts
!     send these back
      nzmaxSed = nzmax
      ibmaxSed = ibmax
      numSedMixedTracers = 20
      numSedBuriedTracers = 2
      globalweathflx = weathflx
      sedsa_loc = sedsa
      sedmask_loc(1:nx_loc)=0.d0
      where (map_sed(1:nx_loc,1) .ne. 0.) sedmask_loc(1:nx_loc) = 1.d0
#endif
      
CSPK  Check
      if (dtts.ne.DeltaT) then
         print*,'ERROR: MOBI_INI: dtts does not match DeltaT'
         stop
      endif
      
      if (debugFlag.eq.1) then
        print*,'Tracer time steps: dtts=',dtts,' c2dtts=',c2dtts
#if defined O_mobi
        print*,' dtnpzd=',dtnpzd
#endif        
      endif     

      END
