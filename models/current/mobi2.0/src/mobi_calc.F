
      SUBROUTINE MOBI_CALC(lSize, nx_loc, 
     &     day_loc, relyr_loc, 
     &     TEMP, SALT, TR_surf_glob,
     &     dz_loc,z,
# if defined O_carbon
#if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
     &     co2atm, 
#endif     
     &     winds_loc,
#endif      
# if defined O_c14ccn_data
     &     dc14ccnnatm, dc14ccnsatm, dc14ccneatm,
#endif      
#if defined O_carbon && defined O_carbon_13_coupled && defined O_TMM_interactive_atmosphere
     &     c13o2atm,
#endif      
#  if defined O_npzd_fe_limitation
     &     fe_dissolved_loc,
#  endif
#if defined O_npzd_iron
     &     fe_adep_loc, fe_detr_flux, fe_hydr_loc, 
#endif
#  if defined O_embm
     &     swr_loc,
#  endif
#  if defined O_ice
#   if !defined O_ice_cpts
     &     aice_loc, hice_loc, hsno_loc,
#   endif
#  endif
     &     emp_loc, emp_glob, 
# if defined O_carbon
     &     gasexfluxloc, totfluxloc, 
# endif     
#if defined O_carbon && defined O_carbon_13_coupled && defined O_TMM_interactive_atmosphere
     &     c13gasexfluxloc,
#endif
     &     debugFlag)

      implicit none

#include "size.h"

      integer lSize, nx_loc
      real day_loc, relyr_loc
      real TEMP(lSize), SALT(lSize), TR_surf_glob(nsrc)
      real dz_loc(lSize), z(km)

# if defined O_carbon
#if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
#   if defined O_carbon_co2_2d
      real co2atm(nx_loc)
#   else
      real co2atm
#   endif      
#endif      
      real winds_loc(nx_loc)
#endif

# if defined O_c14ccn_data
      real dc14ccnnatm, dc14ccnsatm, dc14ccneatm
#endif      

#if defined O_carbon && defined O_carbon_13_coupled && defined O_TMM_interactive_atmosphere
      real c13o2atm
#endif
      
#  if defined O_npzd_fe_limitation
      real fe_dissolved_loc(lSize)
#  endif
#if defined O_npzd_iron
      real fe_adep_loc(nx_loc), fe_detr_flux(lSize), fe_hydr_loc(lSize)
#endif

#  if defined O_embm
	  real swr_loc(nx_loc)
#  endif
#  if defined O_ice
#   if !defined O_ice_cpts
      real aice_loc(nx_loc), hice_loc(nx_loc), hsno_loc(nx_loc)
#   endif
#  endif
      real emp_loc(nx_loc), emp_glob
# if defined O_carbon
      real gasexfluxloc(nx_loc), totfluxloc(nx_loc)
# endif
#if defined O_carbon && defined O_carbon_13_coupled && defined O_TMM_interactive_atmosphere
      real c13gasexfluxloc(nx_loc)
#endif
      integer debugFlag

      integer i, k, j, ii, ip, n
      real tarea, tsflx, rsocn, tmp
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW

      integer joff, is, ie, js, je, itr

#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "coord.h"
#include "csbc.h"
#include "grdvar.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "tmngr.h"
# if defined O_npzd || defined O_carbon_14
#include "npzd.h"
#include "calendar.h"
# endif
# if defined O_save_convection || defined O_carbon_14
#include "diaga.h"
      real rrc14std
# endif
# if defined O_ice
#  if defined O_ice_cpts
#include "cpts.h"
#  endif
#include "ice.h"
# endif
# if defined O_npzd
C     real t_in(km), po4_in(km)
C#  if defined O_npzd_o2
C     real o2_in(km)
C#  endif
C#  if defined O_carbon_13
C     real s_in(km), dic_in(km), alk_in(km), co2_in, dic13_in(km)
C#  endif
C#  if defined O_npzd_nitrogen
C     real no3_in(km), sgb_in(km)
C#   if defined O_npzd_nitrogen_15
C     real din15_in(km)
C#   endif
C#  endif
C#  if defined O_npzd_fe_limitation
C     real felimit_in(km), felimit_D_in(km)
C#  endif
C     real expo, tnpzd(km,ntnpzd)
#  if defined O_embm
#include "atm.h"
#   if defined O_carbon_13
#include "cembm.h"
#   endif
#  endif
C# endif
C# if defined O_carbon_fnpzd
C#include "calendar.h"
# endif

# if defined O_npzd || defined O_carbon_14
      real src(imt,km,jsmw:jemw,nsrc)
      common/mobicomm/src
# endif

#ifdef O_TMM_partial_step_topo
CSPK Recompute local grid-dependent data
      dzt(:) = c0
      dzw(:) = c0
      zt(:) = c0
      zw(:) = c0

      do k=1,km_loc
        dzt(k) = dz_loc(k)
        zt(k) = z(k)      
      enddo

      zw(1) = dz_loc(1)      
      do k=2,km_loc
        zw(k) = zw(k-1) + dz_loc(k)
      enddo
      
      do k=1,km_loc-1
        dzw(k) = zt(k+1) - zt(k)
      enddo
      dzw(0) = zt(1)
      dzw(km_loc) = zw(km_loc) - zt(km_loc)

      do k=1,km_loc
        c2dzt(k) = c2*dzt(k)
        dzt2r(k) = c1/c2dzt(k)
      enddo

      dzwr(km_loc)  = c1/dzw(km_loc)
      dzw2r(km_loc) = p5/dzw(km_loc)

      do k=1,km_loc
        dzwr(k-1)    = c1/dzw(k-1)
        dzw2r(k-1)   = p5/dzw(k-1)
      enddo

      do k=1,km_loc
        dztr(k)  = c1/dzt(k)
      enddo
CCCCC

CSPK from S/R mobi_init
!     calculate sinking speed of detritus divided by grid width
      do k=1,km_loc
!     linear increase wd0-200m with depth
         if (zt(k) .lt. mwz) then
            wd(k) = (wd0+mw*zt(k))/daylen/dzt(k) ! [s-1]
# if defined O_npzd_caco3
            wc(k) = (wc0+mw_c*zt(k))/daylen/dzt(k)    ! [s-1]
# endif            
         else
            wd(k) = (wd0+mw*mwz)/daylen/dzt(k) ! [s-1]
# if defined O_npzd_caco3
            wc(k) = (wc0+mw_c*mwz)/daylen/dzt(k)    ! [s-1]
# endif 
         endif
         rkwz(k) = 1./(kw*dzt(k))
      enddo
      ztt(1)=0.0
      do k=1,km_loc-1
         ztt(k+1)=(-1)*zw(k)
      enddo

CSPK from S/R mobi_init
# if defined O_carbon || defined O_npzd_alk

!---------------------------------------------------------------------
!     calculate variables used in calcite remineralization
!---------------------------------------------------------------------

      rcak(1) = -(exp(-zw(1)/dcaco3)-1.0)/dzt(1)
      rcab(1) = 1./dzt(1)
      do k=2,km_loc
        rcak(k) = -(exp(-zw(k)/dcaco3))/dzt(k)
     &          + (exp(-zw(k-1)/dcaco3))/dzt(k)
        rcab(k) = (exp(-zw(k-1)/dcaco3))/dzt(k)
      enddo
# endif

C end recompute
#endif ! O_TMM_partial_step_topo

CSPK  Set some parameters
	  dayoyr = day_loc ! relative day number referenced to the beginning of the current year
	  relyr = relyr_loc ! number of years (and fractional years) of model

	  j = 1

CSPK Copy T/S from TMM to UVic
CSPK Other tracers should already have been copied by a previous call to S/R mobi_copy_data
      ii=1
      do i=1,nx_loc
        do k=1,kmt(i,j)
	      t(i,k,j,itemp,taum1) = TEMP(ii) ! !degree Celsius
	      t(i,k,j,isalt,taum1) = SALT(ii)       
          ii=ii+1
        enddo
      enddo
C     if (debugFlag.eq.1) then
C       write(*,*)'SPK: ',
C    &   t(1,1,j,itemp,taum1),t(1000,1,j,itemp,taum1),
C    &   t(nx_loc,1,j,itemp,taum1),TEMP(3000),TEMP(lSize)
C     endif  

CSPK set SBC arrays
CSPK note: we set tau, taum1 and taup1 to all point to the same location in S/R mobi_ini
	   do n=1,nt
		 if (trsbcindex(n) .ne. 0) then
		   call set_sbc(0, j, j, 1, nx_loc, trsbcindex(n), n, 
     &                  taup1, .false.)
C	   if (debugFlag.eq.1) then
C		 write(*,*)'SPK2: ',trsbcindex(n),itemp,n,mapt(n)
C	   endif  
     
		 endif
	   enddo                      
       
CSPK The sequence here is based on the time stepping in UVic_ESCM.F

CSPK copy over winds and other data needed subsequently by S/R gasbc and S/R tracer
#ifdef O_carbon
	  sbc(1:nx_loc,j,iws) = winds_loc(1:nx_loc)
#endif

#  if defined O_embm
	  dnswr(1:nx_loc,j) = swr_loc(1:nx_loc)
#  endif

#  if defined O_ice
#   if !defined O_ice_cpts
	  aice(1:nx_loc,j,2) = aice_loc(1:nx_loc)
	  hice(1:nx_loc,j,2) = hice_loc(1:nx_loc)
	  hsno(1:nx_loc,j,2) = hsno_loc(1:nx_loc)
#   endif
#  endif

CSPK set atmospheric CO2 here as the call to co2ccndata in gasbc is to a stub routine where 
C    it would otherwise be read from file
#if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
#   if defined O_carbon_co2_2d
        at(1:nx_loc,j,2,ico2) = co2atm(1:nx_loc)
#   else
        co2ccn = co2atm
#   endif
#endif

#if defined O_carbon && defined O_carbon_13_coupled && defined O_TMM_interactive_atmosphere
        c13ccn = c13o2atm
        dc13ccn = 1000.*(c13ccn/(co2ccn-c13ccn)/rc13std - 1.)        
#endif

CSPK set global averages here if necessary as the calls in gasbc are to stub routines
#  if !defined O_constant_flux_reference
      itr = 0
      do n=3,nt
C        if (itrc(n) .ne. 0) then
         itr = itr + 1
		 gaost(n) = TR_surf_glob(itr)
C        endif
      enddo
# endif

CSPK set atmospheric dc14 here as the call to c14data in gasbc is to a stub routine where 
C    it would otherwise be read from file
# if defined O_c14ccn_data
      dc14ccnn = dc14ccnnatm
      dc14ccne = dc14ccneatm
      dc14ccns = dc14ccnsatm
#endif      

CSPK These arrays are first set to zero in S/R glsbc and then updated in S/R mtlm. They 
C    are used in S/R gasbc but as we're not using the land model they should not be 
C    referenced anywhere. Still, we set them to zero anyway to be safe.
#if defined O_mtlm && defined O_carbon
	  sbc(1:nx_loc,j,inpp) = 0.0
	  sbc(1:nx_loc,j,isr) = 0.0
	  sbc(1:nx_loc,j,iburn) = 0.0
#endif

!-----------------------------------------------------------------------
!       get the atmospheric S.B.C.
!-----------------------------------------------------------------------

        call gasbc (1, nx_loc, 1, jmt)

CSPK S/R gasbc modifies co2ccn. This should not have any effect if co2emit and global mean 
C    dic flux calculated by a call to stub routine areaavg are exactly zero as we want them 
C    to be. But to be safe, we reset its value again here.
#if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
#   if defined O_carbon_co2_2d
        at(1:nx_loc,j,2,ico2) = co2atm(1:nx_loc)
#   else
        co2ccn = co2atm
#   endif
#endif

CSPK S/R gasbc modifies c13ccn. We reset it here. 
#if defined O_carbon && defined O_carbon_13_coupled && defined O_TMM_interactive_atmosphere
        c13ccn = c13o2atm
        dc13ccn = 1000.*(c13ccn/(co2ccn-c13ccn)/rc13std - 1.)        
#endif

CSPK S/R gasbc computes the air-sea gas exchange flux in umol cm-2 s-1. We convert it here 
C    to mol m-2 s-1.
#if defined O_carbon
      gasexfluxloc(1:nx_loc) = sbc(1:nx_loc,j,idicflx)*1.e-2
#endif

CSPK S/R gasbc computes the air-sea gas exchange flux in umol cm-2 s-1. We convert it here 
C    to mol m-2 s-1.
#if defined O_carbon && defined O_carbon_13_coupled && defined O_TMM_interactive_atmosphere
      c13gasexfluxloc(1:nx_loc) = sbc(1:nx_loc,j,idic13flx)*1.e-2
#endif      

#  if defined O_npzd_fe_limitation
      ii=1
      do i=1,nx_loc
		do k=1,kmt(i,j)
		  if (k.le.kmfe) then
			fe_dissolved(i,j,k,1) = fe_dissolved_loc(ii)*1.e9  ! input is in nM; convert here as per code in npzd_src.F where it is otherwise done when read from file
		  endif
		  ii=ii+1  
	  enddo
#  endif

C Note: we do this after call to S/R gasbc because sbc(1,1,idfeadep) is set to zero there (before being read from file)
#if defined O_npzd_iron
      sbc(1:nx_loc,j,idfeadep)=fe_adep_loc(1:nx_loc) ! note input is in molFe/m^2/s; it is converted to mmolFe/m^2/s in tracer.F so we don't do it here
      ii=1
      do i=1,nx_loc
		do k=1,kmt(i,j)
          fe_hydr(i,j,k) = fe_hydr_loc(ii)*1000.0 ! input is molFe/m^3/s; we convert it to mmolFe/m^3/s as per code in npzd_src.F where it is otherwise done when read from file
          ii=ii+1
        enddo  
      enddo  
#endif

C Salt flux for virtual flux calculation in gosbc; emp_loc should be in cm/s
C Note: emp_glob is now not used and global mean freshwater flux is set to zero in gosbc. 
#if defined O_embm
      sbc(1:nx_loc,j,isflx) = emp_loc(1:nx_loc)*socn
#endif
      
C#if defined O_mtlm
C!-----------------------------------------------------------------------
C!       get land S.B.C.s
C!-----------------------------------------------------------------------
C
C       call glsbc (1, imt, 1, jmt)
C
C!----------------------------------------------------------------------
C!       call the land-surface and vegetation  model once for each time
C!       step until one segment of "segtim" days is complete.
C!-----------------------------------------------------------------------
C
C       do loop=1,ntspls
C         call mtlm (1, imt, 1, jmt)
C       enddo
C
C#endif

C
C!-----------------------------------------------------------------------
C!       get ocean S.B.C.s
C!-----------------------------------------------------------------------

      call gosbc (1, nx_loc, 1, jmt)

C
C!-----------------------------------------------------------------------
C!       call the ocean model once for each time step until one
C!       segment of "segtim" days is complete. hold ocean S.B.C. fixed
C!       during each segment and predict average S.B.C. for atmos
C!-----------------------------------------------------------------------
C
C       do loop=1,ntspos
C         call mom

CSPK S/R mom calls S/R setvbc and S/R tracer

       js = 1
       je = 1
       is = 1
       ie = nx_loc
       joff = 0 ! just to be safe; for the standard UVic options joff=0

       call setvbc(joff, js, je, is, ie)
       call tracer(joff, js, je, is, ie)
       
C      if (debugFlag.eq.1) then
C      n=idic
C      print*,'TEST:',idic,isdic,itrc(n),src(1,1,1,itrc(n)),stf(1,1,n),
C    &    sbc(1,1,idicflx),dztr(1),tmask(1,1,1)
C      endif
CSPK Add surface fluxes (calculated in setvbc.F) as a body force term in surface layer. Be very careful of tracer indexing!
      do j=js,je
        do i=is,ie
		  do n=3,nt
C	  if (itrc(n) .ne. 0) then
		   src(i,1,j,itrc(n)) = src(i,1,j,itrc(n))+stf(i,j,n)*dztr(1)
C	  endif
		  enddo
		enddo
	  enddo	  
C      if (debugFlag.eq.1) then
C      n=idic
C      print*,'TEST:',idic,isdic,itrc(n),src(1,1,1,itrc(n)),stf(1,1,n),
C    &    sbc(1,1,idicflx),dztr(1),tmask(1,1,1)
C      endif

CSPK Virtual flux is added above to the air-sea gas exchange flux computed in S/R gasbc. Units 
C    are umol cm-2 s-1. We convert here to mol m-2 s-1.
#if defined O_carbon
      totfluxloc(1:nx_loc) = sbc(1:nx_loc,j,idicflx)*1.e-2
#endif
		
C      if (flag.eq.1) then
C      do n=3,nt
C        if (itrc(n) .ne. 0) then
C         print*,n-2,mapt(n),src(1,1:km_loc,1,itrc(n))
C        endif
C      enddo
C      endif
      
      END
      
