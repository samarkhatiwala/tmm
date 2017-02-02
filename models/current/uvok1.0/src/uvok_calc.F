#ifdef O_TMM
      SUBROUTINE UVOK_CALC(kmt_loc, tlat_loc, day_loc, relyr_loc, 
     &     TEMP, SALT, TR_surf_glob,dz_loc,z,
# if defined O_carbon
#if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
     &     co2atm, 
#endif     
     &     winds_loc,
#endif      
# if defined O_c14ccn_data
     &     dc14ccnnatm, dc14ccnsatm, dc14ccneatm,
#endif      
#  if defined O_npzd_nitrogen_15
     &     sg_bathy_loc,
#  endif
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
     & emp_loc, emp_glob,
# if defined O_carbon
     & gasexfluxloc, totfluxloc, 
# endif
     & debugFlag)

      implicit none

#include "size.h"

      integer kmt_loc
      real tlat_loc, day_loc, relyr_loc
      real TEMP(kmt_loc), SALT(kmt_loc)
      real TR_surf_glob(nsrc)
      real dz_loc(kmt_loc)
      real z(km)

# if defined O_carbon
#if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
      real co2atm
#endif      
      real winds_loc
#endif

# if defined O_c14ccn_data
      real dc14ccnnatm, dc14ccnsatm, dc14ccneatm
#endif      
      
#  if defined O_npzd_nitrogen
!      real sg_bathy_loc(kmt_loc)
#  endif
#  if defined O_npzd_fe_limitation
      real fe_dissolved_loc(kmt_loc) ! only the first min(kmfe,kmt_loc) elements are used
#  endif
#if defined O_npzd_iron
      real fe_adep_loc, fe_detr_flux, fe_hydr_loc(kmt_loc)
#endif

#  if defined O_embm
	  real swr_loc
#  endif
#  if defined O_ice
#   if !defined O_ice_cpts
      real aice_loc, hice_loc, hsno_loc
#   endif
#  endif
      real emp_loc, emp_glob
# if defined O_carbon
      real gasexfluxloc, totfluxloc
# endif 
      integer debugFlag

      integer i, k, j, ip, n
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
      common/uvokcomm/src
# endif

#ifdef O_TMM_partial_step_topo
CSPK Recompute local grid-dependent data
#if defined O_mom && defined O_npzd
      dzt(:) = c0
      dzw(:) = c0
      zt(:) = c0
      zw(:) = c0

      do k=1,kmt_loc
        dzt(k) = dz_loc(k)
        zt(k) = z(k)      
      enddo

      zw(1) = dz_loc(1)      
      do k=2,kmt_loc
        zw(k) = zw(k-1) + dz_loc(k)
      enddo
      
      do k=1,kmt_loc-1
        dzw(k) = zt(k+1) - zt(k)
      enddo
      dzw(0) = zt(1)
      dzw(kmt_loc) = zw(kmt_loc) - zt(kmt_loc)

      do k=1,kmt_loc
        c2dzt(k) = c2*dzt(k)
        dzt2r(k) = c1/c2dzt(k)
      enddo

      dzwr(kmt_loc)  = c1/dzw(kmt_loc)
      dzw2r(kmt_loc) = p5/dzw(kmt_loc)

      do k=1,kmt_loc
        dzwr(k-1)    = c1/dzw(k-1)
        dzw2r(k-1)   = p5/dzw(k-1)
      enddo

      do k=1,kmt_loc
        dztr(k)  = c1/dzt(k)
      enddo
#endif
CCCCC
#endif !O_TMM_partial_step_topo

CSPK  Set index of deepest point and other parameters
	  kmt(1,1) = kmt_loc
	  tlat(1,1) = tlat_loc
	  dayoyr = day_loc ! relative day number referenced to the beginning of the current year
	  relyr = relyr_loc ! number of years (and fractional years) of model

CSPK Copy T/S from TMM to UVOK
CSPK Other tracers should already have been copied by a previous call to S/R uvok_copy_data
	  t(1,1:kmt_loc,1,itemp,taum1) = TEMP(1:kmt_loc) ! !degree Celsius
	  t(1,1:kmt_loc,1,isalt,taum1) = SALT(1:kmt_loc)       

CSPK set SBC arrays: code based on UVic_ESCM.F
CSPK note: we set tau, taum1 and taup1 to all point to the same location in S/R uvok_ini
       i = 1
       j = 1
C       do j=1,jmt
C         do i=1,imt
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
C            if (issdop .ne. 0) sbc(i,j,issdon) = t(i,1,j,idop,taup1)
            if (issphyt .ne. 0) sbc(i,j,issphyt) = t(i,1,j,iphyt,taup1)
            if (isszoop .ne. 0) sbc(i,j,isszoop) = t(i,1,j,izoop,taup1)
            if (issdetr .ne. 0) sbc(i,j,issdetr) = t(i,1,j,idetr,taup1)
#  if defined O_npzd_nitrogen
            if (issno3 .ne. 0) sbc(i,j,issno3) = t(i,1,j,ino3,taup1)
C            if (issdon .ne. 0) sbc(i,j,issdon) = t(i,1,j,idon,taup1)
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
#  if defined O_npzd_iron
       if (issdfe .ne. 0) sbc(i,j,issdfe) = t(i,1,j,idfe,taup1)
       if (issdetrfe .ne. 0) sbc(i,j,issdetrfe) = t(i,1,j,idetrfe,taup1)
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
C         enddo
C       enddo

CSPK The sequence here is based on the time stepping in UVic_ESCM.F

CSPK copy over winds and other data needed subsequently by S/R gasbc and S/R tracer
#ifdef O_carbon
      sbc(i,j,iws) = winds_loc
#endif

#  if defined O_npzd_nitrogen_15
	  do k=1,kmt_loc
	    sg_bathy(1,1,k) = sg_bathy_loc(k)
	  enddo
#  endif

#  if defined O_embm
	  dnswr(1,1) = swr_loc
#  endif

#  if defined O_ice
#   if !defined O_ice_cpts
	  aice(1,1,2) = aice_loc
	  hice(1,1,2) = hice_loc
	  hsno(1,1,2) = hsno_loc
#   endif
#  endif

CSPK set atmospheric CO2 here as the call to co2ccndata in gasbc is to a stub routine where 
C    it would otherwise be read from file
#if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
#   if defined O_carbon_co2_2d
        at(1,1,2,ico2) = co2atm
#   else
        co2ccn = co2atm
#   endif
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
      sbc(i,j,inpp) = 0.0
      sbc(i,j,isr) = 0.0
      sbc(i,j,iburn) = 0.0
#endif

!-----------------------------------------------------------------------
!       get the atmospheric S.B.C.
!-----------------------------------------------------------------------

        call gasbc (1, imt, 1, jmt)

CSPK S/R gasbc modifies co2ccn. This should not have any effect if co2emit and global mean 
C    dic flux calculated by a call to stub routine areaavg are exactly zero as we want them 
C    to be. But to be safe, we reset its value again here.
#if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
#   if defined O_carbon_co2_2d
        at(1,1,2,ico2) = co2atm
#   else
        co2ccn = co2atm
#   endif
#endif

CSPK S/R gasbc computes the air-sea gas exchange flux in umol cm-2 s-1. We convert it here 
C    to mol m-2 s-1.
#if defined O_carbon
      gasexfluxloc = sbc(1,1,idicflx)*1.e-2
#endif

#  if defined O_npzd_fe_limitation
	  do k=1,kmt_loc
	    if (k.le.kmfe) then
	      fe_dissolved(1,1,k,1) = fe_dissolved_loc(k)*1.e9  ! input is in nM; convert here as per code in npzd_src.F where it is otherwise done when read from file
	    endif  
	  enddo
#  endif

C Note: we do this after call to S/R gasbc because sbc(1,1,idfeadep) is set to zero there (before being read from file)
#if defined O_npzd_iron
      sbc(1,1,idfeadep)=fe_adep_loc ! note input is in molFe/m^2/s; it is converted to mmolFe/m^2/s in tracer.F so we don't do it here
	  do k=1,kmt_loc      
        fe_hydr(1,1,k) = fe_hydr_loc(k)*1000.0 ! input is molFe/m^3/s; we convert it to mmolFe/m^3/s as per code in npzd_src.F where it is otherwise done when read from file
      enddo  
#endif

C Salt flux for virtual flux calculation in gosbc; emp_loc should be in cm/s
C Note: emp_glob is now not used and global mean freshwater flux is set to zero in gosbc. 
#if defined O_embm
      sbc(1,1,isflx) = emp_loc*socn
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

      call gosbc (1, imt, 1, jmt)

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
       ie = 1
       joff = 0

       call setvbc(joff, js, je, is, ie)
       call tracer(joff, js, je, is, ie)
       
C      if (debugFlag.eq.1) then
!      n=ialk
!      print*,'TEST:',idic,isdic,itrc(n),src(1,1,1,itrc(n)),stf(1,1,n),
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
      totfluxloc = sbc(1,1,idicflx)*1.e-2
#endif
		
C      if (flag.eq.1) then
C      do n=3,nt
C        if (itrc(n) .ne. 0) then
C         print*,n-2,mapt(n),src(1,1:kmt_loc,1,itrc(n))
C        endif
C      enddo
C      endif
      
      END
#endif      