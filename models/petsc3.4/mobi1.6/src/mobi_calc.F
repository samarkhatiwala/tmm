#include "MOBI_TMM_OPTIONS.h"

      SUBROUTINE MOBI_CALC(kmt_loc, tlat_loc, day_loc, relyr_loc, 
     &     TEMP, SALT, TR_surf_glob,
# if defined O_carbon
     &     co2atm, winds_loc,
#endif      
#  if defined O_npzd_nitrogen
     &     sg_bathy_loc,
#  endif
#  if defined O_npzd_fe_limitation
     &     fe_loc,
#  endif
#  if defined O_embm
     &     swr_loc,
#  endif
#  if defined O_ice
#   if !defined O_ice_cpts
     &     aice_loc, hice_loc, hsno_loc,
#   endif
#  endif
     & emp_loc, emp_glob, debugFlag)

      implicit none

#include "size.h"

      integer kmt_loc
      real tlat_loc, day_loc, relyr_loc
      real TEMP(km), SALT(km)
      real TR_surf_glob(nsrc)
C      real TR_loc(km,nsrc), JTR_loc(km,nsrc)

# if defined O_carbon
      real co2atm
      real winds_loc
#endif      
#  if defined O_npzd_nitrogen
      real sg_bathy_loc(km)
#  endif
#  if defined O_npzd_fe_limitation
      real fe_loc(km) ! only the first three elements are used
#  endif
#  if defined O_embm
	  real swr_loc
#  endif
#  if defined O_ice
#   if !defined O_ice_cpts
      real aice_loc, hice_loc, hsno_loc
#   endif
#  endif
      real emp_loc, emp_glob
      integer debugFlag

      integer i, k, j, ip, n, jrow, js
      real tarea, tsflx, rsocn, tmp
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW

      integer je, limit, joff, is, ie, m, index, itr

#include "param.h"
#include "pconst.h"
#include "stdunits.h"
C#include "accel.h"
#include "coord.h"
C#include "cregin.h"
#include "csbc.h"
C#include "emode.h"
#include "grdvar.h"
C#include "hmixc.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
C#include "switch.h"
C#include "timeavgs.h"
#include "tmngr.h"
C#include "vmixc.h"
# if defined O_npzd || defined O_carbon_14
#include "npzd.h"
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
      real t_in(km), po4_in(km)
#  if defined O_npzd_o2
      real o2_in(km)
#  endif
#  if defined O_carbon_13
      real s_in(km), dic_in(km), alk_in(km), co2_in, dic13_in(km)
#  endif
#  if defined O_npzd_nitrogen
      real no3_in(km), sgb_in(km)
#   if defined O_npzd_nitrogen_15
      real din15_in(km)
#   endif
#  endif
#  if defined O_npzd_fe_limitation
      real felimit_in(km), felimit_D_in(km)
#  endif
      real expo, tnpzd(km,ntnpzd)
#  if defined O_embm
#include "atm.h"
#   if defined O_carbon_13
#include "cembm.h"
#   endif
#  endif
# endif
# if defined O_carbon_fnpzd
#include "calendar.h"
# endif

# if defined O_npzd || defined O_carbon_14
      real ai, hi, hs, rctheta, dayfrac, swr, declin
      real src(imt,km,jsmw:jemw,nsrc)
      common/mobicomm/src
# endif

	  kmt(1,1) = kmt_loc
	  tlat(1,1) = tlat_loc
	  dayoyr = day_loc ! relative day number referenced to the beginning of the current year
	  relyr = relyr_loc ! number of years (and fractional years) of model

C Copy tracers over to UVIC/MOBI arrays
C The sequence of variables is the same as that in mobi_ini.F
C We always pass T/S
	  t(1,1:kmt_loc,1,itemp,taum1) = TEMP(1:kmt_loc) ! !degree Celsius
	  t(1,1:kmt_loc,1,isalt,taum1) = SALT(1:kmt_loc)       
#  if defined O_npzd || defined O_carbon_14	  
C	  itr=0
C      do n=3,nt
C        if (itrc(n) .ne. 0) then
C         itr = itr + 1
C		 t(1,1:kmt_loc,1,n,taum1) = TR_loc(1:kmt_loc,itr)
C        endif
C      enddo
CCSPK check
C      if (itr.ne.nsrc) then
C        print*,'ERROR in mobi_calc: number of tracers incorrect!'
C        stop
C      endif  
#  endif

CSPK set SBC arrays
       i = 1
       j = 1
	   if (isst .ne. 0) sbc(i,j,isst) = t(i,1,j,itemp,taum1)
	   if (isss .ne. 0) sbc(i,j,isss) = t(i,1,j,isalt,taum1)
# if defined O_carbon
            if (issdic .ne. 0) sbc(i,j,issdic) = t(i,1,j,idic,taum1)
#  if defined O_carbon_13
            if (issdic13 .ne. 0) sbc(i,j,issdic13) = 
     &           t(i,1,j,idic13,taum1)
#  endif
#  if defined O_carbon_14
            if (issc14 .ne. 0) sbc(i,j,issc14) = t(i,1,j,ic14,taum1)
#  endif
# endif
# if defined O_npzd_alk
            if (issalk .ne. 0) sbc(i,j,issalk) = t(i,1,j,ialk,taum1)
# endif
# if defined O_npzd_o2
            if (isso2 .ne. 0) sbc(i,j,isso2) = t(i,1,j,io2,taum1)
# endif
# if defined O_npzd
            if (isspo4 .ne. 0) sbc(i,j,isspo4) = t(i,1,j,ipo4,taum1)
            if (issdop .ne. 0) sbc(i,j,issdon) = t(i,1,j,idop,taum1)
            if (issphyt .ne. 0) sbc(i,j,issphyt) = t(i,1,j,iphyt,taum1)
            if (isszoop .ne. 0) sbc(i,j,isszoop) = t(i,1,j,izoop,taum1)
            if (issdetr .ne. 0) sbc(i,j,issdetr) = t(i,1,j,idetr,taum1)
#  if defined O_npzd_nitrogen
            if (issno3 .ne. 0) sbc(i,j,issno3) = t(i,1,j,ino3,taum1)
            if (issdon .ne. 0) sbc(i,j,issdon) = t(i,1,j,idon,taum1)
            if (issdiaz .ne. 0) sbc(i,j,issdiaz) = t(i,1,j,idiaz,taum1)
#   if defined O_npzd_nitrogen_15
            if (issdin15 .ne. 0) sbc(i,j,issdin15) 
     &           = t(i,1,j,idin15,taum1)
            if (issdon15 .ne. 0) sbc(i,j,issdon15) 
     &           = t(i,1,j,idon15,taum1)
            if (issphytn15 .ne. 0) sbc(i,j,issphytn15) 
     &           = t(i,1,j,iphytn15,taum1)
            if (isszoopn15 .ne. 0) sbc(i,j,isszoopn15) 
     &           = t(i,1,j,izoopn15,taum1)
            if (issdetrn15 .ne. 0) sbc(i,j,issdetrn15) 
     &           = t(i,1,j,idetrn15,taum1)
            if (issdiazn15 .ne. 0) sbc(i,j,issdiazn15) 
     &           = t(i,1,j,idiazn15,taum1)
#   endif
#  endif
#  if defined O_carbon_13
            if (issdoc13 .ne. 0) sbc(i,j,issdoc13) = 
     &           t(i,1,j,idoc13,taum1)
            if (issphytc13 .ne. 0) sbc(i,j,issphytc13) = 
     &           t(i,1,j,iphytc13,taum1)
            if (isszoopc13 .ne. 0) sbc(i,j,isszoopc13) = 
     &           t(i,1,j,izoopc13,taum1)
            if (issdetrc13 .ne. 0) sbc(i,j,issdetrc13) = 
     &           t(i,1,j,idetrc13,taum1)
#   if defined O_npzd_nitrogen
            if (issdiazc13 .ne. 0) sbc(i,j,issdiazc13) = 
     &           t(i,1,j,idiazc13,taum1)
#   endif
#  endif
#endif

#ifdef O_carbon
      sbc(i,j,iws) = winds_loc
#endif

#if defined O_mtlm && defined O_carbon
      sbc(i,j,inpp) = 0.0
      sbc(i,j,isr) = 0.0
      sbc(i,j,iburn) = 0.0
#endif
      
CSPK need to set global averages as well
      itr = 0
      do n=3,nt
C        if (itrc(n) .ne. 0) then
         itr = itr + 1
		 gaost(n) = TR_surf_glob(itr)
C        endif
      enddo

C     &  call areaavg (sbc(1,1,issdic), dmsk, gaost(idic))
C     &  call areaavg (sbc(1,1,issdic13), dmsk, gaost(idic13))
C     &  call areaavg (sbc(1,1,issalk), dmsk, gaost(ialk))
C     &  call areaavg (sbc(1,1,isso2), dmsk, gaost(io2))
C     &  call areaavg (sbc(1,1,isspo4), dmsk, gaost(ipo4))
C     &  call areaavg (sbc(1,1,issdop), dmsk, gaost(idop))
C     &  call areaavg (sbc(1,1,issphyt), dmsk, gaost(iphyt))
C     &  call areaavg (sbc(1,1,isszoop), dmsk, gaost(izoop))
C     &  call areaavg (sbc(1,1,issdetr), dmsk, gaost(idetr))
C     &  call areaavg (sbc(1,1,issno3), dmsk, gaost(ino3))
C     &  call areaavg (sbc(1,1,issdon), dmsk, gaost(idon))
C     &  call areaavg (sbc(1,1,issdiaz), dmsk, gaost(idiaz))
C     &  call areaavg (sbc(1,1,issdin15), dmsk, gaost(idin15))
C     &  call areaavg (sbc(1,1,issdon15), dmsk, gaost(idon15))
C     &  call areaavg (sbc(1,1,issphytn15), dmsk, gaost(iphytn15))
C     &  call areaavg (sbc(1,1,isszoopn15), dmsk, gaost(izoopn15))
C     &  call areaavg (sbc(1,1,issdetrn15), dmsk, gaost(idetrn15))
C     &  call areaavg (sbc(1,1,issdiazn15), dmsk, gaost(idiazn15))
C     &  call areaavg (sbc(1,1,issdoc13), dmsk, gaost(idoc13))
C     &  call areaavg (sbc(1,1,issphytc13), dmsk, gaost(iphytc13))
C     &  call areaavg (sbc(1,1,isszoopc13), dmsk, gaost(izoopc13))
C     &  call areaavg (sbc(1,1,issdetrc13), dmsk, gaost(idetrc13))
C     &  call areaavg (sbc(1,1,issdiazc13), dmsk, gaost(idiazc13))
C     &  call areaavg (sbc(1,1,issc14), dmsk, gaost(ic14))
C     &  call areaavg (sbc(1,1,isscfc11), dmsk, gaost(icfc11))
C     &  call areaavg (sbc(1,1,isscfc12), dmsk, gaost(icfc12))

CSPK: code taken from tracer.F
CSPK set some other data
#   if defined O_carbon_co2_2d
        at(1,1,2,ico2) = co2atm
#   else
        co2ccn = co2atm
#   endif

#  if defined O_npzd_fe_limitation
	  do k=1,kmt_loc
	    if (k.le.3) then
	      fe_dissolved(1,1,k,1) = fe_loc(k)
	    endif  
	  enddo
#  endif

#  if defined O_npzd_nitrogen
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

       js = 1
       je = 1
       is = 1
       ie = 1
       joff = 0

!-----------------------------------------------------------------------
!       get the atmospheric S.B.C.
!-----------------------------------------------------------------------
	  call gasbc (1, imt, 1, jmt)

C      if (flag.eq.1) then
C      do j=js,je
C        do i=is,ie
C# if defined O_carbon
C          print*,'idicflx',sbc(i,j,idicflx)
C#  if defined O_carbon_13
C          print*,'idic13flx',sbc(i,j,idic13flx)
C#  endif
C#  if defined O_carbon_14
C          print*,'ic14flx',sbc(i,j,ic14flx)
C#  endif
C# endif
C# if defined O_npzd_alk
C          print*,'ialkflx',sbc(i,j,ialkflx)
C# endif
C# if defined O_npzd_o2
C          print*,'io2flx',sbc(i,j,io2flx)
C# endif
C# if defined O_npzd
C          print*,'ipo4flx',sbc(i,j,ipo4flx)
C          print*,'idopflx',sbc(i,j,idopflx)
C#  if !defined O_npzd_no_vflux
C          print*,'iphytflx',sbc(i,j,iphytflx)
C          print*,'izoopflx',sbc(i,j,izoopflx)
C          print*,'idetrflx',sbc(i,j,idetrflx)
C#  endif
C#  if defined O_npzd_nitrogen
C          print*,'ino3flx',sbc(i,j,ino3flx)
C          print*,'idonflx',sbc(i,j,idonflx)
C#   if !defined O_npzd_no_vflux
C          print*,'idiazflx',sbc(i,j,idiazflx)
C#   endif
C#   if defined O_npzd_nitrogen_15
C          print*,'idin15flx',sbc(i,j,idin15flx)
C          print*,'idon15flx',sbc(i,j,idon15flx)
C#    if !defined O_npzd_no_vflux
C          print*,'iphytn15flx',sbc(i,j,iphytn15flx)
C          print*,'izoopn15flx',sbc(i,j,izoopn15flx)
C          print*,'idetrn15flx',sbc(i,j,idetrn15flx)
C          print*,'idiazn15flx',sbc(i,j,idiazn15flx)
C#    endif
C#   endif          
C#  endif
C#  if defined O_carbon_13
C          print*,'idoc13flx',sbc(i,j,idoc13flx)
C#   if !defined O_npzd_no_vflux
C          print*,'iphytc13flx',sbc(i,j,iphytc13flx)
C          print*,'izoopc13flx',sbc(i,j,izoopc13flx)
C          print*,'idetrc13flx',sbc(i,j,idetrc13flx)
C#    if defined O_npzd_nitrogen
C          print*,'idiazc13flx',sbc(i,j,idiazc13flx)
C#    endif
C#   endif
C#  endif
C# endif
C        enddo
C      enddo
C      endif

CSPK: code taken from gosbc.F
#  if defined O_carbon || defined O_npzd_alk || defined O_npzd_o2 || defined O_npzd || defined O_cfcs_data || defined O_cfcs_data_transient
!-----------------------------------------------------------------------
!     add normalized virtual fluxes to other tracers
!-----------------------------------------------------------------------
      tarea = 0.
      tsflx = emp_glob*rho0*socn 
      rsocn = 1./socn
C      do j=2,jmtm1
C        do i=2,imtm1
C          if (tmsk(i,j) .ge. 0.5) then
C            area = dxt(i)*dyt(j)*cst(j)
C            tarea = tarea + area
C            tsflx = tsflx + sbc(i,j,isflx)*area
C          endif
C        enddo
C      enddo
C      tsflx = tsflx/tarea
C      do j=2,jmtm1
C        do i=2,imtm1
       i = 1
       j = 1
       sbc(i,j,isflx)=emp_loc*rho0*socn       
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
#  endif

CSPK: code taken from tracer.F
# if defined O_npzd
!-----------------------------------------------------------------------
!     ocean biogeochemistry and isotopes
!-----------------------------------------------------------------------
      declin = sin((mod(relyr,1.) - 0.22)*2.*pi)*0.4   ! declination

      nbio = c2dtts/dtnpzd 
      dtbio = c2dtts/nbio
      rdtts = 1./c2dtts
      rnbio = 1./nbio
      do j=js,je
        jrow = j + joff
        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then
#  if defined O_ice
#   if defined O_ice_cpts
            ai = 0.
            hi = 0.
            hs = 0.
            do n=1,ncat
               ai =  ai + A(i,jrow,2,n)
               hi =  hi + heff(i,jrow,2,n)
               hs =  hs + hseff(i,jrow,2,n)
            enddo
#   else
            ai = aice(i,jrow,2)
            hi = hice(i,jrow,2)
            hs = hsno(i,jrow,2)
#   endif
#  else
            ai = 0.
            hi = 0.
            hs = 0.
#  endif
!     calculate day fraction and incoming solar
!     angle of incidence = lat - declin, refraction index = 1.33
            rctheta = max(-1.5, min(1.5, tlat(i,jrow)/radian - declin))
            rctheta = kw/sqrt(1. - (1. - cos(rctheta)**2.)/1.33**2.)
            dayfrac = min( 1., -tan(tlat(i,jrow)/radian)*tan(declin))
            dayfrac = max(1e-12, acos(max(-1., dayfrac))/pi)
#  if defined O_embm
            swr = tap*dnswr(i,jrow)*1e-3*(1. + ai*(exp(-ki*(hi + hs))
     &                                          - 1.))
#  else
            swr = tap*200.
#  endif

C         print*,'MCALC',rctheta,tlat(i,jrow),radian,declin,dayfrac
!        set prognostic variables
            tnpzd(:,1) = t(i,:,j,ipo4,taum1)
            tnpzd(:,2) = t(i,:,j,idop,taum1)
            tnpzd(:,3) = t(i,:,j,iphyt,taum1)
            tnpzd(:,4) = t(i,:,j,izoop,taum1)
            tnpzd(:,5) = t(i,:,j,idetr,taum1)
#  if defined O_npzd_nitrogen
            tnpzd(:,6) = t(i,:,j,idic,taum1)
            tnpzd(:,7) = t(i,:,j,ino3,taum1)
            tnpzd(:,8) = t(i,:,j,idon,taum1)
            tnpzd(:,9) = t(i,:,j,idiaz,taum1)
#   if defined O_npzd_nitrogen_15
            tnpzd(:,10) = t(i,:,j,idin15,taum1)
            tnpzd(:,11) = t(i,:,j,idon15,taum1)
            tnpzd(:,12) = t(i,:,j,iphytn15,taum1)
            tnpzd(:,13) = t(i,:,j,izoopn15,taum1)
            tnpzd(:,14) = t(i,:,j,idetrn15,taum1)
            tnpzd(:,15) = t(i,:,j,idiazn15,taum1)
#   endif
#  endif
#  if defined O_carbon_13
            tnpzd(:,16) = t(i,:,j,idic13,taum1)
            tnpzd(:,17) = t(i,:,j,idoc13,taum1)
            tnpzd(:,18) = t(i,:,j,iphytc13,taum1)
            tnpzd(:,19) = t(i,:,j,izoopc13,taum1)
            tnpzd(:,20) = t(i,:,j,idetrc13,taum1)
#   if defined O_npzd_nitrogen
            tnpzd(:,21) = t(i,:,j,idiazc13,taum1)
#   endif
#  endif
!        set other input variables
            t_in(:) = t(i,:,j,itemp,taum1) !degree Celsius
            po4_in(:) = t(i,:,j,ipo4,taum1)
            o2_in(:) = t(i,:,j,io2,taum1)*1000. ! uM
#  if defined O_carbon_13
            s_in(:) = 1.e3*t(i,:,j,isalt,taum1) + 35.0
            dic_in(:) = t(i,:,j,idic,taum1) !mol/m^3
            alk_in(:)  = t(i,:,j,ialk,taum1) !eq/m^3
#   if defined O_carbon_co2_2d
            co2_in = at(i,j,2,ico2)
#   else
            co2_in = co2ccn
#   endif
            dic13_in(:) = tnpzd(:,16) !mol/m^3
#  endif
#  if defined O_npzd_nitrogen
            no3_in(:) = tnpzd(:,7)
#   if defined O_npzd_nitrogen_15
            din15_in(:) = tnpzd(:,10)
#   endif
            sgb_in(:) = sg_bathy(i,j,:)
#  endif
#  if defined O_npzd_fe_limitation
            do k=1,kmt(i,jrow)
               call fe_limit(i,j,k,felimit_in(k),felimit_D_in(k))
            enddo
#  endif
!  end input
            call mobi (kmt(i,jrow), c2dtts, rctheta, dayfrac, swr, tnpzd
     &,                t_in, po4_in
# if defined O_npzd_o2
     &,                o2_in
# endif
# if defined O_carbon_13
     &,                s_in, dic_in, alk_in, co2_in, dic13_in
# endif
# if defined O_npzd_nitrogen
     &,                no3_in, sgb_in
#  if defined O_npzd_nitrogen_15
     &,                din15_in
#  endif
# endif
#  if defined O_npzd_fe_limitation
     &,                felimit_in, felimit_D_in
# endif
# if defined O_sed
     &,                sedcorgflx, sedcalflx)
# endif
     &,                src(i,:,j,:))

#  if defined O_sed
            if (addflxo .and. eots) then
               sbc(i,jrow,irorg) = sbc(i,jrow,irorg) + sedcorgflx
               sbc(i,jrow,ircal) = sbc(i,jrow,ircal) + sedcalflx
            endif
#  endif
#  if defined O_time_averages && defined O_save_npzd

!-----------------------------------------------------------------------
!             accumulate time averages
!-----------------------------------------------------------------------
            if (timavgperts .and. .not. euler2) then
             expo = rprca
             ta_rprocal(i,jrow) =  ta_rprocal(i,jrow) + rprca
             do k = 1, kmt(i,jrow)
              ta_rnpp(i,k,jrow) = ta_rnpp(i,k,jrow) + rnpp(k)
              ta_rgraz(i,k,jrow) = ta_rgraz(i,k,jrow) + rgraz(k)
              ta_rgraz_Z(i,k,jrow) = ta_rgraz_Z(i,k,jrow) + rgraz_Z(k)
              ta_rgraz_Det(i,k,jrow) = ta_rgraz_Det(i,k,jrow)
     &                                 + rgraz_Det(k)
              ta_rmorp(i,k,jrow) = ta_rmorp(i,k,jrow) + rmorp(k)
              ta_rmorpt(i,k,jrow)= ta_rmorpt(i,k,jrow) + rmorpt(k)
              ta_rmorz(i,k,jrow) = ta_rmorz(i,k,jrow) + rmorz(k)
              ta_rexcr(i,k,jrow) = ta_rexcr(i,k,jrow) + rexcr(k)
#   if defined O_npzd_nitrogen
              ta_rnpp_D(i,k,jrow) = ta_rnpp_D(i,k,jrow) + rnpp_D(k)
              ta_rgraz_D(i,k,jrow) = ta_rgraz_D(i,k,jrow) + rgraz_D(k)
              ta_rmorp_D(i,k,jrow) = ta_rmorp_D(i,k,jrow) + rmorp_D(k)
              ta_rmorpt_D(i,k,jrow) = ta_rmorpt_D(i,k,jrow)
     &                               + rmorpt_D(k)
              ta_rnfix(i,k,jrow) = ta_rnfix(i,k,jrow) + rnfix(k)
              ta_rbdeni(i,k,jrow) = ta_rbdeni(i,k,jrow) + rbdeni(k)
#   endif
#   if defined O_npzd_extra_diagnostics
              ta_ravej(i,k,jrow) = ta_ravej(i,k,jrow) + ravej(k)
              ta_ravej_D(i,k,jrow) = ta_ravej_D(i,k,jrow) + ravej_D(k)
              ta_rgmax(i,k,jrow) = ta_rgmax(i,k,jrow) + rgmax(k)
              ta_rno3P(i,k,jrow) = ta_rno3P(i,k,jrow) + rno3P(k)
              ta_rpo4P(i,k,jrow) = ta_rpo4P(i,k,jrow) + rpo4P(k)
              ta_rpo4_D(i,k,jrow) = ta_rpo4_D(i,k,jrow) + rpo4_D(k)
#   endif
              expo = expo*dztr(k)
              ta_rremi(i,k,jrow) = ta_rremi(i,k,jrow) + rremi(k)
              ta_rexpo(i,k,jrow) = ta_rexpo(i,k,jrow) + rexpo(k)
              expo = expo - rprca*rcak(k)
              ta_rexpocal(i,k,jrow) = ta_rexpocal(i,k,jrow) + expo
              expo = expo*dzt(k)
#   if defined O_npzd_nitrogen && defined O_npzd_o2
              ta_rwcdeni(i,k,jrow) = ta_rwcdeni(i,k,jrow) 
     &           + rwcdeni(k)
#   endif
             enddo
            endif
#  endif
          endif
        enddo
      enddo

# endif
# if defined O_carbon_fnpzd

!-----------------------------------------------------------------------
!    option to write and read fixed biological fluxes for carbon and
!    alkalinity. code assumes that restarts line up with the first
!    record. the repeating time period and restarts should probably
!    line up with the start of the year.
!-----------------------------------------------------------------------

      do j=js,je
        jrow = j + joff

!       adjust counters and read fluxes, if first row of slab
        if (jrow .eq. 2) then
          mfnpzd = mfnpzd + 1
          if (mfnpzd .gt. mxfnpzd) mfnpzd = 1
          if (nfnpzd .le. mxfnpzd) nfnpzd = nfnpzd + 1
!         read fluxes, if enough records written
          if (nfnpzd .gt. mxfnpzd) then
            fname = new_file_name ("O_fnpzd.nc")
            inquire (file=trim(fname), exist=exists)
            if (exists) then
              ib(:) = 1
              ic(:) = 1
              ic(1) = imtm2
              ic(2) = jmtm2
              ic(3) = km
              ib(4) = mfnpzd
              print*, "Reading fnpzd from: ", trim(fname),
     &          " record:", mfnpzd
              call openfile (fname, iou)
#  if defined O_carbon
              call getvara ('O_fnpzd1', iou, imtm2*jmtm2*km, ib, ic
     &,         tmpijk, c1, c0)
              fnpzd(2:imtm1,2:jmtm1,1:km,1)=tmpijk(1:imtm2,1:jmtm2,1:km)
#  endif
#  if defined O_npzd_alk
              call getvara ('O_fnpzd2', iou, imtm2*jmtm2*km, ib, ic
     &,         tmpijk, c1, c0)
              fnpzd(2:imtm1,2:jmtm1,1:km,2)=tmpijk(1:imtm2,1:jmtm2,1:km)
#  endif
            endif
          endif
        endif

        if (nfnpzd .le. mxfnpzd) then

!         store fluxes, if enough records written
          do i=is,ie
            do k=1,km
#  if defined O_carbon
              fnpzd(i,jrow,k,1) = src(i,k,j,isdic)
#  endif
#  if defined O_npzd_alk
              fnpzd(i,jrow,k,2) = src(i,k,j,isalk)
#  endif
            enddo
          enddo
!         write fluxes, if last row of slab
          if (jrow .eq. jmt-1) then
            fname = new_file_name ("O_fnpzd.nc")
            call openfile (fname, iou)
            print*, "Writing fnpzd to: ", trim(fname),
     &          " record:", mfnpzd
            ib(:) = 1
            ic(:) = 1
            ic(1) = imtm2
            ic(2) = jmtm2
            ic(3) = km
            ib(4) = mfnpzd
            time = relyr
            call putvars ('time', iou, mfnpzd, time, c1, c0)
#  if defined O_carbon
            tmpijk(1:imtm2,1:jmtm2,1:km) = fnpzd(2:imtm1,2:jmtm1,1:km,1)
            call putvara ('O_fnpzd1', iou, imtm2*jmtm2*km, ib, ic
     &,       tmpijk, c1, c0)
#  endif
#  if defined O_npzd_alk
            tmpijk(1:imtm2,1:jmtm2,1:km) = fnpzd(2:imtm1,2:jmtm1,1:km,2)
            call putvara ('O_fnpzd2', iou, imtm2*jmtm2*km, ib, ic
     &,       tmpijk, c1, c0)
#  endif
          endif

        else

!         set fluxes to source terms
          do i=is,ie
            do k=1,km
#  if defined O_carbon
              src(i,k,j,isdic) = fnpzd(i,jrow,k,1)
#  endif
#  if defined O_npzd_alk
              src(i,k,j,isalk) = fnpzd(i,jrow,k,2)
#  endif
            enddo
          enddo

        endif

      enddo
# endif
# if defined O_carbon && defined O_carbon_14

!-----------------------------------------------------------------------
!     set source for c14
!-----------------------------------------------------------------------
      do j=js,je
        jrow = j + joff
        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then
            do k=1,kmt(i,jrow)
#  if defined O_npzd
              src(i,k,j,isc14) = src(i,k,j,isdic)*rc14std
     &                         - 3.836e-12*t(i,k,j,ic14,taum1)
#  else
              src(i,k,j,isc14) = - 3.836e-12*t(i,k,j,ic14,taum1)
#  endif
            enddo
          endif
        enddo
      enddo
# endif
# if defined O_sed && !defined O_sed_uncoupled

!-----------------------------------------------------------------------
!     set source terms from sediment model
!-----------------------------------------------------------------------
      do j=js,je
        jrow = j + joff
        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then
            k = kmt(i,jrow)
#  if defined O_carbon && defined O_npzd
            src(i,k,j,isdic) = src(i,k,j,isdic) + sbc(i,j,ibdicfx)
#   if defined O_global_sums
            if (addflxo .and. eots) dtoic = dtoic - sbc(i,j,ibdicfx)
     &        *c2dtts*dtxcel(k)*dxt(i)*dyt(jrow)*cst(jrow)*dzt(k)
#   endif
#  endif
#  if defined O_npzd_alk
            src(i,k,j,isalk) = src(i,k,j,isalk) + sbc(i,j,ibalkfx)
#  endif
          endif
        enddo
      enddo
# endif
C      if (flag.eq.1) then
C      print*,'TEST:',idic,idic13,src(1,1,1,itrc(idic)),src(1,1,1,isdic),sbc(1,1,idicflx),dztr(1)
C      endif
CSPK add surface fluxes as a body force term in surface layer
CSPK code based on setvbc.F
      do j=js,je
        do i=is,ie
# if defined O_carbon
          src(i,1,j,isdic) = src(i,1,j,isdic) 
     &     + sbc(i,j,idicflx)*dztr(1)
#  if defined O_carbon_13
          src(i,1,j,isdic13) = src(i,1,j,isdic13) 
     &     + sbc(i,j,idic13flx)*dztr(1)
#  endif
#  if defined O_carbon_14
          src(i,1,j,isc14) = src(i,1,j,isc14) 
     &     + sbc(i,j,ic14flx)*dztr(1)
#  endif
# endif
# if defined O_npzd_alk
          src(i,1,j,isalk) = src(i,1,j,isalk) 
     &     + sbc(i,j,ialkflx)*dztr(1)
# endif
# if defined O_npzd_o2
          src(i,1,j,iso2) = src(i,1,j,iso2) 
     &     + sbc(i,j,io2flx)*dztr(1)
# endif
# if defined O_npzd
          src(i,1,j,ispo4) = src(i,1,j,ispo4) 
     &     + sbc(i,j,ipo4flx)*dztr(1)
          src(i,1,j,isdop) = src(i,1,j,isdop) 
     &     + sbc(i,j,idopflx)*dztr(1)
#  if !defined O_npzd_no_vflux
          src(i,1,j,isphyt) = src(i,1,j,isphyt) 
     &     + sbc(i,j,iphytflx)*dztr(1)
          src(i,1,j,iszoop)) = src(i,1,j,iszoop) 
     &     + sbc(i,j,izoopflx)*dztr(1)
          src(i,1,j,isdetr) = src(i,1,j,isdetr) 
     &     + sbc(i,j,idetrflx)*dztr(1)
#  endif
#  if defined O_npzd_nitrogen
          src(i,1,j,isno3) = src(i,1,j,isno3) 
     &     + sbc(i,j,ino3flx)*dztr(1)
          src(i,1,j,isdon) = src(i,1,j,isdon) 
     &     + sbc(i,j,idonflx)*dztr(1)
#   if !defined O_npzd_no_vflux
          src(i,1,j,isdiaz) = src(i,1,j,isdiaz) 
     &     + sbc(i,j,idiazflx)*dztr(1)
#   endif
#   if defined O_npzd_nitrogen_15
          src(i,1,j,isdin15) = src(i,1,j,isdin15) 
     &     + sbc(i,j,idin15flx)*dztr(1)
          src(i,1,j,isdon15) = src(i,1,j,isdon15) 
     &     + sbc(i,j,idon15flx)*dztr(1)
#    if !defined O_npzd_no_vflux
          src(i,1,j,isphytn15) = src(i,1,j,isphytn15) 
     &     + sbc(i,j,iphytn15flx)*dztr(1)
          src(i,1,j,iszoopn15) = src(i,1,j,iszoopn15) 
     &     + sbc(i,j,izoopn15flx)*dztr(1)
          src(i,1,j,isdetrn15) = src(i,1,j,isdetrn15)) 
     &     + sbc(i,j,idetrn15flx)*dztr(1)
          src(i,1,j,isdiazn15) = src(i,1,j,isdiazn15) 
     &     + sbc(i,j,idiazn15flx)*dztr(1)
#    endif
#   endif          
#  endif
#  if defined O_carbon_13
          src(i,1,j,isdoc13) = src(i,1,j,isdoc13) 
     &     + sbc(i,j,idoc13flx)*dztr(1)
#   if !defined O_npzd_no_vflux
          src(i,1,j,isphytc13) = src(i,1,j,isphytc13) 
     &     + sbc(i,j,iphytc13flx)*dztr(1)
          src(i,1,j,iszoopc13) = src(i,1,j,iszoopc13) 
     &     + sbc(i,j,izoopc13flx)*dztr(1)
          src(i,1,j,isdetrc13) = src(i,1,j,isdetrc13) 
     &     + sbc(i,j,idetrc13flx)*dztr(1)
#    if defined O_npzd_nitrogen
          src(i,1,j,isdiazc13) = src(i,1,j,isdiazc13) 
     &     + sbc(i,j,idiazc13flx)*dztr(1)
#    endif
#   endif
#  endif
# endif
        enddo
      enddo

C      do j=js,je
C        do i=is,ie
C# if defined O_carbon
C          print*,'idicflx',sbc(i,j,idicflx)
C#  if defined O_carbon_13
C          print*,'idic13flx',sbc(i,j,idic13flx)
C#  endif
C#  if defined O_carbon_14
C          print*,'ic14flx',sbc(i,j,ic14flx)
C#  endif
C# endif
C# if defined O_npzd_alk
C          print*,'ialkflx',sbc(i,j,ialkflx)
C# endif
C# if defined O_npzd_o2
C          print*,'io2flx',sbc(i,j,io2flx)
C# endif
C# if defined O_npzd
C          print*,'ipo4flx',sbc(i,j,ipo4flx)
C          print*,'idopflx',sbc(i,j,idopflx)
C#  if !defined O_npzd_no_vflux
C          print*,'iphytflx',sbc(i,j,iphytflx)
C          print*,'izoopflx',sbc(i,j,izoopflx)
C          print*,'idetrflx',sbc(i,j,idetrflx)
C#  endif
C#  if defined O_npzd_nitrogen
C          print*,'ino3flx',sbc(i,j,ino3flx)
C          print*,'idonflx',sbc(i,j,idonflx)
C#   if !defined O_npzd_no_vflux
C          print*,'idiazflx',sbc(i,j,idiazflx)
C#   endif
C#   if defined O_npzd_nitrogen_15
C          print*,'idin15flx',sbc(i,j,idin15flx)
C          print*,'idon15flx',sbc(i,j,idon15flx)
C#    if !defined O_npzd_no_vflux
C          print*,'iphytn15flx',sbc(i,j,iphytn15flx)
C          print*,'izoopn15flx',sbc(i,j,izoopn15flx)
C          print*,'idetrn15flx',sbc(i,j,idetrn15flx)
C          print*,'idiazn15flx',sbc(i,j,idiazn15flx)
C#    endif
C#   endif          
C#  endif
C#  if defined O_carbon_13
C          print*,'idoc13flx',sbc(i,j,idoc13flx)
C#   if !defined O_npzd_no_vflux
C          print*,'iphytc13flx',sbc(i,j,iphytc13flx)
C          print*,'izoopc13flx',sbc(i,j,izoopc13flx)
C          print*,'idetrc13flx',sbc(i,j,idetrc13flx)
C#    if defined O_npzd_nitrogen
C          print*,'idiazc13flx',sbc(i,j,idiazc13flx)
C#    endif
C#   endif
C#  endif
C# endif
C        enddo
C      enddo

C      if (flag.eq.1) then
C      do n=3,nt
C        if (itrc(n) .ne. 0) then
C         print*,n-2,mapt(n),src(1,1:kmt_loc,1,itrc(n))
C        endif
C      enddo
C      endif
      
      END
      
