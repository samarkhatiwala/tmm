#ifdef O_TMM
#include "MOBI_TMM_OPTIONS.h"
#endif

      subroutine gasbc (is, ie, js, je)

#if defined O_embm
!=======================================================================
!     calculate boundary conditions for the atmospheric model
!=======================================================================

      implicit none

      integer ie, is, je, js, i, iem1, isp1, j, jem1, jsp1, k, n

      real sss, sst, xconv, t_in, s_in, dic_in, ta_in, co2_in
      real atmpres, pH, co2star,  dco2star, pCO2
      real dpco2, CO3, Omega_c, Omega_a, scco2, piston_vel, avgflxc
      real calday, f, sco2, o2sat, o2sato, o2surf, piston_o2, cfc11ccn
      real cfc12ccn, wt, sccfc, piston_cfc, sol_cfc, cfcsat, ao, tarea
      real tdc14ccn, h_r, d, f1, f2, f3, f4, f5, area, C2K, tmp, zero
      real depth
      real ak, aaqg, adicg, r13a, r13dic, batmc13, bdic13

#include "size.h"
#include "npzd.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "coord.h"
#include "csbc.h"
# if defined O_mom
#include "mw.h"
# endif
# if defined O_ice
#  if defined O_ice_cpts
#include "cpts.h"
#  endif
#include "ice.h"
# endif
#include "switch.h"
#include "tmngr.h"
#include "cembm.h"
#include "atm.h"
#include "insolation.h"
#include "calendar.h"
#include "grdvar.h"
#include "levind.h"
#include "solve.h"
# if defined O_mtlm
#include "mtlm.h"
c#  if defined O_mtlm_carbon_13
c      include "mtlmc13.h"
c#  endif
c#  if defined O_mtlm_carbon_14
c      include "mtlmc14.h"
c#  endif
# endif
# if defined O_save_carbon_carbonate_chem
#include "diaga.h"
# endif

# if !defined O_embm_annual
      real cosz(is:ie,js:je)
# endif
      real dmsk(is:ie,js:je)

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

!     xconv is constant to convert piston_vel from cm/hr -> cm/s
!     here it is 100.*a*xconv (100 => m to cm, a=0.337, xconv=1/3.6e+05)
      xconv = 33.7/3.6e+05
      xconv = xconv*0.75
      C2K = 273.15
# if defined O_carbon
      atmpres = 1.0       !atm
      zero = 0.
#  if defined O_carbon_13
!     fractionation factors from Zhang et al. 1995 Geochi.        
!     Cosmochi. Acta
      ak = 0.99915 ! kinetic fractionation 
      aaqg = 0.998764           ! aquatic - gas fractionation

#   if defined O_carbon_13_coupled
      r13a = c13ccn/(co2ccn-c13ccn)
#   else
      r13a = (dc13ccn*0.001 + 1.)*rc13std
#   endif
      batmc13 = ak*aaqg*r13a
#  endif
# endif
# if defined O_mom
#  if !defined O_constant_flux_reference || defined O_cfcs_data || defined O_cfcs_data_transient

      dmsk(:,:) = 1.
      where (kmt(:,:) .eq. 0) dmsk(:,:) = 0.
#  endif
#  if !defined O_constant_flux_reference

!-----------------------------------------------------------------------
!     calculate new global average sea surface flux references
!-----------------------------------------------------------------------
      if (issdic .ne. 0)
     &  call areaavg (sbc(1,1,issdic), dmsk, gaost(idic))
      if (issdic13 .ne. 0)
     &  call areaavg (sbc(1,1,issdic13), dmsk, gaost(idic13))
      if (issalk .ne. 0)
     &  call areaavg (sbc(1,1,issalk), dmsk, gaost(ialk))
      if (isso2 .ne. 0)
     &  call areaavg (sbc(1,1,isso2), dmsk, gaost(io2))
      if (isspo4 .ne. 0)
     &  call areaavg (sbc(1,1,isspo4), dmsk, gaost(ipo4))
      if (issdop .ne. 0)
     &  call areaavg (sbc(1,1,issdop), dmsk, gaost(idop))
#   if !defined O_npzd_no_vflux
      if (issphyt .ne. 0)
     &  call areaavg (sbc(1,1,issphyt), dmsk, gaost(iphyt))
      if (isszoop .ne. 0)
     &  call areaavg (sbc(1,1,isszoop), dmsk, gaost(izoop))
      if (issdetr .ne. 0)
     &  call areaavg (sbc(1,1,issdetr), dmsk, gaost(idetr))
#   endif
      if (issno3 .ne. 0)
     &  call areaavg (sbc(1,1,issno3), dmsk, gaost(ino3))
      if (issdon .ne. 0)
     &  call areaavg (sbc(1,1,issdon), dmsk, gaost(idon))
#   if !defined O_npzd_no_vflux
      if (issdiaz .ne. 0)
     &  call areaavg (sbc(1,1,issdiaz), dmsk, gaost(idiaz))
#   endif
      if (issdin15 .ne. 0) 
     &  call areaavg (sbc(1,1,issdin15), dmsk, gaost(idin15))
      if (issdon15 .ne. 0)
     &  call areaavg (sbc(1,1,issdon15), dmsk, gaost(idon15))
#   if !defined O_npzd_no_vflux
      if (issphytn15 .ne. 0)
     &  call areaavg (sbc(1,1,issphytn15), dmsk, gaost(iphytn15))
      if (isszoopn15 .ne. 0)
     &  call areaavg (sbc(1,1,isszoopn15), dmsk, gaost(izoopn15))
      if (issdetrn15 .ne. 0)
     &  call areaavg (sbc(1,1,issdetrn15), dmsk, gaost(idetrn15))
      if (issdiazn15 .ne. 0)
     &  call areaavg (sbc(1,1,issdiazn15), dmsk, gaost(idiazn15))
#   endif
      if (issdoc13 .ne. 0)
     &  call areaavg (sbc(1,1,issdoc13), dmsk, gaost(idoc13))
#   if !defined O_npzd_no_vflux
      if (issphytc13 .ne. 0)
     &  call areaavg (sbc(1,1,issphytc13), dmsk, gaost(iphytc13))
      if (isszoopc13 .ne. 0)
     &  call areaavg (sbc(1,1,isszoopc13), dmsk, gaost(izoopc13))
      if (issdetrc13 .ne. 0)
     &  call areaavg (sbc(1,1,issdetrc13), dmsk, gaost(idetrc13))
      if (issdiazc13 .ne. 0)
     &  call areaavg (sbc(1,1,issdiazc13), dmsk, gaost(idiazc13))
#   endif
      if (issc14 .ne. 0)
     &  call areaavg (sbc(1,1,issc14), dmsk, gaost(ic14))
#  endif
#  if defined O_cfcs_data || defined O_cfcs_data_transient

!     always update cfc flux references since initial is usually zero
      if (isscfc11 .ne. 0)
     &  call areaavg (sbc(1,1,isscfc11), dmsk, gaost(icfc11))
      if (isscfc12 .ne. 0)
     &  call areaavg (sbc(1,1,isscfc12), dmsk, gaost(icfc12))
#  endif
# endif

!-----------------------------------------------------------------------
!     zero totals for new accumulation
!-----------------------------------------------------------------------
      atatm = 0.
      flux(:,:,:) = 0.
# if defined O_plume
      subflux(:,:,:) = 0.
# endif
# if defined O_convect_brine
      cbf(:,:,:) = 0.
      cba(:,:,:) = 0.
# endif
      sbc(:,:,ihflx) = 0.
      sbc(:,:,isflx) = 0.
      sbc(:,:,iro) = 0.
# if defined O_mtlm
      sbc(:,:,iat) = 0.
      sbc(:,:,irh) = 0.
      sbc(:,:,ipr) = 0.
      sbc(:,:,ips) = 0.
      sbc(:,:,iaws) = 0.
      sbc(:,:,iswr) = 0.
# endif
# if defined O_carbon
      sbc(:,:,idicflx) = 0.
#  if defined O_carbon_13
      sbc(:,:,idic13flx) = 0.
#  endif
#  if defined O_carbon_14
      sbc(:,:,ic14flx) = 0.
#  endif
# endif
# if defined O_npzd_alk
      sbc(:,:,ialkflx) = 0.
# endif
# if defined O_npzd_o2
      sbc(:,:,io2flx) = 0.
# endif
# if defined O_npzd
      sbc(:,:,ipo4flx) = 0.
      sbc(:,:,idopflx) = 0.
#  if !defined O_npzd_no_vflux
      sbc(:,:,iphytflx) = 0.
      sbc(:,:,izoopflx) = 0.
      sbc(:,:,idetrflx) = 0.
#  endif
#  if defined O_npzd_nitrogen
      sbc(:,:,ino3flx) = 0.
      sbc(:,:,idonflx) = 0.
#   if !defined O_npzd_no_vflux
      sbc(:,:,idiazflx) = 0.
#   endif
#   if defined O_npzd_nitrogen_15
      sbc(:,:,idin15flx) = 0.
      sbc(:,:,idon15flx) = 0.
#    if !defined O_npzd_no_vflux
      sbc(:,:,iphytn15flx) = 0.
      sbc(:,:,izoopn15flx) = 0.
      sbc(:,:,idetrn15flx) = 0.
      sbc(:,:,idiazn15flx) = 0.
#    endif
#   endif
#   if defined O_carbon_13
      sbc(:,:,idoc13flx) = 0.
#    if !defined O_npzd_no_vflux
      sbc(:,:,iphytc13flx) = 0.
      sbc(:,:,izoopc13flx) = 0.
      sbc(:,:,idetrc13flx) = 0.
#     if defined O_npzd_nitrogen
      sbc(:,:,idiazc13flx) = 0.
#     endif
#    endif
#   endif
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient
      sbc(:,:,icfc11flx) = 0.
      sbc(:,:,icfc12flx) = 0.
# endif
# if defined O_solar_data || defined O_solar_data_transient

!-----------------------------------------------------------------------
!     set solar constant
!-----------------------------------------------------------------------
      call solardata
# endif
# if !defined O_embm_annual

!-----------------------------------------------------------------------
!     update insolation for the current day
!-----------------------------------------------------------------------
!     subroutine decl is expecting a 365.25 day year
      calday = dayoyr*365.25/yrlen
      call decl (calday, eccen, obliq, mvelp, lambm0, sindec, eccf)
      i = (ie-is+1)*(je-js+1)
!     get average zenith angle
      call zenith (i, c0, daylen, daylen, tlat, tlon, sindec, cosz)
      solins(is:ie,js:je) = solarconst*eccf*cosz(is:ie,js:je)
# endif
# if defined O_volcano_data || defined O_volcano_data_transient

!-----------------------------------------------------------------------
!     set anomalous volcanic forcing (modify solins)
!-----------------------------------------------------------------------
      call volcdata
      dmsk(:,:) = 1.
      call areaavg (solins, dmsk, tmp)
      if (tmp .gt. 0.) solins(:,:) = solins(:,:)*(tmp - volcfor)/tmp
# endif
# if defined O_co2emit_data_transient

!-----------------------------------------------------------------------
!     set co2 emissions
!-----------------------------------------------------------------------
      call co2emitdata
#  if defined O_carbon_co2_2d

!-----------------------------------------------------------------------
!     set co2 emissions distribution
!-----------------------------------------------------------------------
      call co2distdata
#  endif
# endif
# if defined O_co2ccn_data || defined O_co2ccn_data_transient || defined O_co2emit_track_co2

!-----------------------------------------------------------------------
!     set co2 concentration or emissions by tracking average co2
!-----------------------------------------------------------------------
      call co2ccndata
# endif
# if defined O_co2emit_track_sat || defined O_embm_vcs

!-----------------------------------------------------------------------
!     set co2 emissions by tracking average surface air temperature
!-----------------------------------------------------------------------
      call satdata
# endif
# if defined O_carbon_14
#  if defined O_c14ccn_data || defined O_c14ccn_data_transient

!-----------------------------------------------------------------------
!     set c14 concentration
!-----------------------------------------------------------------------
      call c14data
#   if defined O_c14ccn_data
      tdc14ccn = 0.
      tarea = 0.
#   endif
#  endif
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient

!-----------------------------------------------------------------------
!     set CFC concentration
!-----------------------------------------------------------------------
      call cfcdata
# endif
# if defined O_aggfor_data || defined O_aggfor_data_transient

!-----------------------------------------------------------------------
!     set additional greenhouse gas forcing
!-----------------------------------------------------------------------
      call aggdata
# endif

!-----------------------------------------------------------------------
!     update any atmospheric data
!-----------------------------------------------------------------------
      call data (is, ie, js, je)
# if defined O_embm_awind

!-----------------------------------------------------------------------
!     calculate winds with new feedback
!-----------------------------------------------------------------------
      call add_awind (is, ie, js, je)
# endif
# if defined O_embm && defined O_sealev_data_transient

!-----------------------------------------------------------------------
!     set anomalous sea level
!-----------------------------------------------------------------------
#  if defined O_sealev_data_transient &&  defined O_sealev_salinity
      dsealev = sealev
      call sealevdata
      dsealev = sealev - dsealev
#  else
      call sealevdata
#  endif

      area = ocnsa/atmsa
      do j=2,jmtm1
        do i=2,imtm1
          elev_sealev(i,j) = sealev*(1. - area)*tmsk(i,j)
     &                     - sealev*area*(1. - tmsk(i,j))
        enddo
      enddo
# endif

!-----------------------------------------------------------------------
!     calculate freezing point of sea water using UNESCO (1983)
!-----------------------------------------------------------------------

# if defined O_save_carbon_carbonate_chem
      sspH(:,:) = 0.
      ssCO3(:,:) = 0.
      ssOc(:,:) = 0.
      ssOa(:,:) = 0.
      sspCO2(:,:) = 0.
# endif
      do j=jsp1,jem1
        do i=isp1,iem1

          if (tmsk(i,j) .ge. 0.5) then

            sss = 1000.0*sbc(i,j,isss) + 35.0
            frzpt(i,j) = -.0575*sss + 1.71e-3*sss**1.5 - 2.155e-4*sss**2
# if defined O_mom
#  if defined O_carbon || defined O_npzd_o2 || defined O_cfcs_data || defined O_cfcs_data_transient
            sst = sbc(i,j,isst)
!           put reasonable limits on sst and sss for chemistry flux calculations
            sst = min(35.,max(sst,-2.))
            sss = min(45.,max(sss,0.))
#   if defined O_ice
#    if defined O_ice_cpts
            ao = 1.
            do n=1,ncat
              ao =  ao - A(i,j,2,n)
            enddo
#    else
            ao = 1. - aice(i,j,2)
#    endif
#   else
            ao = 1.
#   endif
#  endif
#  if defined O_carbon

!-----------------------------------------------------------------------
!           calculate ocean carbon fluxes
!-----------------------------------------------------------------------
            t_in = sst
            s_in = sss
            dic_in = sbc(i,j,issdic)
#   if defined O_npzd_alk
            ta_in = sbc(i,j,issalk)
#   else
            ta_in = 2.36775*sss/(socn*1000.)
#   endif
#   if defined O_carbon_co2_2d
            co2_in = at(i,j,2,ico2)
#   else
            co2_in = co2ccn
#   endif
            call co2calc_SWS (t_in, s_in, dic_in, ta_in, co2_in, atmpres
     &,                       zero, pH, co2star, dco2star, pCO2, dpco2
     &,                       CO3, Omega_c, Omega_a)
#   if defined O_save_carbon_carbonate_chem
            sspH(i,j) = pH
            ssCO3(i,j) = CO3
            ssOc(i,j) = Omega_c
            ssOa(i,j) = Omega_a
            sspCO2(i,j) = pCO2
#   endif
!           Schmidt number for CO2
            scco2 = 2073.1 - 125.62*sst + 3.6276*sst**2
     &            - 0.043219*sst**3
            piston_vel = ao*xconv*((sbc(i,j,iws)*0.01)**2)
     &                  *((scco2/660.)**(-0.5))
!           dic in umol cm-3 or (mol m-3) => flux in umol cm-2 s-1
            sbc(i,j,idicflx) = piston_vel*dco2star
#   if defined O_carbon_co2_2d
!           convert from umol cm-2 s-1 => g cm-2 s-1
            flux(i,j,ico2) = sbc(i,j,idicflx)*12.e-6 + flux(i,j,ico2)
#   endif
#   if defined O_carbon_13
            adicg = 1.01051 - 1.05e-4*sst ! DIC-gas fractionation

            r13dic = sbc(i,j,issdic13)
     &               / (sbc(i,j,issdic)-sbc(i,j,issdic13))
            r13dic = min(r13dic, 2.*rc13std)
            r13dic = max(r13dic, 0.5*rc13std)
            bdic13 = ak*aaqg*r13dic/adicg

            sbc(i,j,idic13flx) = piston_vel
     &           *((batmc13/(1+batmc13))*(dco2star + co2star)
     &           - (bdic13/(1+bdic13))*co2star)
c            sbc(i,j,idic13flx) = sbc(i,j,idicflx)*rc13std/(1+rc13std)
#   endif
#   if defined O_carbon_14

!-----------------------------------------------------------------------
!           calculate ocean c14 fluxes
!-----------------------------------------------------------------------
#    if defined O_c14ccn_data
            if (tlat(i,j) .gt. 20.) then
              dc14ccn = dc14ccnn
            elseif (tlat(i,j) .lt. -20.) then
              dc14ccn = dc14ccns
            else
              dc14ccn = dc14ccne
            endif
            tarea = tarea + dxt(i)*dyt(j)*cst(j)
            tdc14ccn = tdc14ccn + dc14ccn*dxt(i)*dyt(j)*cst(j)
#    endif

            sbc(i,j,ic14flx) = piston_vel*((dco2star + co2star)
     &                         *(1 + dc14ccn*0.001)*rc14std
     &                       - co2star*sbc(i,j,issc14)/sbc(i,j,issdic))
#   endif
#  endif
#  if defined O_npzd_o2

!-----------------------------------------------------------------------
!           calculate ocean oxygen fluxes
!-----------------------------------------------------------------------
!           Schmidt number for O2
            sco2 = 1638.0 - 81.83*sst + 1.483*sst**2 - 0.008004*sst**3
!           piston velocity for O2
            piston_o2 = ao*xconv*((sbc(i,j,iws)*0.01)**2)
     &                  *(sco2/660.0)**(-0.5)
!           oxygen saturation concentration [mol/m^3]
            f1 = alog((298.15 - sst)/(C2K + sst))
            f2 = f1*f1
            f3 = f2*f1
            f4 = f3*f1
            f5 = f4*f1
            o2sat = exp (2.00907 + 3.22014*f1 + 4.05010*f2
     &             + 4.94457*f3 - 2.56847E-1*f4 + 3.88767*f5
     &             + sss*(-6.24523e-3 - 7.37614e-3*f1 - 1.03410e-2*f2
     &             - 8.17083E-3*f3) - 4.88682E-7*sss*sss)
!           Convert from ml/l to mol/m^3
            o2sat = o2sat/22391.6*1000.0
            sbc(i,j,io2flx) = piston_o2*(o2sat - sbc(i,j,isso2))
#  endif
#  if defined O_cfcs_data || defined O_cfcs_data_transient

!-----------------------------------------------------------------------
!           calculate ocean CFC11 fluxes
!-----------------------------------------------------------------------
            if (tlat(i,j) .gt. 10.) then
              cfc11ccn = cfc11ccnn
            elseif (tlat(i,j) .lt. -10.) then
              cfc11ccn = cfc11ccns
            else
              wt = (tlat(i,j) + 10.)/20.
              cfc11ccn = cfc11ccnn*wt + cfc11ccns*(1. - wt)
            endif
!           Schmidt number for CFC11
            sccfc = 3501.8 -210.31*sst + 6.1851*sst**2 -0.07513*sst**3
!           piston velocity for CFC
            piston_cfc = ao*xconv*((sbc(i,j,iws)*0.01)**2)
     &                 *(sccfc/660.0)**(-0.5)
!           cfc saturation concentration [mol/m^3]
            f1 = (sst + 273.16)*0.01
            d = (0.091459 - 0.0157274*f1)*f1 - 0.142382
            sol_cfc = exp(-229.9261 + 319.6552/f1 + 119.4471*alog(f1)
     $              - 1.39165*f1*f1 + sss*d )
!           conversion from mol/(l * atm) to mol/(m3 * pptv)

            cfcsat =  1.0e-12 *1000.*sol_cfc*cfc11ccn
            sbc(i,j,icfc11flx) = piston_cfc*(cfcsat - sbc(i,j,isscfc11))

!-----------------------------------------------------------------------
!           calculate ocean CFC12 fluxes
!-----------------------------------------------------------------------
            if (tlat(i,j) .gt. 10.) then
              cfc12ccn = cfc12ccnn
            elseif (tlat(i,j) .lt. -10.) then
              cfc12ccn = cfc12ccns
            else
              wt = (tlat(i,j) + 10.)/20.
              cfc12ccn = cfc12ccnn*wt + cfc12ccns*(1. - wt)
            endif
!           Schmidt number for CFC12
            sccfc = 3845.4 -228.95*sst + 6.1908*sst**2 -0.067430*sst**3
 !          piston velocity for CFC12
            piston_cfc = ao*xconv*((sbc(i,j,iws)*0.01)**2)
     &                 *(sccfc/660.0)**(-0.5)
!           cfc saturation concentration [mol/m^3]
            f1 = (sst + 273.16)*0.01
            d = (0.091015 - 0.0153924*f1)*f1 - 0.143566
            sol_cfc = exp(-218.0971 + 298.9702/f1 + 113.8049*alog(f1)
     $              - 1.39165*f1*f1 + sss*d )
!           conversion from mol/(l * atm) to mol/(m3 * pptv)
            cfcsat =  1.0e-12 *1000.*sol_cfc*cfc12ccn
            sbc(i,j,icfc12flx) = piston_cfc*(cfcsat - sbc(i,j,isscfc12))
#  endif
# endif
# if defined O_carbon && defined O_mtlm

          else
!-----------------------------------------------------------------------
!           calculate land carbon fluxes
!-----------------------------------------------------------------------
#  if defined O_carbon_co2_2d
!           convert from kg m-2 s-1 => g cm-2 s-1
            flux(i,j,ico2) = (sbc(i,j,inpp) - sbc(i,j,isr)
     &                     - sbc(i,j,iburn))*0.1 + flux(i,j,ico2)
#  else
!           convert from kg m-2 s-1 => umol cm-2 s-1
            sbc(i,j,idicflx) = (sbc(i,j,inpp) - sbc(i,j,isr)
     &                       - sbc(i,j,iburn))*0.1/12.e-6
#   if defined O_carbon_13
#    if defined O_mtlm_carbon_13
            sbc(i,j,idic13flx) = (sbc(i,j,inpp13) - sbc(i,j,isr13)
     &                       - sbc(i,j,iburn13))*0.1/12.e-6
#    else
            sbc(i,j,idic13flx) = (sbc(i,j,inpp) - sbc(i,j,isr)
     &                       - sbc(i,j,iburn))*0.1/12.e-6
     &                         *rc13std/(1+rc13std)
#    endif
#   endif
#  endif
#  if defined O_carbon_14
#   if defined O_mtlm_carbon_14
            sbc(i,j,ic14flx) = (sbc(i,j,inpp14) - sbc(i,j,isr14)
     &                       - sbc(i,j,iburn14))*0.1/12.e-6
#   else
!           use the carbon flux scaled by rc14std for c14
            sbc(i,j,ic14flx) = (sbc(i,j,inpp) - sbc(i,j,isr)
     &                       - sbc(i,j,iburn))*rc14std*0.1/12.e-6
#   endif
#  endif
# endif
          endif
        enddo
      enddo

# if defined O_carbon
#  if defined O_mtlm && defined O_global_sums
      ntlbc = ntlbc + 1

#  endif
!-----------------------------------------------------------------------
!     set boundary conditions for carbon
!-----------------------------------------------------------------------
      call setbcx (sbc(1,1,idicflx), imt, jmt)
#  if defined O_carbon_13
      call setbcx (sbc(1,1,idic13flx), imt, jmt)
#  endif
#  if defined O_carbon_co2_2d
#   if !defined O_co2ccn_user && !defined O_co2ccn_data && !defined O_co2ccn_data_transient
      flux(:,:,ico2) = flux(:,:,ico2) - co2emit*co2dist(:,:,2)
#   endif
      call setbcx (flux(1,1,ico2), imt, jmt)
#  else
      dmsk(:,:) = 1.
      call areaavg (sbc(1,1,idicflx), dmsk, avgflxc)
#   if !defined O_co2ccn_user && !defined O_co2ccn_data && !defined O_co2ccn_data_transient
      co2ccn = co2ccn + (co2emit - avgflxc*12.e-6)*segtim*daylen*gtoppm
#   endif
#   if defined O_carbon_13
      dmsk(:,:) = 1.
      call areaavg (sbc(1,1,idic13flx), dmsk, avgflxc)
#    if defined O_carbon_13_coupled
      c13ccn = c13ccn + (co2emit - avgflxc*12.e-6)*segtim*daylen*gtoppm
!     calculate dc13ccn from c13ccn and co2ccn
      dc13ccn = 1000.*(c13ccn/(co2ccn-c13ccn)/rc13std - 1.)
#    endif
#   endif
#  endif
      carbemit = carbemit + co2emit*atmsa*segtim*daylen*1e-15
#  if defined O_global_sums
!     convert from g cm-2 s-1 to umol s-1 for conversion later
      dtoic = dtoic - co2emit*atmsa*segtim*daylen/12e-6
#  endif
#  if defined O_carbon_14
      dmsk(:,:) = 1.
      call areaavg (sbc(1,1,ic14flx), dmsk, avgflxc)
#   if defined O_carbon_14_coupled
      c14ccn = c14ccn + (c14prod - avgflxc)*12.e-6*segtim*daylen*gtoppm
!     calculate dc14ccn from c14ccn and co2ccn
      dc14ccn = 1000.*(c14ccn/co2ccn/rc14std - 1.)
#   endif

!-----------------------------------------------------------------------
!     set boundary conditions for c14
!-----------------------------------------------------------------------
      call setbcx (sbc(1,1,ic14flx), imt, jmt)
#   if defined O_c14ccn_data
      if (tarea .gt. 0) dc14ccn = tdc14ccn/tarea
#   endif
#  endif
# endif
# if defined O_npzd_o2

!-----------------------------------------------------------------------
!     set boundary conditions for oxygen
!-----------------------------------------------------------------------
      call setbcx (sbc(1,1,io2flx), imt, jmt)
# endif
# if defined O_cfcs_data || defined O_cfcs_data_transient

!-----------------------------------------------------------------------
!     set boundary conditions for CFC11
!-----------------------------------------------------------------------
      call setbcx (sbc(1,1,icfc11flx), imt, jmt)

!-----------------------------------------------------------------------
!     set boundary conditions for CFC12
!-----------------------------------------------------------------------
      call setbcx (sbc(1,1,icfc12flx), imt, jmt)
# endif

!-----------------------------------------------------------------------
!     calculate CO2 forcing
!-----------------------------------------------------------------------
      call co2forc

!-----------------------------------------------------------------------
!     set flags to calculate new coefficients
!-----------------------------------------------------------------------
      newcoef(:,:) = .true.
# if defined O_crop_data_transient || defined O_pasture_data_transient || defined O_agric_data_transient

!-----------------------------------------------------------------------
!     update boundary conditions over vegetation
!-----------------------------------------------------------------------
      call gvsbc
# endif
# if defined O_time_averages

!-----------------------------------------------------------------------
!     zero time averages if not in an averaging period
!-----------------------------------------------------------------------
      if (.not. timavgperts) call ta_embm_tavg (is, ie, js, je, 0)
# endif
# if defined O_time_step_monitor

!-----------------------------------------------------------------------
!     zero time step integrals if not in an averaging period
!-----------------------------------------------------------------------
      if (.not. tsiperts) call ta_embm_tsi (is, ie, js, je, 0)
# endif
#endif

      return
      end
