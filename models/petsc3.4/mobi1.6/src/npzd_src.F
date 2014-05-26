#ifdef O_TMM
#include "MOBI_TMM_OPTIONS.h"
#endif
#if defined O_mom && defined O_npzd 
      
      subroutine mobii

!     initialize MOBI parameters
      implicit none
#include "size.h"
#include "npzd.h"
#include "calendar.h"
#include "coord.h"
#include "stdunits.h"
#include "scalar.h"
      integer k, ioun
      real alpha, par
# if defined O_npzd_fe_limitation
#include "param.h"
#include "pconst.h"
      integer ib(10), ic(10), iou, m, j, i
      character (120) :: fname, new_file_name
      logical exists
      real c1e9
      real tmpijkm(1:100,1:100,1:3,1:12)
# endif
      namelist /npzd/   alpha, kw, kc, abio, bbio, cbio, k1n, nup
     &,                 gamma1, gbio, epsbio, nuz, nud0, LFe
     &,                 wd0, par, dtnpzd, redctn, ki, redptn, capr
     &,                 dcaco3, nupt0, redotn, jdiar, kzoo, geZ
     &,                 zprefP, zprefZ, zprefDet, zprefD, kfe, kfe_D
     &,                 nupt0_D, sgbdfac, diazntp, nup_D, dfr, dfrt
      
!     set defaults for namelist npzd

      alpha = 0.16  ! Initial slope P-I curve [(W/m^2)^(-1)/day]
      kw = 0.04  ! Light attenuation due to water [1/m]
      kc = 0.047  ! Light atten. by phytoplankton [1/(m*mmol/m^3)]
      ki = 5.0  ! Light attenuation through sea ice & snow
      abio = 0.6  ! a; Maximum growth rate parameter [1/day]
      bbio = 1.066  ! b
      cbio = 1.0  ! c [1/deg_C]
      k1n = 0.7  ! Half saturation constant for N uptake [mmol/m^3]
      nup = 0.03  ! Specific mortality rate (Phytoplankton) [1/day]
      nup_D = 0.0
      nupt0 = 0.015  ! Specific mortality rate (Phytoplankton) [1/day]
      nupt0_D = 0.001
      gamma1 = 0.70  ! gama1; Assimilation efficiency (zpk)
      gbio = 0.38  ! Maximum grazing rate [1/day]
      epsbio = 1.6  ! Prey capture rate [(mmol/m^3)^(-2)day^(-1)]
      nuz = 0.06  ! Quadratic mortality (zpk) [(mmol/m^3)^(-2)day^(-1)]
      nud0 = 0.07  ! remineralization rate [1/day]
      LFe = 1.0  ! Iron limitation
      wd0 = 16.0  ! Sinking speed of detritus at surface [m/day]
      par = 0.43  ! fraction of photosythetically active radiation
      dtnpzd = dtts/4.  ! time step of biology [s]
      redctn = 7.1  ! C/N Redfield ratio (Schneider et al., 2003, GBC)
      redptn = 1./16.  ! P/N Redfield ratio
      capr = 0.022  ! carbonate to carbon production ratio
      dcaco3 = 650000.0  ! remineralisation depth of calcite [cm]
      redotn = 10.6  ! O2/N Redfield ratio (Anderson & Sarmiento, 1994, GBC)
      jdiar = 0.08  ! factor reducing the growth rate of diazotrophs
      kzoo = 0.15  ! half sat. constant for Z grazing [mmol N m^3]
      geZ = 0.6   ! Zooplankton growth efficiency
      sgbdfac = 1.0
      diazntp = 28.
      dfr = 0.08   ! refractory/semi-labile DOM fraction
      dfrt = 0.01
# if defined O_npzd_nitrogen
      zprefP = 0.30       ! Zooplankton preference for P
      zprefZ = 0.30       ! Zooplankton preference for other Z
      zprefDet = 0.30     ! Zooplankton preference for Detritus
      zprefD = 0.10       ! Zooplankton preference for diazotrophs
# else
      zprefP = 0.35       ! Zooplankton preference for P
      zprefZ = 0.35       ! Zooplankton preference for other Z
      zprefDet = 0.30     ! Zooplankton preference for Detritus
# endif
# if defined O_npzd_fe_limitation
      kmfe = 3            ! no. of depth levels used for Fe limitation 
      kfe = 0.12          ! Half saturation constant for Fe limitation
      kfe_D = 0.12        ! Half saturation constant for Diaz Fe limitation
# endif

!     read namelist
      call getunit (ioun, 'control.in', 'f s r')
      read  (ioun, npzd, end=108)
108   continue
      write (stdout, npzd)
      call relunit (ioun)

!     convert units of NPZD parameters to MOM units
      redctn = redctn*1.e-3
      redotn = redotn*1.e-3
      redotp = redotn/redptn
      redctp = redctn/redptn
      redotc = redotn/redctn
      redntp = 1./redptn
      redntc = 1./redctn
      diazptn = 1./diazntp
      k1p   = k1n*redptn
      kw = kw*1.e-2
      kc = kc*1.e-2
      ki = ki*1.e-2
      wd0 = wd0*1.e2
      alpha = alpha/daylen
      abio = abio/daylen
      nup = nup/daylen
      nup_D = nup_D/daylen
      nupt0 = nupt0/daylen
      nupt0_D = nupt0_D/daylen
      gbio = gbio/daylen
      epsbio = epsbio/daylen
      nuz = nuz/daylen
      nud0 = nud0/daylen
      tap = 2.*alpha*par
!     calculate sinking speed of detritus divided by grid width
      do k=1,km
!     linear increase wd0-200m with depth
         if (zt(k) .lt. 100000) then
            wd(k) = (wd0+5.5e-2*zt(k))/daylen/dzt(k) ! [s-1]
         else
            wd(k) = (wd0+5.5e-2*100000)/daylen/dzt(k) ! [s-1]
         endif
         rkwz(k) = 1./(kw*dzt(k))
      enddo
      ztt(1)=0.0
      do k=1,km
         ztt(k+1)=(-1)*zw(k)
      enddo
# if defined O_npzd_fe_limitation
      do k=km,1,-1
         if (zw(k) .gt. 24000) kmfe = k-1
      enddo
#ifndef O_TMM
! read in the dissolved Fe conc from O_dissolved_fe.nc
      fe_dissolved(:,:,:,:) = 0.
      ib(:) = 1
      ic(:) = imtm2
      ic(2) = jmtm2
      ic(3) = kmfe
      ic(4) = 12

      fname = new_file_name ("O_fe_dissolved.nc")
      inquire (file=trim(fname), exist=exists)
      if (exists) then
         c1e9 = 1000000000
         call openfile (trim(fname), iou)
         call getvara ('O_dissolved_fe', iou, ic(1)*ic(2)*ic(3)*ic(4)
     &,                 ib, ic, tmpijkm, c1e9, c0)
         fe_dissolved(2:imtm1,2:jmtm1,1:kmfe,:) = tmpijkm(1:imtm2
     &,                                            1:jmtm2,1:kmfe,:)
!         print*,'fe_dissolved(50,50,1,1)=',fe_dissolved(50,50,1,1)
         do m=1,12
            do k=1,kmfe
               do j=1,jmt
                  fe_dissolved(1,j,k,m) = fe_dissolved(imtm1,j,k,m)
                  fe_dissolved(imt,j,k,m) = fe_dissolved(2,j,k,m)
               enddo
               do i=1,imt
                  fe_dissolved(i,1,k,m) = fe_dissolved(i,2,k,m)
                  fe_dissolved(i,jmt,k,m) = fe_dissolved(2,j,k,m)
               enddo
            enddo
         enddo
      else
         print*,"Warning => Cannot find", trim(fname)
      endif
#endif
# endif

# if defined O_carbon || defined O_npzd_alk

!---------------------------------------------------------------------
!     calculate variables used in calcite remineralization
!---------------------------------------------------------------------

      rcak(1) = -(exp(-zw(1)/dcaco3)-1.0)/dzt(1)
      rcab(1) = -1./dzt(1)
      do k=2,km
        rcak(k) = -(exp(-zw(k)/dcaco3))/dzt(k)
     &          + (exp(-zw(k-1)/dcaco3))/dzt(k)
        rcab(k) = (exp(-zw(k-1)/dcaco3))/dzt(k)
      enddo
# endif
      return
      end
!     END mobii

      subroutine mobi (kmx, twodt, rctheta, dayfrac, swr, tnpzd
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
# if defined O_npzd_fe_limitation
     &,                felimit_in, felimit_D_in
# endif
# if defined O_sed
     &,                sedcorgflx, sedcalflx)
# endif
     &,                src)

!     main driver for Model of Ocean Biogeochemistry and Isotopes (MOBI)
!     input: kmx, twodt, rctheta, dayfrac, swr, tnpzd, t_in, 
!            po4_in, (o2_in, s_in, dic_in, alk_in, co2_in, dic13_in, no3_in,
!            sgb_in, din15_in, felimit_in, felimit_D_in)
!     output: src, (sedcorgflx, sedcalflx)
!             if O_save_npzd is switch on additional output (rnpp etc.)
!                is in npzd.h 
!     copied from tracer.F (Andreas Schmittner, Oct 2, 2013)
!     modified to a vertical column model 
      implicit none
      integer k,kmx
#include "size.h"
#include "coord.h"
#include "grdvar.h"
#include "mw.h"
#include "npzd.h"
      real snpzd(ntnpzd), tnpzd(km,ntnpzd), src(km,nsrc)
      real dic_npzd_sms(km), twodt, rctheta, gl, impo, expo
      real npp, remi, excr, graz, morp, morpt, morz, temp, swr, dayfrac
      real graz_Det, graz_Z, phin, dz, prca, dprca, nud, bct, bctz
      real t_in(km), po4_in(km), sedrr
# if defined O_npzd_extra_diagnostics
      real avej, avej_D, gmax, no3P, po4P, po4_D
# endif
# if defined O_npzd_o2
      real o2_in(km), fo2, so2
# endif
# if defined O_carbon_13
      real s_in(km), dic_in(km), alk_in(km), co2_in, dic13_in(km)
      real rc13impo, rc13expo, ac13_DIC_aq
      real ac13_aq_POC, ac13b, prca13, rdic13, rtdic13
      real pH, co2star, dco2star, pCO2, dpco2, CO3
      real Omega_c, Omega_a, pCO2_old, pCO2_new
      real atmpres, depth
# endif
# if defined O_npzd_nitrogen
      real no3_in(km), sgb_in(km), lno3, sg_bdeni
      real morp_D, lntp
      real npp_D, graz_D, morpt_D, no3flag, wcdeni, bdeni(km), nfix(km)
#  if defined O_npzd_nitrogen_15
      real din15_in(km), din15flag, eps_bdeni
      real rn15impo, rn15expo, uno3, rno3, bwcdeni, bbdeni
#  endif
# endif
# if defined O_npzd_fe_limitation
      real felimit_D_in(km), felimit_in(km)
# endif
# if defined O_sed
      real sedcorgflx, sedcalflx
# endif


      expo = 0.0
      impo = 0.0
      phin = 0.0                ! integrated phytoplankton
      prca = 0.0                ! integrated production of calcite

      sedrr = 0.0
# if defined O_save_npzd
      rsedrr = 0.0
      rprca = 0.0
# endif
# if defined O_npzd_nitrogen_15
      rn15impo = 0.0
      rn15expo = 0.0
# endif
# if defined O_carbon_13
      rc13impo = 0.0
      rc13expo = 0.0
      prca13 = 0.0
# endif

!1111 start of main k-loop
      do k=1,kmx
# if defined O_npzd_nitrogen_15
         rn15impo = rn15expo
# endif
# if defined O_carbon_13
!        c13 biological fractionation 
         atmpres = 1.0          !atm
         depth = zt(k)/100.     !m
         call co2calc_SWS (t_in(k), s_in(k), dic_in(k), alk_in(k)
     &,                    co2_in, atmpres, depth, pH, co2star, dco2star
     &,                    pCO2, dpco2
     &,                    CO3, Omega_c, Omega_a)

         ac13_DIC_aq = -1.0512994e-4*t_in(k)+1.011765
!        Popp et al. (1989) Am. J. Sci.
         ac13_aq_POC = -0.017*log10(co2star*1000.)+1.0034
         ac13b = ac13_aq_POC/ac13_DIC_aq

         rc13impo = rc13expo
# endif
         swr = swr*exp(-kc*phin)

         phin = phin + max(tnpzd(k,3), trcmin)*dzt(k)
# if defined O_npzd_nitrogen
     &                + max(tnpzd(k,9), trcmin)*dzt(k)
# endif
         gl = swr*exp(ztt(k)*rctheta)
         impo = expo*dztr(k)
         bct = bbio**(cbio*t_in(k))
         bctz = (0.5*(tanh(o2_in(k) - 8.)+1))
     &           *bbio**(cbio*min(t_in(k),20.))
# if defined O_npzd_o2
!        decrease remineralisation rate in oxygen minimum zone
         nud = nud0*(0.75+0.25*tanh(o2_in(k) - 3.0))
# else
         nud = nud0
# endif

!-----------------------------------------------------------------------
!        call the npzd ecosystem dynamics model
!-----------------------------------------------------------------------

         call npzd_src (tnpzd(k,:), gl, bct, impo, dzt(k)
     &,                 dayfrac, wd(k), rkwz(k), nud, snpzd, expo
     &,                 graz, morp, morz, graz_Det
     &,                 graz_Z
# if defined O_save_npzd
     &,                 npp, morpt, remi, excr
#  if defined O_npzd_nitrogen
     &,                 npp_D, graz_D, morp_D, morpt_D
     &,                 nfix(k)
#  endif
#  if defined O_npzd_extra_diagnostics
     &,                 avej, avej_D, gmax, no3P
     &,                 po4P, po4_D
#  endif
# endif
# if defined O_npzd_fe_limitation
     &,                 felimit_in(k), felimit_D_in(k)
# endif
     &,                 bctz
# if defined O_npzd_nitrogen_15
     &,                 rn15impo, rn15expo
# endif
# if defined O_carbon_13
     &,                 rc13impo, rc13expo, ac13b
# endif
     &                      )
! These are source/sink terms 
         snpzd = snpzd*rdtts
         expo = expo*rnbio
# if defined O_npzd_nitrogen_15
         rn15expo = rn15expo*rnbio
# endif
# if defined O_carbon_13
         rc13expo = rc13expo*rnbio
# endif
# if defined O_save_npzd
         rexpo(k) = expo
         rgraz(k) = graz*rnbio 
         rgraz_Det(k) = graz_Det*rnbio
         rgraz_Z(k) = graz_Z*rnbio
         rmorp(k) = morp*rnbio
         rmorz(k) = morz*rnbio
         rnpp(k) = npp*rnbio
         rmorpt(k) = morpt*rnbio
         rremi(k) = remi*rnbio
         rexcr(k) = excr*rnbio
#  if defined O_npzd_nitrogen
         rnpp_D(k) = npp_D*rnbio
         rgraz_D(k) = graz_D*rnbio
         rmorp_D(k) = morp_D*rnbio
         rmorpt_D(k) = morpt_D*rnbio
         rnfix(k) = nfix(k)*rnbio
#  endif
#  if defined O_npzd_extra_diagnostics
         ravej(k) = avej*rnbio
         ravej_D(k) = avej_D*rnbio
         rgmax(k) = gmax*rnbio
         rno3P(k) = no3P*rnbio
         rpo4P(k) = po4P*rnbio
         rpo4_D(k) = po4_D*rnbio
#  endif
# endif
!-----------------------------------------------------------------------
!             calculate detritus at the bottom and remineralize
!-----------------------------------------------------------------------

         if (k .eq. kmx) then
            sedrr = expo*dzt(k)
# if defined O_sed
            sedcorgflx = sedcorgflx + expo*redctn*twodt*dzt(k)
# endif
# if defined O_save_npzd
            rremi(k) = rremi(k) + expo
            rsedrr = rsedrr + sedrr
# endif
            snpzd(1) = snpzd(1) + redptn*expo
            snpzd(6) = snpzd(6) + redctn*expo
# if defined O_carbon_13
            snpzd(16) = snpzd(16) + rc13expo*redctn*expo
# endif
# if defined O_npzd_nitrogen
!-----------------------------------------------------------------------
!           benthic denitrification model of Bohlen et al., 2012, GBC (eq. 10)
!           NO3 is removed out of bottom water nitrate.
!           See Somes et al., 2012, BGS for additional details/results
!-----------------------------------------------------------------------
!           limit denitrification as nitrate approaches 0 uM
            no3flag = 0.5+sign(0.5,no3_in(k)-trcmin)
#  if defined O_npzd_nitrogen_15
            din15flag = 0.5+sign(0.5,din15_in(k)-trcmin)
#  endif
            lno3 = 0.5*tanh(no3_in(k)*10 - 5.0)
                
            bdeni(k) = (0.06 + 0.19*0.99**(max(o2_in(k),trcmin)
     &               - max(no3_in(k),trcmin)))*max(expo,trcmin)*redctn
     &                 *1.e3

            bdeni(k) = min(bdeni(k), expo)
            bdeni(k) = max(bdeni(k), 0.)

            bdeni(k) = bdeni(k)*(0.5 + lno3)*no3flag
#  if defined O_npzd_nitrogen_15
     &                *din15flag
#  endif

            snpzd(7) = snpzd(7) + expo - bdeni(k)

#  if defined O_npzd_nitrogen_15
            rno3 = max(din15_in(k),trcmin*rn15std/(1+rn15std))
     &            /max(no3_in(k)-din15_in(k),trcmin*rn15std/(1+rn15std))
            rno3 = min(rno3, 2.*rn15std)
            rno3 = max(rno3, rn15std/2.)
            eps_bdeni = 6.0*exp(-2.5e-6*(zt(k)+dzt(k)/2.))
            bbdeni = rno3 - eps_bdeni*rno3/1000.

            snpzd(10) = snpzd(10) + rn15expo*expo 
     &               - bbdeni/(1+bbdeni)*bdeni(k)
#  endif

            expo = 0.
            sg_bdeni = 0.
         else
            if (sgb_in(k) .gt. 0) then
               sedrr = sgb_in(k)*expo*dzt(k)

               no3flag = 0.5+sign(0.5,no3_in(k)-trcmin)
#  if defined O_npzd_nitrogen_15
               din15flag = 0.5+sign(0.5,din15_in(k)-trcmin)
#  endif
               lno3 = 0.5*tanh(no3_in(k)*10 - 5.0)

               sg_bdeni = (0.06 + 0.19*0.99**(max(o2_in(k),trcmin)
     &                 - max(no3_in(k),trcmin)))
     &                   *max(expo*sgb_in(k),trcmin)*redctn*1.e3

               sg_bdeni = min(sg_bdeni, sgb_in(k)*expo)
               sg_bdeni = max(sg_bdeni, 0.)

               if (zt(k) .le. 200000.) sg_bdeni = sgbdfac*sg_bdeni

               sg_bdeni = sg_bdeni*(0.5 + lno3)*no3flag
#  if defined O_npzd_nitrogen_15
     &                   *din15flag
#  endif

               snpzd(1) = snpzd(1) + sgb_in(k)*expo*redptn
               snpzd(6) = snpzd(6) + sgb_in(k)*expo*redctn
               snpzd(7) = snpzd(7) + sgb_in(k)*expo - sg_bdeni
#  if defined O_carbon_13
               snpzd(16) = snpzd(16) + rc13expo*sgb_in(k)*expo*redctn
#  endif
#  if defined O_npzd_nitrogen_15
               rno3 = max(din15_in(k), trcmin*rn15std/(1+rn15std))
     &            /max(no3_in(k)-din15_in(k),trcmin*rn15std/(1+rn15std))
               rno3 = min(rno3, 2.*rn15std)
               rno3 = max(rno3, rn15std/2.)
               eps_bdeni = 6.0*exp(-2.5e-6*(zt(k)))
               bbdeni = rno3 - eps_bdeni*rno3/1000.

               snpzd(10) = snpzd(10) + rn15expo*sgb_in(k)*expo 
     &                     - bbdeni/(1+bbdeni)*sg_bdeni
#  endif

               expo = expo - sgb_in(k)*expo
#  if defined O_save_npzd
               rremi(k) = rremi(k) + sgb_in(k)*expo
               rsedrr = rsedrr + sedrr
#  endif
            else
               sg_bdeni = 0.
            endif
            bdeni(k) = sg_bdeni
# endif
         endif
# if defined O_npzd_nitrogen && defined O_save_npzd
         rbdeni(k) = bdeni(k)
# endif

!-----------------------------------------------------------------------
!             set source/sink terms
!-----------------------------------------------------------------------

         src(k,ispo4) = snpzd(1)
         src(k,isdop) = snpzd(2)
         src(k,isphyt) = snpzd(3)
         src(k,iszoop) = snpzd(4)
         src(k,isdetr) = snpzd(5)
# if defined O_npzd_nitrogen
         src(k,isdic) = snpzd(6)
         src(k,isno3) = snpzd(7)
         src(k,isdon) = snpzd(8)
         src(k,isdiaz) = snpzd(9)
#  if defined O_npzd_nitrogen_15
         src(k,isdin15) = snpzd(10)
         src(k,isdon15) = snpzd(11)
         src(k,isphytn15) = snpzd(12)
         src(k,iszoopn15) = snpzd(13)
         src(k,isdetrn15) = snpzd(14)
         src(k,isdiazn15) = snpzd(15)
#  endif
# endif
# if defined O_carbon_13
         src(k,isdic13) = snpzd(16)
         src(k,isdoc13) = snpzd(17)
         src(k,isphytc13) = snpzd(18)
         src(k,iszoopc13) = snpzd(19)
         src(k,isdetrc13) = snpzd(20)
#  if defined O_npzd_nitrogen
         src(k,isdiazc13) = snpzd(21)
#  endif
# endif
!        production of calcite
         dprca = (morp+morz+(graz+graz_Z)*(1.-gamma1))
     &              *capr*redctn*rnbio
         prca = prca + dprca*dzt(k)

!        These are sources and sinks of DIC (i.e. remin - pp)
!        all are based on po4 uptake and remineralization
!        dprca is a correction term
# if defined O_carbon
         dic_npzd_sms(k) = snpzd(6)
         src(k,isdic) = (snpzd(6) - dprca)
#  if defined O_carbon_13
         rtdic13 = max(dic13_in(k),trcmin*rc13std/(1+rc13std))
     &             / max(dic_in(k),trcmin)
         rtdic13 = min(rtdic13, 2.*rc13std/(1+rc13std))
         rtdic13 = max(rtdic13, 0.5*rc13std/(1+rc13std))

         prca13 = prca13 + dprca*dzt(k)*rtdic13

         src(k,isdic13) = src(k,isdic13) - rtdic13*dprca
#  endif
# endif
# if defined O_npzd_alk
         src(k,isalk) = (-snpzd(6)*redntc*1.e-3 - 2.*dprca)
# endif
!        calculate total export to get total import for next layer
         expo = expo*dzt(k)

      enddo
!1111 end of main k-loop

# if defined O_save_npzd
      rprca = prca
# endif

# if defined O_npzd_o2
!2222 start of second k-loop
      do k=1,kmx
!        limit oxygen consumption below concentrations of
!        5umol/kg as recommended in OCMIP
         fo2 = 0.5*tanh(o2_in(k)*1000. - 2.5)
!        sink of oxygen
! O2 is needed to generate the equivalent of NO3 from N2 during N2 fixation
! 0.5 H2O + 0.5 N2+1.25O2 -> HNO3
! note that so2 is -dO2/dt
         so2 = dic_npzd_sms(k)*redotc + nfix(k)*rnbio*1.25e-3
c         so2 = dic_npzd_sms(k)*redotc

         src(k,iso2) = -so2*(0.5 + fo2)
#  if defined O_npzd_nitrogen
!        add denitrification as source term for NO3
         no3flag = 0.5+sign(0.5,no3_in(k)-trcmin)
#   if defined O_npzd_nitrogen_15
         din15flag = 0.5+sign(0.5,din15_in(k)-trcmin)
#   endif
         lno3 = 0.5*tanh(no3_in(k) - 2.5)
         lntp = 0.5*tanh(no3_in(k)/(redntp*po4_in(k))*100 - 60.)
!        800 = 0.8*1000 = (elec/mol O2)/(elec/mol NO3)*(mmol/mol)
         wcdeni = 800.*no3flag*so2*(0.5 - fo2)*(0.5 + lno3)
#   if defined O_npzd_nitrogen_15
     &           *din15flag
#   endif
         wcdeni = max(wcdeni, 0.)
         src(k,isno3) = src(k,isno3) - wcdeni
#   if defined O_npzd_nitrogen_15
!        calculate isotope effect of water column denitrification
         uno3 = wcdeni*twodt/no3_in(k)
         uno3 = min(uno3, 0.999)
         uno3 = max(uno3, trcmin)
         rno3 = max(din15_in(k),trcmin*rn15std/(1+rn15std))
     &             / max(no3_in(k)-din15_in(k),
     &                   trcmin*rn15std/(1+rn15std))
         rno3 = min(rno3, 2.*rn15std)
         rno3 = max(rno3, rn15std/2.)
         bwcdeni = rno3 + 25.*(1-uno3)/uno3*log(1-uno3)*rno3/1000.

         src(k,isdin15) = src(k,isdin15) - (bwcdeni/(1+bwcdeni))*wcdeni
#   endif
#   if defined O_deni_alk
! Correct the ALK stoichiometry to account for N2 fixation
         src(k,isalk) = src(k,isalk) + wcdeni*1.e-3
!bdeni
         src(k,isalk) = src(k,isalk) + bdeni(k)*1.e-3
! Now add account for N2 fixation (ALK production is tied to PO4 change and
! thus in the case of N2 fixation is not correct).
         src(k,isalk) = src(k,isalk) - nfix(k)*rnbio*1.e-3
#   endif
#   if defined O_save_npzd
         rwcdeni(k) = wcdeni
#   endif
#  endif
      enddo
!2222 end of second k-loop
# endif

!-----------------------------------------------------------------------
!           remineralize calcite
!-----------------------------------------------------------------------
!3333 start of third k-loop
      do k=1,kmx-1
# if defined O_carbon
         src(k,isdic) = src(k,isdic) + prca*rcak(k)
#  if defined O_carbon_13
         src(k,isdic13) = src(k,isdic13) + prca13*rcak(k)
#  endif
# endif
# if defined O_npzd_alk
         src(k,isalk) = src(k,isalk) + 2.*prca*rcak(k)
# endif
      enddo
!3333 end of third k-loop
# if defined O_sed
      sedcalflx = sedcalflx + (prca*rcab(kmx)
     &                       - prca*rcak(kmx))*twodt*dzt(kmx)
# endif
# if defined O_carbon
      src(kmx,isdic) = src(kmx,isdic) + prca*rcab(kmx)
#  if defined O_carbon_13
      src(kmx,isdic13) = src(kmx,isdic13) + prca13*rcab(kmx)
#  endif
# endif
# if defined O_npzd_alk
      src(kmx,isalk) = src(kmx,isalk) + 2.*prca*rcab(kmx)
# endif
      return
      end
!     END mobi

      subroutine npzd_src (bioin, gl, bct, impo, dzt
     &,                    dayfrac, wwd, rkw, nud, bioout, expoout
     &,                    grazout, morpout, morzout, graz_Det_out
     &,                    graz_Zout
# if defined O_save_npzd
     &,                    nppout, morptout, remiout, excrout
#  if defined O_npzd_nitrogen
     &,                    npp_Dout, graz_Dout, morp_Dout, morpt_Dout
     &,                    nfixout
#  endif
#  if defined O_npzd_extra_diagnostics
     &,                    avej_out, avej_D_out, gmax_out, no3P_out
     &,                    po4P_out, po4_D_out
#  endif
# endif
# if defined O_npzd_fe_limitation
     &,                    felimit, felimit_D
# endif
     &,                    bctz
# if defined O_npzd_nitrogen_15
     &,                    rn15impo, rn15expoout
# endif
# if defined O_carbon_13
     &,                    rc13impo, rc13expoout, ac13b
# endif
     &                     ) 

!=======================================================================
!     computes source terms of the NPZD model
!     initial version of code adapted from Xavier Giraud:
!     Giraud et al. 2000, J Mar Res, 58, 609-630
!     original model reference:
!     Oeschlies and Garcon 1999, Global Biogeochem. Cycles 13, 135-160
!     Schmittner et al. 2005,  Global Biogeochem. Cycles 19, GB3004,
!     doi:10.1029/2004GB002283.
!     Schmittner et al. 2008, Global Biogeochem. Cycles 22, GB1013
!
!     This version was modified by David Keller and corrects the zooplankton
!     grazing formulation.  Note that zooplankton are now allowed to graze
!     on themselves and detritus, in addition to phyt. and diazotrophs.
!     The calculation of light has also been corrected.
!     Keller et al. 2012, Geosci. Model Dev., 5, 1195-1220
!
!     The nitrogen isotope model (nitrogen_15) has been developed by 
!     Chris Somes and is described in
!     Somes et al. 2010, Global Biogeochem. Cycles 24, GB4019
!     Somes et al. 2010, Geophys. Res. Lett. 37, L23605 and
!     Somes et al. 2013, Biogeosc. 10, 5889-5910
!
!     The carbon isotope model (carbon_13) has been written by Andreas 
!     Schmittner and is documented in
!     Schmittner et al. 2013 Biogeosc. 10, 5793-5816
!     Chris Somes modified the code converting from alpha to beta formulation
!
!     Note that nutrient (N) represents phosphate
!
!     input variables:
!
!       bioin(1:5) = N,DON,P,Z,D [mmol m-3]
!       bioin(6)   = nitrate [mmol m-3]
!       bioin(7)   = DON [mmol m-3]
!       bioin(8)   = diazotrophs [mmol m-3]
!
!       gl         = 2.*light at top of grid box
!       nbio       = number of time steps
!       dtbio        = time step [s]
!       bct        = bbio**(cbio*temperature)
!       impo       = import of detritus from above [mmol m-3]
!       dzt        = depth of grid box [cm]
!       dayfrac    = day length (fraction: 0 < dayfrac < 1)
!       wwd        = sinking speed of detritus/dzt
!       rkw        = reciprical of kw*dzt(k)
!       nud        = remineralisation rate of detritus [s-1]
!
!     output variables:
!
!       bioout     = change from bioin [mmol m-3]
!       nppout     = net primary production [mmol m-3]
!       grazout    = grazing [mmol m-3]
!       morpout    = quadratic mortality of phytoplankton [mmol m-3]
!       morptout   = specific mortality of phytoplankton [mmol m-3]
!       morzout    = mortality of zooplankton [mmol m-3]
!       remiout    = remineralisation [mmol m-3]
!       excrout    = excretion [mmol m-3]
!       expoout    = detrital export [mmol m-3]
!       npp_Dout   = NPP of diazotrophs
!       graz_Dout  = grazing of diazotrophs
!       morp_Dout  = mortality of diazotrophs
!       nfixout    = rate of N2 fixation
!       graz_Det_out = grazing of detritus
!       graz_Zout   = grazing on othe zooplankton
!       avej_out    = light-depend phyt. growth rate
!       avej_D_out  = light-depend Diaz growth rate
!       gmax_out    = temp-depend. zoo growth rate
!       no3P_out    = no3 depend. phyt growth rate
!       po4P_out    = po4 depend. phyt growth rate
!       po4_D_out   = po4 depend. Diaz growth rate
!
!       New grazing formulation variables and parameters
!
!       The following terms determine ingestion according to a
!       a Holling II curve (i.e. Michaelis Menten):
!
!       Ingestion = max_graz_rate * (Ft/(Ft + kzoo))
!
!       where Ft is the weighted measure of the total food available
!       and equals the sum of the different prey types times the
!       preference of Z for that type of prey 
!
!       zprefP   = Z preference for P
!       zprefD   = Z preference for Diaz
!       zprefDet = Z preference for detritus
!       zprefZ   = Z preference for other Z
!       kzoo = half saturation coefficienct for Z ingestion mmol N m-3
!       ing_P    = zooplankton ingestion of phytoplankon
!       ing_D    = zooplankton ingestion of diazotrophs
!       ing_Det  = zooplankton ingestion of detritus
!       ing_Z    = zooplankton ingestion of other zooplankton
!       thetaZ   = Michaelis-Menten denominator
!
!       felimit = Fe limitation parameter
!       felmit_D = Fe limitation parameter for diazotrophs
!
!=======================================================================

      implicit none

      integer n

      real gl, f1, biopo4, biophyt, biozoop, biodetr, jmax, u_P, g_P
      real morp, morpt, morz, remi, excr, expo, impo, nppout, grazout
      real morpout, morptout, morzout, remiout, excrout, expoout
      real avej_out, avej_D_out, gmax_out, no3P_out, po4P_out, po4_D_out
      real dzt, po4flag, phytflag, zoopflag, detrflag, wwd, rkw, gd
      real nupt, nud, biono3, u_D,npp_D, npp_Dout, no3flag, biodiaz, bct
      real diazflag, g_D,graz_D, morpt_D, jmax_D, gd_D, avej_D, no3upt_D
      real morpt_Dout, graz_Dout, nfixout, biop2, u1, u2, phi1, phi2
      real avej, graz_Det_out, graz_Zout, thetaZ, ing_P, ing_D, dayfrac
      real ing_Det, ing_Z, g_Z, g_Det, graz_Z, graz_Det, gmax
      real no3P, po4P, po4_D, felimit, bctz, felimit_D
      real excrdiaz, biodop, biodon, dopflag, donflag, recy
      real recy_p, npp, graz, biodin15, biodon15, biophytn15, biozoopn15
      real biodetrn15, biodiazn15, bassim, fcassim, bexcr, fcexcr
      real bnfix, fcnfix, rn15impo, rn15expoout, dig, dig_P, dig_Z
      real dig_Det, dig_D, excr_P, excr_Z, excr_Det, excr_D, sf, sf_P
      real sf_Z, sf_Det, sf_D, nr_excr_D, uno3, rno3, rzoop
      real rtdin15, rtphytn15, rtzoopn15, rtdetrn15, rtdiazn15, rtdon15
      real din15flag
      real don15flag, phytn15flag, zoopn15flag, detrn15flag, diazn15flag
      real limP, limP_D, uptdop_D, dopupt_D, morp_D, nupt_D, morp_Dout
      real udon, rdon, brecy, fcrecy, biodic, uptdop, dopupt
      real biodic13, biophytc13, biozoopc13, biodetrc13, biodiazc13
      real dic13flag, phytc13flag, zoopc13flag, detrc13flag, diazc13flag
      real biodoc13, doc13flag, biodoc
      real rdic13, rtphytc13, rtzoopc13, rtdetrc13, rtdiazc13, ac13b
      real rc13impo, rc13expoout, dicflag, bc13npp, rtdic13, fcnpp
      real rtdoc13

#include "size.h"
#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "calendar.h"
#include "npzd.h"

      real bioin(ntnpzd), bioout(ntnpzd)

      bioout(:) = 0.0
      biopo4 = bioin(1)
      biodop = bioin(2)
      biophyt = bioin(3)
      biozoop = bioin(4)
      biodetr = bioin(5)
# if defined O_npzd_nitrogen
      biodic = bioin(6)
      biono3 = bioin(7)
      biodon = bioin(8)
      biodiaz = bioin(9)
#  if defined O_npzd_nitrogen_15
      biodin15 = bioin(10)
      biodon15 = bioin(11)
      biophytn15 = bioin(12)
      biozoopn15 = bioin(13)
      biodetrn15 = bioin(14)
      biodiazn15 = bioin(15)
#  endif
# endif
# if defined O_carbon_13
      biodic13 = bioin(16)
      biodoc13 = bioin(17)
      biophytc13 = bioin(18)
      biozoopc13 = bioin(19)
      biodetrc13 = bioin(20)
#  if defined O_npzd_nitrogen
      biodiazc13 = bioin(21)
#  endif
# endif
!       flags prevent negative values by setting outgoing fluxes to
!       zero if tracers are lower than trcmin
      po4flag = 0.5 + sign(0.5,biopo4 - trcmin)
      dopflag = 0.5 + sign(0.5,biodop - trcmin)
      phytflag = 0.5 + sign(0.5,biophyt - trcmin)
      zoopflag = 0.5 + sign(0.5,biozoop - trcmin)
      detrflag = 0.5 + sign(0.5,biodetr - trcmin)
# if defined O_npzd_nitrogen
      no3flag = 0.5 + sign(0.5,biono3 - trcmin)
      donflag = 0.5 + sign(0.5,biodon - trcmin)
      diazflag = 0.5 + sign(0.5,biodiaz - trcmin)
#  if defined O_npzd_nitrogen_15
      din15flag = 0.5 + sign(0.5, biodin15 - trcmin)
      don15flag = 0.5 + sign(0.5, biodon15 - trcmin)
      phytn15flag = 0.5 + sign(0.5, biophytn15 - trcmin)
      zoopn15flag = 0.5 + sign(0.5, biozoopn15 - trcmin)
      detrn15flag = 0.5 + sign(0.5, biodetrn15 - trcmin)
      diazn15flag = 0.5 + sign(0.5, biodiazn15 - trcmin)
#  endif
# endif
# if defined O_carbon_13
      dic13flag = 0.5 + sign(0.5,biodic13 - trcmin)
      doc13flag = 0.5 + sign(0.5,biodoc13 - trcmin)
      phytc13flag = 0.5 + sign(0.5,biophytc13 - trcmin)
      zoopc13flag = 0.5 + sign(0.5,biozoopc13 - trcmin)
      detrc13flag = 0.5 + sign(0.5,biodetrc13 - trcmin)
      diazc13flag = 0.5 + sign(0.5,biodiazc13 - trcmin)
# endif

!     limit tracers to positive values
      bioin(:) = max(bioin(:), trcmin)
      biopo4 = max(biopo4, trcmin)
      biodop = max(biodop, trcmin)
      biophyt = max(biophyt, trcmin)
      biozoop = max(biozoop, trcmin)
      biodetr = max(biodetr, trcmin)
# if defined O_npzd_nitrogen
      biodic = max(biodic, trcmin)
      biono3 = max(biono3, trcmin)
      biodon = max(biodon, trcmin)
      biodiaz = max(biodiaz, trcmin)
#  if defined O_npzd_nitrogen_15
      biodin15 = max(biodin15, trcmin)
      biodon15 = max(biodon15, trcmin)
      biophytn15 = max(biophytn15, trcmin)
      biozoopn15 = max(biozoopn15, trcmin)
      biodetrn15 = max(biodetrn15, trcmin)
      biodiazn15 = max(biodiazn15, trcmin)
#  endif
# endif
# if defined O_carbon_13
      biodic13 = max(biodic13, trcmin)
      biodoc13 = max(biodoc13, trcmin)
      biophytc13 = max(biophytc13, trcmin)
      biozoopc13 = max(biozoopc13, trcmin)
      biodetrc13 = max(biodetrc13, trcmin)
      biodiazc13 = max(biodiazc13, trcmin)
# endif
!     photosynthesis after Evans & Parslow (1985)
!     notation as in JGOFS report No. 23 p. 6
      f1 = exp((-kw - kc*(biophyt+biodiaz))*dzt)
# if defined O_npzd_fe_limitation
! In the following "felimit" is determined by an iron mask and
! is used to limit phytoplankton growth in HNLC regions
      jmax = abio*bct*felimit
# else
      jmax = abio*bct
# endif
      gd = jmax*dayfrac
      u1 = max(gl/gd,1.e-6)
      u2 = u1*f1
!     for the following approximation ensure that u1 < 20
      phi1 = log(u1+sqrt(1.+u1**2.))-(sqrt(1.+u1**2.)-1.)/u1
      phi2 = log(u2+sqrt(1.+u2**2.))-(sqrt(1.+u2**2.)-1.)/u2
      avej = gd*(phi1 - phi2)/((kw+kc*(biophyt+biodiaz))*dzt)
! Make the max grazing rate a function of temperature
      gmax = gbio*bctz
! Note bctz, which sets an upper limit on the effects of temp on the
! grazing rate, is set in tracers.F
# if defined O_npzd_nitrogen
#  if defined O_npzd_fe_limitation
      jmax_D = max(0.,abio*bct*felimit_D)*jdiar
#  else
      jmax_D = max(0.,abio*bct)*jdiar
#  endif
      gd_D = max(1.e-14,jmax_D*dayfrac)
      u1 = max(gl/gd_D,1.e-6)
      u2 = u1*f1
!     for the following approximation ensure that u1 < 20
      phi1 = log(u1+sqrt(1.+u1**2.))-(sqrt(1.+u1**2.)-1.)/u1
      phi2 = log(u2+sqrt(1.+u2**2.))-(sqrt(1.+u2**2.)-1.)/u2
      avej_D = gd_D*(phi1 - phi2)/((kw+kc*(biophyt+biodiaz))*dzt)
# endif
# if defined O_npzd_nitrogen
!     check grazing preferences = 1 for N case
        IF ((zprefP + zprefDet + zprefZ + zprefD).ne.1) THEN
           zprefP = 0.3
           zprefZ = 0.3
           zprefDet = 0.3
           zprefD = 0.1
        END IF
# else
!      check that grazing preferences = 1 for no N case
        IF ((zprefP + zprefDet + zprefZ).ne.1) THEN
           zprefP = 0.35
           zprefZ = 0.35
           zprefDet = 0.30
        END IF
# endif
      nupt = nupt0*bct
      nupt_D = nupt0_D*bct

      expoout = 0.0
      grazout = 0.0
      morpout = 0.0
      morzout = 0.0
      graz_Det_out = 0.0
      graz_Zout = 0.0
# if defined O_npzd_nitrogen_15
      rn15expoout = 0.0
# endif
# if defined O_carbon_13
      rc13expoout = 0.0
# endif
# if defined O_save_npzd
      nppout = 0.0
      morptout = 0.0
      remiout = 0.0
      excrout = 0.0
#  if defined O_npzd_nitrogen
      npp_Dout = 0.0
      graz_Dout = 0.0
      morpt_Dout = 0.0
      morp_Dout = 0.0
      nfixout = 0.0
#  endif
#  if defined O_npzd_extra_diagnostics
      avej_out = 0.0
      avej_D_out = 0.0
      gmax_out = 0.0
      no3P_out = 0.0
      po4P_out = 0.0
      po4_D_out = 0.0
#  endif
# endif

      do n=1,nbio

!       growth rate of phytoplankton
        if (0.4*biodop/(k1p + biodop) .gt. biopo4/(k1p + biopo4)) then
           uptdop = 1.
           limP = biodop/(k1p + biodop)
c Chris: is this correct "_D" ?
c           u_P = min(avej_D, jmax_D*limP*0.4)
           u_P = min(avej, jmax*limP*0.4)
        else
           uptdop = 0.
           limP = biopo4/(k1p + biopo4)
           u_P = min(avej, jmax*limP)
        endif
        po4P = jmax*biopo4/(k1p + biopo4) 

# if defined O_npzd_nitrogen
!       nitrate limitation
        u_P = min(u_P, jmax*biono3/(k1n + biono3))
        no3P = jmax*biono3/(k1n + biono3)
!       growth rate of diazotrophs smaller than other phytoplankton and
!       not nitrate limited
        if (0.4*biodop/(k1p + biodop) .gt. biopo4/(k1p + biopo4)) then
           uptdop_D = 1.
           limP_D = biodop/(k1p + biodop)
           u_D = min(avej_D, jmax_D*limp_D*0.4)
        else
           uptdop_D = 0.
           limP_D = biopo4/(k1p + biopo4)
           u_D = min(avej_D, jmax_D*limp_D)
        endif
        po4_D = jmax_D*biopo4/(k1p + biopo4)
!       Set the grazing coefficients for the N case
        thetaZ = zprefP*biophyt + zprefDet*biodetr + zprefZ*biozoop
     &       + zprefD*biodiaz + kzoo
        ing_P = zprefP/thetaZ
        ing_Det = zprefDet/thetaZ
        ing_Z = zprefZ/thetaZ
        ing_D = zprefD/thetaZ
# else
!       If "else" then set the grazing coefficients for the no N case
!       note kzoo is in terms of N so convert to P
        thetaZ = zprefP*biophyt + zprefDet*biodetr + zprefZ*biozoop 
     &       +kzoo*redptn
        ing_P = zprefP/thetaZ
        ing_Det = zprefDet/thetaZ
        ing_Z = zprefZ/thetaZ
# endif
        npp = u_P*biophyt
        dopupt = npp*uptdop
# if defined O_npzd_nitrogen
        npp_D = max(0.,u_D*biodiaz)
!       grazing on diazotrophs
        g_D = gmax*ing_D*biodiaz
        graz_D = g_D*biozoop
        morpt_D = nupt_D*biodiaz ! linear mortality
        morp_D = nup_D*biodiaz
        no3upt_D = biono3/(k1n + biono3)*npp_D ! nitrate uptake
        dopupt_D = npp_D*uptdop_D
# endif
!       grazing on P
        g_P = gmax*ing_P*biophyt
        graz = g_P*biozoop
!       grazing on Z
        g_Z = gmax*ing_Z*biozoop
        graz_Z = g_Z*biozoop
!       grazing on Detritus
        g_Det = gmax*ing_Det*biodetr
        graz_Det = g_Det*biozoop
!
        morp = nup*biophyt
        morpt = nupt*biophyt
        recy = nupt/3000.*biodon
        recy_p = nupt/1000.*biodop
        morz = nuz*biozoop*biozoop
        remi = nud*bct*biodetr
        expo = wwd*biodetr

        graz = graz*phytflag*zoopflag
# if defined O_npzd_nitrogen_15
     &         *phytn15flag*zoopn15flag
# endif
        graz_Z = graz_Z*zoopflag
# if defined O_npzd_nitrogen_15
     &          *zoopn15flag
# endif
        graz_Det = graz_Det*detrflag*zoopflag
# if defined O_npzd_nitrogen_15
     &            *detrn15flag*zoopn15flag
# endif
        morp = morp*phytflag
# if defined O_npzd_nitrogen_15
     &        *phytn15flag
# endif
        morpt = morpt*phytflag
# if defined O_npzd_nitrogen_15
     &         *phytn15flag
# endif
        morz = morz*zoopflag
# if defined O_npzd_nitrogen_15
     &        *zoopn15flag
# endif
        remi = remi*detrflag
# if defined O_npzd_nitrogen_15
     &        *detrn15flag
# endif
        expo = expo*detrflag
# if defined O_npzd_nitrogen_15
     &        *detrn15flag
# endif
        recy_p = recy_p*dopflag
# if defined O_npzd_nitrogen
        if (uptdop .eq. 1) then
           npp = npp*dopflag*no3flag
        else
           npp = npp*po4flag*no3flag
        endif
#  if defined O_npzd_nitrogen_15
        npp = npp*din15flag
#  endif
        if (uptdop_D .eq. 1) then
           npp_D = npp_D*dopflag
        else
           npp_D = npp_D*po4flag
        endif
        graz_D = graz_D*diazflag*zoopflag
#  if defined O_npzd_nitrogen_15
     &          *diazn15flag*zoopn15flag
#  endif
        morpt_D = morpt_D*diazflag
#  if defined O_npzd_nitrogen_15
     &           *diazn15flag
#  endif
        morp_D = morp_D*diazflag
#  if defined O_npzd_nitrogen_15
     &          *diazn15flag
#  endif
        no3upt_D = no3upt_D*no3flag
#  if defined O_npzd_nitrogen_15
     &            *din15flag
#  endif
        recy = recy*donflag
#  if defined O_npzd_nitrogen_15
     &        *don15flag
#  endif
# else
        npp = npp*po4flag
# endif

!   Excretion is the difference between
!   the zooplankton assimilation and growth efficiencies
# if defined O_npzd_nitrogen
c        excr_P = (gamma1-geZ)*graz
c        excr_Z = (gamma1-geZ)*graz_Z
c        excr_Det = (gamma1-geZ)*graz_Det
c        excr_D = (gamma1-geZ)*graz_D

        dig_P = gamma1*graz
        dig_Z = gamma1*graz_Z
        dig_Det = gamma1*graz_Det
        dig_D = gamma1*graz_D*(redntp/diazntp)
        dig = dig_P + dig_Z + dig_Det + dig_D

        excr_P = gamma1*(1-geZ)*graz
        excr_Z = gamma1*(1-geZ)*graz_Z
        excr_Det = gamma1*(1-geZ)*graz_Det
        excr_D = gamma1*(1-geZ)*graz_D*(redntp/diazntp)
        excr = excr_P + excr_Z + excr_Det + excr_D
        nr_excr_D = gamma1*graz_D*(1-(redntp/diazntp)) 
     &       + (1-gamma1)*graz_D*(1-(redntp/diazntp))

        sf_P = (1-gamma1)*graz
        sf_Z = (1-gamma1)*graz_Z
        sf_Det = (1-gamma1)*graz_Det
        sf_D = (1-gamma1)*graz_D*(redntp/diazntp)
        sf = sf_P + sf_Z + sf_Det + sf_D
# else
c        excr = (gamma1-geZ)*(graz+graz_Z+graz_Det)
# endif
# if defined O_npzd_nitrogen_15
!       calculate isotope parameters
!       See Somes et al., 2010, GBC for details/results
        uno3 = npp*dtbio/biono3
        uno3 = min(uno3, 0.999)
        uno3 = max(uno3, trcmin)
        rno3 = biodin15/(biono3-biodin15)
        rno3 = min(rno3, 2*rn15std)
        rno3 = max(rno3, rn15std/2.)
        bassim = rno3 + 6.*(1-uno3)/uno3*log(1-uno3)*rno3/1000.
        fcassim = bassim/(1+bassim)

        udon = recy*dtbio/biodon
        udon = min(udon, 0.999)
        udon = max(udon, trcmin)
        rdon = biodon15/(biodon-biodon15)
        rdon = min(rdon, 2*rn15std)
        rdon = max(rdon, rn15std/2.)
        brecy = rdon + 1.3*(1-udon)/udon*log(1-udon)*rdon/1000.
        fcrecy = brecy/(1+brecy)

        rzoop = biozoopn15/(biozoop-biozoopn15)
        rzoop = min(rzoop, 2.*rn15std)
        rzoop = max(rzoop, rn15std/2.)
        bexcr = rzoop - 4.*rzoop/1000.
        fcexcr = bexcr/(1+bexcr)

        bnfix = rn15std - 1.*rn15std/1000.
        fcnfix = bnfix/(1+bnfix)

        rtdin15 = biodin15/biono3
        rtdin15 = min(rtdin15, 2*rn15std/(1+rn15std))
        rtdin15 = max(rtdin15, rn15std/(1+rn15std)/2.)

        rtdon15 = biodon15/biodon
        rtdon15 = min(rtdon15, 2*rn15std/(1+rn15std))
        rtdon15 = max(rtdon15, rn15std/(1+rn15std)/2.)

        rtphytn15 = biophytn15/biophyt
        rtphytn15 = min(rtphytn15, 2.*rn15std/(1+rn15std))
        rtphytn15 = max(rtphytn15, rn15std/(1+rn15std)/2.)

        rtzoopn15 = biozoopn15/biozoop
        rtzoopn15 = min(rtzoopn15, 2.*rn15std/(1+rn15std))
        rtzoopn15 = max(rtzoopn15, rn15std/(1+rn15std)/2.)

        rtdetrn15 = biodetrn15/biodetr
        rtdetrn15 = min(rtdetrn15, 2.*rn15std/(1+rn15std))
        rtdetrn15 = max(rtdetrn15, rn15std/(1+rn15std)/2.)

        rtdiazn15 = biodiazn15/biodiaz
        rtdiazn15 = min(rtdiazn15, 2.*rn15std/(1+rn15std))
        rtdiazn15 = max(rtdiazn15, rn15std/(1+rn15std)/2.)
# endif
# if defined O_carbon_13
        rdic13 = biodic13/(biodic-biodic13)
        rdic13 = min(rdic13, 2.*rc13std)
        rdic13 = max(rdic13, 0.5*rc13std)
        bc13npp = ac13b*rdic13
        fcnpp = bc13npp/(1+bc13npp)
        
        rtdic13 = biodic13/biodic
        rtdic13 = min(rtdic13, 2*rc13std/(1+rc13std))
        rtdic13 = max(rtdic13, 0.5*rc13std/(1+rc13std))

        rtdoc13 = biodoc13/(biodon*redctn)
        rtdoc13 = min(rtdoc13, 2*rc13std/(1+rc13std))
        rtdoc13 = max(rtdoc13, 0.5*rc13std/(1+rc13std))

        rtphytc13 = biophytc13/(biophyt*redctn)
        rtphytc13 = min(rtphytc13, 2.*rc13std/(1+rc13std))
        rtphytc13 = max(rtphytc13, 0.5*rc13std/(1+rc13std))

        rtzoopc13 = biozoopc13/(biozoop*redctn)
        rtzoopc13 = min(rtzoopc13, 2.*rc13std/(1+rc13std))
        rtzoopc13 = max(rtzoopc13, 0.5*rc13std/(1+rc13std))

        rtdetrc13 = biodetrc13/(biodetr*redctn)
        rtdetrc13 = min(rtdetrc13, 2.*rc13std/(1+rc13std))
        rtdetrc13 = max(rtdetrc13, 0.5*rc13std/(1+rc13std))

        rtdiazc13 = biodiazc13/(biodiaz*redctn)
        rtdiazc13 = min(rtdiazc13, 2.*rc13std/(1+rc13std))
        rtdiazc13 = max(rtdiazc13, 0.5*rc13std/(1+rc13std))
# endif
# if defined O_npzd_nitrogen
!       nutrients equation
        biopo4 = biopo4 + dtbio*(redptn*((1.-dfrt)*morpt + excr + remi 
     &       - (npp-dopupt)) + diazptn*(morpt_D - (npp_D-dopupt_D)) 
     &       + recy_p)
!       DOP equation
        biodop = biodop + dtbio*(redptn*dfr*morp + redptn*dfrt*morpt 
     &       - redptn*dopupt - diazptn*dopupt_D - recy_p)
!       phytoplankton equation
        biophyt = biophyt + dtbio*(npp - morp - graz - morpt)
!       zooplankton equation
        biozoop = biozoop + dtbio*(dig - morz - graz_Z - excr)
!       detritus equation
        biodetr = biodetr + dtbio*((1.-dfr)*morp + sf + morz - remi 
     &       - graz_Det - expo + impo + morp_D*(redntp/diazntp))
!       DIC equation
        biodic = biodic + dtbio*redctn*((1.-dfrt)*morpt + excr + remi 
     &       - npp + morpt_D - npp_D + recy + nr_excr_D
     &       + morp_D*(1-(redntp/diazntp)))
!       nitrate (NO3) equation
        biono3 = biono3 + dtbio*((1.-dfrt)*morpt + excr + morpt_D  
     &       + morp_D*(1-(redntp/diazntp)) + nr_excr_D + remi + recy 
     &       - npp - no3upt_D)
!       DON equation
        biodon = biodon + dtbio*(dfr*morp + dfrt*morpt - recy)
!       diazotroph equation
        biodiaz = biodiaz + dtbio*(npp_D - morp_D - morpt_D - graz_D)
# else
!       nutrients equation
        biopo4 = biopo4 + dtbio*redptn*(remi + excr - npp + morpt)
!       phytoplankton equation
        biophyt = biophyt + dtbio*(npp - morp - graz - morpt)
!       zooplankton equation
        biozoop = biozoop + dtbio*(geZ*(graz + graz_Det + graz_Z) - morz
     &       - graz_Z)
!       detritus equation
        biodetr = biodetr + dtbio*((1.-gamma1)*(graz + graz_Det 
     &       + graz_Z) + morp + morz - remi - graz_Det - expo + impo)
# endif
# if defined O_npzd_nitrogen_15
!       isotope equations
        biodin15 = biodin15 + dtbio*(rtphytn15*(1.-dfrt)*morpt 
     &       + fcexcr*excr + rtdiazn15*morpt_D + rtdiazn15*nr_excr_D 
     &       + rtdiazn15*morp_D*(1-(redntp/diazntp)) + rtdetrn15*remi 
     &       + fcrecy*recy - fcassim*npp - fcassim*no3upt_D)

        biodon15 = biodon15 + dtbio*(dfr*rtphytn15*morp 
     &       + dfrt*rtphytn15*morpt - fcrecy*recy)

        biophytn15 = biophytn15 + dtbio*(fcassim*npp - rtphytn15*morp 
     &       - rtphytn15*graz - rtphytn15*morpt)

        biozoopn15 = biozoopn15 + dtbio*(rtphytn15*dig_P 
     &       + rtzoopn15*dig_Z + rtdetrn15*dig_Det + rtdiazn15*dig_D 
     &       - rtzoopn15*morz - rtzoopn15*graz_Z - fcexcr*excr)

        biodetrn15 = biodetrn15 + dtbio*(rtphytn15*(1.-dfr)*morp 
     &       + rtphytn15*sf_P + rtzoopn15*sf_Z + rtdetrn15*sf_Det 
     &       + rtdiazn15*sf_D + rtzoopn15*morz - rtdetrn15*remi 
     &       - rtdetrn15*graz_Det - rtdetrn15*expo + rn15impo*impo 
     &       + rtdiazn15*morp_D*(redntp/diazntp))

        biodiazn15 = biodiazn15 + dtbio*(fcnfix*(npp_D-no3upt_D)
     &       + fcassim*no3upt_D - rtdiazn15*morp_D - rtdiazn15*graz_D 
     &       - rtdiazn15*morpt_D)
# endif
# if defined O_carbon_13
! !!!!!!!!!!!!!!!!!! isotope equations
        biodic13 = biodic13 + dtbio*redctn*(rtphytc13*(1.-dfrt)*morpt 
     &       + rtzoopc13*excr + rtdiazc13*morpt_D + rtdiazc13*nr_excr_D 
     &       + rtdiazc13*morp_D*(1-(redntp/diazntp)) + rtdetrc13*remi 
     &       + rtdoc13*recy - fcnpp*npp - fcnpp*npp_D)

        biodoc13 = biodoc13 + dtbio*redctn*(dfr*rtphytc13*morp 
     &       + rtphytc13*dfrt*morpt - rtdoc13*recy)

        biophytc13 = biophytc13 + dtbio*redctn*(fcnpp*npp 
     &       - rtphytc13*morp - rtphytc13*graz - rtphytc13*morpt)

        biozoopc13 = biozoopc13 + dtbio*redctn*(rtphytc13*dig_P 
     &       + rtzoopc13*dig_Z + rtdetrc13*dig_Det + rtdiazc13*dig_D 
     &       - rtzoopc13*morz - rtzoopc13*graz_Z - rtzoopc13*excr)

        biodetrc13 = biodetrc13 + dtbio*redctn*(rtphytc13*(1.-dfr)*morp 
     &       + rtphytc13*sf_P + rtzoopc13*sf_Z + rtdetrc13*sf_Det 
     &       + rtdiazc13*sf_D + rtzoopc13*morz - rtdetrc13*remi 
     &       - rtdetrc13*graz_Det - rtdetrc13*expo + rc13impo*impo 
     &       + rtdiazc13*morp_D*(redntp/diazntp))

        biodiazc13 = biodiazc13 + dtbio*redctn*(fcnpp*npp_D 
     &       - rtdiazc13*morp_D - rtdiazc13*graz_D - rtdiazc13*morpt_D)
# endif
        expoout = expoout + expo
# if defined O_npzd_nitrogen_15
        rn15expoout = rn15expoout + rtdetrn15
# endif
# if defined O_carbon_13
        rc13expoout = rc13expoout + rtdetrc13
# endif
        grazout = grazout + graz
        morpout = morpout + morp
        morzout = morzout + morz
        graz_Det_out = graz_Det_out + graz_Det
        graz_Zout = graz_Zout + graz_Z
# if defined O_save_npzd
        nppout = nppout + npp
        morptout = morptout + morpt
        remiout = remiout + remi
        excrout = excrout + excr
#  if defined O_npzd_nitrogen
        npp_Dout = npp_Dout + npp_D
        graz_Dout = graz_Dout + graz_D
        morpt_Dout = morpt_Dout + morpt_D
        morp_Dout = morp_Dout + morp_D
        nfixout = nfixout + npp_D - no3upt_D
#  endif
#  if defined O_npzd_extra_diagnostics
        avej_out = avej_out + avej
        avej_D_out = avej_D_out + avej_D
        gmax_out = gmax_out + gmax
        no3P_out = no3P_out + no3P
        po4P_out = po4P_out + po4P
        po4_D_out = po4_D_out + po4_D
#  endif
# endif
!     flags prevent negative values by setting outgoing fluxes to
!     zero if tracers are lower than trcmin
        if (po4flag .eq. 1) po4flag = 0.5 + sign(0.5,biopo4 - trcmin)
        if (dopflag .eq. 1) dopflag = 0.5 + sign(0.5,biodop - trcmin)
        if (phytflag .eq. 1) phytflag = 0.5 + sign(0.5,biophyt - trcmin)
        if (zoopflag .eq. 1) zoopflag = 0.5 + sign(0.5,biozoop - trcmin)
        if (detrflag .eq. 1) detrflag = 0.5 + sign(0.5,biodetr - trcmin)
# if defined O_npzd_nitrogen
        if (no3flag .eq. 1) no3flag = 0.5 + sign(0.5,biono3 - trcmin)
        if (donflag .eq. 1) donflag = 0.5 + sign(0.5,biodon - trcmin)
        if (diazflag .eq. 1) diazflag = 0.5 + sign(0.5,biodiaz - trcmin)
#  if defined O_npzd_nitrogen_15
        if (din15flag .eq. 1) 
     &     din15flag = 0.5 + sign(0.5, biodin15 - trcmin)
        if (don15flag .eq. 1)
     &     don15flag = 0.5 + sign(0.5, biodon15 - trcmin)
        if (phytn15flag .eq. 1) 
     &     phytn15flag = 0.5 + sign(0.5, biophytn15 - trcmin)
        if (zoopn15flag .eq. 1)
     &     zoopn15flag = 0.5 + sign(0.5, biozoopn15 - trcmin)
        if (detrn15flag .eq. 1)
     &     detrn15flag = 0.5 + sign(0.5, biodetrn15 - trcmin)
        if (diazn15flag .eq. 1)
     &     diazn15flag = 0.5 + sign(0.5, biodiazn15 - trcmin)
#  endif
# endif
# if defined O_carbon_13
        if (dic13flag .eq. 1)
     &       dic13flag = 0.5 + sign(0.5,biodic13 - trcmin)
        if (doc13flag .eq. 1)
     &       doc13flag = 0.5 + sign(0.5,biodoc13 - trcmin)
        if (phytc13flag .eq. 1)
     &       phytc13flag = 0.5 + sign(0.5,biophytc13 - trcmin)
        if (zoopc13flag .eq. 1)
     &       zoopc13flag = 0.5 + sign(0.5,biozoopc13 - trcmin)
        if (detrc13flag .eq. 1)
     &       detrc13flag = 0.5 + sign(0.5,biodetrc13 - trcmin)
#  if defined O_npzd_nitrogen
        if (diazc13flag .eq. 1)
     &       diazc13flag = 0.5 + sign(0.5,biodiazc13 - trcmin)
#  endif
# endif
      enddo

      bioout(1) = biopo4 - bioin(1)
      bioout(2) = biodop - bioin(2)
      bioout(3) = biophyt - bioin(3)
      bioout(4) = biozoop - bioin(4)
      bioout(5) = biodetr - bioin(5)
# if defined O_npzd_nitrogen
      bioout(6) = biodic - bioin(6)
      bioout(7) = biono3 - bioin(7)
      bioout(8) = biodon - bioin(8)
      bioout(9) = biodiaz - bioin(9)
#  if defined O_npzd_nitrogen_15
      bioout(10) = biodin15 - bioin(10)
      bioout(11) = biodon15 - bioin(11)
      bioout(12) = biophytn15 - bioin(12)
      bioout(13) = biozoopn15 - bioin(13)
      bioout(14) = biodetrn15 - bioin(14)
      bioout(15) = biodiazn15 - bioin(15)
#  endif
# endif
# if defined O_carbon_13
      bioout(16) = biodic13 - bioin(16)
      bioout(17) = biodoc13 - bioin(17)
      bioout(18) = biophytc13 - bioin(18)
      bioout(19) = biozoopc13 - bioin(19)
      bioout(20) = biodetrc13 - bioin(20)
      bioout(21) = biodiazc13 - bioin(21)
# endif
#endif
      return
      end

#if defined O_npzd_fe_limitation
      subroutine fe_limit(i,j,k,felimit,felimit_D)
      implicit none
      integer i,j,k,m,fe_n,fe_jlo,fe_m,fe_k
      parameter (fe_n = 14)
      real fe_x(fe_n), fe_y(fe_n), fe_dy, fe_conc
      real felimit, felimit_D
#include "size.h"
#include "npzd.h"
#include "tmngr.h"
      if (k.le.kmfe) then
#ifndef O_TMM      
! create x (time) and y (data) arrays for fe interpolation
! first zero them out just in case
         fe_y(:) = 0
         fe_x(:) = 0
! write in values for days 0 and 365 for interp. bounds
         fe_y(1) = fe_dissolved(i,j,k,1)
         fe_x(1) = 0
         fe_y(14) = fe_dissolved(i,j,k,12)
         fe_x(14) = 365
! now create the rest of the array based on the BLING
! data which is from day 16 of each month
         do m=2,13
            fe_y(m) = fe_dissolved(i,j,k,m-1)
            if (m.eq.2) then
               fe_x(2) = 16
            else if (m.eq.3) then
               fe_x(3) = 47
            else if (m.eq.4) then
               fe_x(4) = 75
            else if (m.eq.5) then
               fe_x(5) = 106
            else if (m.eq.6) then
               fe_x(6) = 136
            else if (m.eq.7)  then
               fe_x(7) = 167
            else if (m.eq.8) then
               fe_x(8) = 197
            else if (m.eq.9) then
               fe_x(9) = 228
            else if (m.eq.10) then
               fe_x(10) = 259
            else if (m.eq.11) then
               fe_x(11) = 289
            else if (m.eq.12) then
               fe_x(12) = 320
            else if (m.eq.13) then
               fe_x(13) = 350
            endif
         enddo
         fe_jlo = 2
         fe_m = 4
! find fe data day index
         call hunt (fe_x,fe_n,dayoyr,fe_jlo)
! initialize the fe data array at the right day 
         fe_k = min(max(fe_jlo-(fe_m-1)/2,1),fe_n+1-fe_m)
! interpolate the fe data, note this does not use the whole array
         call polint (fe_x(fe_k),fe_y(fe_k)
     &,                         fe_m,dayoyr,fe_conc,fe_dy)
! calculate the fe limitation term                
#else
         fe_conc = fe_dissolved(i,j,k,1) ! this is the time-interpolated value
#endif
         felimit = fe_conc/(kfe + fe_conc)
         felimit_D = fe_conc/(kfe_D + fe_conc)
      else
         felimit = 1
         felimit_D = 1
      end if
!     print*,'felimit=',felimit
      return
      end

      subroutine polint (xa,ya,n,x,y,dy)
!
      implicit none
!
      integer n,i,m,ns, nmax
      parameter (nmax = 10) ! largest anticipated value of n
      real dy,x,y,xa(n),ya(n)
      real den,dif,dift,ho,hp,w,c(nmax),d(nmax)
!
! Given arrays xa and ya, each of length n, and a given value x, this routine
! returns a value y, and an error estimate dy. If P(x) is the polynomial of
! degree N-1 such that P(xai) = yai, i=1,....,n then the returned value y=P(x)
!
!
      ns=1
      dif=abs(x-xa(1))
!
      do i=1,n
         dift = abs(x-xa(i))
         IF (dift.lt.dif) THEN
            ns = i
            dif = dift
         END IF
!
       c(i) = ya(i)
       d(i) = ya(i)
!
      enddo
!
      y = ya (ns)        ! this is the initial approximation to y
      ns = ns-1
!
      do m=1,n-1
         do i=1,n-m
            ho = xa(i) - x
            hp = xa (i+m) -x
            w = c(i+1) - d(i)
            den = ho - hp
            IF (den.eq.0.) PAUSE 'failure in polint'
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
         enddo
         IF (2*ns.lt.n-m) THEN
            dy = c(ns+1)
         ELSE
            dy = d(ns)
            ns = ns - 1
         END IF
         y = y + dy
      enddo
      return
      end
!
      subroutine hunt (xa,n,x,jlo)
!
      integer jlo,n,inc,jhi,jm
      real x, xa(n)
      logical ascnd
!
! Given an array xa(1:n), and given a value x, returns a value jlo
! such that x is between xa(jlo) and xa(jlo +1). xa(1:n) must be
! monotonic, either increasing or decreasing. jlo=0 or jlo=n is
! returned to indicate that x is out of range. jlo on input is taken
! as the initial guess for jlo on output
!
      ascnd = xa(n).ge.xa(1)    ! True if ascending order of table, false otherwise
      IF (jlo.le.0.or.jlo.gt.n) THEN ! Input guess not useful. Go to bisection
         jlo = 0
         jhi = n+1
         goto 3
      END IF
      inc = 1                   ! Set the hunting increment
      IF (x.ge.xa(jlo).eqv.ascnd) THEN ! Hunting up:
 1       jhi = jlo + inc
         IF (jhi.gt.n) THEN     ! Done hunting, since off end of table
            jhi = n + 1
         ELSE IF (x.ge.xa(jhi).eqv.ascnd) THEN ! Not done hunting
            jlo = jhi
            inc = inc + inc     ! so double the increment
      goto 1              ! and try again
      END IF                    ! Done hunting, value bracketed.
      ELSE                      ! Hunt down:
         jhi = jlo
 2       jlo = jhi-inc
         IF (jlo.lt.1) THEN     !Done hunting, since off end of table
            jlo = 0
         ELSE IF (x.lt.xa(jlo).eqv.ascnd) THEN ! Not done hunting
            jhi = jlo
            inc =  inc + inc    ! so double the increment
            goto 2              ! and try again
         END IF                 ! Done hunting, value bracketed
      END IF                    ! Hunt is done, so begin the final bisection phase:
!
 3    IF (jhi-jlo.eq.1) THEN
         IF (x.eq.xa(n)) jlo = n-1
         IF (x.eq.xa(1)) jlo = 1
         RETURN
      END IF
      jm = (jhi + jlo)/2
      IF (x.ge.xa(jm).eqv.ascnd) THEN
         jlo = jm
      ELSE
         jhi = jm
      END IF
!
      GOTO 3
!
      end
#endif
