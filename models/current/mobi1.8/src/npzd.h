!====================== include file "npzd.h" =========================

!   variables for npzd model

!   ntnpzd   = number of npzd tracers
!   nbio     = number of npzd timesteps per ocean timestep
!   trcmin   = minimum tracer for flux calculations
!   tap      = 2*alpha*par with
!   alpha    = initial slope P-I curve [(W/m^2)^(-1)/day] and
!   par      = fraction of photosythetically active radiation
!   kw       = light attenuation due to water [1/m]
!   kc       = light attenuation by phytoplankton [1/(m*mmol m-3)]
!   ki       = light attenuation through sea ice & snow
!   abio     = maximum growth rate parameter [1/day]
!   bbio     = b
!   cbio     = [1/deg_C]
!   k1n      = half saturation constant for N uptake [mmol m-3]
!   nup      = specific mortality rate (Phytoplankton) [day-1]
!   gamma1   = assimilation efficiency (zpk)
!   gbio     = maximum grazing rate at 0 deg C [day-1]
!   nuz      = quadratic mortality (zpk)
!   nud0     = remineralization rate [day-1]
!   LFe      = Iron limitation
!   wd       = sinking speed of detritus [m day-1]
!   ztt      = depth to top of grid cell [cm]
!   rkwz     = reciprical of light attenuation times grid depth
!   dtnpzd   = time step of biology
!   capr     = carbonate to carbon production ratio
!   dcaco3   = remineralisation depth of calcite [cm]
!   rcak     = array used in calculating calcite remineralization
!   rcab     = array used in calculating bottom calcite remineralization
!   nupt0    = specific mortality rate (Phytoplankton) [1/day]
!   wd0      = sinking speed of detritus at surface [m/day]
!   k1p      = half saturation constant for P uptake
!   jdiar    = factor reducing the growth rate of diazotrophs
!   redctn   = C/N Redfield ratio (includes mol to mmol conversion)
!   redctp   = C/P Redfield ratio (includes mol to mmol conversion)
!   redptn   = P/N Redfield ratio
!   redntp   = N/P Redfield ratio
!   redotn   = O/N Redfield ratio (includes mol to mmol conversion)
!   redotp   = O/P Redfield ratio (includes mol to mmol conversion)
!   rnbio    = reciprical of nbio
!   rdtts    = reciprical of dtts [s-1]
!   dtbio    = npzd time step [s]
!   rnpp     = rate of net primary production [nmol cm-3 s-1]
!   rgraz    = rate of grazing [nmol cm-3 s-1]
!   rmorp    = rate of mortality of phytoplankton [nmol cm-3 s-1]
!   rmorz    = rate of mortality of zooplankton [nmol cm-3 s-1]
!   rremi    = rate of remineralization [nmol cm-3 s-1]
!   rexcr    = rate of excretion [nmol cm-3 s-1]
!   rexpo    = rate of export through the bottom [nmol cm-3 s-1]
!   rnpp_D   = npp for diazotraphs [nmol cm-3 s-1]
!   rgraz_D  = rgraz for diazotraphs [nmol cm-3 s-1]
!   rmorpt_D = rmorp for diazotraphs [nmol cm-3 s-1]
!   rnfix    = rate of nitrogen fixation [nmol cm-3 s-1]
!   rdeni    = rate of denitrification [nmol cm-3 s-1]
!   kzoo     = half saturation constant for Z grazing
!   zprefP   = Z preference for grazing on P
!   zprefD   = Z preference for grazing on Diaz
!   zprefZ   = Z preference for grazing on other Z
!   zprefDet = Z preference for grazing on Detritus
!   rgraz_Det = rate of grazing on Detritus [nmol cm-3 s-1]
!   rgraz_Z   = rate of grazing on other Zooplankton [nmol cm-3 s-1]
!   geZ      = growth efficiency of zooplankton
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! New diagnostic output
!   ravej      = light-dependant growth rate of P
!   ravej_D    = light-dependant growth rate of Diaz
!   rgmax      = temp-dependant growth rate of zoo
!   rno3P      = nitrate-dependant growth rate of P
!   rpo4P       = phosphate-dependant growth rate of P
!   rpo4_D     = phosphate-dependant growth rate of D
!
!   fe_dissolved = dissolved iron concentration
!   kfe = Fe limitation half saturation parameter
!   kfe_D = Fe limitation half sat. param. for diaz.
!
!++++ Climate engineering ++++++++++++++
!   kpipe = ocean pipe coordinate
!   kpipe_fe = ocean pipe coordinate for fe limitation
!

      integer ntnpzd, nbio, kmfe

      parameter (ntnpzd = 5
#if defined O_npzd_nitrogen
     &                  + 4
# if defined O_npzd_nitrogen_15
     &                  + 6
# endif
#endif
#if defined O_carbon_13
     &                  + 5
# if defined O_npzd_nitrogen
     &                  + 1
# endif
#endif
     &                     )
      common /npzd_i/ nbio

      real trcmin
      parameter (trcmin=5e-12)

#if defined O_npzd_nitrogen_15
      real rn15std

      parameter (rn15std=0.0036765)
#endif
#if defined O_carbon_13
      real rc13std
      parameter (rc13std=0.0112372)
#endif
#if defined O_carbon_14
!     rc14std   = standard c14/c12 ratio
      real rc14std
      parameter (rc14std=1.176e-12)
#endif

# if defined O_npzd
      real tap, kw, kc, ki, abio, bbio, cbio, k1n, nup, gamma1, gbio
      real epsbio, nuz, nud0, LFe, wd, ztt, rkwz, dtnpzd
      real capr, dcaco3, rcak, rcab, nupt0, wd0, k1p, jdiar, redctn
      real redctp, redptn, redntp, redotn, redotp, rnbio, rdtts, dtbio
      real rnpp, rgraz, rmorp, rmorpt, rmorz, rremi, rexcr, rexpo
      real rnpp_D, rgraz_D, rmorpt_D, rnfix, kzoo, zprefP, rmorp_D
      real zprefD, zprefZ, zprefDet, rgraz_Det, rgraz_Z, geZ, kfe
      real ravej, ravej_D, rgmax, rno3P, rpo4P, rpo4_D, kfe_D, kpipe
      real kpipe_fe, rwcdeni, rbdeni, rsedrr, sgbdfac, nupt0_D
      real diazntp, diazptn, nup_D, dfr, redotc, redntc, dfrt
      real redptc, rprca

      common /npzd_r/ tap, kw, kc, ki, abio, bbio, cbio, k1n, nup
      common /npzd_r/ gamma1, gbio, epsbio, nuz, nud0, LFe, dfr
      common /npzd_r/ wd(km), ztt(km), rkwz(km), dtnpzd, capr
      common /npzd_r/ dcaco3, rcak(km), rcab(km), nupt0, wd0, k1p
      common /npzd_r/ jdiar, redctn, redctp, redptn, redntp, redotn
      common /npzd_r/ redotp, rnbio, rdtts, dtbio, geZ
      common /npzd_r/ kzoo, zprefP, zprefD, zprefZ, zprefDet, kfe, kfe_D
      common /npzd_r/ sgbdfac, nupt0_D, diazntp, diazptn, nup_D
      common /npzd_r/ redotc, redntc, dfrt, redptc
# if defined O_npzd_fe_limitation
      real fe_dissolved
      common /fe_dissolved_r/ fe_dissolved(imt,jmt,km,12)
      common /fe_dissolved_i/ kmfe      
# endif
# if defined O_save_npzd
      common /npzd_r/ rnpp(kpzd), rgraz(kpzd), rmorp(kpzd), rmorpt(kpzd)
      common /npzd_r/ rmorz(kpzd), rexcr(kpzd), rremi(km), rexpo(km)
      common /npzd_r/ rgraz_Det(kpzd), rgraz_Z(kpzd), rsedrr, rprca
#  if defined O_npzd_nitrogen
      common /npzd_r/ rnpp_D(kpzd), rgraz_D(kpzd), rmorp_D(kpzd)	 
      common /npzd_r/ rmorpt_D(kpzd), rnfix(kpzd), rwcdeni(km)
      common /npzd_r/ rbdeni(km)
#   if defined O_npzd_extra_diagnostics
      common /npzd_r/ ravej(kpzd), ravej_D(kpzd), rgmax(kpzd), rno3P(kpzd)
      common /npzd_r/ rpo4P(kpzd), rpo4_D(kpzd)
#   endif
#  endif
# endif
#endif
