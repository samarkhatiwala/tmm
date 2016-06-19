!======================== include file "ice.h" =========================

!     arrays for the ice model
!     hice     = ice thickness (tau-1, tau, tau+1)
!     aice     = ice area (tau-1, tau, tau+1)
!     tice     = ice surface temperature
!     hsno     = effective thickness of snow
!     psno     = precipitation as snow (+ accumulation, - melt)
!     frzpt    = freezing temperature of sea water
!   cpts ice model ("ice_cpts")
!     A        = area
!     heff     = effective thickness of ice
!     Ts       = ice surface temperature
!     groice   = thermo ice thickness change from diagnostic
!     E        = enthalpy of each ice layer
!     strength = pressure (from Rothrock), must have ncat > 3
!     hseff    = effective thickness of snow
!   ice dynamics ("ice_evp")
!     uice     = u ice velocity
!     vice     = v ice velocity
!     pice     = pressure due to internal stress
!     xint     = x component of ice interaction
!     yint     = y component of ice interaction
!     del      = delta for ice dynamics models
!     eI       = divergence of ice velocity
!     nseg     = number of domains needed for ice calculations
!     jsi      = start of limited domain for ice calculations
!     jei      = end of limited domain for ice calculations
!   brine convection ("convect_brine")
!     cbf      = flux from rejected brine over ice categories
!     cba      = area of brine rejection for each categories
!     cba0     = area of the cell without brine rejection
!   land ice ("landice_data")
!     hicel    = land ice thickness
!     aicel    = land ice area

      integer nseg, jsi, jei

      real hice, aice, tice, hsno, psno, frzpt
      real A, heff, Ts, groice, E, hseff, strength
      real uice, vice, pice, xint, yint, del, eI
      real cbf, cba, cba0, hicel, aicel

      common /ice_r/ hice(imt,jmt,2), aice(imt,jmt,2), tice(imt,jmt)
      common /ice_r/ hsno(imt,jmt,2), psno(imt,jmt), frzpt(imt,jmt)
#if defined O_ice_cpts && defined O_ice
      common /ice_r/ A(imt,jmt,itme,ncat), heff(imt,jmt,itme,ncat)
      common /ice_r/ Ts(imt,jmt,itme,ncat), groice(imt,jmt,ncat)
      common /ice_r/ E(imt,jmt,itme,ntilay), hseff(imt,jmt,itme,ncat)
#else
      integer ncat
      parameter (ncat=1)
#endif
#if defined O_ice_cpts && defined O_ice_cpts_roth_press && defined O_ice
      common /ice_r/ strength(imt,jmt,itme)
#endif
#if defined O_ice_evp && defined O_ice
      common /ice_r/ uice(imt,jmt), vice(imt,jmt), pice(imt,jmt)
      common /ice_r/ xint(imt,jmt), yint(imt,jmt)
      common /ice_r/ del(imt,jmt), eI(imt,jmt)
      common /ice_i/ nseg, jsi(jmt), jei(jmt)
#endif
#if defined O_convect_brine
      common /ice_r/ cbf(imt,jmt,0:ncat), cba(imt,jmt,0:ncat)
      common /ice_r/ cba0(imt,jmt)
# endif
#if defined O_landice_data || defined O_landice_data_transient
      common /ice_r/ hicel(imt,jmt,3), aicel(imt,jmt,3)
#endif
#if defined O_time_averages
!   time averaged arrays
!     ta_hsno     = time average snow thickness
!     ta_hice     = time average ice thickness
!     ta_aice     = time average ice area
!     ta_tice     = time average ice surface temperature
!   ice dynamics ("ice_evp")
!     ta_uice     = time average u ice velocity
!     ta_vice     = time average v ice velocity
!     ta_pice     = time average pressure
!     ta_xint     = time average x component of ice interaction
!     ta_yint     = time average y component of ice interaction
!   land ice ("landice_data")
!     ta_hicel    = time average land ice thickness
!     ta_aicel    = time average land ice area

      real ta_hice, ta_aice, ta_tice, ta_hsno, ta_heff, ta_A, ta_Ts
      real ta_hseff, ta_uice, ta_vice, ta_pice, ta_xint, ta_yint
      real ta_aicel, ta_hicel

# if defined O_ice
      common /ice_r/ ta_hice(imt,jmt), ta_aice(imt,jmt)
      common /ice_r/ ta_tice(imt,jmt), ta_hsno(imt,jmt)
# endif
# if defined O_ice_cpts && defined O_ice
      common /ice_r/ ta_heff(imt,jmt,ncat), ta_A(imt,jmt,ncat)
      common /ice_r/ ta_Ts(imt,jmt,ncat), ta_hseff(imt,jmt,ncat)
# endif
# if defined O_ice_evp && defined O_ice
      common /ice_r/ ta_uice(imt,jmt), ta_vice(imt,jmt)
      common /ice_r/ ta_pice(imt,jmt), ta_xint(imt,jmt)
      common /ice_r/ ta_yint(imt,jmt)

# endif
# if defined O_landice_data
      common /ice_r/ ta_hicel(imt,jmt), ta_aicel(imt,jmt)
# endif
#endif
