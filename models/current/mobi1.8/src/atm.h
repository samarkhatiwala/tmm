!======================== include file "atm.h" =========================

!     arrays for the energy-moisture balance model

!     note: units for heat flux are: g s-3 (or mW m-2)
!           units for fresh water flux are in g cm-2 s-1
!           downward is into the top surface (ocean, ice or land)
!           upward is into the bottom of the atmosphere
!           outward is out of the top of the atmosphere
!           inward is into the top of the atmosphere

!     at          = tracers (previous and most recent)
!     elev        = elevations (cm)
!     flux        = flux accumulator (1=heat, 2=water, 3=carbon, ...)
!     rh          = relative humidity (%/100)
!     tmsk        = tracer grid ocean mask (0 = land, 1 = ocean)
!     umsk        = velocity grid ocean mask (0 = land, 1 = ocean)
!     dn          = northward atmospheric tracer diffusivity
!     de          = eastward atmospheric tracer diffusivity
!     fcor        = Coriolis factor
!   heat fluxes
!     solins      = incoming short wave radiation (g s-3)
!     dnswr       = downward surface shortwave flux absorbed (g s-3)
!     uplwr       = upward surface longwave flux (g s-3)
!     upsens      = upward surface sensible heat flux (g s-3)
!     upltnt      = upward surface latent heat flux (g s-3)
!     outlwr      = outgoing atmosphere longwave flux (g s-3)
!   fresh water fluxes
!     precip      = precipitation (g cm-2 s-1)
!     evap        = evaporation (g cm-2 s-1)
!     disch       = discharge from continents (g cm-2 s-1)
!     vflux       = normalized virtual flux
!   land model
!     soilm       = soil moisture, as depth in bucket (g cm -2)
!     runoff      = water runoff from continents (g cm-2 s-1)
!     surf        = land surface temperature (C)
!   annual average temperature ("embm_awind" and "embm_adiff")
!     rtbar       = running average atmospheric temperature (C)
!     atbar       = average temperature accumulator (C)
!     tbar        = climatological average atmospheric temperature (C)
!   wind feedback ("embm_awind")
!     awx         = anomalous x component of wind (cm s-1)
!     awy         = anomalous y component of wind (cm s-1)
!     apress      = anomalous pressure (g cm-2 s-2)
!   sulphates ("sulphate_data_transient")
!     sulph       = sulphate forcing
!   distributed CO2 emissions ("co2emit_data" and "carbon_co2_2d")
!     co2dist     = distributed co2 forcing
!   sea level ("sealev" or "sealev_data")
!     elev_sealev = elevation anomaly due to sea level changes
!   flux adjustment ("save_flxadj")
!     flxadj      = flux adjustment
!   co2 emissions from tracking co2 ("co2emit_track_co2")
!     track_co2   = array of global average co2 values
!   co2 emissions from tracking sat("co2emit_track_sat" or "embm_vcs")
!     track_sat   = array of global average sat values
!   indicies for atmospheric tracer array
!     isat        = index for surface air temperature
!     ishum       = index for surface specific humidity
!     ico2        = index for atmospheric co2
!     mapat       = map for atmospheric tracer names
!   mtrack        = maximum size for tracking arrays

      character(10) :: mapat(nat)
      common /embm_c/ mapat

      integer isat, ishum, ico2, mtrack
      parameter (mtrack=3650)
      common /embm_i/ isat, ishum, ico2

      real at, elev, flux, rh, tmsk, umsk, dn, de, fcor, solins, dnswr
      real uplwr, upsens, upltnt, outlwr, precip, evap, disch, vflux
      real soilm, runoff, surf, rtbar, atbar, tbar, awx, awy, apress
      real sulph, co2dist, elev_sealev, flxadj, track_co2, track_sat

      common /embm_r/ at(imt,jmt,2,nat), elev(imt,jmt)
#if defined O_ice_evp || defined O_embm_awind
      common /embm_r/ flux(imt,jmt,nat+2), rh(imt,jmt)
#else
      common /embm_r/ flux(imt,jmt,nat), rh(imt,jmt)
#endif
      common /embm_r/ tmsk(imt,jmt), umsk(imt,jmt), dn(imt,jmt,nat)
      common /embm_r/ de(imt,jmt,nat), fcor(imt,jmt), solins(imt,jmt)
      common /embm_r/ dnswr(imt,jmt), uplwr(imt,jmt), upsens(imt,jmt)
      common /embm_r/ upltnt(imt,jmt), outlwr(imt,jmt), precip(imt,jmt)
      common /embm_r/ evap(imt,jmt), disch(imt,jmt), vflux(imt,jmt)
      common /embm_r/ soilm(imt,jmt,2), runoff(imt,jmt)
      common /embm_r/ surf(imt,jmt)
#if defined O_embm_awind || defined O_embm_adiff
      common /embm_r/ rtbar(imt,jmt), atbar(imt,jmt), tbar(imt,jmt)
#endif
#if defined O_embm_awind
      common /embm_r/ awx(imt,jmt), awy(imt,jmt),apress(imt,jmt)
#endif
#if defined O_sulphate_data || defined O_sulphate_data_transient
      common /embm_r/ sulph(imt,jmt,3)
#endif
#if defined O_carbon_co2_2d
# if defined O_co2emit_data || O_co2emit_data_transient
      common /embm_r/ co2dist(imt,jmt,3)
# endif
#endif
#if defined O_sealev || defined O_sealev_data
      common /embm_r/ elev_sealev(imt,jmt)
#endif
#if defined O_save_flxadj
      common /embm_r/ flxadj(imt,jmt,2)
#endif
#if defined O_co2emit_track_co2 || defined O_co2emit_track_co2_transient
      common /embm_r/ track_co2(mtrack)
#endif
#if defined O_co2emit_track_sat || defined O_co2emit_track_sat_transient || defined O_embm_vcs
      common /embm_r/ track_sat(mtrack)
#endif

#if defined O_time_averages
!   time averaged arrays
!     ta_at         = time averaged atmospheric tracers
!     ta_solins     = time averaged incoming shortwave flux
!     ta_arswr      = time averaged atmospheric reflected shortwave flux
!     ta_dnswr      = time averaged downward shortwave flux
!     ta_absin      = time averaged absorbed shortwave coming in
!     ta_absout     = time averaged absorbed shortwave going out
!     ta_uplwr      = time averaged upward longwave flux
!     ta_upsens     = time averaged upward sensible heat flux
!     ta_outlwr     = time averaged outgoing longwave flux
!     ta_precip     = time averaged precipitation
!     ta_psno       = time averaged precipitation as snow
!     ta_evap       = time averaged evaporation
!     ta_disch      = time averaged discharge
!     ta_vflux      = time averaged virtual flux
!     ta_ws         = time averaged surface wind speed
!     ta_runoff     = time averaged runoff
!     ta_soilm      = time averaged soil moisture
!     ta_surf       = time averaged land surface temperature
!     ta_wx         = time averaged x component of wind
!     ta_wy         = time averaged y component of wind
!   windstress feedback ("embm_awind")
!     ta_awx        = time averaged x component of wind
!     ta_awy        = time averaged y component of wind
!     ta_rtbar      = time averaged running average temperature
!     ta_apress     = time averaged anomalous pressure
!   flux adjustment ("save_flxadj")
!     ta_flxadj     = time averaged flux adjustment
!   diffusivities ("save_embm_diff")
!     ta_dn         = time averaged dn
!     ta_de         = time averaged de
!   flux of co2 ("carbon_co2_2d")
!     ta_flux_co2   = time averaged flux of co2 (g/cm2/s)
!   flux of co2 emissions ("co2emit_data" and "carbon_co2_2d")
!     ta_co2emit    = time averaged emissions (g/cm2/s)
!   sulphates ("sulphate_data_transient")
!     ta_sulph      = time averaged sulph

      real ta_at, ta_solins, ta_arswr, ta_dnswr, ta_absin, ta_absout
      real ta_uplwr, ta_upsens, ta_outlwr, ta_precip, ta_psno
      real ta_evap, ta_disch, ta_vflux, ta_ws, ta_soilm, ta_runoff
      real ta_surf, ta_wx, ta_wy, ta_awx, ta_awy, ta_rtbar, ta_apress
      real ta_flxadj, ta_dn, ta_de, ta_flux_co2, ta_co2emit, ta_sulph

      common /embm_r/ ta_at(imt,jmt,nat), ta_solins(imt,jmt)
      common /embm_r/ ta_arswr(imt,jmt), ta_dnswr(imt,jmt)
      common /embm_r/ ta_absin(imt,jmt), ta_absout(imt,jmt)
      common /embm_r/ ta_uplwr(imt,jmt), ta_upsens(imt,jmt)
      common /embm_r/ ta_outlwr(imt,jmt), ta_precip(imt,jmt)
      common /embm_r/ ta_psno(imt,jmt), ta_evap(imt,jmt)
      common /embm_r/ ta_disch(imt,jmt), ta_vflux(imt,jmt)
      common /embm_r/ ta_ws(imt,jmt)

      common /embm_r/ ta_runoff(imt,jmt), ta_soilm(imt,jmt)
      common /embm_r/ ta_surf(imt,jmt), ta_wx(imt,jmt,nat)
      common /embm_r/ ta_wy(imt,jmt,nat)
# if defined O_embm_awind
      common /embm_r/ ta_awx(imt,jmt), ta_awy(imt,jmt)
      common /embm_r/ ta_rtbar(imt,jmt), ta_apress(imt,jmt)
# endif
# if defined O_save_flxadj
      common /embm_r/ ta_flxadj(imt,jmt,2)
# endif
# if defined O_save_embm_diff
      common /embm_r/ ta_dn(imt,jmt,nat), ta_de(imt,jmt,nat)
# endif
# if defined O_carbon_co2_2d
      common /embm_r/ ta_flux_co2(imt,jmt)
# endif
# if defined O_carbon_co2_2d
#  if defined O_co2emit_data || defined O_co2emit_data_transient
      common /embm_r/ ta_co2emit(imt,jmt)
#  endif
# endif
# if defined O_sulphate_data || defined O_sulphate_data_transient
      common /embm_r/ ta_sulph(imt,jmt)
# endif
#endif
