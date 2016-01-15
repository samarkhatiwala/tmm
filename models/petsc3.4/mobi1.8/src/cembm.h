!======================= include file "cembm.h" ========================

!     parameters for use in the energy balance model (also see atm.h)

!     namix        = time steps between mixing (set in atmos.in)
!     lf           = time step flag (1=>leapfrog, 2=>forward)
!     niats        = number of ice advection sub time steps
!     nivts        = time steps between recalculating ice velocities
!     nivc         = time step counter for nivts
!     ns           = number of subcycles for advection and diffusion
!     pyear        = default paleo calendar year (-/+ = BC/AD)
!     dts          = time step (2*dtatm=>leapfrog, dtatm=>forward)
!     co2ccn       = atmospheric CO2 concentration (ppmv)
!     co2emit      = atmospheric CO2 emissions flux (g cm-2 s-1)
!     co2emit_fuel = emissions flux from fossil fuels (g cm-2 s-1)
!     co2emit_land = emissions flux from land use change (g cm-2 s-1)
!     anthro       = radiative forcing by atmospheric CO2
!     co2for       = atmospheric CO2 forcing term (in units of heat flux)
!     c14ccn       = atmospheric C14 concentration (ppmv)
!     dc14ccn      = atmospheric dC14 concentration (permil)
!     dc14ccnn     = northern hemisphere atmospheric dC14 concentration
!     dc14ccne     = equatorial atmospheric dC14 concentration
!     dc14ccns     = southern hemisphere atmospheric dC14 concentration
!     c14flx       = atmospheric C14 flux (umol cm-2 s-1)
!     c14prod      = atmospheric C14 production (umol cm-2 s-1)
!     c13ccn       = atmospheric C13 concentration (ppmv)
!     dc13ccn      = atmospheric dC13 concentration (permil)
!     cfc11ccnn    = northern hemisphere atmospheric CFC11 concentration
!     cfc11ccns    = southern hemisphere atmospheric CFC11 concentration
!     cfc12ccnn    = northern hemisphere atmospheric CFC12 concentration
!     cfc12ccns    = southern hemisphere atmospheric CFC12 concentration
!     scatter      = proportion of solar scattered by the atmosphere
!     solarconst   = solar constant (g s-3)
!     cssh         = constant used in calculation of ssh (g g-1)
!     cdatm        = drag coefficient (dimensionless)
!     cpatm        = atmospheric heat capacity (cm**2 s-2 K-1)
!     sht          = scale height for temperature
!     shq          = scale height for specific humidity
!     shc          = scale height for carbon
!     rhoatm       = density of air at sea surface (g cm-3)
!     esatm        = atmosphere emissivity times Stefan's constant
!     rhoocn       = representative sea surface density
!     esocn        = ocean emissivity times Stefan's constant
!     vlocn        = latent heat of vapourization of water
!     cdice        = drag coefficient (dimensionless)
!     dampice      = time scale for freezing first layer under ice (days)
!     rhoice       = ice density (g cm-3)
!     rhosno       = snow density (g cm-3)
!     esice        = ice emissivity times Stefan's constant
!     slice        = latent heat of sublimation of ice
!     flice        = latent heat of fusion of ice (cm2 s-2)
!     condice      = ice conductivity (g*cm s-3 K-1)
!     tsno         = air temperature for accumulating snow
!     hsno_max     = maximum snow depth
!     totaltime    = total time for long term averages
!     rlapse       = lapse rate (K cm-1)
!     soilmax      = soil water field capacity (cm)
!     eslnd        = land emissivity time Stefan's constant
!     pass         = atmospheric transmission coefficient (%/100)
!     ice_calb     = ice coalbedo (%/100)
!     sno_calb     = snow coalbedo (%/100)
!     pcfactor     = precip - cloud correlation factor ( %/100)
!     rf1          = factor used in calculating lapse rate reduction
!     rf2          = factor used in calculating lapse rate reduction
!     co2_yr       = co2 forcing year (-/+ = BC/AD)
!     agric_yr     = agricutural forcing year (-/+ = BC/AD)
!     landice_yr   = ice sheet forcing year (-/+ = BC/AD)
!     solar_yr     = solar forcing year (-/+ = BC/AD)
!     orbit_yr     = orbital forcing year (-/+ = BC/AD)
!     volcano_yr   = volcanic forcing year (-/+ = BC/AD)
!     sulph_yr     = sulphate forcing year (-/+ = BC/AD)
!     aggfor_yr    = additional greenhouse gas forcing year (-/+ = BC/AD)
!     cfcs_yr      = cfc forcing year (-/+ = BC/AD)
!     c14_yr       = c14 forcing year (-/+ = BC/AD)
!     ice_yr       = redundant (same as landice_yr, kept for compatibility)
!     dalt_v       = dalton number over vegetation (no vegetation model)
!     dalt_o       = dalton number over ocean
!     dalt_i       = dalton number over ice
!     rhmax        = maximum relative humidity
!     volcfor      = anomalous volcanic forcing (g/s**3)
!     aggfor       = anomalous additional greenhouse gas forcing (g s-3)
!     aggfor_os    = additional greenhouse gas forcing offset (g s-3)
!     atmsa        = atmospheric surface area (cm 2)
!     ocnsa        = exposed ocean surface area (cm 2)
!     sealev       = anomalous sea level change (cm)
!     dsealev      = anomalous sea level change in segtim (cm)
!     sealev_yr    = sea level forcing year (-/+ = BC/AD)
!     itrack_co2   = index for averaging co2 array
!     ntrack_co2   = number of points for average co2
!     itrack_sat   = index for averaging sat array
!     ntrack_sat   = number of points for average sat
!     vcsref       = climate sensitivity reference temperature (C)
!     vcsfac       = climate sensitivity factor (mW/m2/C)
!     vcsyri       = year to start changing the climate sensitivity (year)
!     gtoppm       = conversion from g cm-2 to ppmv
!     carbemit     = accumulated co2 emissions (Pg)
!     adiff        = anomolous diffusion factor (percent/100/C)
!     dtbar        = global average sat anomoly (C)

      integer namix, lf, niats, nivts, nivc, ns, itrack_co2
      integer ntrack_co2, itrack_sat, ntrack_sat

      common /cembm_i/ namix, lf, niats, nivts, nivc, ns, itrack_co2
      common /cembm_i/ ntrack_co2, itrack_sat, ntrack_sat

      real pyear, dts, co2ccn, co2emit, co2emit_fuel, co2emit_land
      real anthro, co2for, c14ccn, dc14ccn, dc14ccnn, dc14ccne, dc14ccns
      real c14prod, c14flx, c13ccn, dc13ccn
      real cfc11ccnn, cfc11ccns, cfc12ccnn, cfc12ccns, scatter
      real solarconst, cssh, cdatm, cpatm, sht, shq, shc, rhoatm, esatm
      real rhoocn, esocn, vlocn, cdice, dampice, rhoice, rhosno, esice
      real slice, flice, condice, tsno, hsno_max, totaltime, rlapse
      real soilmax, eslnd, pass, ice_calb, sno_calb, pcfactor, rf1, rf2
      real co2_yr, agric_yr, landice_yr, solar_yr, orbit_yr, volcano_yr
      real sulph_yr, aggfor_yr, cfcs_yr, c14_yr, ice_yr, dalt_v, dalt_o
      real dalt_i, rhmax, volcfor, aggfor, aggfor_os, atmsa, ocnsa
      real sealev, dsealev, sealev_yr, vcsref, vcsfac, vcsyri, gtoppm
      real carbemit, adiff, dtbar

      common /cembm_r/ pyear, dts, co2ccn, co2emit, co2emit_fuel
      common /cembm_r/ co2emit_land, anthro, co2for, c14ccn, dc14ccn
      common /cembm_r/ dc14ccnn, dc14ccne, dc14ccns, c14prod, c14flx
      common /cembm_r/ c13ccn, dc13ccn, cfc11ccnn
      common /cembm_r/ cfc11ccns, cfc12ccnn, cfc12ccns, cssh, cdatm
      common /cembm_r/ cpatm, sht, shq, shc, rhoatm, esatm, rhoocn
      common /cembm_r/ scatter, solarconst, esocn, dampice, rhoice
      common /cembm_r/ rhosno, esice, slice, flice, condice, vlocn
      common /cembm_r/ cdice, tsno, hsno_max, totaltime, rlapse, soilmax
      common /cembm_r/ eslnd, pass, ice_calb, sno_calb, pcfactor, rf1
      common /cembm_r/ rf2, co2_yr, agric_yr, landice_yr, solar_yr
      common /cembm_r/ orbit_yr, volcano_yr, sulph_yr, aggfor_yr
      common /cembm_r/ cfcs_yr, c14_yr, ice_yr, dalt_v, dalt_o, dalt_i
      common /cembm_r/ rhmax, volcfor, aggfor, aggfor_os, atmsa, ocnsa
      common /cembm_r/ sealev, dsealev, sealev_yr, vcsref, vcsfac
      common /cembm_r/ vcsyri, gtoppm, carbemit, adiff, dtbar

!     ntatsa        = time step counter for time averaging
!     ntatia        = number of time averaged time step integrals
!     tai_sat       = average integrated elevated surface air temperature
!     tai_shum      = average integrated surface specific humidity
!     tai_precip    = average integrated precipitation
!     tai_evap      = average integrated evaporation
!     tai_ohice     = total integrated sea ice volume
!     tai_oaice     = total integrated sea ice area
!     tai_hsno      = total integrated snow volume
!     tai_lhice     = total integrated land ice volume
!     tai_laice     = total integrated land ice area
!     tai_co2ccn    = average integrated CO2 concentration
!     tai_co2emit   = average integrated CO2 emissions
!     tai_dc14ccn   = average integrated dC14 concentration
!     tai_ac14flx    = average dC14 flux (land & ocean)
!     tai_dc13ccn   = average integrated dC13 concentration
!     tai_cfc11ccn  = average integrated CFC11 concentration
!     tai_cfc12ccn  = average integrated CFC12 concentration
!     tai_maxit     = average maximum iterations for atmospheric solver
!     tai_nsat      = average northern hemisphere surface air temperature
!     tai_ssat      = average southern hemisphere surface air temperature
!     tai_nshum     = average northern hemisphere surface specific humidity
!     tai_sshum     = average southern hemisphere surface specific humidity
!     tai_nprecip   = average northern hemisphere precipitation
!     tai_sprecip   = average southern hemisphere precipitation
!     tai_nevap     = average northern hemisphere evaporation
!     tai_sevap     = average southern hemisphere evaporation
!     tai_nohice    = total northern hemisphere sea ice volume
!     tai_sohice    = total southern hemisphere sea ice volume
!     tai_noaice    = total northern hemisphere sea ice area
!     tai_soaice    = total southern hemisphere sea ice area
!     tai_nhsno     = total northern hemisphere snow volume
!     tai_shsno     = total southern hemisphere snow volume
!     tai_nlhice    = total northern hemisphere land ice volume
!     tai_slhice    = total southern hemisphere land ice volume
!     tai_nlaice    = total northern hemisphere land ice area
!     tai_slaice    = total southern hemisphere land ice area
!     tai_lsat      = average sat over land
!     tai_osat      = average sat over ocean
!     tai_lprecip   = average precip over land
!     tai_oprecip   = average precip over ocean
!     tai_levap     = average evap over land
!     tai_oevap     = average evap over ocean
!     tai_solins    = average incoming solar
!     tai_upsens    = average surface sensible heat
!     tai_uplwr     = average surface upward longwave radiation
!     tai_outlwr    = average outgoing longwave radiation
!     tai_dnswr     = average downward (absorbed) shortwave
!     tai_absswr    = net absorbed shortwave radiation
!     tai_netrad    = net radiation at the top of the atmosphere
!     tai_palb      = average planetary albedo
!     tai_aalb      = average atmospheric albedo
!     tai_salb      = average surface albedo
!     tai_lsalb     = average surface albedo over land
!     tai_osalb     = average surface albedo over ocean
!     tai_sst       = average sea surface temperature
!     tai_sss       = average sea surface salinity
!     tai_ssdic     = average sea surface dic
!     tai_ssc14     = average sea surface c14
!     tai_ssalk     = average sea surface alkalinity
!     tai_sso2      = average sea surface o2
!     tai_sspo4     = average sea surface po4
!     tai_ssdop     = average sea surface dop
!     tai_ssno3     = average sea surface no3
!     tai_ssdon     = average sea surface don
!     tai_sscfc11   = average sea surface cfc11
!     tai_sscfc12   = average sea surface cfc12
!     tai_sulph     = average tropospheric sulphate forcing
!     tai_volc      = average volcanic forcing
!     tai_agg       = average additional greenhouse gas forcing
!     tai_catm      = average total carbon in atmosphere
!     tai_carbemit  = average taccumulated co2 emissions

      integer ntatsa, ntatia

      common /cembm_i/ ntatsa, ntatia

      real tai_sat, tai_shum, tai_precip, tai_evap, tai_ohice, tai_oaice
      real tai_hsno, tai_lhice, tai_laice, tai_co2ccn, tai_co2emit
      real tai_dc14ccn, tai_cfc11ccn, tai_cfc12ccn, tai_maxit, tai_nsat
      real tai_dc13ccn, tai_ac14flx
      real tai_ssat, tai_nshum, tai_sshum, tai_nprecip, tai_sprecip
      real tai_nevap, tai_sevap, tai_nohice, tai_sohice, tai_noaice
      real tai_soaice, tai_nhsno, tai_shsno, tai_nlhice, tai_slhice
      real tai_nlaice, tai_slaice, tai_lsat, tai_osat, tai_lprecip
      real tai_oprecip, tai_levap, tai_oevap, tai_solins, tai_upsens
      real tai_uplwr, tai_outlwr, tai_dnswr, tai_absswr, tai_netrad
      real tai_palb, tai_aalb, tai_salb, tai_lsalb, tai_osalb, tai_sst
      real tai_sss, tai_ssdic, tai_ssc14, tai_ssalk, tai_sso2, tai_sspo4
      real tai_ssno3, tai_sscfc11, tai_sscfc12, tai_sulph, tai_volc
      real tai_agg, tai_catm, tai_carbemit, tai_ssdop, tai_ssdon
      real tai_ssdin15, tai_ssdon15, tai_ssdic13, tai_ssdoc13

      common /cembm_r/ tai_sat, tai_shum, tai_precip, tai_evap
      common /cembm_r/ tai_ohice, tai_oaice, tai_hsno, tai_lhice
      common /cembm_r/ tai_laice, tai_co2ccn, tai_co2emit, tai_dc14ccn
      common /cembm_r/ tai_dc13ccn, tai_ac14flx
      common /cembm_r/ tai_cfc11ccn, tai_cfc12ccn, tai_maxit, tai_nsat
      common /cembm_r/ tai_ssat, tai_nshum, tai_sshum, tai_nprecip
      common /cembm_r/ tai_sprecip, tai_nevap, tai_sevap, tai_nohice
      common /cembm_r/ tai_sohice, tai_noaice, tai_soaice, tai_nhsno
      common /cembm_r/ tai_shsno, tai_nlhice, tai_slhice, tai_nlaice
      common /cembm_r/ tai_slaice, tai_lsat, tai_osat, tai_lprecip
      common /cembm_r/ tai_oprecip, tai_levap, tai_oevap, tai_solins
      common /cembm_r/ tai_upsens, tai_uplwr, tai_outlwr, tai_dnswr
      common /cembm_r/ tai_absswr, tai_netrad, tai_palb, tai_aalb
      common /cembm_r/ tai_salb, tai_lsalb, tai_osalb, tai_sst
      common /cembm_r/ tai_sss, tai_ssdic, tai_ssc14, tai_ssalk
      common /cembm_r/ tai_sso2, tai_sspo4, tai_ssno3, tai_sscfc11
      common /cembm_r/ tai_sscfc12, tai_sulph, tai_volc, tai_agg
      common /cembm_r/ tai_catm, tai_carbemit, tai_ssdop, tai_ssdon
      common /cembm_r/ tai_ssdin15, tai_ssdon15, tai_ssdic13
      common /cembm_r/ tai_ssdoc13
