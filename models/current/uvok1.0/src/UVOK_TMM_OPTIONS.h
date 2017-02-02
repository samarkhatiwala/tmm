#define O_TMM
#define O_co2ccn_user
#undef O_TMM_interactive_atmosphere
#undef O_TMM_partial_step_topo
#define O_even_fluxes
#undef O_read_my_kmt 
#define O_read_my_grid 
#define O_cyclic
#undef O_time_averages
#undef O_time_step_monitor
#define O_sbc_in_memory
#define O_fourfil
#define O_constant_flux_reference

#define O_embm
#define O_embm_mgrid
#define O_embm_annual
#define O_embm_adiff
#undef O_embm_awind

#define O_ice
#define O_ice_evp  
#define O_ice_fourfil
#undef O_correct_ice_to_ocean_stress

#undef O_mtlm
#undef O_mtlm_segday

#define O_mom
#define O_ramdrive 
#define O_conjugate_gradient 
#define O_sf_5_point
#define O_stream_function
#define O_consthmix 
#define O_constvmix 
#define O_fullconvect
#define O_save_convection
#define O_stability_tests
#undef O_gyre_components 
#define O_term_balances
#define O_energy_analysis 
#define O_meridional_overturning
#define O_tracer_averages
#define O_gent_mcwilliams
#define O_isopycmix 
#define O_fct
#define O_npzd
#define O_npzd_alk
#define O_npzd_nitrogen
#define O_npzd_no_vflux
#define O_npzd_o2
#define O_save_npzd
#define O_tidal_kv
#define O_highmix_Southern_Ocean
#define O_anisotropic_viscosity
#define O_npzd_fe_limitation
#define O_zoop_graz_upper_temp_limit

#define O_carbon
#define O_carbon_14
#define O_save_carbon_carbonate_chem
#define O_save_carbon_totals

#undef O_co2ccn_data
#undef O_agric_data
#undef O_landice_data
#undef O_solar_data

#define O_tai_otsf
#define O_tai_ns
#define O_tai_lo
#define O_tai_slh
#define O_tai_rad
#define O_units_temperature_Celsius
#define O_units_time_years
#define O_save_time_relyear0


#undef O_volcano_data
#undef O_volcano_data_transient
#undef O_co2emit_data_transient
#undef O_co2ccn_data_transient
#undef O_co2emit_track_co2
#undef O_co2emit_track_sat
#undef O_embm_vcs
#undef O_c14ccn_data
#undef O_c14ccn_data_transient
#undef O_aggfor_data
#undef O_aggfor_data_transient
#undef O_sealev_data_transient
#undef O_ice_cpts
#undef O_crop_data_transient
#undef O_pasture_data_transient
#undef O_agric_data_transient
#undef O_carbon_fnpzd


