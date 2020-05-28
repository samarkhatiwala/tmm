!======================= include file "size.h" =========================

!-----------------------------------------------------------------------
!     USER INPUT:
!-----------------------------------------------------------------------

!     imt    = number of grid points in the longitudinal direction
!              (calculated points are from 2 through imt-1. End points
!               are boundaries)
!     jmt    = number of grid points (latitude rows) in the latitudinal
!              direction (calculated points are from 2 through jmt-1.
!              End points are boundaries)
!     km     = number of grid points in the vertical direction
!              (calculated points are from 1 through km)
!     nt     = number of tracers (temperature, salinity, ...)
!     nsrc   = number of tracer with sources
!     kpzd   = depth for limited npzd model
!     jmz    = size for "unrotated" zonal averages
!     jmzm1  = jmz minus one
!     mnisle = maximum number of islands (unconnected land masses)
!     maxipp = maximum number of all island perimeter points
!-----------------------------------------------------------------------

      integer imt, jmt, km, nt, nsrc, kpzd, nat, jmz, jmzm1, mnisle
      integer maxipp, jmw, jsmw, jemw
#ifndef O_TMM
      parameter (imt=  102, jmt=  102, km= 19)
#else
      parameter (imt=  7000, jmt=  1, km= 19)
#endif
      parameter (nt=2 ! temp, salt
#if defined O_matrix
     $             +km
#endif      
#if defined O_carbon
     $             +1 ! dic
# if defined O_carbon_13
     $             +1 ! dic13
# endif
# if defined O_carbon_14
     $             +1 ! dic14
# endif
#endif
#if defined O_cfcs_data || defined O_cfcs_data_transient
     $             +2 !
#endif
#if defined O_npzd_alk
     $             +1 ! alk
#endif
#if defined O_npzd_o2
     $             +1 ! 
#endif
#if defined O_npzd
     $             +4 ! po4, phyt, zoop, detr
# if defined O_npzd_caco3
     $             +2 ! cocc, caco3
#  if defined O_kk_ballast
     $             +1 !
#  endif
# endif
# if defined O_npzd_nitrogen
     $             +4 ! no3, diaz, dop, don
#  if defined O_npzd_nitrogen_15
     $             +6 ! din15, phytn15, zoopn15, detrn15, diazn15, don15
#   if defined O_npzd_caco3
     $             +1 ! coccn15
#   endif		 
#  endif
# endif
# if defined O_carbon_13
     $             +3 ! phytc13, zoopc13, detrc13
#  if defined O_npzd_caco3
     $             +2 ! coccc13, caco3c13
#  endif		 
#  if defined O_npzd_nitrogen
     $             +2 ! diazc13, doc13
#  endif
# endif
# if defined O_npzd_iron
     $               +2 ! dfe, detrfe
# endif
#endif
     $               )
      parameter (nsrc=0
#if defined O_carbon
     $               +1 ! dic
# if defined O_carbon_13
     $               +1 ! dic13
# endif
# if defined O_carbon_14
     $               +1 ! dic14
# endif
#endif
#if defined O_npzd_alk
     $               +1 ! alk
#endif
#if defined O_npzd_o2
     $               +1 ! o2
#endif
#if defined O_npzd
     $               +4 ! po4, phyt, zoop, detr
# if defined O_kk_ballast
     $               +1 
# endif
# if defined O_npzd_caco3
     $               +2 ! cocc, caco3
# endif
# if defined O_npzd_nitrogen
     $               +4 ! no3, diaz, dop, don
#  if defined O_npzd_nitrogen_15
     $               +6 ! din15, phytn15, zoopn15, detrn15, diazn15, don15
#   if defined O_npzd_caco3
     $               +1 ! coccn15
#   endif		 
#  endif
# endif
# if defined O_carbon_13
     $               +3 ! phytc13, zoopc13, detrc13
#  if defined O_npzd_caco3
     $               +2 ! coccc13, caco3c13
#  endif
#  if defined O_npzd_nitrogen
     $               +2 ! diazc13, doc13
#  endif
# endif
# if defined O_npzd_iron
     $               +2 ! dfe, detrfe
# endif
#endif
     $                 )
      parameter (kpzd=km)

      parameter (nat=2
#if defined O_carbon && defined O_carbon_co2_2d
     $              +1
#endif
     $, jmz=jmt, jmzm1=jmz-1)
      parameter (mnisle=50, maxipp=5000)

#if !defined O_min_window
      parameter (jmw=jmt)
#else

!     for UNI-TASKING: "jmw" is set to the minimum for each option class
!     "jmw" may be increased up to "jmt"

# if defined O_fourth_order_tracer_advection || defined O_fct || defined O_quicker || defined O_pressure_gradient_average || defined O_biharmonic
#  if defined O_pressure_gradient_average
#   if defined O_biharmonic  || defined O_fourth_order_tracer_advection || defined O_fct || defined O_quicker
      parameter (jmw=5)
#   else
      parameter (jmw=4)
#   endif
#  else
      parameter (jmw=4)
#  endif
# else
      parameter (jmw=3)
# endif
#endif

!-----------------------------------------------------------------------
!     set first and last calculated row within the MW. other rows
!     are used as buffers
!-----------------------------------------------------------------------

!     jsmw   = 1st calculated row within the MW
!     jemw   = last calculated row within the MW
#ifndef O_TMM
      parameter (jsmw=2, jemw=jmw-1)
#else
      parameter (jsmw=1, jemw=1)
#endif
! Moses-Triffid land model

! POINTS = Maximum number of points in grid.
! STEPSM = Maximum number of timesteps in a day.
! klmax = maximum ocean depth levels over which the land model can exist

      integer POINTS, STEPSM, klmax
      parameter (POINTS=14300, STEPSM=24, klmax=0)

! NNVG  = Number of non-vegetation surface types.
! NPFT  = Number of plant functional types.
! NTYPE = Number of surface types.
! SOIL  = Index of the surface type 'Soil'
! Land surface types :
!     1 - Broadleaf Tree
!     2 - Needleleaf Tree
!     3 - C3 Grass
!     4 - C4 Grass
!     5 - Shrub
!     6 - Soil

      integer NNVG, NPFT, NTYPE, SOIL
      parameter (NNVG=4, NPFT=5, NTYPE=6, SOIL=6)
