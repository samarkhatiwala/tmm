!====================== include file "scalar.h" ========================

!     various scalar quantities:

!     dtts   = time step for density & tracers (in seconds)
!     dtuv   = time step for baroclinic velocity (in seconds)
!     dtsf   = time step for barotropic velocity (in seconds)
!     c2dtts = 2*dtts
!     c2dtuv = 2*dtuv
!     c2dtsf = 2*dtsf
!     acor   = (>0, 0) = (implicit, explicit) treatment of Coriolis
!               term for internal and external modes.
!     rho0   = mean density for Bousinessq approximation
!     rho0r  = 1/rho0
!     omega  = earth`s rotation rate (radians/sec)
!     radius = earth`s radius (cm)
!     grav   = earth`s gravitational acceleration (cm/sec**2)
!     cdbot  = bottom drag coefficient
!     ncon   = number of  passes through convective code in tracer
!     gcor   = time centring for Coriolis term

!     taux0  = constant zonal windstress (dynes/cm**2) for idealized
!              equatorial studies
!     tauy0  = constant meridional windstress (dynes/cm**2) for
!              idealized equatorial studies

      integer ncon
      common /scalar_i/ ncon

      real dtts, dtuv, dtsf, c2dtts, c2dtuv, c2dtsf, acor, rho0
      real rho0r, omega, radius, grav, cdbot, gcor, taux0, tauy0
      common /scalar_r/ dtts, dtuv, dtsf, c2dtts, c2dtuv, c2dtsf, acor
      common /scalar_r/ rho0, rho0r, omega, radius, grav, cdbot, gcor
      common /scalar_r/ taux0, tauy0

!     various non dimensional quantities:

!     radian = degrees per radian
!     pi     = something good to eat

      real radian, pi
      common /ndcon/ radian, pi
