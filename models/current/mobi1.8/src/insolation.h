!===================== include file "insolation.h" =====================

!     parameters used in calculating solar insolation

!     eccen  = orbital eccentricity
!     obliq  = obliquity (degrees)
!     mvelp  = moving vernal equinox longitude of perihelion (degrees)
!     lambm0 = Mean long of perihelion at vernal equinox (radians)
!     sindec = sine of solar declination angle in rad
!     eccf   = Earth-sun distance factor (ie. (1/r)**2)

      real eccen, obliq , mvelp, lambm0, sindec, eccf

      common /insolation_r/ eccen, obliq , mvelp, lambm0, sindec, eccf
