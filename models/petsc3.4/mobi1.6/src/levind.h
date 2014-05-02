!====================== include file "levind.h" ========================

!     vertical level indicators which define model geometry & bottom
!     topography:

!     kmt = number of vertical boxes over "t" points
!     kmu = number of vertical boxes over "u,v" points

      integer kmt, kmu
      real sg_bathy

      common /levind/ kmt(imt,jmt), kmu(imt,jmt)
      common /levind/ sg_bathy(imt,jmt,km)
