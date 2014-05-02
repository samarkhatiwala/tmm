!======================= include file "param.h" ========================

!     main parameter file which sets ocean characteristics:

!     nvar   = number of prognostic variables
!     lseg   = maximum number of longitudinal stream function segments
!     nlatpr = maximum number of latitudes for matrix printouts
!              on diagnostic time steps
!     nhreg  = number of regions in the horizontal used for averaging
!              tracers.
!     nvreg  = number of regions in the vertical used for term balance
!              calculations. note "nvreg" is not used for tracer
!              averages
!     numreg = total number of regions ( = product of nhreg & nvreg)
!              used for term balance calculations

!     nvarbh = number of prognostic variables using biharmonic mixing

!     ncrows = number of calculated rows within the MW.
!              (the remaining rows are buffer rows).


      integer lseg, nlatpr, nhreg, nvreg, numreg, nvar, nvarbh
      integer imtm1, kmm1, imtp1, imtm2, jmtp1, jmtm1, jmtm2, jscan
      integer kmp1, kmp2, imtkm, nwds, nkflds, nslab, ntmin2, ncrows

      parameter (lseg=5, nlatpr=10)
      parameter (nhreg=3, nvreg=1, numreg=nhreg*nvreg)
      parameter (nvar=nt+2)

#if defined O_biharmonic
      parameter (nvarbh=nt+2)
#endif

      parameter (imtm1=imt-1, kmm1=km-1)
      parameter (imtp1=imt+1, imtm2=imt-2
     &,          jmtp1=jmt+1, jmtm1=jmt-1, jmtm2=jmt-2
#if defined O_symmetry
     &,          jscan=jmtm2+1
#else
     &,          jscan=jmtm2
#endif
     &,          kmp1=km+1, kmp2=km+2
     &,          imtkm=imt*km, nwds=imt*jmt, nkflds=2
     &,          nslab=imt*nvar*km, ntmin2=nt+1/nt)

#if defined O_fourth_order_tracer_advection || defined O_fct || defined O_quicker || defined O_pressure_gradient_average || defined O_biharmonic
# if defined O_pressure_gradient_average
#  if defined O_biharmonic || defined O_fourth_order_tracer_advection || defined O_fct || defined O_quicker
      parameter (ncrows = jmw - 4 + 2*(jmw/jmt))
#  else
      parameter (ncrows = jmw - 3 + jmw/jmt)
#  endif
# else
      parameter (ncrows = jmw - 3 + jmw/jmt)
# endif
#else
      parameter (ncrows = jmw - 2)
#endif
