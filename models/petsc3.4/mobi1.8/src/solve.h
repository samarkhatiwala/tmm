!======================== include file "solve.h" ========================

!     variables needed for solving atmospheric advection and diffusion

!     newcoef = logical flag for calculating new coefficients

      logical newcoef(2,nat)
      common /solver_l/ newcoef

#if defined O_embm_explicit || defined O_embm_explicit_q
!     filter  = filter for east-west diffusion and advection
!     dce     = eastern diffusion coefficient with grid terms
!     dcn     = northern diffusion coefficient with grid terms
!     dfs     = holder for previous northern coefficient
!     ace     = eastern advective coefficient with grid terms
!     acn     = northern advective coefficient with grid terms
!     afs     = holder for previous northern coefficient

      real filter, dce, dcn, dfs
      common /solver_r/ filter(jmt), dce(imt,jmt,nat)
      common /solver_r/ dcn(imt,jmt,nat), dfs(imt)
# if defined O_embm_adv_q
      real ace, acn, afs
      common /solver_r/ ace(imt,jmt,nat), acn(imt,jmt,nat), afs(imt)
# endif
#endif
#if !defined O_embm_explicit || defined O_embm_explicit_q

      integer iimtm2, jjmtm2, nord, nelm
# if defined O_embm_solve2x
      parameter (iimtm2 = imtm2/2)
# else
      parameter (iimtm2 = imtm2)
# endif
# if defined O_embm_solve2y
      parameter (jjmtm2 = jmtm2/2)
# else
      parameter (jjmtm2 = jmtm2)
# endif
      parameter (nord = iimtm2*jjmtm2)

!     itin    = requested maximum iterations
!     itout   = actual iterations
!     newcoef = logical flag for calculating new coefficients
!     bv      = right hand side vector (b)
!     xv      = left hand side vector (x)
!     epsin   = requested maximum error
!     epsout  = actual error

      integer itin(nat), itout(nat)
      common /solver_i/ itin, itout

# if defined O_embm_mgrid
!     for mgrid routine storage is by compass coefficient
!     ap, an, as, ae, aw = centre, north, south, east, and west coef
!     ie. ap*xp = an*xn + as*xs + ae*xe + aw*xw + bp

      integer levelin, levelout
      common /solver_i/ levelin, levelout

      real(kind=8) bv, xv, epsin, epsout, an, as, ae, aw, ap
      common /solver_r/ bv(nord), xv(nord), epsin(nat), epsout(nat)
      common /solver_r/ an(nord,2,nat), as(nord,2,nat), ae(nord,2,nat)
      common /solver_r/ aw(nord,2,nat), ap(nord,2,nat)
# elif defined O_embm_slap
!     for slap routines storage is by row and column index
!     ia     = row index
!     ja     = column index
!     nelm   = number of nonzero elements in A
!     niaux  = size of integer work space
!     iaux   = integer work space
!     nraux  = size of real work space
!     ar     = coefficient matrix (stored form of A)
!     raux   = real work space

      integer nraux, niaux
      parameter (nelm = 3*nord + 2*iimtm2*(jjmtm2 - 1))
      parameter (nraux = nelm + 16*nord + 132)
      parameter (niaux = nelm + 3*nord + 33)

      integer ia, ja, iaux
      common /solver_i/ ia(nelm), ja(nelm), iaux(niaux)

      real bv, xv, epsin, epsout, ar, raux
      common /solver_r/ bv(nord), xv(nord), epsin(nat), epsout(nat)
      common /solver_r/  ar(nelm,2,nat), raux(nraux)
# elif defined O_embm_adi
!     for ADI routine, storage is by compass coefficient
!     (an, ans, as) = north-central-south  coefs.
!     (ae, aew, aw) = east-central-west coefs.
!     ie. (ans+aew)*xc = an*xn + as*xs + ae*xe + aw*xw + bc

      real bv, xv, epsin, epsout, an, ans, as, ae, aew, aw
      common /solver_r/ bv(nord), xv(nord), epsin(nat), epsout(nat)
      common /solver_r/ an(nord,2,nat), ans(nord,2,nat), as(nord,2,nat)
      common /solver_r/ ae(nord,2,nat), aew(nord,2,nat), aw(nord,2,nat)
# elif defined O_embm_essl
!     for ESSL DSRIS routine storage is by rows
!     ia     = index of ar for the first entry of a row
!     ja     = column index
!     iparm  = integer solver parameters
!     nelm   = number of nonzero elements in A
!     naux1  = size of work space 1 (may change if solver type changed)
!     naux2  = size of work space 2 (may change if solver type changed)
!     ar     = coefficient matrix (stored form of A)
!     aux1   = work space 1
!     aux2   = work space 2
!     rparm  = real solver parameters

      integer naux1, naux2
      parameter (nelm = 3*nord + 2*iimtm2*(jjmtm2 - 1))
      parameter (naux1 = 3*nelm + 13*nord + 60)
      parameter (naux2 = 7*nord)

      integer ia, ja, iparm
      common /solve_i/ ia(nord+1), ja(nelm), iparm(6)

      real(kind=8) bv, xv, epsin, epsout, ar, aux1, aux2, rparm
      common /solver_r/ bv(nord), xv(nord), epsin(nat), epsout(nat)
      common /solve_r/ ar(nelm,2,nat), aux1(naux1,2,nat)
      common /solve_r/ aux2(naux2), rparm(3)
# elif defined O_embm_sparskit
!     for storage is by rows
!     ia     = index of ar for the first entry of a row
!     ja     = column index
!     ipar   = integer solver parameters
!     nelm   = number of nonzero elements in A
!     nwork  = size of work space (may change if solver type changed)
!     ar     = coefficient matrix (stored form of A)
!     fpar   = real solver parameters
!     work   = work space

      integer nwork
      parameter (nelm = 3*nord + 2*iimtm2*(jjmtm2 - 1))
      parameter (nwork = 3*nelm + 13*nord + 60)

      integer ia, ipar, ja, jau, ju
      common /solve_i/ ia(nord+1), ipar(16), ja(nelm), jau(nelm)
      common /solve_i/ ju(nelm)

      real(kind=8) bv, xv, epsin, epsout, ar, aur, fpar, work
      common /solver_r/ bv(nord), xv(nord), epsin(nat), epsout(nat)
      common /solve_r/ ar(nelm,2,nat), aur(nelm,2,nat), fpar(16)
      common /solve_r/ work(nwork,2,nat)
# endif
#endif

!     grid terms for the atmospheric solver

      real dwgrd, degrd, azgrd, dsgrd, dngrd, asgrd, angrd
      common /solve_r/ dwgrd(2:imtm1), degrd(2:imtm1), azgrd(2:imtm1)
      common /solve_r/ dsgrd(2:jmtm1), dngrd(2:jmtm1), asgrd(2:jmtm1)
      common /solve_r/ angrd(2:jmtm1)
#if defined O_embm_solve2x
      real wti, xgrd
      common /solve_r/ wti(imt), xgrd(imt)
#endif
#if defined O_embm_solve2y
      real wtj, ygrd
      common /solve_r/ wtj(jmt), ygrd(jmt)
#endif
#if defined O_embm_solve2x || defined O_embm_solve2y

!     grid ratio for coarse grid atmospheric solver

      real gr
      common /solve_r/ gr(imt,jmt)
#endif
