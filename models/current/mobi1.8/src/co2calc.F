      subroutine co2calc_SWS (t, s, dic_in, ta_in, co2_in, atmpres
     &,                       depth, ph, co2star, dco2star, pCO2, dpco2
     &,                       CO3, Omega_c, Omega_a)

!-------------------------------------------------------------------------

! Modified from co2calc.f (RCS version 1.8, OCMIP-2)
! - by A. Mouchet, 2004:

! NOW; "All" constants are given on seawater H scale (hSWS)
! - HOWEVER, dissociation constants for S and F are on the 'free' H scale
!            (this is necessary to work with hSWS)
! - this routine corrects for inconsistencies in the previous version.

! - Other changes:
!   * use ta_iter_SWS instead of ta_iter_1;
!   * hSWS replaces htotal;
!   * moved upward the calculation of bt, st and ft
!     (needed in evaluation of kb);
!   * added various comments and references.

! subroutine CO2CALC_SWS

! PURPOSE
!        Calculate delta co2* from total alkalinity and total CO2 at
! temperature (t), salinity (s) and "atmpres" atmosphere total pressure.

! USAGE
!       call co2calc_SWS(t,s,dic_in,ta_in,pt_in,sit_in
!    &                  ,phlo,phhi,ph,co2_in,atmpres
!    &                  ,co2star,dco2star,pCO2,dpco2)

! INPUT
!        dic_in  = total inorganic carbon (mol/m^3)
!                  where 1 T = 1 metric ton = 1000 kg
!        ta_in   = total alkalinity (eq/m^3)
!        pt_in   = inorganic phosphate (mol/m^3)
!        sit_in  = inorganic silicate (mol/m^3)
!        t       = temperature (degrees C)
!        s       = salinity (PSU)
!        phlo    = lower limit of pH range
!        phhi    = upper limit of pH range
!        co2_in  = atmospheric mole fraction CO2 in dry air (ppmv)
!        atmpres = atmospheric pressure in atmospheres (1 atm==1013.25mbar)
!        depth   = ocean depth (m)

!       Note: arguments dic_in, ta_in, pt_in, sit_in, and co2_in are
!             used to initialize variables dic, ta, pt, sit, and co2.
!             * Variables dic, ta, pt, and sit are in the common block
!               "species".
!             * Variable co2 is a local variable.
!             * Variables with "_in" suffix have different units
!               than those without.
!             pressure dependency added (Millero 1995).

! OUTPUT
!        co2star  = CO2*water (mol/m^3)
!        dco2star = delta CO2 (mol/m^3)
!        pco2     = oceanic pCO2 (ppmv)
!        dpco2    = Delta pCO2, i.e, pCO2ocn - pCO2atm (ppmv)

! IMPORTANT: Some words about units - (JCO, 4/4/1999)
!     - Models carry tracers in mol/m^3 (on a per volume basis)
!     - Conversely, this routine, which was written by observationalists
!       (C. Sabine and R. Key), passes input arguments in umol/kg
!       (i.e., on a per mass basis)
!     - I have changed things slightly so that input arguments are in mol/m^3,
!     - Thus, all input concentrations (dic_in, ta_in, pt_in, and st_in)
!       should be given in mol/m^3; output arguments "co2star" and "dco2star"
!       are likewise in mol/m^3.

! FILES and PROGRAMS NEEDED
!        drtsafe
!        ta_iter_SWS 

!--------------------------------------------------------------------------

      implicit none

      real permil, pt, pt_in, sit, sit_in, ta, ta_in, dic, dic_in
      real permeg, co2, co2_in, tk, t, tk100, tk1002, dlogtk, s
      real sqrtis, s2, sqrts, s15, scl, bt, st, ft, ff, x1, phhi
      real x2, phlo, xacc, hSWS, drtsafe, hSWS2, co2star, co2starair
      real atmpres, dco2star, ph, pco2, dpco2, invtk, is, is2
      real k0, k1, k12, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi, C2K
      real Kspc, Kspa, ca, co3, omega_c, omega_a, depth, pres
      real t2, dvc, dva, dk, pitkR, p2itkR, b_x, delta_x, FugFac, rt_x

      common /const/ k0, k1, k12, k2, kw, kb, ks, kf, k1p, k2p, k3p
      common /const/ ksi, ff
      common /species/ bt, st, ft, sit, pt, dic, ta

!     ---------------------------------------------------------------------
!     Change units from the input of mol/m^3 -> mol/kg:
!     (1 mol/m^3)  x (1 m^3/1024.5 kg)
!     where the ocean's mean surface density is 1024.5 kg/m^3
!     Note: mol/kg are actually what the body of this routine uses
!     for calculations.
!     ---------------------------------------------------------------------
      permil = 1.0 / 1024.5
!     To convert input in mol/m^3 -> mol/kg
      pt = pt_in*permil
      sit = sit_in*permil
      ta = ta_in*permil
      dic = dic_in*permil
      C2K = 273.15
      pres = depth*0.1  ! [decibars ~ meters of depth to bars]

!     ---------------------------------------------------------------------
!     Change units from uatm to atm. That is, atm is what the body of
!     this routine uses for calculations.
!     ---------------------------------------------------------------------
      permeg = 1.e-6
!     To convert input in uatm -> atm
      co2 = co2_in*permeg
!     ---------------------------------------------------------------------

!*************************************************************************
! Calculate all constants needed to convert between various measured
! carbon species. References for each equation are noted in the code.
! Version 2 of "Handbook of Methods for the Analysis of the Various
! Parameters of the Carbon Dioxide System in Seawater", DOE, 1994
! (SOP No. 3, p25-26).

! Derive simple terms used more than once

      tk = C2K + t
      tk100 = tk/100.0
      tk1002 = tk100*tk100
      invtk = 1.0/tk
      dlogtk = log(tk)
      is = 19.924*s/(1000.-1.005*s)
      is2 = is*is
      sqrtis = sqrt(is)
      s2 = s*s
      t2 = t*t
      sqrts = sqrt(s)
      s15 = s**1.5
      scl = s/1.80655
      pitkR = pres/tk/83.15 ! % 83.15 = the gas constant R in cm3 bar k-1
      p2itkR = pres*pitkR

!------------------------------------------------------------------------
! Calculate concentrations for borate, sulfate, and fluoride

! Uppstrom (1974)
      bt = 0.000232 * scl/10.811
! Morris & Riley (1966)
      st = 0.14 * scl/96.062
! Riley (1965)
      ft = 0.000067 * scl/18.9984

!------------------------------------------------------------------------
! f = k0(1-pH2O)*correction term for non-ideality
! Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)

      ff = exp(-162.8301 + 218.2968/tk100  + 90.9241*log(tk100)
     &   - 1.47696*tk1002 + s*(.025695 - .025225*tk100
     &   + 0.0049867*tk1002))

! K0 from Weiss 1974
      k0 = exp(93.4517/tk100 - 60.2409 + 23.3585 * log(tk100)
     &   + s*(.023517 - 0.023656 * tk100 + 0.0047036 * tk1002))
!
! Now calculate FugFac according to Weiss (1974) Marine Chemistry
! delta_x and b_x are in cm3/mol
      rt_x = 83.1451*tk !rt_x = the gas constant times the temperature
      delta_x = (57.7 - 0.118*tk)
      b_x = -1636.75 + 12.0408*tk - 0.0327957*tk*tk
      b_x = b_x + 3.16528*1e-5*tk*tk*tk
      FugFac = exp((b_x+2*delta_x)*1/rt_x) ! Note that 1 is the atmospheric pressure in bars
!
# if defined k0_as_solubility
      k0 = exp(93.4517/tk100 - 60.2409 + 23.3585 * log(tk100)
     &   + s*(.023517 - 0.023656 * tk100 + 0.0047036 * tk1002))
# endif
!------------------------------------------------------------------------
! k1 = [H][HCO3]/[H2CO3]
! k2 = [H][CO3]/[HCO3]     on hSWS
! Millero p.664 (1995) using Mehrbach et al. data on SEAWATER scale
! (Original reference: Dickson and Millero, DSR, 1987)

      k1 = 10**(-1.*(3670.7*invtk - 62.008 + 9.7944*dlogtk
     &   - 0.0118*s + 0.000116*s2))
# if defined O_save_carbon_carbonate_depth || O_carbon_13
     &     *exp((25.5 - 0.1271*t)*pitkR
     &   + 0.5*(-3.08e-3 + 8.77e-5*t)*p2itkR)
# endif

      k2 = 10**(-1*(1394.7*invtk + 4.777 - 0.0184*s + 0.000118*s2))
# if defined O_save_carbon_carbonate_depth || O_carbon_13
     &     *exp((15.82 + 0.0219*t)*pitkR
     &   + 0.5*(1.13e-3 - 1.475e-4*t)*p2itkR)
# endif

!------------------------------------------------------------------------
! k1p = [H][H2PO4]/[H3PO4] on hSWS
! Millero p.670 (1995)

      k1p = exp(-4576.752*invtk + 115.540 - 18.453*dlogtk
     &    + (-106.736*invtk + 0.69171)*sqrts
     &    + (-0.65643*invtk - 0.01844)*s)
# if defined O_save_carbon_carbonate_depth || O_carbon_13
     &      *exp((14.51 - 0.1211*t + 3.21e-4*t2)*pitkR
     &    + 0.5*(-2.67e-3 + 4.27e-5*t)*p2itkR)
# endif

!------------------------------------------------------------------------
! k2p = [H][HPO4]/[H2PO4] on hSWS
! Millero p.670 (1995)

      k2p = exp(-8814.715*invtk + 172.1033 - 27.927*dlogtk
     &    + (-160.340*invtk + 1.3566)*sqrts
     &    + (0.37335*invtk - 0.05778)*s)
# if defined O_save_carbon_carbonate_depth || O_carbon_13
     &      *exp((23.12 - 0.1758*t + 2.647e-3*t2)*pitkR
     &    + 0.5*(-5.15e-3 + 9.0e-5*t)*p2itkR)
# endif

!------------------------------------------------------------------------
! k3p = [H][PO4]/[HPO4] on hSWS
! Millero p.670 (1995)

      k3p = exp(-3070.75*invtk - 18.126
     &    + (17.27039*invtk + 2.81197)*sqrts
     &    + (-44.99486*invtk - 0.09984)*s)
# if defined O_save_carbon_carbonate_depth || O_carbon_13
     &      *exp((26.57 - 0.202*t + 3.042e-3*t2)*pitkR
     &    + 0.5*(-4.08e-3 + 7.14e-5*t)*p2itkR)
# endif

!------------------------------------------------------------------------
! ksi = [H][SiO(OH)3]/[Si(OH)4] on hSWS
! Millero p.671 (1995) using data from Yao and Millero (1995)
! change to (mol/ kg soln)
! depth dependancy assumed to be the same as boric acid
! typo in Millero 1994 corrected in sign of 0.1622

       ksi = exp(-8904.2*invtk + 117.400 - 19.334*dlogtk
     &     + (-458.79*invtk + 3.5913)*sqrtis
     &     + (188.74*invtk - 1.5998)*is
     &     + (-12.1652*invtk + 0.07871)*is2
     &     + log(1.0-0.001005*s))
# if defined O_save_carbon_carbonate_depth || O_carbon_13
     &       *exp((29.48 - 0.1622*t - 2.608e-3*t2)*pitkR
     &     + 0.5*(-2.84e-3)*p2itkR)
# endif

!------------------------------------------------------------------------
! kw = [H][OH] on hSWS
! Millero p.670 (1995) using composite data
! pressure dependancy in Millero 1994 corrected for sea water from
! Millero 1983

       kw = exp(-13847.26*invtk + 148.9802 - 23.6521*dlogtk
     &    + (118.67*invtk - 5.977 + 1.0495*dlogtk)*sqrts - 0.01615*s)
# if defined O_save_carbon_carbonate_depth || O_carbon_13
     &      *exp((20.02 - 0.1119*t + 1.409e-3*t2)*pitkR
     &    + 0.5*(-5.13e-3 + 7.94e-5*t)*p2itkR)
# endif

!------------------------------------------------------------------------
! ks = [H][SO4]/[HSO4] on free H scale
! Dickson (1990, J. chem. Thermodynamics 22, 113)
! change to (mol/ kg soln)

      ks = exp(-4276.1*invtk + 141.328 - 23.093*dlogtk
     &   + (-13856*invtk + 324.57 - 47.986*dlogtk)*sqrtis
     &   + (35474*invtk - 771.54 + 114.723*dlogtk)*is
     &   - 2698*invtk*is**1.5 + 1776*invtk*is2
     &   + log(1.0 - 0.001005*s))
# if defined O_save_carbon_carbonate_depth || O_carbon_13
     &     *exp((18.03 - .0466*t - 3.16e-4*t2)*pitkR
     &   + 0.5*(-4.53e-3 + 9.0e-5*t)*p2itkR)
# endif

!------------------------------------------------------------------------
! kf = [H][F]/[HF] on free H scale
! Dickson and Riley (1979)
! change to (mol/ kg soln)

      kf = exp(1590.2*invtk - 12.641 + 1.525*sqrtis
     &   + log(1.0 - 0.001005*s))
# if defined O_save_carbon_carbonate_depth || O_carbon_13
     &     *exp((9.78 + 9.0e-3*t + 9.42e-4*t2)*pitkR
     &   + 0.5*(-3.91e-3 + 5.4e-5*t)*p2itkR)
# endif

!------------------------------------------------------------------------
! kb = [H][BO2]/[HBO2] on hSWS
! Dickson p.673 (1990)
! change from htotal to hSWS
! typo in Millero 1994 corrected in sign of 0.1622

      kb = exp((-8966.90 - 2890.53*sqrts - 77.942*s
     &   + 1.728*s15 - 0.0996*s2)*invtk
     &   + (148.0248 + 137.1942*sqrts + 1.62142*s)
     &   + (-24.4344 - 25.085*sqrts - 0.2474*s)*dlogtk
     &   + 0.053105*sqrts*tk
     &   + log((1+(st/ks)+(ft/kf))/(1+(st/ks))))
# if defined O_save_carbon_carbonate_depth || O_carbon_13
     &     *exp((29.48 - 0.1622*t - 2.608e-3*t2)*pitkR
     &   + 0.5*(-2.84e-3)*p2itkR)
# endif

!*************************************************************************

! Calculate [H+] SWS when DIC and TA are known at T, S and 1 atm.
! The solution converges to err of xacc. The solution must be within
! the range x1 to x2.

! If DIC and TA are known then either a root finding or iterative method
! must be used to calculate hSWS. In this case we use the Newton-Raphson
! "safe" method taken from "Numerical Recipes" (function "rtsafe.f" with
! error trapping removed).

! As currently set, this procedure iterates about 12 times. The x1 and x2
! values set below will accomodate ANY oceanographic values. If an initial
! guess of the pH is known, then the number of iterations can be reduced to
! about 5 by narrowing the gap between x1 and x2. It is recommended that
! the first few time steps be run with x1 and x2 set as below. After that,
! set x1 and x2 to the previous value of the pH +/- ~0.5. The current
! setting of xacc will result in co2star accurate to 3 significant figures
! (xx.y). Making xacc bigger will result in faster convergence also, but this
! is not recommended (xacc of 10**-9 drops precision to 2 significant figures).

! Hardwire constants (Andreas Schmittner Oct 3, 2013)
      phhi = 6.
      phlo = 10.
      sit_in = 7.6875e-03       !mol/m^3
      pt_in = 0.5125e-3         !mol/m^3

! Parentheses added around negative exponents (Keith Lindsay)

      x1 = 10.0**(-phhi)
      x2 = 10.0**(-phlo)
      xacc = 1.e-10
      hSWS = drtsafe(x1, x2, xacc)

! Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2,
! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)

      hSWS2 = hSWS*hSWS
      co2star = dic*hSWS2/(hSWS2 + k1*hSWS + k1*k2)
      co2starair = co2*ff*atmpres
      dco2star = co2starair - co2star
      ph = -log10(hSWS)

! Add two output arguments for storing pCO2
! Should we be using K0 or ff for the solubility here?
# if defined k0_as_solubility
      pCO2 = co2star/k0
# else
      pCO2 = co2star/(k0*FugFac)
c      pCO2 = co2star/ff ! if not O_pCO2_bug_fix
# endif
      dpCO2 = pCO2 - co2starair
c      dpCO2 = pCO2 - co2*atmpres ! if not O_pCO2_bug_fix

# if defined O_save_carbon_carbonate_chem
! Calculate [CO3] in mol/kg from co2star and hSWS
      CO3 = k1*k2*co2star/hSWS2

! Solubility products of calcite and aragonite at sea level
! (Ksp0) as a function of temperature  and salinity. Based on
! Mucci, 1983 as written by Sarmiento and Gruber, 2006 p366

! Ksp0 for calcite:
      Kspc = exp(-395.8293 + (6537.773/tk) + 71.595*alog(tk)
     &     - 0.17959*tk + (-1.78938 + (410.64/tk)
     &     + 0.0065453*tk)*s**0.5 - 0.17755*s + 0.0094979*s15)

! Ksp0 for aragonite:
      Kspa = exp(-395.9180 + (6685.079/tk) + 71.595*log(tk)
     &     - 0.17959*tk + (-0.157481 + (202.938/tk)
     &     + 0.0039780*tk)*s**0.5 - 0.23067*s + 0.0136808*s15)

# endif
# if defined O_save_carbon_carbonate_depth
!------------------------------------------------------------------------
! Pressure dependent solubility products (KspP) based on Millero, 1983 as
! written by Sarmiento and Gruber, 2006 p366

! Molal Volume for calcite and aragonite (temperature in Celsius)
      DVc = -65.28 + 0.397*t - 0.005155*t**2
     &    + (19.816 - 0.0441*t - 0.00017*t**2)*(s/35.)**(0.5)

      DVa = -65.50 + 0.397*t - 0.005155*t**2
     &    + (19.82 - 0.0441*t - 0.00017*t**2)*(s/35.)**(0.5)

! Compressability for calcite and aragonite (temperature in Celsius)
      DK = 0.01847 + 0.0001956*t - 0.000002212*t**2
     &   + (-0.03217 - 0.0000711*t + 0.000002212)*(s/35.)**(0.5)

! KspP for calcite, pressure (pres) in bars
      Kspc = Kspc*exp(-DVc*pitkR + 0.5*DK*p2itkR)

! KspP for aragonite, pressure (pres) in bars
      Kspa = Kspa*exp(-DVa*pitkR + 0.5*DK*p2itkR)

# endif
# if defined O_save_carbon_carbonate_chem
! Calculate saturation state (Omega) using a fix value
! for [Ca] in mol/kg
      Ca = 10.28E-3

! calcite (Omgc)
      Omega_c = Ca*CO3/Kspc

! aragonite (Omga)
      Omega_a = Ca*CO3/Kspa

# endif
! Convert units of output arguments
! Note: co2star and dco2star are calculated in mol/kg within this routine
! Thus Convert now from mol/kg -> mol/m^3
      co2star  = co2star/permil
      dco2star = dco2star/permil
      CO3 = CO3/permil

! Note: pCO2 and dpCO2 are calculated in atm above.
! Thus convert now to uatm
      pCO2 = pCO2/permeg
      dpCO2 = dpCO2/permeg

      return
      end

      real function drtsafe(x1, x2, xacc)

      implicit none

      integer maxit, j
      real x1, fl, df, x2, fh, xl, xh, swap, dxold, dx, f, temp, xacc

!     File taken from Numerical Recipes. Modified  R.M.Key 4/94

      maxit = 100
      call ta_iter_SWS (x1, fl, df)
      call ta_iter_SWS (x2, fh, df)
      if (fl .lt. 0.0) then
        xl = x1
        xh = x2
      else
        xh = x1
        xl = x2
        swap = fl
        fl = fh
        fh = swap
      endif
      drtsafe = 0.5*(x1 + x2)
      dxold = abs(x2 - x1)
      dx = dxold
      call ta_iter_SWS (drtsafe, f, df)
      do J=1,maxit
        if (((drtsafe - xh)*df - f)*((drtsafe - xl)*df - f) .ge. 0.
     &    .or. abs(2.0*f) .gt. abs(dxold*df)) then
          dxold=  dx
          dx = 0.5*(xh-xl)
          drtsafe = xl + dx
          if (xl .eq. drtsafe) return
        else
          dxold = dx
          dx = f/df
          temp = drtsafe
          drtsafe = drtsafe - dx
          if (temp .eq. drtsafe) return
        endif
        if (abs(dx) .lt. xacc) return
        call ta_iter_SWS (drtsafe, f, df)
        if (f .lt. 0.0) then
          xl = drtsafe
          fl = f
        else
          xh = drtsafe
          fh = f
        endif
      enddo

      return
      end

      subroutine ta_iter_SWS (x, fn, df)

      implicit none

      real x2, x, x3, c, st, ft, a, a2, da, b, b2, db, fn, dic, bt
      real pt, sit, ta, df, ff, k12p, k123p
      real k0, k1, k12, k2, kw, kb, ks, kf, k1p, k2p, k3p, ksi

      common /const/ k0, k1, k12, k2, kw, kb, ks, kf, k1p, k2p, k3p
      common /const/ ksi, ff
      common /species/ bt, st, ft, sit, pt, dic, ta

! Modified from ta_iter_1.f (RCS version 1.2, OCMIP-2)
! - by A. Mouchet, 2004:
! Fixed Problems w/ version of ta_iter_1.f used in OCMIP-2 (vers. 1.2)
!  1) fixed errors in signs, parenthesis and coefficient c in derivative
!  2) changed from Total to Seawater Scale
!     * c defined for seawater H scale;
!     * fn and df adapted to KF on free H scale
!     * comments have been adapted

! This routine expresses TA as a function of DIC, hSWS and constants.
! It also calculates the derivative of this function with respect to
! hSWS. It is used in the iterative solution for hSWS. In the call
! "x" is the input value for hSWS, "fn" is the calculated value for TA
! and "df" is the value for dTA/dhSWS

      x2=x*x
      x3=x2*x
      k12 = k1*k2
      k12p = k1p*k2p
      k123p = k12p*k3p
      c = 1.0 + st/ks + ft/kf
      a = x3 + k1p*x2 + k12p*x + k123p
      a2=a*a
      da = 3.0*x2 + 2.0*k1p*x + k12p
      b = x2 + k1*x + k12
      b2=b*b
      db = 2.0*x + k1

!     fn = hco3+co3+borate+oh+hpo4+2*po4+silicate-hfree-hso4-hf-h3po4-ta
!===========================================================================

      fn = k1*x*dic/b +
     &     2.0*dic*k12/b +
     &     bt/(1.0 + x/kb) +
     &     kw/x +
     &     pt*k12p*x/a +
     &     2.0*pt*k123p/a +
     &     sit/(1.0 + x/ksi) -
     &     x/c -
     &     st/(1.0 + ks/(x/c)) -
     &     ft/(1.0 + kf/(x/c)) -
     &     pt*x3/a -
     &     ta

!     df = dfn/dx

      df = ((k1*dic*b) - k1*x*dic*db)/b2 -
     &     2.0*dic*k12*db/b2 -
     &     bt/kb/(1.0+x/kb)**2 -
     &     kw/x2 +
     &     (pt*k12p*(a - x*da))/a2 -
     &     2.0*pt*k123p*da/a2 -
     &     sit/ksi/(1.0+x/ksi)**2 -
     &     1.0/c -
     &     st *(1.0 + ks/(x/c))**(-2) *(ks*c/x2) -
     &     ft*(1.0 + kf/(x/c))**(-2) *(kf*c/x2) -
     &     pt*x2*(3.0*a-x*da)/a2

      return
      end
