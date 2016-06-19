      real*8 FUNCTION sw_dens0(S,T)

C %           Copyright (C) CSIRO, Phil Morgan 1992
C %
C % USAGE:  dens0 = sw_dens0(S,T)
C %
C % DESCRIPTION:
C %    Density of Sea Water at atmospheric pressure using
C %    UNESCO 1983 (EOS 1980) polynomial.
C %
C % INPUT:  (all must have same dimensions)
C %   S = salinity    [psu      (PSS-78)]
C %   T = temperature [degree C (ITS-90)]
C %
C % OUTPUT:
C %   dens0 = density  [kg/m^3] of salt water with properties S,T,
C %           P=0 (0 db gauge pressure)
C %
C % AUTHOR:  Phil Morgan 92-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)

C % REFERENCES:
C %     Unesco 1983. Algorithms for computation of fundamental properties of
C %     seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
C %
C %     Millero, F.J. and  Poisson, A.
C %     International one-atmosphere equation of state of seawater.
C %     Deep-Sea Res. 1981. Vol28A(6) pp625-629.
C %=========================================================================
C 
C % Modifications
C % 03-12-12. Lindsay Pender, Converted to ITS-90.
C 
C % CALLER: general purpose, sw_dens.m
C % CALLEE: sw_smow.m

      IMPLICIT NONE

      real*8 S,T

      real*8 b0,b1,b2,b3,b4,c0,c1,c2,d0
      parameter (b0=8.24493d-1,b1=-4.0899d-3,b2=7.6438d-5)
      parameter (b3=-8.2467d-7,b4=5.3875d-9)
      parameter (c0=-5.72466d-3,c1=+1.0227d-4,c2=-1.6546d-6)
      parameter (d0=4.8314d-4)
      real*8 sw_smow
      real*8 T68

      T68 = T*1.00024
      
C       %     UNESCO 1983 eqn(13) p17.
      
      sw_dens0 = sw_smow(T) 
     &           + (b0 + (b1 + (b2 + (b3 + b4*T68)*T68)*T68)*T68)*S
     &           + (c0 + (c1 + c2*T68)*T68)*S*sqrt(S) + d0*(S**2)

      return
      end

      real*8 FUNCTION sw_dens(S,T,P)

C % SW_DENS    Density of sea water
C %=========================================================================
C % SW_DENS  $Id: sw_dens.m,v 1.1 2003/12/12 04:23:22 pen078 Exp $
C %          Copyright (C) CSIRO, Phil Morgan 1992.
C %
C % USAGE:  dens = sw_dens(S,T,P)
C %
C % DESCRIPTION:
C %    Density of Sea Water using UNESCO 1983 (EOS 80) polynomial.
C %
C % INPUT:  (all must have same dimensions)
C %   S = salinity    [psu      (PSS-78)]
C %   T = temperature [degree C (ITS-90)]
C %   P = pressure    [db]
C %       (P may have dims 1x1, mx1, 1xn or mxn for S(mxn) )
C %
C % OUTPUT:
C %   dens = density  [kg/m^3]
C %
C % AUTHOR:  Phil Morgan 92-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)
C %
C % REFERENCES:
C %    Fofonoff, P. and Millard, R.C. Jr
C %    Unesco 1983. Algorithms for computation of fundamental properties of
C %    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
C %
C %    Millero, F.J., Chen, C.T., Bradshaw, A., and Schleicher, K.
C %    " A new high pressure equation of state for seawater"
C %    Deap-Sea Research., 1980, Vol27A, pp255-264.
C %=========================================================================
C 
C % Modifications
C % 99-06-25. Lindsay Pender, Fixed transpose of row vectors.
C % 03-12-12. Lindsay Pender, Converted to ITS-90.

      IMPLICIT NONE

      real*8 S,T,P

      real*8 sw_dens0,sw_seck
      real*8 densP0,K,Ploc

      densP0 = sw_dens0(S,T)
      K = sw_seck(S,T,P)
      Ploc = P/10.d0 ! convert from db to atm pressure units
      sw_dens = densP0/(1-Ploc/K)

      return
      end

      real*8 FUNCTION sw_smow(T)

C % SW_SMOW    Denisty of standard mean ocean water (pure water)
C %=========================================================================
C % SW_SMOW  $Id: sw_smow.m,v 1.1 2003/12/12 04:23:22 pen078 Exp $
C %          Copyright (C) CSIRO, Phil Morgan 1992.
C %
C % USAGE:  dens = sw_smow(T)
C %
C % DESCRIPTION:
C %    Denisty of Standard Mean Ocean Water (Pure Water) using EOS 1980.
C %
C % INPUT:
C %   T = temperature [degree C (ITS-90)]
C %
C % OUTPUT:
C %   dens = density  [kg/m^3]
C %
C % AUTHOR:  Phil Morgan 92-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)
C %
C % DISCLAIMER:
C %   This software is provided "as is" without warranty of any kind.
C %   See the file sw_copy.m for conditions of use and licence.
C %
C % REFERENCES:
C %     Unesco 1983. Algorithms for computation of fundamental properties of
C %     seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
C %     UNESCO 1983 p17  Eqn(14)
C %
C %     Millero, F.J & Poisson, A.
C %     INternational one-atmosphere equation of state for seawater.
C %     Deep-Sea Research Vol28A No.6. 1981 625-629.    Eqn (6)
C %=========================================================================
C 
C % Modifications
C % 99-06-25. Lindsay Pender, Fixed transpose of row vectors.
C % 03-12-12. Lindsay Pender, Converted to ITS-90.

      IMPLICIT NONE

      real*8 T

      real*8 a0,a1,a2,a3,a4,a5      
      parameter (a0=999.842594d0,a1=6.793952d-2,a2=-9.095290d-3)
      parameter (a3=1.001685d-4,a4=-1.120083d-6,a5=6.536332d-9)
      real*8 T68
C %----------------------
C % DEFINE CONSTANTS
C %----------------------


      T68 = T*1.00024d0
      sw_smow = a0 + (a1 + (a2 + (a3 + (a4 + a5*T68)*T68)*T68)*T68)*T68

      return
      end

      real*8 FUNCTION sw_seck(S,T,P)

C % SW_SECK    Secant bulk modulus (K) of sea water
C %=========================================================================
C % SW_SECK  $Id: sw_seck.m,v 1.1 2003/12/12 04:23:22 pen078 Exp $
C %          Copyright (C) CSIRO, Phil Morgan 1992.
C %
C % USAGE:  dens = sw_seck(S,T,P)
C %
C % DESCRIPTION:
C %    Secant Bulk Modulus (K) of Sea Water using Equation of state 1980.
C %    UNESCO polynomial implementation.
C %
C % INPUT:  (all must have same dimensions)
C %   S = salinity    [psu      (PSS-78) ]
C %   T = temperature [degree C (ITS-90)]
C %   P = pressure    [db]
C %       (alternatively, may have dimensions 1*1 or 1*n where n is columns in S)
C %
C % OUTPUT:
C %   K = Secant Bulk Modulus  [bars]
C %
C % AUTHOR:  Phil Morgan 92-11-05, Lindsay Pender (Lindsay.Pender@csiro.au)
C %
C % DISCLAIMER:
C %   This software is provided "as is" without warranty of any kind.
C %   See the file sw_copy.m for conditions of use and licence.
C %
C % REFERENCES:
C %    Fofonoff, P. and Millard, R.C. Jr
C %    Unesco 1983. Algorithms for computation of fundamental properties of
C %    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
C %    Eqn.(15) p.18
C %
C %    Millero, F.J. and  Poisson, A.
C %    International one-atmosphere equation of state of seawater.
C %    Deep-Sea Res. 1981. Vol28A(6) pp625-629.
C %=========================================================================
C 
C % Modifications
C % 99-06-25. Lindsay Pender, Fixed transpose of row vectors.
C % 03-12-12. Lindsay Pender, Converted to ITS-90.

      IMPLICIT NONE

      real*8 S,T,P

      real*8 h0,h1,h2,h3
      parameter (h3=-5.77905d-7,h2=+1.16092d-4,h1=+1.43713d-3,
     &           h0=+3.239908d0)

      real*8 k0,k1,k2
      parameter (k2=5.2787d-8,k1=-6.12293d-6,k0=+8.50935d-5)

      real*8 e0,e1,e2,e3,e4
      parameter (e4=-5.155288d-5,e3=+1.360477d-2,e2=-2.327105d0,
     &          e1=+148.4206d0,e0=19652.21d0)

      real*8 j0      
      parameter (j0=1.91075d-4)

      real*8 i0,i1,i2
      parameter (i2=-1.6078d-6,i1=-1.0981d-5,i0=2.2838d-3)
      
      real*8 m0,m1,m2
      parameter (m2=9.1697d-10,m1=+2.0816d-8,m0=-9.9348d-7)
      
      real*8 f0,f1,f2,f3
      parameter (f3=-6.1670d-5,f2=+1.09987d-2,f1=-0.603459d0,
     &           f0=+54.6746d0)
     
      real*8 g0,g1,g2
      parameter (g2=-5.3009d-4,g1=+1.6483d-2,g0=+7.944d-2)
            
      real*8 T68,AW,BW,KW,SR,A,B,K00,Ploc
      
C       %--------------------------------------------------------------------
C       % COMPUTE COMPRESSION TERMS
C       %--------------------------------------------------------------------
      Ploc = P/10.d0 !  %convert from db to atmospheric pressure units
      T68 = T*1.00024d0
      
C       % Pure water terms of the secant bulk modulus at atmos pressure.
C       % UNESCO eqn 19 p 18
      
      AW  = h0 + (h1 + (h2 + h3*T68)*T68)*T68
      
      BW  = k0 + (k1 + k2*T68)*T68
            
      KW  = e0 + (e1 + (e2 + (e3 + e4*T68)*T68)*T68)*T68   ! eqn 19
      
C       %--------------------------------------------------------------------
C       % SEA WATER TERMS OF SECANT BULK MODULUS AT ATMOS PRESSURE.
C       %--------------------------------------------------------------------
      
      SR = sqrt(S)
      
      A  = AW + (i0 + (i1 + i2*T68)*T68 + j0*SR)*S
            
      B = BW + (m0 + (m1 + m2*T68)*T68)*S ! eqn 18
      
      K00 = KW + (  f0 + (f1 + (f2 + f3*T68)*T68)*T68 
     &     +   (g0 + (g1 + g2*T68)*T68)*SR)*S ! eqn 16
      
      sw_seck = K00 + (A + B*Ploc)*Ploc ! eqn 15

      return
      end
