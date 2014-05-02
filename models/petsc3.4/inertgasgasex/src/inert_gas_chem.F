      real*8 FUNCTION gas_mol_fract(gasID)

      IMPLICIT NONE

      integer gasID
      
      if (gasID.eq.1) then !, 'Ne')
        gas_mol_fract = 0.00001818d0
      elseif (gasID.eq.2) then !'Ar')
        gas_mol_fract = 0.00934d0
      elseif (gasID.eq.3) then !'Kr')
        gas_mol_fract = 0.00000114d0
      elseif (gasID.eq.4) then !'Xe')
        gas_mol_fract = 9.d-8
      elseif (gasID.eq.5) then !'N2')
        gas_mol_fract = 0.780840d0
      elseif (gasID.eq.6) then !'Ar36')
        gas_mol_fract = 0.00934d0*0.003651267d0
      elseif (gasID.eq.7) then ! 'O2')
        gas_mol_fract = 0.209460d0
      else
        stop 'Gas name must be Ne, Ar, Kr, Xe, N2, O2 or Ar36'
      endif

      return
      end

      real*8 FUNCTION solfun(S,T,gasID)

      IMPLICIT NONE

      real*8 S,T
      integer gasID

      real*8 Nesol,Arsol,Krsol,Xesol,N2sol,Ar36sol,O2sol
      
      if (gasID.eq.1) then !, 'Ne')
        solfun = Nesol(S,T)
      elseif (gasID.eq.2) then !'Ar')
        solfun = Arsol(S,T)
      elseif (gasID.eq.3) then !'Kr')
        solfun = Krsol(S,T)
      elseif (gasID.eq.4) then !'Xe')
        solfun = Xesol(S,T)
      elseif (gasID.eq.5) then !'N2')
        solfun = N2sol(S,T)
      elseif (gasID.eq.6) then !'Ar36')
        solfun = Ar36sol(S,T)
      elseif (gasID.eq.7) then ! 'O2')
        solfun = O2sol(S,T)
      else
        stop 'Gas name must be Ne, Ar, Kr, Xe, N2, O2 or Ar36'
      endif

      return
      end

      real*8 FUNCTION Nesol(S,T)

      IMPLICIT NONE

      real*8 S,T

      real*8 temp_S,conc
      real*8 A0,A1,A2,B0,B1
      
C % convert T to scaled temperature
      temp_S = log((298.15d0 - T)/(273.15d0 + T))

C % constants from Table 4 of Hamme and Emerson 2004
      A0 = 2.18156d0
      A1 = 1.29108d0
      A2 = 2.12504d0
      B0 = -5.94737d-3
      B1 = -5.13896d-3

C % Eqn (1) of Hamme and Emerson 2004
      conc = exp(A0 + A1*temp_S + A2*(temp_S**2) + S*(B0 + B1*temp_S))

C % Convert from nmol/kg to umol/kg
      Nesol = conc/1000.d0

      return
      end

      real*8 FUNCTION Arsol(S,T)

      IMPLICIT NONE

      real*8 S,T

      real*8 temp_S
      real*8 A0,A1,A2,A3,B0,B1,B2
      
C % convert T to scaled temperature
      temp_S = log((298.15d0 - T)/(273.15d0 + T))

C % constants from Table 4 of Hamme and Emerson 2004
      A0 = 2.79150d0
      A1 = 3.17609d0
      A2 = 4.13116d0
      A3 = 4.90379d0
      B0 = -6.96233d-3
      B1 = -7.66670d-3
      B2 = -1.16888d-2

C % Eqn (1) of Hamme and Emerson 2004
      Arsol = exp(A0 + A1*temp_S + A2*(temp_S**2) + A3*(temp_S**3) + 
     &        S*(B0 + B1*temp_S + B2*(temp_S**2)))

      return
      end

      real*8 FUNCTION Krsol(S,T)

      IMPLICIT NONE

      real*8 S,T

      real*8 temp_abs,conc
      real*8 A1,A2,A3,A4,B1,B2,B3

C % convert T to scaled temperature
      temp_abs = T + 273.15d0

C % constants from Table 2 Weiss and Kyser for mL/kg
      A1 = -112.6840d0
      A2 = 153.5817d0
      A3 = 74.4690d0
      A4 = -10.0189d0
      B1 = -0.011213d0
      B2 = -0.001844d0
      B3 = 0.0011201d0

C % Eqn (7) of Weiss and Kyser
      conc = exp(A1 + A2*100.d0/temp_abs + A3*log(temp_abs/100.d0) 
     &       + A4*temp_abs/100.d0 + S*(B1 + B2*temp_abs/100.d0 
     &       + B3*(temp_abs/100.d0)**2))

C % Convert concentration from mL/kg to umol/kg
C % Molar volume at STP is calculated from Dymond and Smith (1980) "The virial coefficients of pure gases and mixtures", Clarendon Press, Oxford.
      Krsol = conc/22.3511d-3
      
      return
      end

      real*8 FUNCTION Xesol(S,T)

      IMPLICIT NONE

      real*8 S,T

      real*8 temp_S
      real*8 A0,A1,A2,B0,B1

C % convert T to scaled temperature
      temp_S = log((298.15d0 - T)/(273.15d0 + T))

C % constants from fit procedure of Hamme and Emerson 2004 to Wood and Caputi data
      A0 = -7.48679d0
      A1 = 5.08637d0
      A2 = 4.22243d0
      B0 = -8.15683d-3
      B1 = -1.20129d-2

C % Eqn (1) of Hamme and Emerson 2004
      Xesol = exp(A0 + A1*temp_S + A2*(temp_S**2) + S*(B0 + B1*temp_S))      

      return
      end

      real*8 FUNCTION N2sol(S,T)

      IMPLICIT NONE

      real*8 S,T

      real*8 temp_S
      real*8 A0,A1,A2,A3,B0,B1,B2

C % convert T to scaled temperature
      temp_S = log((298.15d0 - T)/(273.15d0 + T))

C % constants from Table 4 of Hamme and Emerson 2004
      A0 = 6.42931d0
      A1 = 2.92704d0
      A2 = 4.32531d0
      A3 = 4.69149d0
      B0 = -7.44129d-3
      B1 = -8.02566d-3
      B2 = -1.46775d-2

C % Eqn (1) of Hamme and Emerson 2004
      N2sol = exp(A0 + A1*temp_S + A2*(temp_S**2) + A3*(temp_S**3) 
     &        + S*(B0 + B1*temp_S + B2*(temp_S**2)))

      return
      end

      real*8 FUNCTION Ar36sol(S,T)

      IMPLICIT NONE

      real*8 S,T

      real*8 p1,p2
      parameter (p1=-6.95652173913754e-06,p2=1.00122391304348)
      real*8 alpha_fit,X36,XAr,X40
      real*8 gas_mol_fract,Arsol

C % alpha experimentally determined
C       t_alpha = [2 25];
C       alpha = [1.00121 1.00105];
C 
C       alpha_fit = polyval(polyfit(t_alpha,alpha,1),T);

      alpha_fit=p1*T + p2

      X36 = gas_mol_fract(6)
      XAr = gas_mol_fract(2)
      X40 = 0.996349
      Ar36sol = Arsol(S,T)/((alpha_fit*(XAr*X40)/X36)+1.d0)

      return
      end

      real*8 FUNCTION O2sol(S,T)

      IMPLICIT NONE

      real*8 S,T

      real*8 temp_S
      real*8 A0,A1,A2,A3,A4,A5,B0,B1,B2,B3,C0

C % convert T to scaled temperature
      temp_S = log((298.15d0 - T)/(273.15d0 + T))

C % constants from Table 1 of Garcia & Gordon for the fit to Benson and Krause (1984)
      A0 = 5.80871d0 
      A1 = 3.20291
      A2 = 4.17887
      A3 = 5.10006d0
      A4 = -9.86643d-2
      A5 = 3.80369d0
      B0 = -7.01577d-3
      B1 = -7.70028d-3
      B2 = -1.13864d-2
      B3 = -9.51519d-3
      C0 = -2.75915d-7

C % Corrected Eqn (8) of Garcia and Gordon 1992
      O2sol = exp(A0 + A1*temp_S + A2*(temp_S**2) + A3*(temp_S**3) 
     &        + A4*(temp_S**4) + A5*(temp_S**5) + S*(B0 + B1*temp_S 
     &        + B2*(temp_S**2) + B3*(temp_S**3)) + C0*(S**2))

      return
      end

      real*8 FUNCTION moldiff(S,T,gasID)

      IMPLICIT NONE

      real*8 S,T
      integer gasID
      
      real*8 alphaD

      if (gasID.eq.1) then !, 'Ne')
        moldiff=((2.33d0 + (.0550d0*T)) + ((5.21d-4*T)*T))*1.d-9
      elseif (gasID.eq.2) then !'Ar')
        moldiff=((1.58d0 + (0.0045d0*T)) + ((0.00125d0*T)*T))*1.d-9
      elseif (gasID.eq.3) then !'Kr')
        moldiff=((0.87d0 + (.0287d0*T)) + ((4.25d-4*T)*T))*1.d-9
      elseif (gasID.eq.4) then !'Xe')
        moldiff=((0.636d0 + (.0266d0*T)) + ((3.00d-4*T)*T))*1.d-9
      elseif (gasID.eq.5) then !'N2')
        moldiff=((1.227d0 + (0.0579d0*T)) + ((4.71d-4*T)*T))*1.d-9
      elseif (gasID.eq.6) then !'Ar36')
        moldiff=((1.58d0 + (0.0045d0*T)) + ((0.00125d0*T)*T))*1.d-9
C     % adjust diffusivity by fractionation factor alpha
        alphaD = 0.995d0
        moldiff = moldiff/(alphaD**2)
      elseif (gasID.eq.7) then ! 'O2')
        moldiff=((1.45d0 + (0.018d0*T)) + ((0.001d0*T)*T))*1.d-9
      else
        stop 'Gas name must be Ne, Ar, Kr, Xe, N2, O2 or Ar36'
      endif

      return
      end

