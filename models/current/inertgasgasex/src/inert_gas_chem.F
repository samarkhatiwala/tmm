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
      real*8 A1,A2,A3,A4,B1,B2,B3,C1,MV
      
C % convert T to scaled temperature
      temp_S = (T+273.15d0)/100.d0

C % constants from Table 3 of Jenkins et al. 2019
      A1 = -2.74133d2
      A2 = 3.52620d2
      A3 = 2.26968d2
      A4 = -3.71339d1
      B1 = -6.38600d-2
      B2 = 3.53260d-2
      B3 = -5.32580d-3
      C1 = 1.28233d-5
      MV = 1.00045d0

C % Eqn (7) of Jenkins et al. 2019 (in mol/kg)
      conc = exp(A1 + A2/temp_S + A3*log(temp_S) + A4*temp_S + 
     &        S*(B1 + B2*temp_S + B3*(temp_S**2)) + C1*(S**2))/MV

C % Convert from mol/kg to umol/kg
      Nesol = conc*1000000.d0

      return
      end

      real*8 FUNCTION Arsol(S,T)

      IMPLICIT NONE

      real*8 S,T

      real*8 temp_S,conc
      real*8 A1,A2,A3,A4,B1,B2,B3,C1,MV
      
C % convert T to scaled temperature
      temp_S = (T+273.15d0)/100.d0

C % constants from Table 3 of Jenkins et al. 2019
      A1 = -2.27461d2
      A2 = 3.05435d2
      A3 = 1.80528d2
      A4 = -2.79945d1
      B1 = -6.69420d-2
      B2 = 3.72010d-2
      B3 = -5.63640d-3
      C1 = -5.30325d-6
      MV = 9.99038d-1

C % Eqn (7) of Jenkins et al. 2019 (in mol/kg)
      conc = exp(A1 + A2/temp_S + A3*log(temp_S) + A4*temp_S + 
     &        S*(B1 + B2*temp_S + B3*(temp_S**2)) + C1*(S**2))/MV

C % Convert from mol/kg to umol/kg
      Arsol = conc*1000000.d0

      return
      end

      real*8 FUNCTION Krsol(S,T)

      IMPLICIT NONE

      real*8 S,T

      real*8 temp_S,conc
      real*8 A1,A2,A3,A4,B1,B2,B3,C1,MV
      
C % convert T to scaled temperature
      temp_S = (T+273.15d0)/100.d0

C % constants from Table 3 of Jenkins et al. 2019
      A1 = -1.22469d2
      A2 = 1.53565d2
      A3 = 7.01969d1
      A4 = -8.52524d0
      B1 = -4.95220d-2
      B2 = 2.44340d-2
      B3 = -3.39680d-3
      C1 = 4.19208d-6
      MV = 9.97232d-1

C % Eqn (7) of Jenkins et al. 2019 (in mol/kg)
      conc = exp(A1 + A2/temp_S + A3*log(temp_S) + A4*temp_S + 
     &        S*(B1 + B2*temp_S + B3*(temp_S**2)) + C1*(S**2))/MV

C % Convert from mol/kg to umol/kg
      Krsol = conc*1000000.d0

      return
      end

      real*8 FUNCTION Xesol(S,T)

      IMPLICIT NONE

      real*8 S,T

      real*8 temp_S,conc
      real*8 A1,A2,A3,A4,B1,B2,B3,C1,MV
      
C % convert T to scaled temperature
      temp_S = (T+273.15d0)/100.d0

C % constants from Table 3 of Jenkins et al. 2019
      A1 = -2.24510d2
      A2 = 2.92823d2
      A3 = 1.57613d2
      A4 = -2.26690d1
      B1 = -8.49150d-2
      B2 = 4.79960d-2
      B3 = -7.35950d-3
      C1 = 6.69292d-6
      MV = 9.93098d-1

C % Eqn (7) of Jenkins et al. 2019 (in mol/kg)
      conc = exp(A1 + A2/temp_S + A3*log(temp_S) + A4*temp_S + 
     &        S*(B1 + B2*temp_S + B3*(temp_S**2)) + C1*(S**2))/MV

C % Convert from mol/kg to umol/kg
      Xesol = conc*1000000.d0

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
      parameter (p1=-6.3492708401d-6,p2=1.0011708d0)
      real*8 alpha_fit,X36,XAr,X40
      real*8 gas_mol_fract,Arsol

C % alpha experimentally determined (Seltzer et al., 2019)

      alpha_fit=p1*T + p2

      X36 = gas_mol_fract(6)
      XAr = gas_mol_fract(2)
      X40 = 0.996349d0
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
      A1 = 3.20291d0
      A2 = 4.17887d0
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

      real*8 R
      parameter (R = 8.314510d0)
      
      real*8 alphaD, A, Ea

      if (gasID.eq.1) then !, 'Ne')
        A = 1.6080d-6
        Ea = 14840.d0
      elseif (gasID.eq.2) then !'Ar')
        A = 2.227d-6
        Ea = 16680.d0
      elseif (gasID.eq.3) then !'Kr')
        A = 6.3930d-6
        Ea = 20200.d0
      elseif (gasID.eq.4) then !'Xe')
        A = 9.0070d-6
        Ea = 21610.d0
      elseif (gasID.eq.5) then !'N2')
        A = 3.4120d-6
        Ea = 18500.d0
      elseif (gasID.eq.6) then !'Ar36')
C     % adjust diffusivity by fractionation factor alpha
        alphaD = 0.9963d0
        A = 2.227d-6*(1.d0/(alphaD**2))
        Ea = 16680.d0
C     % adjust diffusivity by fractionation factor alpha
        alphaD = 0.9963d0
        moldiff = moldiff/(alphaD**2)
      elseif (gasID.eq.7) then ! 'O2')
        A = 4.286d-6
        Ea = 18700.d0
      elseif (gasID.eq.8) then ! 'He')
        A = 0.8180d-6
        Ea = 11700.d0        
      else
        stop 'Gas name must be Ne, Ar, Kr, Xe, N2, O2, Ar36 or He'
      endif

      moldiff = (1.d0 - 0.049d0*S/35.5d0)*A*exp(-Ea/(R*(T + 273.16d0)))

      return
      end
