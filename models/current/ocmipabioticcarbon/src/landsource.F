      SUBROUTINE LANDSOURCE(C,pCO2a,D,deltaTsg,Fl,S)

C     == input_constants ==
C SPK: UNITS:
C  Inputs: 
C    C in PgC
C    pCO2a in ppm
C    D in PgC/y
C  Outputs:
C    Fl in PgC/y
C    S in PgC/y
C      Same Values and Units as in Eliseev and Mokhov, 2007
C      Cini=kgC or ppmv
C      kd=adimensional
C      kc=kM=ppmv
C      deltaTo= degrees K
C      Al=Ap=Ar=As=1/year
C      Q10p=Q10r=Q10s=adimensional
C      uc=kgC/ppmv
C      ut=kgC/K

C     == Definition of parameters ==
	  REAL*8 KC, KM, KD, DELTAO, AL, AR, AP, AS, Q10P, Q10R, Q10S, UC
	  REAL*8 UT
	  PARAMETER(KC=29.d0, KM=150.d0, KD=0.27d0, DELTATO=10.d0)
	  PARAMETER(AL=1.d0/(11.d0), AR=AL, AP=1.d0/(3.4d0))
	  PARAMETER(AS=1.d0/(30.d0))
	  PARAMETER(Q10P=1.5d0, Q10R=2.15d0, Q10S=2.4d0)
	  PARAMETER(UC=1.31d12, UT=0.33d14)

C     == Definition of passed variables == 
C     Cv, Cs are Carbon terrestrial vegetation, soil initial stocks
C     Cvs is Carbon steady state living biomass, corrected for agriculture harvesting
C     pCO2a is atmostpheric CO2
      REAL*8 C(3), pCO2a, D, deltaTsg, Fl, S(3)
      REAL*8 P, Rp, L, Rs, NPP

C     define gf=fertility factor (adimensional)
      REAL*8 gF
       
      Cv    =  C(1)
      Cs    =  C(2)
      Cvs   =  C(3)
      
      
      IF (pCO2a.LT.kc) THEN
       gf=0.d0
      ELSE
       gf=(pCO2a-kc)/(kM+pCO2a-kc)
      ENDIF

C     == Landsource calculations == 
C     P = C prod rate of photosynthesis
C     Rp = Autotrophic (biota) resp rate
C     L = Littefall
C     Rs = Heterotrophic (soil) resp rate
C     NPP = Terrestrial vegetation net primary production
C     Land Flux
      
      P=Ap*gf*Cvs*(Q10p)**(deltaTsg/deltaTo)
      Rp=Ar*Cv*(Q10r)**(deltaTsg/deltaTo)
      L=Al*Cv
      Rs=As*Cs*(Q10s)**(deltaTsg/deltaTo)
      NPP=P-Rp
      Fl=NPP-Rs
      
C      write(*,*) Fl

C     Sv= S(1), Ss= S(2) are Carbon terrestrial vegetation, soil calculated stocks
C     Svs= S(3) is Carbon steady state living biomass, corrected for agriculture harvesting
      S(1) =  NPP-L-D
      S(2) =  L-Rs
      S(3) = -kd*D

      END