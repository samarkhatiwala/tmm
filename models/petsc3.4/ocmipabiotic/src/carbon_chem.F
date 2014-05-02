CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CC  New efficient pCO2 solver, Mick Follows         CC
CC                             Taka Ito             CC
CC                             Stephanie Dutkiewicz CC
CC  20 April 2003                                   CC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CStartOfInterFace
      SUBROUTINE CALC_CO2_APPROX(
     I                       t,s,dic_in,pt_in,sit_in,ta_in,
     I                       k1local,k2local,
     I                       k1plocal,k2plocal,k3plocal,
     I                       kslocal,kblocal,kwlocal,
     I                       ksilocal,kflocal,
     I                       fflocal,btlocal,stlocal,ftlocal,
     U                       pHlocal,co2s)
C     /==========================================================\
C     | SUBROUTINE CALC_PCO2_APPROX                              |
C     \==========================================================/
      IMPLICIT NONE

C     == GLobal variables ==
C     == Routine arguments ==
C       diclocal = total inorganic carbon (mol/m^3)
C             where 1 T = 1 metric ton = 1000 kg
C       ta  = total alkalinity (eq/m^3)
C       pt  = inorganic phosphate (mol/^3)
C       sit = inorganic silicate (mol/^3)
C       t   = temperature (degrees C)
C       s   = salinity (PSU)
      real*8 t,s,pt_in,sit_in,ta_in,dic_in,pHlocal,co2s
      real*8 fflocal,btlocal,stlocal,ftlocal,k1local,k2local,k1plocal
      real*8 k2plocal,k3plocal,kslocal,kblocal,kwlocal,ksilocal,kflocal
CEndOfInterface

C     == Local variables ==
      real*8 permil, permeg
      parameter (permil=1.d0/1024.5d0, permeg=1.d-6)
      real*8 pt,ta,sit,diclocal
      real*8 phguess,cag,bohg,hguess,stuff,gamm,hnew,h3po4g,h2po4g
      real*8 hpo4g,po4g,siooh3g
      real*8 co3local

c ---------------------------------------------------------------------
C Change units from the input of mol/m^3 -> mol/kg:
c (1 mol/m^3)  x (1 m^3/1024.5 kg)
c where the ocean's mean surface density is 1024.5 kg/m^3
c Note: mol/kg are actually what the body of this routine uses
c for calculations.  Units are reconverted back to mol/m^3 at the
c end of this routine.
c To convert input in mol/m^3 -> mol/kg
      pt=pt_in*permil
      sit=sit_in*permil
      ta=ta_in*permil
      diclocal=dic_in*permil
c ---------------------------------------------------------------------
c set first guess and brackets for [H+] solvers
c first guess (for newton-raphson)
      phguess = phlocal
cmick - new approx method
cmick - make estimate of htotal (hydrogen ion conc) using
cmick   appromate estimate of CA, carbonate alkalinity
      hguess = 10.d0**(-phguess)
cmick - first estimate borate contribution using guess for [H+]
      bohg = btlocal*kblocal/(hguess+kblocal)

cmick - first estimate of contribution from phosphate
cmick based on Dickson and Goyet
      stuff = hguess*hguess*hguess
     &           + (k1plocal*hguess*hguess)
     &           + (k1plocal*k2plocal*hguess)
     &           + (k1plocal*k2plocal*k3plocal)
      h3po4g = (pt*hguess*hguess*hguess) / stuff
      h2po4g = (pt*k1plocal*hguess*hguess) / stuff
      hpo4g  = (pt*k1plocal*k2plocal*hguess) / stuff
      po4g   = (pt*k1plocal*k2plocal*k3plocal) / stuff

cmick - estimate contribution from silicate
cmick based on Dickson and Goyet
      siooh3g = sit*ksilocal / (ksilocal + hguess)

cmick - now estimate carbonate alkalinity
      cag = ta - bohg - (kwlocal/hguess) + hguess
     &           - hpo4g - 2.d0*po4g + h3po4g
     &           - siooh3g

cmick - now evaluate better guess of hydrogen ion conc
cmick   htotal = [H+], hydrogen ion conc
      gamm  = diclocal/cag
      stuff = (1.d0-gamm)*(1.d0-gamm)*k1local*k1local
     &          - 4.d0*k1local*k2local*(1.d0-2.d0*gamm)
      hnew  = 0.5d0*( (gamm-1.d0)*k1local + sqrt(stuff) )
cmick - now determine [CO2*]
      co2s  = diclocal/
     &   (1.d0 + (k1local/hnew) + (k1local*k2local/(hnew*hnew)))
cmick - return update pH to main routine
      phlocal = -log10(hnew)

c NOW EVALUATE CO32-, carbonate ion concentration
c used in determination of calcite compensation depth
c Karsten Friis & Mick - Sep 2004
c       co3local = k1local*k2local*diclocal / 
c    &         (hnew*hnew + k1local*hnew + k1local*k2local)

c ---------------------------------------------------------------
C SPK: NOTE on units
C fflocal is the solubility (computed in S/R CARBON_COEFFS) in mol/(kg*atm)
C To convert to mol/(m^3*uatm), multiply ff by 1e-6*1024.5, i.e.
C sol = ff*permeg*rho0. 
C Equilibrium [CO2]aq in mol/m^3 = sol*pCO2_atm*atmpres, where
C            pCO2_atm = atmospheric mole fraction CO2 in dry air (ppmv)
C            atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)
c surface pCO2 (following Dickson and Goyet, DOE...)
C      pCO2surfloc = co2s/fflocal/permeg ! pCO2 in uatm
C      write(*,*)co2s/fflocal/permeg,fflocal   
      co2s = co2s/permil ! convert from mol/kg -> mol/m^3
      return
      end

c=================================================================
c *******************************************************************
c=================================================================
CStartOfInterFace
      SUBROUTINE CARBON_COEFFS(temp,salt,ff,ak0,ak1,ak2,akb,ak1p,ak2p,
     &                         ak3p,aksi,akw,aks,akf,bt,st,ft,N)
C     /==========================================================\
C     | SUBROUTINE CARBON_COEFFS                                 |
C     | determine coefficients for surface carbon chemistry      |
C     | adapted from OCMIP2:  SUBROUTINE CO2CALC                 |
C     | mick follows, oct 1999                                   |
c     | minor changes to tidy, swd aug 2002                      |
C     \==========================================================/
C INPUT
C       diclocal = total inorganic carbon (mol/m^3) 
C             where 1 T = 1 metric ton = 1000 kg
C       ta  = total alkalinity (eq/m^3) 
C       pt  = inorganic phosphate (mol/^3) 
C       sit = inorganic silicate (mol/^3) 
C       t   = temperature (degrees C)
C       s   = salinity (PSU)
C OUTPUT
C IMPORTANT: Some words about units - (JCO, 4/4/1999)
c     - Models carry tracers in mol/m^3 (on a per volume basis)
c     - Conversely, this routine, which was written by observationalists 
c       (C. Sabine and R. Key), passes input arguments in umol/kg  
c       (i.e., on a per mass basis)
c     - I have changed things slightly so that input arguments are in mol/m^3,
c     - Thus, all input concentrations (diclocal, ta, pt, and st) should be 
c       given in mol/m^3; output arguments "co2star" and "dco2star"  
c       are likewise be in mol/m^3.
C--------------------------------------------------------------------------
        IMPLICIT NONE
C     == GLobal variables ==
C     == Routine arguments ==
      integer N
      real*8 temp(N),salt(N)
      real*8 ff(N),ak0(N),ak1(N),ak2(N),akb(N),ak1p(N),ak2p(N),ak3p(N)
      real*8 aksi(N),akw(N),aks(N),akf(N),bt(N),st(N),ft(N)
CEndOfInterface

C LOCAL VARIABLES
      real*8 t,s,ta,pt,sit,tk,tk100,tk1002,dlogtk,sqrtis,sqrts,s15,scl
      real*8 x1,x2,s2,xacc,invtk,is,is2
      integer i
CEOP

      do i=1,N
        t = temp(i)
        s = salt(i)
C terms used more than once
        tk = 273.15d0 + t
        tk100 = tk/100.0d0
        tk1002=tk100*tk100
        invtk=1.0d0/tk
        dlogtk=log(tk)
        is=19.924d0*s/(1000.d0-1.005d0*s)
        is2=is*is
        sqrtis=sqrt(is)
        s2=s*s
        sqrts=sqrt(s)
        s15=s**1.5d0
        scl=s/1.80655d0

C------------------------------------------------------------------------
C Calculate concentrations for borate, sulfate, and fluoride
C Uppstrom (1974)
        bt(i) = 0.000232d0 * scl/10.811d0
C Morris & Riley (1966)
        st(i) = 0.14d0 * scl/96.062d0
C Riley (1965)
        ft(i) = 0.000067d0 * scl/18.9984d0
C------------------------------------------------------------------------
        
C------------------------------------------------------------------------
C f = k0(1-pH2O)*correction term for non-ideality
C Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
        ff(i) = exp(-162.8301d0 + 218.2968d0/tk100  +
     &          90.9241d0*log(tk100) - 1.47696d0*tk1002 +
     &          s * (.025695d0 - .025225d0*tk100 + 
     &          0.0049867d0*tk1002))
C------------------------------------------------------------------------
C K0 from Weiss 1974
        ak0(i) = exp(93.4517d0/tk100 - 60.2409d0 + 
     &        23.3585d0 * log(tk100) +
     &        s * (0.023517d0 - 0.023656d0*tk100 + 
     &        0.0047036d0*tk1002))
C------------------------------------------------------------------------
C k1 = [H][HCO3]/[H2CO3]
C k2 = [H][CO3]/[HCO3]     on hSWS
C Millero p.664 (1995) using Mehrbach et al. data on SEAWATER scale 
C (Original reference: Dickson and Millero, DSR, 1987)
        ak1(i)=10**(-1.d0*(3670.7d0*invtk - 
     &          62.008d0 + 9.7944d0*dlogtk -
     &          0.0118d0*s + 0.000116d0*s2))
        ak2(i)=10**(-1.d0*(1394.7d0*invtk + 4.777d0 - 
     &          0.0184d0*s + 0.000118d0*s2))
C------------------------------------------------------------------------
C k1p = [H][H2PO4]/[H3PO4] on hSWS
C Millero p.670 (1995)
	    ak1p(i) = exp(-4576.752d0*invtk + 115.540d0 - 
     &          18.453d0*dlogtk + 
     &		    (-106.736d0*invtk + 0.69171d0)*sqrts +
     &		    (-0.65643d0*invtk - 0.01844d0)*s)
C------------------------------------------------------------------------
C k2p = [H][HPO4]/[H2PO4] on hSWS
C Millero p.670 (1995)
        ak2p(i) = exp(-8814.715d0*invtk + 172.1033d0 - 
     &          27.927d0*dlogtk +
     &		    (-160.340d0*invtk + 1.3566d0)*sqrts +
     &		    (0.37335d0*invtk - 0.05778d0)*s)
C------------------------------------------------------------------------
C k3p = [H][PO4]/[HPO4] on hSWS
C Millero p.670 (1995)
	    ak3p(i) = exp(-3070.75d0*invtk - 18.126d0 + 
     &		    (17.27039d0*invtk + 2.81197d0) *
     &		    sqrts + (-44.99486d0*invtk - 0.09984d0) * s)
C------------------------------------------------------------------------
C ksi = [H][SiO(OH)3]/[Si(OH)4] on hSWS
C Millero p.671 (1995) using data from Yao and Millero (1995)
C change to (mol/ kg soln)
        aksi(i) = exp(-8904.2d0*invtk + 117.400d0 - 
     &          19.334d0*dlogtk +
     &		    (-458.79d0*invtk + 3.5913d0) * sqrtis +
     &		    (188.74d0*invtk - 1.5998d0) * is +
     &		    (-12.1652d0*invtk + 0.07871d0) * is2 +
     &		    log(1.d0-0.001005d0*s))
C------------------------------------------------------------------------
C kw = [H][OH] on hSWS
C Millero p.670 (1995) using composite data
        akw(i) = exp(-13847.26d0*invtk + 148.9802d0 - 
     &          23.6521d0*dlogtk +
     &		    (118.67d0*invtk - 5.977d0 + 1.0495d0 * dlogtk) *
     &		    sqrts - 0.01615d0 * s)
C------------------------------------------------------------------------
C ks = [H][SO4]/[HSO4] on free H scale
C Dickson (1990, J. chem. Thermodynamics 22, 113)
C change to (mol/ kg soln)
        aks(i)=exp(-4276.1d0*invtk + 141.328d0 - 
     &          23.093d0*dlogtk +
     &		    (-13856.d0*invtk + 324.57d0 - 47.986d0*dlogtk)*sqrtis +
     &   		(35474.d0*invtk - 771.54d0 + 114.723d0*dlogtk)*is -
     &  		2698.d0*invtk*is**1.5d0 + 1776.d0*invtk*is2 +
     &		    log(1.d0 - 0.001005d0*s))
C------------------------------------------------------------------------
C kf = [H][F]/[HF] on free H scale
C Dickson and Riley (1979)
C change to (mol/ kg soln)
        akf(i)=exp(1590.2d0*invtk - 12.641d0 + 1.525d0*sqrtis +
     &		    log(1.d0 - 0.001005d0*s)) 
C------------------------------------------------------------------------
C kb = [H][BO2]/[HBO2] on hSWS
C Dickson p.673 (1990)
C change from htotal to hSWS
        akb(i)=exp( (-8966.90d0 - 2890.53d0*sqrts - 77.942d0*s +
     &		    1.728d0*s15 - 0.0996d0*s2)*invtk +
     &  		(148.0248d0 + 137.1942d0*sqrts + 1.62142d0*s) +
     &  		(-24.4344d0 - 25.085d0*sqrts - 0.2474d0*s) *
     &  		dlogtk + 0.053105d0*sqrts*tk +
     &          log((1.d0+(st(i)/aks(i))+(ft(i)/akf(i))) 
     &          /(1.d0+(st(i)/aks(i)))) )
      enddo

      return
      end

c=================================================================
CStartOfInterFace
      SUBROUTINE CALC_CO2SOLUBILITY(temp,salt,sol,N)
C     /==========================================================\
C     | SUBROUTINE CO2SOLUBILITY                                 |
C     \==========================================================/
C INPUT
C OUTPUT
C--------------------------------------------------------------------------
      IMPLICIT NONE
C     == GLobal variables ==
      integer N
      real*8 temp(N),salt(N)
      real*8 sol(N)
CEndOfInterface
      
C LOCAL VARIABLES
      real*8 rho0, permil, permeg
      parameter (rho0=1024.5d0, permil=1.d0/rho0, permeg=1.d-6)      
      real*8 t,s,tk,tk100,tk1002,ff
      integer i
CEOP

      do i=1,N
        t = temp(i)
        s = salt(i)
        tk = 273.15d0 + t
        tk100 = tk/100.0d0
        tk1002=tk100*tk100

C------------------------------------------------------------------------
C f = k0(1-pH2O)*correction term for non-ideality
C Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
        ff = exp(-162.8301d0 + 218.2968d0/tk100  +
     &        90.9241d0*log(tk100) - 1.47696d0*tk1002 +
     &        s * (.025695d0 - .025225d0*tk100 + 
     &        0.0049867d0*tk1002))
C------------------------------------------------------------------------

C ff is the solubility in mol/(kg*atm)
C To convert to mol/(m^3*uatm), multiply ff by 1e-6*1024.5, i.e.
C sol = ff*permeg*rho0. 
C Equilibrium [CO2]aq in mol/m^3 = sol*pCO2_atm*atmpres, where
C            pCO2_atm = atmospheric mole fraction CO2 in dry air at 1 atm total pres (ppmv)
C            atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)
        sol(i) = ff*permeg*rho0 ! solubility in mol/(m^3*uatm)
      ENDDO
      
      END
      
c=================================================================
CStartOfInterFace
      SUBROUTINE CO2_NEWTON(dic_in,ta_in,pt_in,sit_in,phlo,phhi,
     &                      pH,co2star)
C     /==========================================================\
C     | SUBROUTINE CO2_NEWTON                                 |
C     \==========================================================/
C INPUT
C OUTPUT
C--------------------------------------------------------------------------
      IMPLICIT NONE
C     == GLobal variables ==
C     == Routine arguments ==
      real*8 dic_in,ta_in,pt_in,sit_in,phhi,phlo
      real*8 co2star,pH
CEndOfInterface

C LOCAL VARIABLES
      real*8 rho0, permil, permeg
      parameter (rho0=1024.5d0, permil=1.d0/rho0, permeg=1.d-6)
      real*8 x1,x2,xacc,hSWS,hSWS2
      real*8 k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi
      real*8 bt,st,ft
      real*8 sit,pt,dic,ta
	  common/carbconst/k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi
      common/species/bt,st,ft
      common/chemvars/sit,pt,dic,ta      
      real*8 drtsafe
      external ta_iter_SWS
CEOP

c       To convert input in mol/m^3 -> mol/kg 
      pt=pt_in*permil
      sit=sit_in*permil
      ta=ta_in*permil
      dic=dic_in*permil

      x1 = 10.d0**(-phhi)
      x2 = 10.d0**(-phlo)
      xacc = 1.d-12
      hSWS = drtsafe(ta_iter_SWS,x1,x2,xacc)
C
C Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2, 
C ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
C
	  hSWS2=hSWS*hSWS
	  co2star=dic*hSWS2/(hSWS2 + k1*hSWS + k1*k2)
C      co2starair=xco2*ff*atmpres
C      dco2star=co2starair-co2star
	  ph=-log10(hSWS)

c     Note: co2star is calculated in mol/kg within this routine 
c     Thus Convert now from mol/kg -> mol/m^3
      co2star = co2star/permil

      END

c=================================================================
CStartOfInterFace
      SUBROUTINE CO2_NEWTON_VECTORIZED(dic_in,ta_in,pt_in,sit_in,
     &                    ak0,ak1,ak2,akb,ak1p,ak2p,ak3p,aksi,
     &                    akw,aks,akf,abt,ast,aft,phlo,phhi,
     &                    pH,co2star,dco2ddic,N)
C     /==========================================================\
C     | SUBROUTINE CO2_NEWTON_VECTORIZED                         |
C     \==========================================================/
C INPUT
C OUTPUT
C--------------------------------------------------------------------------
      IMPLICIT NONE
C     == GLobal variables ==
C     == Routine arguments ==
      integer N
      real*8 dic_in(N),ta_in(N),pt_in(N),sit_in(N)
      real*8 ak0(N),ak1(N),ak2(N),akb(N),ak1p(N),ak2p(N),ak3p(N),
     &       aksi(N),akw(N),aks(N),akf(N),abt(N),ast(N),aft(N)      
      real*8 phhi,phlo     
      real*8 pH(N),co2star(N),dco2ddic(N)
CEndOfInterface

C LOCAL VARIABLES
      real*8 rho0, permil, permeg
      parameter (rho0=1024.5d0, permil=1.d0/rho0, permeg=1.d-6)
      integer i
      real*8 x1,x2,xacc,hSWS,hSWS2
      real*8 k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi
      real*8 bt,st,ft
      real*8 sit,pt,dic,ta
	  common/carbconst/k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi
      common/species/bt,st,ft
      common/chemvars/sit,pt,dic,ta      
      real*8 drtsafe
      external ta_iter_SWS
CEOP

      do i=1,N
C       Set common block variables      
c       To convert input in mol/m^3 -> mol/kg 
		pt=pt_in(i)*permil
		sit=sit_in(i)*permil
		ta=ta_in(i)*permil
		dic=dic_in(i)*permil
        k1=ak1(i)
        k2=ak2(i)
        kw=akw(i)
        kb=akb(i)
        ks=aks(i)
        kf=akf(i)
        k1p=ak1p(i)
        k2p=ak2p(i)
        k3p=ak3p(i)
        ksi=aksi(i)
        bt=abt(i)
        st=ast(i)
        ft=aft(i)
        
		x1 = 10.d0**(-phhi)
		x2 = 10.d0**(-phlo)
		xacc = 1.d-12
		hSWS = drtsafe(ta_iter_SWS,x1,x2,xacc)
C
C Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2, 
C ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
C
		hSWS2=hSWS*hSWS
		co2star(i)=dic*hSWS2/(hSWS2 + k1*hSWS + k1*k2)
C      co2starair=xco2*ff*atmpres
C      dco2star=co2starair-co2star
		ph(i)=-log10(hSWS)

c     Note: co2star is calculated in mol/kg within this routine 
c     Thus Convert now from mol/kg -> mol/m^3
		co2star(i) = co2star(i)/permil
        dco2ddic(i) = (hSWS/k2)/(dic_in(i)/co2star(i)) ! approximate d(co2star)/d(dic)
      enddo
      
      END

      
c=================================================================
CStartOfInterFace
      SUBROUTINE TA_ITER_SWS(x,fn,df)
C     /==========================================================\
C     | SUBROUTINE TA_ITER_SWS                                 |
C     \==========================================================/
C INPUT
C OUTPUT
C--------------------------------------------------------------------------
      IMPLICIT NONE
C     == GLobal variables ==
      real*8 k0,k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi
      real*8 bt,st,ft
      real*8 sit,pt,dic,ta
	  common/carbconst/k1,k2,kw,kb,ks,kf,k1p,k2p,k3p,ksi
      common/species/bt,st,ft
      common/chemvars/sit,pt,dic,ta      
C     == Routine arguments ==
      real*8 x,fn,df
CEndOfInterface

C LOCAL VARIABLES
      real*8 k12,k12p,k123p,x2,x3,a,c,da,b,a2,b2,db
CEOP

C
C Modified from ta_iter_1.f (RCS version 1.2, OCMIP-2)
C - by A. Mouchet, 2004:
C Fixed Problems w/ version of ta_iter_1.f used in OCMIP-2 (vers. 1.2)
C  1) fixed errors in signs, parenthesis and coefficient c in derivative
C  2) changed from Total to Seawater Scale 
C     * c defined for seawater H scale; 
C     * fn and df adapted to KF on free H scale
C     * comments have been adapted
C

C
C This routine expresses TA as a function of DIC, hSWS and constants.
C It also calculates the derivative of this function with respect to 
C hSWS. It is used in the iterative solution for hSWS. In the call
C "x" is the input value for hSWS, "fn" is the calculated value for TA
C and "df" is the value for dTA/dhSWS
C
	  x2=x*x
	  x3=x2*x
	  k12 = k1*k2
	  k12p = k1p*k2p
	  k123p = k12p*k3p
	  c = 1.d0 + st/ks + ft/kf
	  a = x3 + k1p*x2 + k12p*x + k123p
	  a2=a*a
	  da = 3.d0*x2 + 2.d0*k1p*x + k12p
	  b = x2 + k1*x + k12
	  b2=b*b
	  db = 2.d0*x + k1
C
C	fn = hco3+co3+borate+oh+hpo4+2*po4+silicate-hfree-hso4-hf-h3po4-ta
C===========================================================================
C
      fn = k1*x*dic/b +
     &	     2.d0*dic*k12/b +
     &	     bt/(1.d0 + x/kb) +
     &	     kw/x +
     &	     pt*k12p*x/a +
     &	     2.d0*pt*k123p/a +
     &	     sit/(1.d0 + x/ksi) -
     &	     x/c -
     &	     st/(1.d0 + ks/(x/c)) -
     &	     ft/(1.d0 + kf/(x/c)) -
     &	     pt*x3/a -
     &	     ta
C
C	df = dfn/dx
C
      df = ((k1*dic*b) - k1*x*dic*db)/b2 -
     &	     2.d0*dic*k12*db/b2 -
     &	     bt/kb/(1.d0+x/kb)**2.d0 -
     &	     kw/x2 +
     &	     (pt*k12p*(a - x*da))/a2 -
     &	     2.d0*pt*k123p*da/a2 -
     &	     sit/ksi/(1.d0+x/ksi)**2.d0 -
     &	     1.d0/c -
     &       st *(1.d0 + ks/(x/c))**(-2.d0)*(ks*c/x2) -
     &	     ft*(1.d0 + kf/(x/c))**(-2.d0)*(kf*c/x2) -
     &	     pt*x2*(3.d0*a-x*da)/a2
C
	  return
	  end

      real*8 FUNCTION DRTSAFE(FUNCD,X1,X2,XACC)
C
C	File taken from Numerical Recipes. Modified  R.M.Key 4/94
C
      IMPLICIT NONE

      INTEGER MAXIT
      REAL*8 x1,x2,xacc
      EXTERNAL funcd
      PARAMETER (MAXIT=100)
      INTEGER j
      REAL*8 df,dx,dxold,f,fh,fl,temp,xh,xl,swap
      
      CALL FUNCD(X1,FL,DF)
      CALL FUNCD(X2,FH,DF)
      IF(FL .LT. 0.d0) THEN
        XL=X1
        XH=X2
      ELSE
        XH=X1
        XL=X2
        SWAP=FL
        FL=FH
        FH=SWAP
      END IF
      DRTSAFE=0.5d0*(X1+X2)
      DXOLD=ABS(X2-X1)
      DX=DXOLD
      CALL FUNCD(DRTSAFE,F,DF)
      DO 100, J=1,MAXIT
        IF(((DRTSAFE-XH)*DF-F)*((DRTSAFE-XL)*DF-F) .GE. 0.d0 .OR.
     &	      ABS(2.d0*F) .GT. ABS(DXOLD*DF)) THEN
          DXOLD=DX
          DX=0.5d0*(XH-XL)
          DRTSAFE=XL+DX
          IF(XL .EQ. DRTSAFE)RETURN
        ELSE
          DXOLD=DX
          DX=F/DF
          TEMP=DRTSAFE
          DRTSAFE=DRTSAFE-DX
          IF(TEMP .EQ. DRTSAFE)RETURN
	END IF
        IF(ABS(DX) .LT. XACC)RETURN
        CALL FUNCD(DRTSAFE,F,DF)
        IF(F .LT. 0.d0) THEN
          XL=DRTSAFE
          FL=F
        ELSE
          XH=DRTSAFE
          FH=F
        END IF
  100  CONTINUE
      RETURN
      END
      