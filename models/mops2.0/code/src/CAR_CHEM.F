C$Header: /Users/ikriest/CVS/mops/CAR_CHEM.F,v 1.2 2016/06/03 09:28:59 ikriest Exp $
C$Name: mops-2_0 $

C Basis for all routines: MIT source code for OCMIP abiotic model,
C modified by Samar Khatiwala, then by Iris Kriest (Aug-Nov 2010)

C NOTE: MIT model units are [mol/m3].
C       Therefore convert DIC, PO4, Silicate and TA input by 
C       factor 1d-3 (1/convert_mol_to_mmol) before using them as input.
C       Therefore vonvert change of CO2 due to air-sea flux by 
C       factor 1d+3 (convert_mol_to_mmol) before using it as output.
C       Leave everything else as before.
C       Take care when later calculating carbonate dissolution.
C       See CO2_SURFFORCING, CAR_INI and CAR_PARAMS.

C THIS FILE CONTAINS:	
C SR car_coffs         - calculate carbonate coefficients
C SR co2_surfac        - calculate surface CO2
C SR co2_surfforcing   - calculate O2 air-sea gas exchange
 
C CALCULATE THE COEFFICIENTS OF THE CARBONATE SYSTEM OVER N LAYERS
C CALCULATE THE COEFFICIENTS OF THE CARBONATE SYSTEM OVER N LAYERS
C CALCULATE THE COEFFICIENTS OF THE CARBONATE SYSTEM OVER N LAYERS
C CALCULATE THE COEFFICIENTS OF THE CARBONATE SYSTEM OVER N LAYERS
C CALCULATE THE COEFFICIENTS OF THE CARBONATE SYSTEM OVER N LAYERS

C PREPROCESSOR OPTIONS:
 
C CALLED BY:	mops_biogeochem_ini
C		bgc_model

C CALLS:	

C INPUT/ARGUMENT LIST:	
C ttemp		temperature, n layers
C stemp		salinity, n layers
C n             number of layers
C 
C OUTPUT/COMMON BLOCK: 
C ak0(n)	k0, Weiss 1974
C ak1(n)	k1 = [H][HCO3]/[H2CO3]
C ak2(n)	k2 = [H][CO3]/[HCO3]
C ak1p(n)	k1p = [H][H2PO4]/[H3PO4]
C ak2p(n)	k2p = [H][HPO4]/[H2PO4]
C ak3p(n)	k3p = [H][PO4]/[HPO4]
C aks(n)	ks = [H][SO4]/[HSO4]
C akb(n)	kb = [H][BO2]/[HBO2]
C akw(n)	kw = [H][OH]
C aksi(n)	ksi = [H][SiO(OH)3]/[Si(OH)4]
C akf(n)	kf = [H][F]/[HF]
C ff(n)		f = k0(1-pH2O)
C bt(n)		borate conc.
C st(n)		sulfate conc.
C ft(n)		fluoride conc.

      subroutine car_coeffs(ttemp,stemp,n)

      implicit none
      
      integer n, i
      real*8 ttemp(n),stemp(n)

#include "CAR_PARAMS.h"

! local coefficients

      real*8 t,tk,tk100,tk1002,invtk,dlogtk
      real*8 is,is2,sqrtis,s,s2,sqrts,s15,scl

      do i=1,n

        t = ttemp(i)
        s = stemp(i)

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

C Calculate concentrations for borate, sulfate, and fluoride
C Uppstrom (1974), Morris & Riley (1966), Riley (1965)
        bt(i) = 0.000232d0 * scl/10.811d0
        st(i) = 0.14d0 * scl/96.062d0
        ft(i) = 0.000067d0 * scl/18.9984d0
        

C f = k0(1-pH2O)*correction term for non-ideality
C Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
        ff(i) = exp(-162.8301d0 + 218.2968d0/tk100  +
     &          90.9241d0*log(tk100) - 1.47696d0*tk1002 +
     &          s * (.025695d0 - .025225d0*tk100 + 
     &          0.0049867d0*tk1002))


C K0 from Weiss 1974
        ak0(i) = exp(93.4517d0/tk100 - 60.2409d0 + 
     &        23.3585d0 * log(tk100) +
     &        s * (0.023517d0 - 0.023656d0*tk100 + 
     &        0.0047036d0*tk1002))


C k1 = [H][HCO3]/[H2CO3]
C k2 = [H][CO3]/[HCO3]     on hSWS
C Millero p.664 (1995) using Mehrbach et al. data on SEAWATER scale 
C (Original reference: Dickson and Millero, DSR, 1987)
        ak1(i)=10**(-1.d0*(3670.7d0*invtk - 
     &          62.008d0 + 9.7944d0*dlogtk -
     &          0.0118d0*s + 0.000116d0*s2))
        ak2(i)=10**(-1.d0*(1394.7d0*invtk + 4.777d0 - 
     &          0.0184d0*s + 0.000118d0*s2))


C k1p = [H][H2PO4]/[H3PO4] on hSWS
C Millero p.670 (1995)
	    ak1p(i) = exp(-4576.752d0*invtk + 115.540d0 - 
     &          18.453d0*dlogtk + 
     &		    (-106.736d0*invtk + 0.69171d0)*sqrts +
     &		    (-0.65643d0*invtk - 0.01844d0)*s)


C k2p = [H][HPO4]/[H2PO4] on hSWS
C Millero p.670 (1995)
        ak2p(i) = exp(-8814.715d0*invtk + 172.1033d0 - 
     &          27.927d0*dlogtk +
     &		    (-160.340d0*invtk + 1.3566d0)*sqrts +
     &		    (0.37335d0*invtk - 0.05778d0)*s)


C k3p = [H][PO4]/[HPO4] on hSWS
C Millero p.670 (1995)
	    ak3p(i) = exp(-3070.75d0*invtk - 18.126d0 + 
     &		    (17.27039d0*invtk + 2.81197d0) *
     &		    sqrts + (-44.99486d0*invtk - 0.09984d0) * s)


C ksi = [H][SiO(OH)3]/[Si(OH)4] on hSWS
C Millero p.671 (1995) using data from Yao and Millero (1995)
C change to (mol/ kg soln)
        aksi(i) = exp(-8904.2d0*invtk + 117.400d0 - 
     &          19.334d0*dlogtk +
     &		    (-458.79d0*invtk + 3.5913d0) * sqrtis +
     &		    (188.74d0*invtk - 1.5998d0) * is +
     &		    (-12.1652d0*invtk + 0.07871d0) * is2 +
     &		    log(1.d0-0.001005d0*s))


C kw = [H][OH] on hSWS
C Millero p.670 (1995) using composite data
        akw(i) = exp(-13847.26d0*invtk + 148.9802d0 - 
     &          23.6521d0*dlogtk +
     &		    (118.67d0*invtk - 5.977d0 + 1.0495d0 * dlogtk) *
     &		    sqrts - 0.01615d0 * s)


C ks = [H][SO4]/[HSO4] on free H scale
C Dickson (1990, J. chem. Thermodynamics 22, 113)
C change to (mol/ kg soln)
        aks(i)=exp(-4276.1d0*invtk + 141.328d0 - 
     &          23.093d0*dlogtk +
     &		    (-13856.d0*invtk + 324.57d0 - 47.986d0*dlogtk)*sqrtis +
     &   		(35474.d0*invtk - 771.54d0 + 114.723d0*dlogtk)*is -
     &  		2698.d0*invtk*is**1.5d0 + 1776.d0*invtk*is2 +
     &		    log(1.d0 - 0.001005d0*s))


C kf = [H][F]/[HF] on free H scale
C Dickson and Riley (1979)
C change to (mol/ kg soln)
        akf(i)=exp(1590.2d0*invtk - 12.641d0 + 1.525d0*sqrtis +
     &		    log(1.d0 - 0.001005d0*s)) 


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

C ----------------------------------------------------------------------

C CALCULATE THE AIR-SEA GAS EXCHANGE OF CO2
C CALCULATE THE AIR-SEA GAS EXCHANGE OF CO2
C CALCULATE THE AIR-SEA GAS EXCHANGE OF CO2
C CALCULATE THE AIR-SEA GAS EXCHANGE OF CO2
C CALCULATE THE AIR-SEA GAS EXCHANGE OF CO2

C PREPROCESSOR OPTIONS:
 
C CALLED BY:	bgc_model

C CALLS:	co2_surface

C INPUT/ARGUMENT LIST:	
C vgas660	exchange coefficient, depends on wind
C atmosp0       atmopheric pressure
C surf_dic      surface DIC [mmol C/m3]
C surf_pho      surface PO4 [mmol P/m3]
C surf_alk      surface alkalinity [mmol eq/m3]
C surf_sil      surface silicate [mmol Si/m3]
C ttemp		surface temperature
C
C INPUT/COMMON BLOCK: 
C scar1		coefficient for Schmidt number, from BGC_INI.F
C scar2		coefficient for Schmidt number, from BGC_INI.F
C scar3		coefficient for Schmidt number, from BGC_INI.F
C scar4		coefficient for Schmidt number, from BGC_INI.F
C atmosp	atmospheric pressure, from BGC_INI.F
C pCO2atm	atmospheric pCO2
C
C OUTPUT/ARGUMENT LIST: 
C co2ex		exchange rate for air-sea gas exchange [mmol C/m2/d]
C co2emp	exchange rate for E-P [mmol C/m2/d]

      subroutine co2_surfforcing(vgas660,atmosp,
     &           surf_dic,surf_pho,surf_alk,surf_sil,ttemp,co2ex,co2emp)

      implicit none

      real*8 surf_dic,surf_pho,surf_alk,surf_sil
      real*8 vgas660,atmosp,ttemp,co2ex,co2emp

#include "CAR_PARAMS.h"

! local coefficients

      real*8 SchmidtNoCO2,kwexch,co2starair,co2star,co2sol
      real*8 sdic,spho,ssil,salk   

      SchmidtNoCO2 = scar1 - scar2*ttemp + scar3*ttemp*ttemp 
     &	 - scar4*ttemp*ttemp*ttemp
	
      KWexch = vgas660/sqrt(SchmidtNoCO2/660.0d0)

c calculate co2star = pCO2 in the surface water 
c calculate co2sol = solubility of CO2 in water


      sdic = surf_dic/convert_mol_to_mmol
      spho = surf_pho/convert_mol_to_mmol
      ssil = surf_sil/convert_mol_to_mmol
      salk = surf_alk/convert_mol_to_mmol

      call co2_surface(sdic,spho,ssil,salk,co2star,co2sol)
      
c sol is solubility of CO2 in mol/(m3*uatm)
C equilibrium [CO2]aq in mol/m^3 = sol*pCO2_atm*atmpres, where
C pCO2_atm = atmospheric mole fraction CO2 in dry air at 1 atm total pres (ppmv)
C atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)

      co2starair = co2sol*pco2atm*atmosp ! equilibrium CO2aq in mol/m^3
      co2ex=-KWexch*(co2star - co2starair)*convert_mol_to_mmol

C E-P: emp in m/s
      co2emp = dicgave*emp*86400.0d0
      
      return
      end


C ----------------------------------------------------------------------

C IK CALCULATE THE SURFACE CO2 FOR AIR SEA GAS EXCHANGE
C IK CALCULATE THE SURFACE CO2 FOR AIR SEA GAS EXCHANGE
C IK CALCULATE THE SURFACE CO2 FOR AIR SEA GAS EXCHANGE
C IK CALCULATE THE SURFACE CO2 FOR AIR SEA GAS EXCHANGE
C IK CALCULATE THE SURFACE CO2 FOR AIR SEA GAS EXCHANGE

C PREPROCESSOR OPTIONS:
 
C CALLED BY:	mops_biogeochem_ini
C		co2_surfforcing

C CALLS:	

C INPUT/ARGUMENT LIST:	
C sdic		surface DIC [mol C/m3]
C spho		surface phosphate [mol P/m3]
C ssil		surface silicate [mol Si/m3]
C sta		surface alkalinity [mol eq/m3]
C
C INPUT/COMMON BLOCK: 
C sph		surface pH, first from BGC_INI, then this routine
C ak1p(1)	k1p = [H][H2PO4]/[H3PO4], from CAR_COFFS
C ak2p(1)	k2p = [H][HPO4]/[H2PO4], from CAR_COFFS
C ak3p(1)	k3p = [H][PO4]/[HPO4], from CAR_COFFS
C akb(1)	kb = [H][BO2]/[HBO2], from CAR_COFFS
C aksi(1)	ksi = [H][SiO(OH)3]/[Si(OH)4], from CAR_COFFS
C ff(1)		f = k0(1-pH2O), from CAR_COFFS
C
C OUTPUT/ARGUMENT LIST: 
C co2s 		surface DIC saturation [mol C/m3]
C sol		solubility [mol C/m3]

      subroutine co2_surface(sdic,spho,ssil,sta,co2s,sol)

      implicit none
    
#include "CAR_PARAMS.h"

      real*8 sdic,spho,ssil,sta,co2s,sol

! local coefficients
      
      real*8 pHlocal,btlocal,kblocal,k1local,k2local,
     &       k1plocal,k2plocal,k3plocal,ksilocal,kwlocal,fflocal
      real*8 pt,sit,ta,dic
      real*8 phguess,hguess,bohg,stuff,h3po4g,h2po4g,hpo4g,po4g,
     &       siooh3g,cag,gamm,hnew

      pHlocal = sph
      btlocal = bt(1)
      kblocal = akb(1)
      k1local = ak1(1)
      k2local = ak2(1)
      k1plocal = ak1p(1)
      k2plocal = ak2p(1)
      k3plocal = ak3p(1)
      ksilocal = aksi(1)
      kwlocal = akw(1)
      fflocal = ff(1)

c change units from the input of mol/m^3 -> mol/kg:
c (1 mol/m^3)  x (1 m^3/1024.5 kg)
      pt=spho*permil
      sit=ssil*permil
      ta=sta*permil
      dic=sdic*permil

c first guess for ph, hydrogen ions, borate, phosphate, silicate
      phguess = phlocal
      hguess = 10.d0**(-phguess)
      bohg = btlocal*kblocal/(hguess+kblocal)
      stuff = hguess*hguess*hguess
     &           + (k1plocal*hguess*hguess)
     &           + (k1plocal*k2plocal*hguess)
     &           + (k1plocal*k2plocal*k3plocal)
      h3po4g = (pt*hguess*hguess*hguess) / stuff
      h2po4g = (pt*k1plocal*hguess*hguess) / stuff
      hpo4g  = (pt*k1plocal*k2plocal*hguess) / stuff
      po4g   = (pt*k1plocal*k2plocal*k3plocal) / stuff
      siooh3g = sit*ksilocal / (ksilocal + hguess)

c estimate carbonate alkalinity
      cag = ta - bohg - (kwlocal/hguess) + hguess
     &           - hpo4g - 2.d0*po4g + h3po4g
     &           - siooh3g

c second guess of hydrogen ions
      gamm  = dic/cag
      stuff = (1.d0-gamm)*(1.d0-gamm)*k1local*k1local
     &          - 4.d0*k1local*k2local*(1.d0-2.d0*gamm)
      hnew  = 0.5d0*( (gamm-1.d0)*k1local + sqrt(stuff) )

c co2*
      co2s  = dic/
     &   (1.d0 + (k1local/hnew) + (k1local*k2local/(hnew*hnew)))
c pH
      sph = -log10(hnew)

c co2* converted from mol/kg to mol/m3 
      co2s = co2s/permil

c fflocal is the solubility (computed in car_coeffs) in mol/(kg*atm)
C To convert to mol/(m^3*uatm), multiply ff by 1e-6*1024.5, i.e.

c solubility of CO2 in mol/(m3*uatm)
      sol=fflocal*permeg*rho0 

C equilibrium [CO2]aq in mol/m^3 = sol*pCO2_atm*atmpres, where
C pCO2_atm = atmospheric mole fraction CO2 in dry air at 1 atm total pres (ppmv)
C atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)

      return
      end

C ----------------------------------------------------------------------

C HERE MAY LATER FOLLOW THE CALCULATION OF CARBONATE UNDERSATURATION
C HERE MAY LATER FOLLOW THE CALCULATION OF CARBONATE UNDERSATURATION
C HERE MAY LATER FOLLOW THE CALCULATION OF CARBONATE UNDERSATURATION
C HERE MAY LATER FOLLOW THE CALCULATION OF CARBONATE UNDERSATURATION
C HERE MAY LATER FOLLOW THE CALCULATION OF CARBONATE UNDERSATURATION

C INPUT/ARGUMENT LIST:	

C INPUT/COMMON BLOCK: 

C OUTPUT/ARGUMENT LIST: 
C dicc		DIC changes due to carbonate dissolution [mmol C/m3/d]
C alkc		Alk changes due to carbonate dissolution [mmol C/m3/d]

c
c      subroutine car_chemistry(...,dicc,alkc,n)
c
c      implicit none
c
c      integer n,i
c      real*8 dicc(n),alkc(n)
c
c#include "CAR_PARAMS.h"
c
! local coefficients
c
c      real*8 
c
c      do i=1,n
c
c      enddo
c      
c      return
c      end

