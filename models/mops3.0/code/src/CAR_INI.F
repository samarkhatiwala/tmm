C$Header: /Users/ikriest/CVS/mops/CAR_INI.F,v 1.2 2016/06/03 09:28:59 ikriest Exp $
C$Name: mops-2_0 $

C Basis for all routines: MIT source code for OCMIP abiotic model,
C modified by Samar Khatiwala, then by Iris Kriest (Aug-Sep. 2010)

C NOTE: MIT model units are [mol/m3].
C       Therefore convert DIC, PO4, Silicate and TA input by 
C       factor 1d-3 (1/convert_mol_to_mmol) before using them as input.
C       Therefore vonvert change of CO2 due to air-sea flux by 
C       factor 1d+3 (convert_mol_to_mmol) before using it as output.
C       Leave everything else as before.
C       Take care when later calculating carbonate dissolution.
C       See CO2_SURFFORCING, CAR_INI and CAR_PARAMS.

C THIS FILE CONTAINS:
C SR car_ini	- initialize carbon model
 
C INITIALIZE THE COEFFICIENTS OF THE CARBONATE SYSTEM OVER N LAYERS
C INITIALIZE THE COEFFICIENTS OF THE CARBONATE SYSTEM OVER N LAYERS
C INITIALIZE THE COEFFICIENTS OF THE CARBONATE SYSTEM OVER N LAYERS
C INITIALIZE THE COEFFICIENTS OF THE CARBONATE SYSTEM OVER N LAYERS
C INITIALIZE THE COEFFICIENTS OF THE CARBONATE SYSTEM OVER N LAYERS

C PREPROCESSOR OPTIONS:
 
C CALLED BY:	mops_biogeochem_ini

C CALLS:	car_coeffs
C		co2_surface

C INPUT/ARGUMENT LIST:	
C ttemp		temperature, n layers
C stemp		salinity, n layers
C n             number of layers
C surf_dic	surface DIC [mmol C/m3]
C surf_pho	surface PO4 [mmol P/m3]
C surf_alk	surface alkalinity [mmol eq/m3]
C surf_sil	surface silicate [mmol Si/m3]

C OUTPUT/COMMON BLOCK
C pCO2atm 	atmospheric pCO2 
C atmosp 	atmopheric pressure
C convert_ik0	convert from mol/m3 to mmol/m3
C rho0		density of seawater to convert from m3 to kg
C permil	1/rho0
C permeg	convert to 1 per million
C phlo		lower pH boundary
C phhi		upper pH boundary 
C sph		initial estimate of pH
C scar1 	coefficient for Schmidt number
C scar2 	coefficient for Schmidt number
C scar3 	coefficient for Schmidt number
C scar4 	coefficient for Schmidt number
C ocmip_alkfac	factor to calculate alkalinit from salt
C ocmip_silfac	fixed silicate concentration

      subroutine car_ini(ttemp,stemp,n,
     &           surf_dic,surf_pho,surf_alk,surf_sil)
       
      integer n,i
      real*8 ttemp(n),stemp(n)
      real*8 sdic,spho,ssil,salk
      real*8 surf_dic,surf_pho,surf_alk,surf_sil
      real*8 sdummy,co2star,co2sol
      
#include "CAR_PARAMS.h"

C pCO2atm is atmospheric mole fraction of CO2 in dry air (ppmv)]
C Multiply with total atmospheric pressure (atmosp [atm], see above)
C to get local partial pressure of CO2 (in uatm)
C This will be done in CO2_SURFFORCING

      pCO2atm = 277.9514d0 !preindustrial value
	
      convert_mol_to_mmol=1000.0d0
      rho0=1024.5d0
      permil=1.d0/rho0
      permeg=1.d-6
	
      phlo=6.d0
      phhi=9.d0
      sph = 8.0d0

      scar1 = 2073.1d0
      scar2 = 125.62d0 
      scar3 = 3.6276d0
      scar4 = 0.043219d0
      
      call car_coeffs(ttemp,stemp,n)
      
      DO i=n+1,car_ktotal
        bt(i)  = 0.0d0
        st(i)  = 0.0d0
        ft(i)  = 0.0d0
        ff(i)  = 0.0d0
        ak0(i) = 0.0d0
        ak1(i) = 0.0d0
        ak2(i) = 0.0d0
        ak1p(i)= 0.0d0
        ak2p(i)= 0.0d0
        ak3p(i)= 0.0d0
        aksi(i)= 0.0d0
        akw(i) = 0.0d0
        aks(i) = 0.0d0
        akf(i) = 0.0d0
        akb(i) = 0.0d0
      ENDDO

      sdic = surf_dic/convert_mol_to_mmol
      spho = surf_pho/convert_mol_to_mmol
      ssil = surf_sil/convert_mol_to_mmol
      salk = surf_alk/convert_mol_to_mmol

      do i=1,10
        call co2_surface(sdic,spho,ssil,salk,co2star,co2sol)
      enddo
      
      return
      end

