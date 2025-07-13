import numpy as np
from numba import njit

car_ktotal=100

convert_mol_to_mmol=1000.0
rho0=1024.5
permil=1.0/rho0
permeg=1.e-6
 
phlo=6.0
phhi=9.0

scar1 = 2073.1
scar2 = 125.62 
scar3 = 3.6276
scar4 = 0.043219

# C Basis for all routines: MIT source code for OCMIP abiotic model,
# C modified by Samar Khatiwala, then by Iris Kriest (Aug-Sep. 2010)
# 
# C NOTE: MIT model units are [mol/m3].
# C       Therefore convert DIC, PO4, Silicate and TA input by 
# C       factor 1d-3 (1/convert_mol_to_mmol) before using them as input.
# C       Therefore vonvert change of CO2 due to air-sea flux by 
# C       factor 1d+3 (convert_mol_to_mmol) before using it as output.
# C       Leave everything else as before.
# C       Take care when later calculating carbonate dissolution.
# C       See CO2_SURFFORCING, CAR_INI and CAR_PARAMS.
# C pCO2atm is atmospheric mole fraction of CO2 in dry air (ppmv)]
# C Multiply with total atmospheric pressure (atmosp [atm], see above)
# C to get local partial pressure of CO2 (in uatm)
# C This will be done in CO2_SURFFORCING

@njit
def car_ini(ttemp,stemp,surf_dic,surf_pho,surf_alk,surf_sil):

# C INPUT/ARGUMENT LIST: 
# C ttemp  temperature, n layers
# C stemp  salinity, n layers
# C n             number of layers
# C surf_dic surface DIC [mmol C/m3]
# C surf_pho surface PO4 [mmol P/m3]
# C surf_alk surface alkalinity [mmol eq/m3]
# C surf_sil surface silicate [mmol Si/m3]

  sdic = surf_dic/convert_mol_to_mmol
  spho = surf_pho/convert_mol_to_mmol
  ssil = surf_sil/convert_mol_to_mmol
  salk = surf_alk/convert_mol_to_mmol
      
  sph=8.0
  for i in range(10):
    co2star,co2sol,sph=__co2_surface(ttemp,stemp,sdic,spho,ssil,salk,sphin=sph)

  return sph
      
@njit
def __co2_surface(t,s,sdic,spho,ssil,sta,sphin=8.0):
#SPK This function is different from the others as the units are in mol/m^3 whereas 
# the others are mmol/m^3. This function is only called from car_ini and co2_surfforcing, 
# which convert from model units of mmol/m^3 to mol/m^3 before calling this routine. 
# As it is never called directly from MOPS I'm making this function private 
# and not exposing it to a MOPS user.

# C INPUT/ARGUMENT LIST: 
# C sdic  surface DIC [mol C/m3]
# C spho  surface phosphate [mol P/m3]
# C ssil  surface silicate [mol Si/m3]
# C sta  surface alkalinity [mol eq/m3]
# C
# C OUTPUT/ARGUMENT LIST: 
# C co2s   surface DIC saturation [mol C/m3]
# C sol  solubility [mol C/m3]
# C pH

  tk = 273.15 + t
  tk100 = tk/100.0
  tk1002=tk100*tk100
  invtk=1.0/tk
  dlogtk=np.log(tk)

  is1=19.924*s/(1000.0-1.005*s)
  is2=is1*is1
  sqrtis=np.sqrt(is1)
  s2=s*s
  sqrts=np.sqrt(s)
  s15=s**1.5
  scl=s/1.80655

# C Calculate concentrations for borate, sulfate, and fluoride
# C Uppstrom (1974), Morris & Riley (1966), Riley (1965)
  bt = 0.000232 * scl/10.811
  st = 0.14 * scl/96.062
  ft = 0.000067 * scl/18.9984
  

# C f = k0(1-pH2O)*correction term for non-ideality
# C Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values)
  ff = np.exp(-162.8301 + 218.2968/tk100  +
     90.9241*np.log(tk100) - 1.47696*tk1002 +
     s * (.025695 - .025225*tk100 + 
     0.0049867*tk1002))


# C K0 from Weiss 1974
  ak0 = np.exp(93.4517/tk100 - 60.2409 + 
     23.3585 * np.log(tk100) +
     s * (0.023517 - 0.023656*tk100 + 
     0.0047036*tk1002))


# C k1 = [H][HCO3]/[H2CO3]
# C k2 = [H][CO3]/[HCO3]     on hSWS
# C Millero p.664 (1995) using Mehrbach et al. data on SEAWATER scale 
# C (Original reference: Dickson and Millero, DSR, 1987)
  ak1=10**(-1.0*(3670.7*invtk - 
     62.008 + 9.7944*dlogtk -
     0.0118*s + 0.000116*s2))
  ak2=10**(-1.0*(1394.7*invtk + 4.777 - 
     0.0184*s + 0.000118*s2))

# C k1p = [H][H2PO4]/[H3PO4] on hSWS
# C Millero p.670 (1995)
  ak1p = np.exp(-4576.752*invtk + 115.540 - 
     18.453*dlogtk + 
     (-106.736*invtk + 0.69171)*sqrts +
     (-0.65643*invtk - 0.01844)*s)

# C k2p = [H][HPO4]/[H2PO4] on hSWS
# C Millero p.670 (1995)
  ak2p = np.exp(-8814.715*invtk + 172.1033 - 
     27.927*dlogtk +
     (-160.340*invtk + 1.3566)*sqrts +
     (0.37335*invtk - 0.05778)*s)

# C k3p = [H][PO4]/[HPO4] on hSWS
# C Millero p.670 (1995)
  ak3p = np.exp(-3070.75*invtk - 18.126 + 
     (17.27039*invtk + 2.81197) *
     sqrts + (-44.99486*invtk - 0.09984) * s)

# C ksi = [H][SiO(OH)3]/[Si(OH)4] on hSWS
# C Millero p.671 (1995) using data from Yao and Millero (1995)
# C change to (mol/ kg soln)
  aksi = np.exp(-8904.2*invtk + 117.400 - 
     19.334*dlogtk +
     (-458.79*invtk + 3.5913) * sqrtis +
     (188.74*invtk - 1.5998) * is1 +
     (-12.1652*invtk + 0.07871) * is2 +
     np.log(1.0-0.001005*s))

# C kw = [H][OH] on hSWS
# C Millero p.670 (1995) using composite data
  akw = np.exp(-13847.26*invtk + 148.9802 - 
     23.6521*dlogtk +
     (118.67*invtk - 5.977 + 1.0495 * dlogtk) *
     sqrts - 0.01615 * s)

# C ks = [H][SO4]/[HSO4] on free H scale
# C Dickson (1990, J. chem. Thermodynamics 22, 113)
# C change to (mol/ kg soln)
  aks=np.exp(-4276.1*invtk + 141.328 - 
     23.093*dlogtk +
     (-13856.0*invtk + 324.57 - 47.986*dlogtk)*sqrtis +
     (35474.0*invtk - 771.54 + 114.723*dlogtk)*is1 -
     2698.0*invtk*is1**1.5 + 1776.0*invtk*is2 +
     np.log(1.0 - 0.001005*s))

# C kf = [H][F]/[HF] on free H scale
# C Dickson and Riley (1979)
# C change to (mol/ kg soln)
  akf=np.exp(1590.2*invtk - 12.641 + 1.525*sqrtis +
     np.log(1.0 - 0.001005*s)) 

# C kb = [H][BO2]/[HBO2] on hSWS
# C Dickson p.673 (1990)
# C change from htotal to hSWS
  akb=np.exp( (-8966.90 - 2890.53*sqrts - 77.942*s +
     1.728*s15 - 0.0996*s2)*invtk +
     (148.0248 + 137.1942*sqrts + 1.62142*s) +
     (-24.4344 - 25.085*sqrts - 0.2474*s) *
     dlogtk + 0.053105*sqrts*tk +
     np.log((1.0+(st/aks)+(ft/akf)) 
     /(1.0+(st/aks))) )

#SPK: sph needs to be passed in and out of this as it is in a common block
  pHlocal = sphin

# c change units from the input of mol/m^3 -> mol/kg:
# c (1 mol/m^3)  x (1 m^3/1024.5 kg)
  pt=spho*permil
  sit=ssil*permil
  ta=sta*permil
  dic=sdic*permil

# c first guess for ph, hydrogen ions, borate, phosphate, silicate
  phguess = pHlocal
  hguess = 10.0**(-phguess)
  bohg = bt*akb/(hguess+akb)
  stuff = hguess*hguess*hguess + \
    (ak1p*hguess*hguess) + \
    (ak1p*ak2p*hguess) + \
    (ak1p*ak2p*ak3p)
  h3po4g = (pt*hguess*hguess*hguess) / stuff
  h2po4g = (pt*ak1p*hguess*hguess) / stuff
  hpo4g  = (pt*ak1p*ak2p*hguess) / stuff
  po4g   = (pt*ak1p*ak2p*ak3p) / stuff
  siooh3g = sit*aksi / (aksi + hguess)

# c estimate carbonate alkalinity
  cag = ta - bohg - (akw/hguess) + hguess \
    - hpo4g - 2.0*po4g + h3po4g \
    - siooh3g

# c second guess of hydrogen ions
  gamm  = dic/cag
  stuff = (1.0-gamm)*(1.0-gamm)*ak1*ak1 \
    - 4.0*ak1*ak2*(1.0-2.0*gamm)
  hnew  = 0.50*( (gamm-1.0)*ak1 + np.sqrt(stuff) )

# c co2*
  co2s  = dic/ \
    (1.0 + (ak1/hnew) + (ak1*ak2/(hnew*hnew)))
# c pH
  sph = -np.log10(hnew)

# c co2* converted from mol/kg to mol/m3 
  co2s = co2s/permil

# c ff is the solubility (computed in car_coeffs) in mol/(kg*atm)
# C To convert to mol/(m^3*uatm), multiply ff by 1e-6*1024.5, i.e.
# 
# c solubility of CO2 in mol/(m3*uatm)
  sol=ff*permeg*rho0 

# C equilibrium [CO2]aq in mol/m^3 = sol*pCO2_atm*atmpres, where
# C pCO2_atm = atmospheric mole fraction CO2 in dry air at 1 atm total pres (ppmv)
# C atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)

  return co2s, sol, sph
  
@njit
def co2_surfforcing(ttemp,stemp,vgas660,atmosp,surf_dic,surf_pho,surf_alk,surf_sil,pCO2atm,emp,dicgave,alkgave,sphin=8.0):

# C INPUT/ARGUMENT LIST: 
# C vgas660 exchange coefficient, depends on wind
# C atmosp0       atmopheric pressure
# C surf_dic      surface DIC [mmol C/m3]
# C surf_pho      surface PO4 [mmol P/m3]
# C surf_alk      surface alkalinity [mmol eq/m3]
# C surf_sil      surface silicate [mmol Si/m3]
# C ttemp  surface temperature
# C
# C OUTPUT/ARGUMENT LIST: 
# C co2ex  exchange rate for air-sea gas exchange [mmol C/m2/d]
# C co2emp exchange rate for E-P [mmol C/m2/d]
# C alkemp exchange rate for E-P [mmol Alk/m2/d]

  SchmidtNoCO2 = scar1 - scar2*ttemp + scar3*ttemp*ttemp \
    - scar4*ttemp*ttemp*ttemp

  KWexch = vgas660/np.sqrt(SchmidtNoCO2/660.0)

# c calculate co2star = pCO2 in the surface water 
# c calculate co2sol = solubility of CO2 in water

  sdic = surf_dic/convert_mol_to_mmol
  spho = surf_pho/convert_mol_to_mmol
  ssil = surf_sil/convert_mol_to_mmol
  salk = surf_alk/convert_mol_to_mmol

#       call co2_surface(sdic,spho,ssil,salk,co2star,co2sol)
  co2star, co2sol, sphout = __co2_surface(ttemp,stemp,sdic,spho,ssil,salk,sphin=sphin)
  
# c sol is solubility of CO2 in mol/(m3*uatm)
# C equilibrium [CO2]aq in mol/m^3 = sol*pCO2_atm*atmpres, where
# C pCO2_atm = atmospheric mole fraction CO2 in dry air at 1 atm total pres (ppmv)
# C atmpres= atmospheric pressure in atmospheres (1 atm==1013.25mbar)

  co2starair = co2sol*pCO2atm*atmosp # equilibrium CO2aq in mol/m^3
  co2ex=-KWexch*(co2star - co2starair)*convert_mol_to_mmol

#SPK: need to pass emp, dicgave and alkgave to this routine as they are in common blocks
# C E-P: emp in m/s
  co2emp = dicgave*emp*86400.0
  alkemp = alkgave*emp*86400.0
  
  return sphout, co2ex, co2emp, alkemp