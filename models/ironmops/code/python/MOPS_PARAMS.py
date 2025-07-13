import numpy as np

ppmToPgC=2.1324

# C NUMERICAL ISSUES RELATED TO BIOGEOCHEMISTRY
vsafe = 1.0e-6

ipo4=0    #PO4
idop=1    #DOP
ioxy=2    #Oxygen
iphy=3    #Phyto-P
izoo=4    #Zoo-P
idet=5    #Detritus-P
idin=6    #DIN
idfe=7    #dissolved Fe
ipfe=8    #particulate Fe
idic=9
ialk=10

# C BGC PARAMETERS FOR SURFACE OCEAN BIOLOGY      
# 
# ! Stoichiometry
rcp=117.0 #redfield ratio C:P
rnp=16.0 #redfield ratio N:P
ro2ut = 169.2545586809128568828164418391679646447301  # L4-SO of Kriest et al. (2023) 
subox = 1.0 #no oxic degradation below this level
subdin = 16.0 #no denitrification below this level

# N2-Fixatioon
# Factors tf2, tf1 and tf0 are a polynomial (2nd order) 
# approximation to the functional relationship by Breitbarth et al. (2007),
# for temperature dependence of Trichodesmium growth, 
# their eq. (2), assuming that their powers relate to "e" (10^), and
# that the last but one sign has to be reversed.
# The relation will be scaled to their max. growth rate.
# Note that the second order approx. is basically similar to their
# function 2 for T-dependent nitrogen fixation multiplied by 4 
# (2 [N atoms per mole] * 12 [light hrs per day]/6 [C-atoms per N-atoms])
tf2 = -0.0042
tf1 = 0.2253
tf0 = -2.7819
tff = 0.2395  
nfix = 0.002272073044 #[mmol N/m3/d]; if PO4>0 and DIN = 0 and T=26.82 this corresponds 
# 	                  !to an integral over 100 m of ~ 200 umol N/m2/d, a 
# 		              !rather high value acc. to Mahaffey et al., 2005
# 		              !alternatively, 2 umol/m3/d = 0.08 nmol/L/hr,
# 		              !a rather high value acc. to Mulholland (2007)
# 
# Phytoplankton parameters
TempB = 15.65 #ref. temperature for T-dependent growth
ACmuphy = 0.6 #max. growth rate [1/day]
ACik = 22.13898097109098993262588095376486307941377 # L4-SO of Kriest et al. (2023) 
ACkpo4 = 0.4995 #half-saturation constant for PO4 uptake [mmol P/m3]
ACkchl  = 0.03*rnp #att. of Phy [1/(m*mmol P/m3)]
ACkw = 0.04 #attenuation of water [1/m]
AClambda = 0.03      #exudation rate [1/day]
AComni = 0.0      #density dep. loss rate [m3/(mmol P * day)]
plambda = 0.01 #phytoplankton mortality

# Zooplankton parameters
ACmuzoo=2.384719471821353179189981186247848654602421 # L4-SO of Kriest et al. (2023) 
ACkphy=np.sqrt(ACmuzoo/1.0)/rnp #zoo half-saturation constant [mmol P]
AClambdaz=0.03 #zooplankton excretion [1/d]
AComniz=4.548 #zooplankton mortality [m3/(mmol P * day)]
ACeff=0.75 #assimilation efficiency
zlambda=0.01 #zooplankton mortality

# DOPparameters
graztodop = 0.00 # fraction grazing that goes into DOP
dlambda  = 0.0005133333049196715389730860613828195004870736 # base value of L4-SO in Kriest et al. (2023) 

# A minimum value for the degradation of PHY and DOP
alimit = 1.0e-3

# Detritus parameters
detlambda = 0.05 #detritus remineralisation rate [1/d]    
detmartin = 1.022589344764692176467310580356695481896168 # L4-SO of Kriest et al. (2023)

# Parameters specific to MOPS: (Kriest and Oschlies 2013, 2015)
burdige_fac = 1.6828 # factor for sediment burial (see Kriest and Oschlies, 2013)
burdige_exp = 0.799  # exponent for sediment burial (see Kriest and Oschlies, 2013)
ACkbaco2 = 1.145532  #Half sat.-constant for oxic degradation (see Kriest and Oschlies, 2015)
ACkbacdin = 23.083559 #Half sat.-constant for suboxic degradation (see Kriest and Oschlies, 2015)

# Parameters for the Fe cycle: all have been converted to deal with Fe in umol/m3 and P in mmol/m3
rfep = 1.06           # Fe:P uptake ratio and internal ratio in organics [umol Fe/mmol P], Nickelsen et al. (2015), conversion: 66.25 [umol Fe/mol N]*16/1000 
kfe  = 0.04           # half-saturation constant for dFe uptake by phytoplankton [umol Fe/m3], Nickelsen et al. (2015) 
kfeorg = 2.3806      # scavenging rate dependent on organic matter [1/((mmol P/m3)**0.58) 1/d], Somes et al. (2021), , Atm+SedHigh_LigVar, conversion: 2.9 [1/(gC/m3)^0.58 1/d] / [1000/(117*12.011)]^0.58
kfeleq = 10.0**2.5  # Fe-ligand stability constant [1/(umol ligand/m3)], Nickelsen et al. (2015), conversion: 10^11.5 [1/(mol lig/L)] * 1e-9 [mol/umol]
kfepre = 0.85         # inorganic scavenging rate [1/(umol Fe/m3) 1/d], Somes et al. (2021), Atm+SedHigh_LigVar, conversion: 850 [1/(mmol Fe/m3) 1/d] / 1000 
fealpha = 0.015       # factor for AOU dependence of ligands [(umol ligand/m3)/(mmol O2/m3)], Somes et al. (2021)
febeta  = 1.92981     # factor for DOP dependence of ligands [(umol ligand/m3)/(mmol DOP/m3)], Somes et al. (2021), conversion: 0.21 [(umol ligand/m3)/(mmol DON/m3)^0.8] * [16 mmol DON / mmol DOP]^0.8
fesedmax = 100.0      # maximum Fe release from sediment [umol Fe/m2/d], Somes et al. (2021)

# ifdef CARBON
frac_caco3 = 0.032 
length_caco3 = 4289.4 
# endif

# Derived parameters
rhno3ut = 0.8*ro2ut - rnp # -HNO3:P ratio for denitrification

detwa = detlambda/detmartin #w = a*z+b for aphotic layers
detwb = 0.0 #offset for detritus sinking


#SPK
#ifdef BUDGET
# ! Stoichiometry required for budget calculations        
# ! H-content (H atoms per P atom) of organic matter follow Anderson stoichiometry for H-Atoms
# ! O-content (O atoms per P atom) of organic matter from H, C, N and P content
# ! H20 produced per mole Porg remineralized after Eqn. 10 of Paulmier et al. (2009) 
# ! H20 produced per mole Porg denitrified after Eqn. 18 of Paulmier et al. (2009) 
rhp   = 175.0 #Anderson, 1995
rop  = 2.0*rcp+0.5*rhp+2.5*rnp+2.5-2.0*ro2ut
orh2o = 0.5*rhp-0.5*rnp-1.5 # Eqn. 10 of Paulmier et al., 2009
srh2o =0.4*rcp+0.6*rhp-0.2*rop-1.0 # Eqn. 18 of Paulmier et al., 2009
#endif

#ifdef CARBON

# ! forc OCMIP constant alkalinity and silicate

ocmip_alkfac = 2310.0*1.0245/34.88
ocmip_silfac = 7.7 
      
#endif

ACik=16.14490902507452672366705659356966862105764
ACmuzoo=2.854995027239606325769952221982350692996988
graztodop=0.
dlambda=0.000626466572010149868724405503505811565467809
kfe=0.2489925986042371879654756111621694003588345
rfep=0.7655590068661256620847349596559183737554122
kfeorg=0.6675512738731916035760799443821156273770612
kfepre=0.04130206314545665129210182814345486690399412

