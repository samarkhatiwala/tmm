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
ro2ut = 194.7075436894855925140657149086109711788595  # from optimisation  O3
subox = 1.0 #no oxic degradation below this level; default value of Kriest et al. (2020)
subdin = 15.77129497301908116636132151100468945514876 # no denitrification below this level, optimal value of MIT28* of Kriest et al. (2020) 

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
nfix = 0.001194890759152962841758761653333104080232374 #[mmol N/m3/d]; optimal value of MIT28* of Kriest et al. (2020)

# Phytoplankton parameters
TempB = 15.65 #ref. temperature for T-dependent growth
ACmuphy = 0.6 #max. growth rate [1/day]
ACik = 11.70844406776519093082294942220755729067605 # Light half-saturation constant [W/m2]; from optimisation O3
ACkpo4 = 0.4995 #half-saturation constant for PO4 uptake [mmol P/m3]
ACkchl  = 0.03*rnp #att. of Phy [1/(m*mmol P/m3)]
ACkw = 0.04 #attenuation of water [1/m]
AClambda = 0.03      #exudation rate [1/day]
AComni = 0.0      #density dep. loss rate [m3/(mmol P * day)]
plambda = 0.01 #phytoplankton mortality

# Zooplankton parameters
ACmuzoo= 3.000 # max. grazing rate [1/d]; from optimisation O3
ACkphy=np.sqrt(ACmuzoo/1.0)/rnp #zoo half-saturation constant [mmol P]
AClambdaz=0.03 #zooplankton excretion [1/d]
AComniz=4.548 #zooplankton mortality [m3/(mmol P * day)]
ACeff=0.75 #assimilation efficiency
zlambda=0.01 #zooplankton mortality

# DOPparameters
graztodop = 0.0003079415719245463227006268081518969292886823 # fraction grazing that goes into DOP; from optimisation O3
dlambda  = 0.000662827652140595182064180623891097576816378 #DOP remineralization rate [1/day]  (SLOW recycling); from optimisation O3

# A minimum value for the degradation of PHY and DOP
alimit = 1.0e-3

# Detritus parameters
detlambda = 0.05 #detritus remineralisation rate [1/d]    
detmartin = 1.0 # from optimisation O3

# Parameters specific to MOPS: (Kriest and Oschlies 2013, 2015)
burdige_fac = 1.6828 # factor for sediment burial (see Kriest and Oschlies, 2013)
burdige_exp = 0.799  # exponent for sediment burial (see Kriest and Oschlies, 2013)
ACkbaco2 = 1.0  #Half sat.-constant for oxic degradation, optimal value of MIT28* of Kriest et al. (2020)
ACkbacdin = 31.95521294352302090965856073978557105874643 #Half sat.-constant for suboxic degradation, optimal value of MIT28* of Kriest et al. (2020)

# For Fe: Parameters all parameters converted to deal with Fe in umol/m3 and P in mmol/m3
rfep = 1.06              # Fe:P uptake ratio and internal ratio in organics [umol Fe/mmol P], Nickelsen et al. (2015), conversion: 66.25 [umol Fe/mol N]*16/1000 
kfemin  = 0.04           # minimum half-saturation constant for dFe uptake by phytoplankton [umol Fe/m3] at low phytoplankton concentrations, Nickelsen et al. (2015) 
kfemax  = 0.40           # maximim half-saturation constant for dFe uptake by phytoplankton [umol Fe/m3] at high phytoplankton concentrations, Nickelsen et al. (2015) 
pmax   = 0.01            # phytoplankton concentration at which kfe starts to increase [mmol P/m3], 0.15 [mmol N/m3] in Nickelsen et al. (2015) 
kfeleq = 10.0**2.5     # Fe-ligand stability constant [1/(umol ligand/m3)], Nickelsen et al. (2015), conversion: 10^11.5 [1/(mol lig/L)] * 1e-9 [mol/umol]
fdoclig = 10.0           # DOC-dependent concentration of ligands as in as in Tagliabue & Voelker (2011; 0.09 1/(mmol C/m3) 1/d) and Aumont et al. (2015) 
flig = 0.60              # constant concentration of total ligands, as in Tagliabue & Voelker (2011; they have 0.4) and Aumont et al. (2015; they have 0.6) 
lambdafemin = 0.00010    # independent dFe scavenging rate [1/d], loosely following Aumont et al. (2015; they have 3e-5)
lambdafeorg = 0.50       # POC-dependent dFe scavenging rate [1/(mmol P/m3) 1/d], loosely following Aumont et al. (2015; they have 0.005 1/d 1/umol/L but also a large aggregation rate of Fecoll (half of Fe-Fe'))
lambdafepre = 0.0002       # quadratic inorganic scavenging rate [1/(umol Fe'/m3) 1/d], Somes et al. (2021), Atm+SedHigh_LigVar, conversion: 850 [1/(mmol Fe/m3) 1/d] / 1000 
fesedmax = 0.0          # maximum Fe release from sediment [umol Fe/m2/d], Somes et al. (2021; have 100)

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

# Third version for Volkmar's optimization tests, 16 October 2015
# rfep      = parvec(1)
# fdoclig   = parvec(2)
# lambdafeorg = parvec(3)
# lambdafepre = parvec(4)

# ACik=16.14490902507452672366705659356966862105764
# ACmuzoo=2.854995027239606325769952221982350692996988
# graztodop=0.
# dlambda=0.000626466572010149868724405503505811565467809
# kfe=0.2489925986042371879654756111621694003588345
# rfep=0.7655590068661256620847349596559183737554122
# kfeorg=0.6675512738731916035760799443821156273770612
# kfepre=0.04130206314545665129210182814345486690399412

