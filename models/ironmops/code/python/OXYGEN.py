import numpy as np
from numba import jit

# C AIR-SEA GAS EXCHANGE OF O2, taken from MIT model
# 
# C define Schmidt no. coefficients for O2
# C based on Keeling et al [GBC, 12, 141, (1998)]
sox1 = 1638.0
sox2 = -81.83
sox3 = 1.483
sox4 = -0.008004
# 
# C coefficients for determining saturation O2
oA0=  2.00907
oA1=  3.22014
oA2=  4.05010
oA3=  4.94457
oA4= -2.56847e-1
oA5=  3.88767
oB0= -6.24523e-3
oB1= -7.37614e-3
oB2= -1.03410e-2
oB3= -8.17083e-3
oC0= -4.88682e-7

@jit
def o2saturation(ttemp,stemp):

# C INPUT/ARGUMENT LIST: 
# C ttemp  temperature
# C stemp  salinity

  aTT  = 298.15 -ttemp
  aTK  = 273.15 +ttemp
  aTS  = np.log(aTT/aTK)
  aTS2 = aTS*aTS
  aTS3 = aTS2*aTS
  aTS4 = aTS3*aTS
  aTS5 = aTS4*aTS
  oCnew  = oA0+oA1*aTS+oA2*aTS2+oA3*aTS3+oA4*aTS4+oA5*aTS5 \
    + stemp*(oB0+oB1*aTS+oB2*aTS2+oB3*aTS3)+oC0*(stemp*stemp)
  o2s = np.exp(oCnew)

# ! Convert from ml/l to mmol/m^3
# ! Note: o2 in mit is in mol/m3; I use mmol/m3, thus coonvert with 1d3
  o2sat = o2s/22391.6 * 1.e3*1.e3

  return o2sat

@jit
def o2_surfforcing(vgas660,atmosp,ttemp,stemp,soxy,ssat):

# C CALCULATE THE AIR-SEA GAS EXCHANGE OF O2      
# C INPUT/ARGUMENT LIST: 
# C vgas660 exchange coefficient, depends on wind
# C atmosp0       atmopheric pressure
# C ttemp  surface temperature
# C stemp  surface salinity
# C soxy  surface oxygen
# C ssat     oxygen saturation
# 
# C OUTPUT/ARGUMENT LIST: 
# C o2ex  exchange rate [mmol O2/m2/d]

  SchmidtNoO2=sox1+sox2*ttemp+sox3*ttemp*ttemp  \
    + sox4*ttemp*ttemp*ttemp

  Kwexch = vgas660/np.sqrt(SchmidtNoO2/660.0)

# ! Determine flux, inc. correction for local atmos surface pressure
  o2ex = Kwexch*(atmosp*ssat-soxy)
  
  return o2ex
  
