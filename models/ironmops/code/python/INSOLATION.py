import numpy as np
from numba import jit

# This is based on the insolation routine in MITgcm

solar = 1360.   #solar constant
albedo = 0.6    #planetary albedo

@jit
def insolation(Time,Yc,daysperyear,sfac,stau):

# C !INPUT PARAMETERS: ===================================================
# C time                 :: current time
# C Yc
# C daysperyear
# 
# C !OUPUT PARAMETERS: ===================================================
# C sfac,stau

#SPK This routine expects array arguments. The output arguments should also be predefined 
#    numpy arrays of the same size as Yc. No checking is done. 
  N = len(Yc)
  
# c find day (****NOTE for year starting in winter*****)
  dayfrac=(
    np.mod(Time,daysperyear*86400.0) 
    /(daysperyear*86400.0)          #fraction of year
    )
  yday = 2.0*3.1416*dayfrac                         #convert to radians
  delta = (
    (0.006918 - (0.399912*np.cos(yday))    # cosine zenith angle
    +(0.070257*np.sin(yday))           #(paltridge+platt)
    -(0.006758*np.cos(2.0*yday)) 
    +(0.000907*np.sin(2.0*yday)) 
    -(0.002697*np.cos(3.0*yday)) 
    +(0.001480*np.sin(3.0*yday)) )
    )
  for j in range(N):
# c latitude in radians
    lat=Yc[j]/180.0*3.1416
    sun1 = -np.sin(delta)/np.cos(delta) * np.sin(lat)/np.cos(lat)
    if (sun1 <= -0.999): sun1=-0.999
    if (sun1 >= 0.999): sun1= 0.999
    dayhrs = np.abs(np.arccos(sun1))
    cosz = (
      ( np.sin(delta)*np.sin(lat)+            #average zenith angle
      (np.cos(delta)*np.cos(lat)*np.sin(dayhrs)/dayhrs) )
      )
    if (cosz <= 0.005): cosz=0.005
    frac = dayhrs/3.1416               #fraction of daylight in day
#  c daily average photosynthetically active solar radiation just below surface
    fluxi = solar*(1.0-albedo)*cosz*frac
#  c
#  c convert to sfac
    if (fluxi > 0.0): sfac[j]=fluxi
#  c very large for polar night
    if (fluxi < 0.00001): sfac[j]=0.00001
#  c daylength; ensure that it lies between 0 and 1 (may be slightly
#  c out of this range in high latitudes)
    fracmin = np.fmin(frac,1.0)
    stau[j] = np.fmax(fracmin,0.0)
