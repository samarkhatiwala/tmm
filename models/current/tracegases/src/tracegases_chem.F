      real*8 FUNCTION N2O_schmidt_number(pt)

C     N2O schmidt number from Wanninkhof (Limnol. Oceanogr. Methods, 2014, Table 1)

      IMPLICIT NONE
      
      real*8 pt ! temperature, deg C

      real*8 A1, A2, A3, A4, A5
      parameter (A1=2356.2d0, A2=-166.38d0, A3=6.3952d0, 
     &           A4=-0.13422d0, A5=0.0011506d0)

      N2O_schmidt_number = 
     &    A1 + A2*pt + A3*(pt**2) + A4*(pt**3) + A5*(pt**4)

      return
      
      END

      real*8 FUNCTION CFC11_schmidt_number(pt)

C     CFC-11 schmidt number from Wanninkhof (Limnol. Oceanogr. Methods, 2014, Table 1)

      IMPLICIT NONE
      
      real*8 pt ! temperature, deg C

      real*8 A1, A2, A3, A4, A5
      parameter (A1=3579.2d0, A2=-222.63d0, A3=7.5749d0, 
     &           A4=-0.14595d0, A5=0.0011874d0)

      CFC11_schmidt_number = 
     &    A1 + A2*pt + A3*(pt**2) + A4*(pt**3) + A5*(pt**4)

      END

      real*8 FUNCTION CFC12_schmidt_number(pt)

C     CFC-12 schmidt number from Wanninkhof (Limnol. Oceanogr. Methods, 2014, Table 1)

      IMPLICIT NONE
      
      real*8 pt ! temperature, deg C

      real*8 A1, A2, A3, A4, A5
      parameter (A1=3828.1, A2=-249.86d0, A3=8.7603d0, 
     &           A4=-0.1716d0, A5=0.001408d0)

      CFC12_schmidt_number = 
     &    A1 + A2*pt + A3*(pt**2) + A4*(pt**3) + A5*(pt**4)

      END     

      real*8 FUNCTION SF6_schmidt_number(pt)

C     SF6 schmidt number from Wanninkhof (Limnol. Oceanogr. Methods, 2014, Table 1)

      IMPLICIT NONE
      
      real*8 pt ! temperature, deg C

      real*8 A1, A2, A3, A4, A5
      parameter (A1=3177.5, A2=-200.57d0, A3=6.8865d0, 
     &           A4=-0.13335d0, A5=0.0010877d0)

      SF6_schmidt_number = 
     &    A1 + A2*pt + A3*(pt**2) + A4*(pt**3) + A5*(pt**4)

      END     

      real*8 FUNCTION N2O_solubility(pt,sal)

C     N2O solubility from eq. 15 of Orr et al. (2016, GMDD). This is the solubility 
C     function 'phi0' defined by Orr et al. It includes the effects of fugacity and 
C     humidity and computes the saturation concentration of the gas in mol/L/atm at 
C     a total atmospheric pressure of 1 atm. The actual saturation concentration is 
C     given by: phi0*Pa*x, where Pa is the total atmospheric pressure in atm and x 
C     is the mixing ratio of the gas. Values below are from Table 2 of Orr et al 
C     (taken from Weiss and Price (1980)).
C
C     pt: temperature, deg C
C     sal: salinity, psu
C     output in mol/L/atm

      IMPLICIT NONE
      
      real*8 pt, sal

      real*8 A1, A2, A3, A4, B1, B2, B3
      parameter (A1=-165.8806d0, A2=222.8743d0, A3=92.0792d0,
     &           A4=-1.48425d0, B1=-0.056235d0, B2=0.031619d0,
     &           B3=-0.0048472d0)
C     A1=-168.2459d0, A2=226.0894d0, A3=93.2817d0, 
C    &           A4=-1.48693d0,B1=-0.060361d0, B2=0.033765d0, 
C    &           B3=-0.0051862d0)

      real*8 tk, tk100, lnphi0

      tk100 = (273.15d0 + pt)*0.01d0

      lnphi0 = A1 + A2/tk100 + A3*log(tk100) + A4*(tk100**2) + 
     &      sal*(B1 + B2*tk100 + B3*(tk100**2))

      N2O_solubility=exp(lnphi0)

      return
      
      END

      real*8 FUNCTION CFC11_solubility(pt,sal)

C     CFC11 solubility from eq. 15 of Orr et al. (2016, GMDD). This is the solubility 
C     function 'phi0' defined by Orr et al. It includes the effects of fugacity and 
C     humidity and computes the saturation concentration of the gas in mol/L/atm at 
C     a total atmospheric pressure of 1 atm. The actual saturation concentration is 
C     given by: phi0*Pa*x, where Pa is the total atmospheric pressure in atm and x 
C     is the mixing ratio of the gas. Values below are from Table 2 of Orr et al 
C     (taken from Warner and Weiss (1985)).
C
C     pt: temperature, deg C
C     sal: salinity, psu
C     output in mol/L/atm

      IMPLICIT NONE
      
      real*8 pt, sal

      real*8 A1, A2, A3, A4, B1, B2, B3
      parameter (A1=-229.9261d0, A2=319.6552d0, A3=119.4471d0,
     &           A4=-1.39165d0, B1=-0.142382d0, B2=0.091459d0,
     &           B3=-0.0157274d0)
C     parameter (A1=-232.04110d0, A2=322.5546d0, A3=120.4956d0, 
C    &           A4=-1.39165d0,B1=-0.146531d0, B2=0.093621d0, 
C    &           B3=-0.0160693d0)

      real*8 tk, tk100, lnphi0

      tk100 = (273.15d0 + pt)*0.01d0

      lnphi0 = A1 + A2/tk100 + A3*log(tk100) + A4*(tk100**2) + 
     &      sal*(B1 + B2*tk100 + B3*(tk100**2))

      CFC11_solubility=exp(lnphi0)

      return
      
      END

      real*8 FUNCTION CFC12_solubility(pt,sal)

C     CFC12 solubility from eq. 15 of Orr et al. (2016, GMDD). This is the solubility 
C     function 'phi0' defined by Orr et al. It includes the effects of fugacity and 
C     humidity and computes the saturation concentration of the gas in mol/L/atm at 
C     a total atmospheric pressure of 1 atm. The actual saturation concentration is 
C     given by: phi0*Pa*x, where Pa is the total atmospheric pressure in atm and x 
C     is the mixing ratio of the gas. Values below are from Table 2 of Orr et al 
C     (taken from Warner and Weiss (1985)).
C
C     pt: temperature, deg C
C     sal: salinity, psu
C     output in mol/L/atm

      IMPLICIT NONE
      
      real*8 pt, sal

      real*8 A1, A2, A3, A4, B1, B2, B3
      parameter (A1=-218.0971d0, A2=298.9702d0, A3=113.8049d0,
     &           A4=-1.39165d0, B1=-0.143566d0, B2=0.091015d0,
     &           B3=-0.0153924d0)
C     parameter (A1=-220.2120d0, A2=301.8695d0, A3=114.8533d0, 
C    &           A4=-1.39165d0,B1=-0.147718d0, B2=0.093175d0, 
C    &           B3=-0.0157340d0)

      real*8 tk, tk100, lnphi0

      tk100 = (273.15d0 + pt)*0.01d0

      lnphi0 = A1 + A2/tk100 + A3*log(tk100) + A4*(tk100**2) + 
     &      sal*(B1 + B2*tk100 + B3*(tk100**2))

      CFC12_solubility=exp(lnphi0)

      return
      
      END     

      real*8 FUNCTION SF6_solubility(pt,sal)

C     SF6 solubility from eq. 15 of Orr et al. (2016, GMDD). This is the solubility 
C     function 'phi0' defined by Orr et al. It includes the effects of fugacity and 
C     humidity and computes the saturation concentration of the gas in mol/L/atm at 
C     a total atmospheric pressure of 1 atm. The actual saturation concentration is 
C     given by: phi0*Pa*x, where Pa is the total atmospheric pressure in atm and x 
C     is the mixing ratio of the gas. Values below are from Table 2 of Orr et al 
C     (taken from Bullister et al. (2002)).
C
C     pt: temperature, deg C
C     sal: salinity, psu
C     output in mol/L/atm

      IMPLICIT NONE
      
      real*8 pt, sal

      real*8 A1, A2, A3, A4, B1, B2, B3
      parameter (A1=-80.0343d0, A2=117.232d0, A3=29.5817d0,
     &           A4=0.0d0, B1=0.0335183d0, B2=-0.0373942d0,
     &           B3=0.00774862d0)
 
      real*8 tk, tk100, lnphi0

      tk100 = (273.15d0 + pt)*0.01d0

      lnphi0 = A1 + A2/tk100 + A3*log(tk100) + A4*(tk100**2) + 
     &      sal*(B1 + B2*tk100 + B3*(tk100**2))

      SF6_solubility=exp(lnphi0)

      return
      
      END     
