C$Header: /Users/ikriest/CVS/mops/insolation.F,v 1.2 2016/06/03 09:28:59 ikriest Exp $
C$Name: mops-2_0 $

CBOP
C !ROUTINE: INSOL

C !INTERFACE: ==========================================================
      SUBROUTINE insolation(N,Time,Yc,daysperyear,sfac,stau)

C !DESCRIPTION:
C find shortwave radiation as function of date and latitude
c based on paltridge and parson
c modified by SPK from MITGCM

C !USES: ===============================================================
      IMPLICIT NONE
C     === Global variables ===

C !INPUT PARAMETERS: ===================================================
C time                 :: current time
       integer N
       real*8 time,Yc(N)
       real*8 daysperyear

C !OUPUT PARAMETERS: ===================================================
       real*8 sfac(N),stau(N)

C !LOCAL VARIABLES: ====================================================
       real*8  solar, albedo, par
       real*8  dayfrac, yday, delta
       real*8  lat, sun1, dayhrs
       real*8  cosz, frac, fluxi, fracmin
       integer j
CEOP

c
      solar = 1360.   !solar constant
      albedo = 0.6    !planetary albedo
C      par = 0.4       !photosynthetically reactive frac

c
c find day (****NOTE for year starting in winter*****)
        dayfrac=mod(Time,daysperyear*86400.d0)
     &                    /(daysperyear*86400.d 0)          !fraction of year
        yday = 2.0*3.1416*dayfrac                         !convert to radians
        delta = (0.006918 - (0.399912*cos(yday))      !cosine zenith angle
     &          +(0.070257*sin(yday))                 !(paltridge+platt)
     &          -(0.006758*cos(2.0*yday))
     &          +(0.000907*sin(2.0*yday))
     &          -(0.002697*cos(3.0*yday))
     &          +(0.001480*sin(3.0*yday)) )
       do j=1,N
c latitude in radians
          lat=YC(j)/180.d0*3.1416
          sun1 = -sin(delta)/cos(delta) * sin(lat)/cos(lat)
          if (sun1.le.-0.999) sun1=-0.999
          if (sun1.ge. 0.999) sun1= 0.999
          dayhrs = abs(acos(sun1))
          cosz = ( sin(delta)*sin(lat)+              !average zenith angle
     &            (cos(delta)*cos(lat)*sin(dayhrs)/dayhrs) )
          if (cosz.le.0.005) cosz=0.005
          frac = dayhrs/3.1416               !fraction of daylight in day
c daily average photosynthetically active solar radiation just below surface
          fluxi = solar*(1.0-albedo)*cosz*frac
c
c convert to sfac
          if (fluxi.gt.0.0) sfac(j)=fluxi
c very large for polar night
          if (fluxi.lt.0.00001) sfac(j)=0.00001
c daylength; ensure that it lies between 0 and 1 (may be slightly
c out of this range in high latitudes)
          fracmin = MIN(frac,1.0)
          stau(j) = MAX(fracmin,0.0)

       enddo !j
      return
      end
c
C==========================================================================
