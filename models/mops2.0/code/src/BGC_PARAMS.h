C$Header: /Users/ikriest/CVS/mops/BGC_PARAMS.h,v 1.2 2016/06/03 09:28:59 ikriest Exp $
C$Name: mops-2_0 $

C EVERYTHING RELATED TO THE MAIN BGC TRACER ARRAY, ITS INDICES AND TYPES

c a dummy for the nominal number of vertical layers
      INTEGER bgc_ktotal

      PARAMETER(bgc_ktotal=100)

c the total number of bgc tracers
      INTEGER bgc_ntracer
      
c the indices of tracers
      INTEGER ipo4,idop,iphy,izoo,ioxy,idet,idin

      PARAMETER(ipo4=1,    !PO4
     &          idop=2,    !DOP
     &          ioxy=3,    !Oxygen
     &          iphy=4,    !Phyto-P
     &          izoo=5,    !Zoo-P
     &          idet=6,    !Detritus-P
     &          idin=7)    !DIN

#ifndef CARBON
      PARAMETER(bgc_ntracer=7)
#else
      PARAMETER(bgc_ntracer=9)
      INTEGER idic,ialk
      PARAMETER(idic=8,ialk=9)
c connect between carbon exchange and P-based BGC
      real*8 ocmip_alkfac,ocmip_silfac
      COMMON/CO2SURFACE/ocmip_alkfac,ocmip_silfac
#endif

c the tracer field
      REAL*8 bgc_tracer

      COMMON/BGC/bgc_tracer(bgc_ktotal,bgc_ntracer)

C EVERYTHING RELATED TO THE BIOGEOCHEMISTRY EVALUATION

c the flux attenuation curve
      real*8 wdet(bgc_ktotal)

c the biogeochemistry constants
      real*8 rcp,rnp,ro2ut,subox,tempB,
     &       acmuphy,acik,ackpo4,ackchl,ackw,aclambda,acomni,
     &       ACMuzoo,ACkphy,AClambdaz,AComniz,ACeff,
     &       graztodop,dlambda,plambda,zlambda,detlambda,
     &       detmartin,detwa,detwb,detwmin,
     &       vsafe,alimit

c the air-sea gas exchange constants
      real*8 sox1,sox2,sox3,sox4,oA0,oA1,oA2,oA3,oA4,oA5,
     &       oB0,oB1,oB2,oB3,oC0

      COMMON/BGCZ/wdet
      COMMON/BGCPARAMS/rcp,rnp,ro2ut,subox,tempB,
     &       acmuphy,acik,ackpo4,ackchl,ackw,aclambda,acomni,
     &       ACMuzoo,ACkphy,AClambdaz,AComniz,ACeff,
     &       graztodop,dlambda,plambda,zlambda,detlambda,
     &       detmartin,detwa,detwb,detwmin,
     &       vsafe,alimit

      COMMON/O2SURFACE/sox1,sox2,sox3,sox4,oA0,oA1,oA2,oA3,oA4,oA5,
     &       oB0,oB1,oB2,oB3,oC0

c sediment burial and O2 sensitivity of OM degradation
      real*8 burdige_fac,burdige_exp,flux_bury,ACkbaco2
      COMMON/BGCSEDPARAMS/burdige_fac,burdige_exp,flux_bury,ACkbaco2

c parameters related to N-Fixation and denitrification
      real*8 tf2,tf1,tf0,tff,nfix,subdin,rhno3ut,ACkbacdin
      COMMON/BGCNPARAMS/tf2,tf1,tf0,tff,nfix,subdin,rhno3ut,ACkbacdin
