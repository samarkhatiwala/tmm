C EVERYTHING RELATED TO THE MAIN BGC TRACER ARRAY, ITS INDICES AND TYPES

c a dummy for the nominal number of vertical layers
      INTEGER bgc_ktotal

      PARAMETER(bgc_ktotal=100)

c the total number of bgc tracers
      INTEGER bgc_ntracer
      
c the indices of tracers
      INTEGER ipo4,idop,iphy,izoo,ioxy,idet

      PARAMETER(ipo4=1,    !PO4
     &          idop=2,    !DOP
     &          ioxy=3,    !Oxygen
     &          iphy=4,    !Phyto-P
     &          izoo=5,    !Zoo-P
     &          idet=6)    !Detritus-P

#ifndef CARBON
      PARAMETER(bgc_ntracer=6)
#else
      PARAMETER(bgc_ntracer=8)
      INTEGER idic,ialk
      PARAMETER(idic=7,ialk=8)
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
     &       detwa,detwb,detwmin,
     &       vsafe,alimit

c the air-sea gas exchange constants
      real*8 sox1,sox2,sox3,sox4,oA0,oA1,oA2,oA3,oA4,oA5,
     &       oB0,oB1,oB2,oB3,oC0

      COMMON/BGCZ/wdet
      COMMON/BGCPARAMS/rcp,rnp,ro2ut,subox,tempB,
     &       acmuphy,acik,ackpo4,ackchl,ackw,aclambda,acomni,
     &       ACMuzoo,ACkphy,AClambdaz,AComniz,ACeff,
     &       graztodop,dlambda,plambda,zlambda,detlambda,
     &       detwa,detwb,detwmin,
     &       vsafe,alimit

      COMMON/O2SURFACE/sox1,sox2,sox3,sox4,oA0,oA1,oA2,oA3,oA4,oA5,
     &       oB0,oB1,oB2,oB3,oC0

      
