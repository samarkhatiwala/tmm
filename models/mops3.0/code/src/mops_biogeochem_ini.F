C$Header: /Users/ikriest/CVS/mops/mops_biogeochem_ini.F,v 1.2 2016/06/03 09:28:59 ikriest Exp $
C$Name: mops-2_0 $

CBOP

C !INTERFACE: ==========================================================
      SUBROUTINE MOPS_BIOGEOCHEM_INI(Nrloc,DeltaT,
#ifdef CARBON
     &                          localph,
#endif      
     &                          thetaloc,saltloc,dzloc,
     &                          drFloc,nzmax,nzeuph,
     &                          numBiogeochemStepsPerOceanStep,
     &                          setDefaults)

C !DESCRIPTION:

C !USES: ===============================================================
      IMPLICIT NONE

#include "BGC_PARAMS.h"
#include "BGC_CONTROL.h"

#ifdef CARBON
#include "CAR_PARAMS.h"
#endif

C !INPUT PARAMETERS: ===================================================
      INTEGER Nrloc
      real*8 DeltaT
#ifdef CARBON
      real*8 localph
#endif      
      real*8 thetaloc(Nrloc),saltloc(Nrloc)
      real*8 dzloc(Nrloc),drFloc(bgc_ktotal)
      integer numBiogeochemStepsPerOceanStep,nzmax,nzeuph
      logical setDefaults

C !LOCAL VARIABLES: ====================================================
      INTEGER k
      real*8 bgc_zu(bgc_ktotal+1),bgc_dz(bgc_ktotal)

#ifdef CARBON
      real*8 sdic,spho,salk,ssil
#endif

C Biogeochemical time step length in days
      bgc_timesteps = numBiogeochemStepsPerOceanStep
      bgc_dt=DeltaT/(86400.0d0*dble(numBiogeochemStepsPerOceanStep))

C Depth of layers for biogeochemistry

      bgc_kmax = nzmax
      bgc_keuph = nzeuph

      do k=1,bgc_ktotal
        bgc_zu(k)=0.d0
      enddo
      	
      do k=2,bgc_kmax+1
        bgc_zu(k)=bgc_zu(k-1)+drFloc(k-1)
      enddo        

C Initialize the tracers
      DO k=1,Nrloc
        bgc_dz(k)=dzloc(k)
      ENDDO
  
      DO k=Nrloc+1,bgc_ktotal
        bgc_dz(k)=0.d0
      ENDDO

      CALL BGC_INI(bgc_zu,setDefaults)

#ifdef CARBON

! forc OCMIP constant alkalinity and silicate

      ocmip_alkfac = 2310.0d0*1.0245d0/34.88d0
      ocmip_silfac = 7.7d0 
      
! Global mean silicate for surface layer

      ssil=ocmip_silfac

! Surface total alkalinity follows the OCMIP protocol:
! ta = 2310*s/sbar where sbar is the model's annual mean surface salinity      

      salk = ocmip_alkfac*saltloc(1)
              
! Initialize carbonate system.
! The constants and initial coefficients will be written into a common block.
      sdic = bgc_tracer(1,idic)
      spho = bgc_tracer(1,ipo4)

      CALL CAR_INI(thetaloc,saltloc,Nrloc,
     &             sdic,spho,salk,ssil)
      
      localph=sph
      	
#endif	

     
      RETURN
      END

