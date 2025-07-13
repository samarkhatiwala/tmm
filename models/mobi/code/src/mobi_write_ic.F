#include "CPP_OPTIONS_CHECK.h"

      program tmmmobiwrite

      implicit none

#include "size.h"
#include "mw.h"
#if defined O_mobi
#include "mobi.h"
#endif
#if defined O_PaTh
#include "protac_thor.h"
#endif
#include "coord.h"
C#include "grdvar.h"
C#include "stdunits.h"

      integer nx
      parameter (nx=1)
      integer nz_loc_arr(nx)
      real tlat_loc(nx), area_loc(nx)
      real S_surf_glob, TR_surf_glob(nt-2)
      real DeltaT
      real drF(km), z(km)
      integer k, dum, iUnit, n
      character(120) :: fname
#if defined O_carbon
      real co2atm
# if defined O_carbon_13
      real dc13atm
# endif     
# if defined O_carbon_14     
      real dc14atm
# endif
#endif ! O_carbon
#if defined O_mobi
      real sg_bathy_loc(km)
# if defined O_mobi_iron
      real fe_hydr_loc(km)
# endif
# if defined O_mobi_silicon
      real si_hydr_loc(km)
# endif     
#endif
#if defined O_PaTh 
      real lith_loc(km), wLith_loc(km)
# if !defined O_mobi
      real pom_loc(km), caco3_loc(km), opal_loc(km)
      real wPOM_loc(km), wCaCO3_loc(km), wOpal_loc(km)
# endif     
#endif
#if defined O_sed  
      integer numOceanStepsPerSedStep
      integer nzmaxSed, ibmaxSed, numSedMixedTracers
      integer numSedBuriedTracers
      real globalweathflx
      real sedsa_loc, sedmask_loc(nx)
#endif

      
      tlat_loc(:)=0.0
      nz_loc_arr(:)=km
      TR_surf_glob(:)=1.0
      DeltaT=28800.0

#if defined O_mobi
      sg_bathy_loc(:)=1.0      
#endif

C     UVic/MOBI expects drF and zt to be in centimeter
      call findunit(iUnit)
      open(UNIT=iUnit, FILE="drF.bin", STATUS="OLD", ACCESS="STREAM")
      read(iUnit)dum
      do k=1,km
        read(iUnit)drF(k)
      enddo  
      close(iUnit)      

      call findunit(iUnit)
      open(UNIT=iUnit, FILE="zt.bin", STATUS="OLD", ACCESS="STREAM")
      read(iUnit)dum
      do k=1,km
        read(iUnit)z(k)
      enddo  
      close(iUnit)      
      
      call mobi_ini(nt-2, nx*km, nx, km, nz_loc_arr,
     &              z, drF, DeltaT, tlat_loc, area_loc,
     &              S_surf_glob, TR_surf_glob,
#if defined O_carbon
     &                    co2atm, 
# if defined O_carbon_13
     &                    dc13atm, 
# endif     
# if defined O_carbon_14     
     &                    dc14atm,      
# endif
#endif ! O_carbon
#if defined O_mobi
     &                    sg_bathy_loc,
# if defined O_mobi_iron
     &                    fe_hydr_loc, 
# endif
# if defined O_mobi_silicon
     &                    si_hydr_loc, 
# endif     
#endif
#if defined O_PaTh 
     &                    lith_loc, wLith_loc, 
# if !defined O_mobi
     &                    pom_loc, caco3_loc, opal_loc, 
     &                    wPOM_loc, wCaCO3_loc, wOpal_loc, 
# endif     
#endif
#if defined O_sed
     &                    numOceanStepsPerSedStep,
     &                    nzmaxSed, ibmaxSed, numSedMixedTracers, 
     &                    numSedBuriedTracers, globalweathflx,
     &                    sedsa_loc, sedmask_loc,
#endif ! O_sed   
     &              1)
     
      do n=3,nt
        call findunit(iUnit)
        fname=trim(mapt(n)) // '.dat'
        open(UNIT=iUnit, FILE=fname, STATUS="UNKNOWN", ACCESS="STREAM",
     &       FORM="FORMATTED")
        do k=1,km
C          We convert zt back to meter        
           write(iUnit,'(F16.3,E22.14)')zt(k)/100.0,t(1,k,1,n,taup1)
        enddo
        close(iUnit)
      enddo

	  call findunit(iUnit)
	  open(UNIT=iUnit, FILE="MOBI_tracer_names.txt", STATUS="UNKNOWN",
     &     ACCESS="STREAM", FORM="FORMATTED")
      do n=3,nt
		 write(iUnit,'(A)')trim(mapt(n))
	  enddo
	  close(iUnit)

#if defined O_PaTh && defined O_mobi
      do k=1,km
        wPOM(k)=(wd(k)*1.e-2)*dzt(k) ! sinking speed in m/s
        wCaCO3(k)=(wc(k)*1.e-2)*dzt(k) ! sinking speed in m/s
        wOpal(k)=(wo(k)*1.e-2)*dzt(k) ! sinking speed in m/s
      enddo  

	  write(*,*)'Sinking speed profile [m/s]:'
	  write(*,*)'k  wPOM   wCaCO3 wOpal'
      do k=1,km
        write(*,*)'k=',k,wPOM(k),wCaCO3(k),wOpal(k)
      enddo  

      call write_r8_field(km,1,1,wPOM,'wPOM.bin')
      call write_r8_field(km,1,1,wCaCO3,'wCaCO3.bin')
      call write_r8_field(km,1,1,wOpal,'wOpal.bin')
      
#endif

      end
      