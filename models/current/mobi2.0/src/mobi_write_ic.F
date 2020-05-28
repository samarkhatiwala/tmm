      program tmmmobiwrite

      implicit none

#include "size.h"
#include "mw.h"
      integer nx
      parameter (nx=1)
      integer nz_loc_arr(nx)
      real tlat_loc(nx), area_loc(nx), sg_bathy_loc(km)
      real S_surf_glob, TR_surf_glob(nt-2)
      real DeltaT
      real zt(km), drF(km)
      integer k, dum, iUnit, n
      character(120) :: fname
      real co2atm, dc13atm, dc14atm      
#ifdef O_sed  
      integer numOceanStepsPerSedStep
      integer nzmaxSed, ibmaxSed, numSedMixedTracers
      integer numSedBuriedTracers
      real globalweathflx
      real sedsa_loc, sedmask_loc(nx)
#endif  
      
      tlat_loc(:)=0.0
      sg_bathy_loc(:)=1.0      
      nz_loc_arr(:)=km
      TR_surf_glob(:)=1.0
      DeltaT=28800.0

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
        read(iUnit)zt(k)
      enddo  
      close(iUnit)      
      
      call mobi_ini(nt-2, nx*km, nx, km, nz_loc_arr,
     &              zt, drF, DeltaT, tlat_loc, area_loc,
     &              sg_bathy_loc,
     &              S_surf_glob, TR_surf_glob,
     &              co2atm, dc13atm, dc14atm,
#ifdef O_sed  
     &              numOceanStepsPerSedStep,
     &              nzmaxSed, ibmaxSed, numSedMixedTracers, 
     &              numSedBuriedTracers, globalweathflx,
     &              sedsa_loc, sedmask_loc,
#endif       
     &              1)
     
      do n=3,nt
        call findunit(iUnit)
        fname=trim(mapt(n)) // '.dat'
        open(UNIT=iUnit, FILE=fname, STATUS="UNKNOWN", ACCESS="STREAM",
     &       FORM="FORMATTED")
        do k=1,km
C          We convert zt back to meter        
           write(iUnit,'(F16.3,E22.14)')zt(k)/100,t(1,k,1,n,taup1)
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

      end
      