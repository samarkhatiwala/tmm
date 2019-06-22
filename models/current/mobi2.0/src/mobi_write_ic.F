      program tmmmobiwrite

      implicit none

#include "size.h"
#include "mw.h"
      integer nx
      parameter (nx=1)
      integer nz_loc_arr(nx)
      real tlat_loc(nx), sg_bathy_loc(km)
      real S_surf_glob, TR_surf_glob(nt-2)
      real DeltaT
      real zt(km), drF(km)
      integer k, dum, iUnit, n
      real buf(km,2)
      character(120) :: fname
      
      tlat_loc(:)=0.0
      sg_bathy_loc(:)=1.0      
      nz_loc_arr(:)=km
      TR_surf_glob(:)=1.0
      DeltaT=28800.0

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
     &              zt, drF, DeltaT, tlat_loc,
     &              sg_bathy_loc,
     &              S_surf_glob, TR_surf_glob, 
     &              1)
     
C     do k=1,km
C       write(*,*)k,drF(k),zt(k)
C     enddo  
      
C     buf(:,1)=zt(:)
      do n=3,nt
        call findunit(iUnit)
        fname=trim(mapt(n)) // '.dat'
        open(UNIT=iUnit, FILE=fname, STATUS="UNKNOWN", ACCESS="STREAM",
     &       FORM="FORMATTED")
        do k=1,km
C          buf(k,2)=t(1,k,1,n,taup1)
           write(iUnit,'(F16.3,E22.14)')zt(k)/100,t(1,k,1,n,taup1)
        enddo
        close(iUnit)
C       call write_r8_field(km,2,1,buf,trim(mapt(n))//'.dat')      
      enddo

	  call findunit(iUnit)
	  open(UNIT=iUnit, FILE="MOBI_tracer_names.txt", STATUS="UNKNOWN",
     &     ACCESS="STREAM", FORM="FORMATTED")
      do n=3,nt
		 write(iUnit,'(A)')trim(mapt(n))
	  enddo
	  close(iUnit)

      end
      