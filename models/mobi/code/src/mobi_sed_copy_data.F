#include "CPP_OPTIONS_CHECK.h"

      SUBROUTINE MOBI_SED_COPY_DATA(nx_loc, sedmixed, sedburied, dir)
      
      implicit none

#include "size.h"
#include "sed.h"

      integer nx_loc, dir
      real sedmixed(nzmax,20,nx_loc)
      real sedburied(ibmax,2,nx_loc)

      integer is, ix

      if (dir.eq.1) then
!      write(*,*)'SPK1: copying from TMM to SED'
		do is=1,ipsed
		  ix = imap(is)

C         copy TMM tracer field to SED tracer array        
		  carb(:,1,is) = sedmixed(:,1,ix)
		  carb(:,2,is) = sedmixed(:,2,ix)
		  carb(:,3,is) = sedmixed(:,3,ix)
		  dcpls(:,1,is) = sedmixed(:,4,ix)
		  dcpls(:,2,is) = sedmixed(:,5,ix)
		  dcpls(:,3,is) = sedmixed(:,6,ix)
		  dcmin(:,1,is) = sedmixed(:,7,ix)
		  dcmin(:,2,is) = sedmixed(:,8,ix)
		  dcmin(:,3,is) = sedmixed(:,9,ix)
C
		  pore(:,is) = sedmixed(:,10,ix)
		  form(:,is) = sedmixed(:,11,ix)
		  o2(:,is) = sedmixed(:,12,ix)
		  orggg(:,is) = sedmixed(:,13,ix)
		  orgml(:,is) = sedmixed(:,14,ix)
		  calgg(:,is) = sedmixed(:,15,ix)
		  calml(:,is) = sedmixed(:,16,ix)
		  dopls(:,is) = sedmixed(:,17,ix)
		  domin(:,is) = sedmixed(:,18,ix)
		  dbpls(:,is) = sedmixed(:,19,ix)
		  dbmin(:,is) = sedmixed(:,20,ix)
C
C         write(*,*)'SPK-copy:',nzmax,ibmax,ipmax
		  buried_mass(:,is) = sedburied(:,1,ix)
		  buried_calfrac(:,is) = sedburied(:,2,ix)
		enddo
!       write(*,*)'SPK2: finished copying from TMM to SED'		
      elseif (dir.eq.2) then
C       copy SED tracer field to TMM tracer array
!      write(*,*)'SPK3: copying from SED to TMM'
		do is=1,ipsed
		  ix = imap(is)

		  sedmixed(:,1,ix) = carb(:,1,is)
		  sedmixed(:,2,ix) = carb(:,2,is)
		  sedmixed(:,3,ix) = carb(:,3,is)
		  sedmixed(:,4,ix) = dcpls(:,1,is)
		  sedmixed(:,5,ix) = dcpls(:,2,is)
		  sedmixed(:,6,ix) = dcpls(:,3,is)
		  sedmixed(:,7,ix) = dcmin(:,1,is)
		  sedmixed(:,8,ix) = dcmin(:,2,is)
		  sedmixed(:,9,ix) = dcmin(:,3,is)
C
		  sedmixed(:,10,ix) = pore(:,is)
		  sedmixed(:,11,ix) = form(:,is)
		  sedmixed(:,12,ix) = o2(:,is)
		  sedmixed(:,13,ix) = orggg(:,is)
		  sedmixed(:,14,ix) = orgml(:,is) 
		  sedmixed(:,15,ix) = calgg(:,is) 
		  sedmixed(:,16,ix) = calml(:,is) 
		  sedmixed(:,17,ix) = dopls(:,is) 
		  sedmixed(:,18,ix) = domin(:,is) 
		  sedmixed(:,19,ix) = dbpls(:,is) 
		  sedmixed(:,20,ix) = dbmin(:,is) 
C
		  sedburied(:,1,ix) = buried_mass(:,is) 
		  sedburied(:,2,ix) = buried_calfrac(:,is)
        enddo     
      else
        stop 'ERROR in MOBI_SED_COPY_DATA: unknown copy direction!'
      endif

      END
      