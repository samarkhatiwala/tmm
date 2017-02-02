#ifdef O_TMM
      SUBROUTINE UVOK_COPY_DATA(kmt_loc, itr, arr, dir)
      
      implicit none

#include "size.h"

      integer kmt_loc, itr, dir
      real arr(kmt_loc)

      integer n

#include "pconst.h"
#include "mw.h"

      real src(imt,km,jsmw:jemw,nsrc)
      common/uvokcomm/src

C      print*,'copy:',kmt_loc,itr,arr(:)
	  n=itr+3  ! itr starts at 0 and first BGC tracer number is 3
      if (dir.eq.1) then
CSPK copy TMM tracer field to UVOK array
		t(1,1:kmt_loc,1,n,taum1) = arr(1:kmt_loc)
      elseif (dir.eq.2) then
CSPK copy UVOK tracer tendency field to TMM array
		arr(:) = c0
		arr(1:kmt_loc) = src(1,1:kmt_loc,1,itrc(n))
      else
        stop 'ERROR in UVOK_COPY_DATA: unknown copy direction!'
      endif

      END
#endif      