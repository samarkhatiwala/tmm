
      SUBROUTINE MOBI_COPY_DATA(lSize, nx_loc, itr, arr, dir)
      
      implicit none

#include "size.h"

      integer lSize, nx_loc, itr, dir
      real arr(lSize)

      integer i, j, k, ii, n

#include "pconst.h"
#include "levind.h"
#include "mw.h"

      real src(imt,km,jsmw:jemw,nsrc)
      common/mobicomm/src

      j=1
      
	  n=itr+3  ! itr starts at 0 and first BGC tracer number is 3
      if (dir.eq.1) then
CSPK copy TMM tracer field to MOBI array
		ii=1
		do i=1,nx_loc
		  do k=1,kmt(i,j)
			t(i,k,j,n,taum1) = arr(ii)
			ii=ii+1
		  enddo
		enddo
      elseif (dir.eq.2) then
CSPK copy MOBI tracer tendency field to TMM array
		arr(:) = c0
		ii=1
		do i=1,nx_loc
		  do k=1,kmt(i,j)
		    arr(ii) = src(i,k,j,itrc(n))
			ii=ii+1
		  enddo
		enddo
      else
        stop 'ERROR in MOBI_COPY_DATA: unknown copy direction!'
      endif

      END
      