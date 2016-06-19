#include "DIC_OPTIONS.h"

      SUBROUTINE MITGCHEM_COPY_DATA(nrloc, itr, ptr, jtr, dt, dir)
      
      implicit none

#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "PTRACERS_SIZE.h"
#include "PTRACERS_PARAMS.h"
#include "PTRACERS_FIELDS.h"
#include "DIC_VARS.h"

      integer nrloc, itr, dir
      real*8 ptr(nrloc), jtr(nrloc), dt

      integer k

C     Note: itr starts at 0
      if (dir.eq.1) then
C       copy TMM tracer field to BGC tracer array
        do k=1,nrloc
		  pTracer(1,1,k,1,1,itr+1) = ptr(k)
		enddo
		do k=nrloc+1,Nr
		  pTracer(1,1,k,1,1,itr+1) = 0.d0
		enddo
      elseif (dir.eq.2) then
C       copy BGC tracer tendency field to TMM tendency array
        do k=1,nrloc
          jtr(k)=(pTracer(1,1,k,1,1,itr+1)-ptr(k))/PTRACERS_dTLev(k)
        enddo  
      elseif (dir.eq.3) then
C     copy BGC tracer field to TMM tracer array
        do k=1,nrloc
		  ptr(k)=pTracer(1,1,k,1,1,itr+1)
		enddo  
      else
        stop 'ERROR in MITGCHEM_COPY_DATA: unknown copy option!'
      endif

      END
      