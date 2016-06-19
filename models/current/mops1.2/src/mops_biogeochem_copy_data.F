C$Header: /Users/ikriest/CVS/mops/mops_biogeochem_copy_data.F,v 1.1.1.1 2015/06/03 17:02:09 ikriest Exp $
C$Name: mops-1_2 $

      SUBROUTINE MOPS_BIOGEOCHEM_COPY_DATA(nrloc,itr,ptr,jtr,deltaT,dir) 
      
      implicit none

#include "BGC_PARAMS.h"

      integer nrloc, itr, dir
      real*8 ptr(nrloc), jtr(nrloc), deltaT

      integer k

C     Note: itr starts at 0
      if (dir.eq.1) then
C       copy TMM tracer field to BGC tracer array
        do k=1,nrloc
		  bgc_tracer(k,itr+1) = ptr(k)
		enddo  
      elseif (dir.eq.2) then
C       copy BGC tracer tendency field to TMM tendency array
        do k=1,nrloc
          jtr(k)=(bgc_tracer(k,itr+1)-ptr(k))/DeltaT
        enddo  
      elseif (dir.eq.3) then
C     copy BGC tracer field to TMM tracer array
        do k=1,nrloc
		  ptr(k)=bgc_tracer(k,itr+1)
		enddo  
      else
        stop 'ERROR in MOPS_BIOGEOCHEM_COPY_DATA: unknown copy option!'
      endif

      END
