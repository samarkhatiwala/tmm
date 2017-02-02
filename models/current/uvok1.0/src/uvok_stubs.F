#ifdef O_TMM
      subroutine areaavg (data, dmsk, avg)

!=======================================================================
C stub routine
!=======================================================================

      end
      
      subroutine setbcx (a, imt, jmtorkm)

!=======================================================================
C stub routine
!=======================================================================

      dimension a(imt,jmtorkm)

      end

      subroutine data (is, ie, js, je)

!=======================================================================
C stub routine
!=======================================================================

      integer is, ie, js, je

      end

      subroutine co2forc

!=======================================================================
C stub routine
!=======================================================================

      end

      subroutine c14data

!=======================================================================
C stub routine
!=======================================================================

      end

      subroutine co2ccndata

!=======================================================================
C stub routine
!=======================================================================

      end

      subroutine defvar (name, ncid, nd, id, rmin, rmax, axis                                                                        &
     &,                  type, lname, sname, units)
!=======================================================================
!     define data

!     input:
!       name  = name of variable to be defined
!       ncid  = unit
!       nd    = number dimensions of data
!       id    = data id
!       rmin  = minimum range (default real)
!       rmax  = maximum range (default real)
!       axis  = axis type
!       type  = data type (D=double,F=float,I=integer,Tn=char*n)
!       lname = long name
!       sname = standard name
!       units = data units
!=======================================================================

      use uvok_diags_mod

      implicit none

      character(*), intent(in) :: name, axis, lname, sname, type, units

      integer :: nd, id(nd), ncid

      real, intent(in) :: rmax, rmin

      id(1) = id(1) + 1 ! increment diagnostic index

	  if (ncid.eq.1) then
		diag2dFileNames(id(1)) = name // '.bin'
	  elseif (ncid.eq.2) then
		diag3dFileNames(id(1)) = name // '.petsc'   
	  endif       

      if (writeFlag.eq.1) then
	   if (ncid.eq.1) then
		write(*,*)'2-d diagnostic #',id(1),' is: ',name,' ',lname,' ',
     &     units
        write(*,*)'  will be written to: ',diag2dFileNames(id(1))
        write(diagsLogFileUnit,'(A,A,A,A,A,A,A)')trim(name),",",
     &   trim(lname),",",trim(units),",",
     &   trim(diag2dFileNames(id(1)))
	   elseif (ncid.eq.2) then
		write(*,*)'3-d diagnostic #',id(1),' is: ',name,' ',lname,' ',
     &     units
        write(*,*)'  will be written to: ',diag3dFileNames(id(1))
        write(diagsLogFileUnit,'(A,A,A,A,A,A,A)')trim(name),",",
     &   trim(lname),",",trim(units),",",
     &   trim(diag3dFileNames(id(1)))
	   endif       
      endif
      
      return
      end

      subroutine putvaramsk (name, ncid, ln, is, ic, din, dm, s, o)
!=======================================================================
!     write data

!     input:
!       name = name of variable to be written
!       ncid = iou unit
!       ln   = length of data
!       is   = starting point for write in each dimension
!       ic   = count (or length) for write in each dimension
!       din  = data to be written (default real)
!       dm   = data mask
!       s    = data scalar
!       o    = data offset
!=======================================================================

      use uvok_diags_mod
      
      implicit none

      character(*), intent(in) :: name

      integer, intent(in) :: ic(10), is(10), ln, ncid

      real, intent(inout) :: din(ln)
      real, intent(in) :: dm(ln), o, s
      real rs

      integer idn, kstart, kend, imode

      idn=is(1) ! diagnostic index
      imode=ic(1) ! mode (1=accumulate, 2=average)
      
      if (idn.le.0) then
        return
      endif
      
	  kstart=is(2)
	  kend=is(3)

      if (imode.eq.1) then ! accumulate
		if (ncid.eq.1) then ! 2-d diagnostic
          diags2d(kstart:kend,idn)=diags2d(kstart:kend,idn)+din(1:ln)
		elseif (ncid.eq.2) then ! 3-d diagnostic
		  diags3d(kstart:kend,idn)=diags3d(kstart:kend,idn)+din(1:ln)
		endif
      elseif (imode.eq.2) then ! time average (dm is the averaging factor) and scale
		if (ncid.eq.1) then ! 2-d diagnostic
          diags2d(kstart:kend,idn)=diags2d(kstart:kend,idn)*dm(1:ln)
		elseif (ncid.eq.2) then ! 3-d diagnostic
		  diags3d(kstart:kend,idn)=diags3d(kstart:kend,idn)*dm(1:ln)
		endif
		rs = 0.0
		if (s .ne. 0.) rs = 1.0/s
C	do i=1,ln
		if (ncid.eq.1) then ! 2-d diagnostic
		  diags2d(kstart:kend,idn)=(diags2d(kstart:kend,idn)-o)*rs
		elseif (ncid.eq.2) then ! 3-d diagnostic       
		  diags3d(kstart:kend,idn)=(diags3d(kstart:kend,idn)-o)*rs
		endif  
C         din(i) = (din(i) - o)*rs
C	enddo
      endif
      
      return
      end

#endif
