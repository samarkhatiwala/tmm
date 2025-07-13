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

      use mobi_diags_mod

      implicit none

      character(*), intent(in) :: name, axis, lname, sname, type, units

      integer :: nd, id(nd), ncid

      real, intent(in) :: rmax, rmin

      id(1) = id(1) + 1 ! increment diagnostic index

	  if (ncid.eq.1) then
		diag2dNames(id(1)) = name
		diag2dFileNames(id(1)) = name // '.bin'
	  elseif (ncid.eq.2) then
		diag3dNames(id(1)) = name
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
!       ncid(5) = (1) diag mode (1=accumulate, 2=average)
!                 (2) diag type (1=2d, 2=3d)
!                 (3)-(5) dimensions of input data (din and dm)
!       ln   = not used
!       is   = not used
!       ic   = not used
!       din  = data to be written (default real)
!       dm   = data mask
!       s    = data scalar
!       o    = data offset
!=======================================================================

      use mobi_diags_mod
      
      implicit none

      character(*), intent(in) :: name

      integer, intent(in) :: ic(10), is(10), ln, ncid(5)

      real, intent(inout) :: din(ncid(3),ncid(4),ncid(5))
      real, intent(in) :: dm(ncid(3),ncid(4),ncid(5)), o, s
      real rs

      integer idn, i, j, k, ii, n, diagMode, diagType

      diagMode=ncid(1) ! 1=accumulate, 2=average
      diagType=ncid(2) ! 1=2d, 2=3d
      
C     Find the diagnostic index corresponding to name
      if (diagType.eq.1) then ! 2-d diagnostic
		idn=-1
		n=1
		do while ((idn .lt. 0) .and. (n.le.num2dDiags))
		  if (trim(diag2dNames(n)).eq.trim(name)) idn=n
		  n=n+1
		enddo
   	  elseif (diagType.eq.2) then ! 3-d diagnostic
		idn=-1
		n=1
		do while ((idn .lt. 0) .and. (n.le.num3dDiags))
		  if (trim(diag3dNames(n)).eq.trim(name)) idn=n
		  n=n+1		  
		enddo      
      endif
      
      if (idn.le.0) then
        return
      endif

      j=1             
      if (diagMode.eq.1) then ! accumulate
		if (diagType.eq.1) then ! 2-d diagnostic
          diags2d(:,idn)=diags2d(:,idn)+din(1:numProfiles,j,1)
		elseif (diagType.eq.2) then ! 3-d diagnostic
		  ii=1
		  do i=1,numProfiles
		    do k=1,kmtdiags(i,j)
		      diags3d(ii,idn)=diags3d(ii,idn)+din(i,j,k)
		      ii=ii+1
		    enddo
		  enddo    
		endif
      elseif (diagMode.eq.2) then ! time average (dm is the averaging factor) and scale
		if (diagType.eq.1) then ! 2-d diagnostic
          diags2d(:,idn)=diags2d(:,idn)*dm(1:numProfiles,j,1)
		elseif (diagType.eq.2) then ! 3-d diagnostic
		  ii=1
		  do i=1,numProfiles
		    do k=1,kmtdiags(i,j)		
		      diags3d(ii,idn)=diags3d(ii,idn)*dm(i,j,k)
		      ii=ii+1
		    enddo
		  enddo    
		endif
		rs = 0.0
		if (s .ne. 0.) rs = 1.0/s
C	do i=1,ln
		if (diagType.eq.1) then ! 2-d diagnostic
		  diags2d(:,idn)=(diags2d(:,idn)-o)*rs
		elseif (diagType.eq.2) then ! 3-d diagnostic       
		  diags3d(:,idn)=(diags3d(:,idn)-o)*rs
		endif  
C         din(i) = (din(i) - o)*rs
C	enddo
      endif
      
      return
      end

