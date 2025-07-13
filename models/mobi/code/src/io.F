C $Header: /Users/spk/CVSREP/ICM_Src/io.F,v 1.1.1.1 2004/06/25 14:29:50 spk Exp $
C $Name:  $
C================================================
C This subroutine reads an integer*4 1d record of length n
C================================================
      SUBROUTINE READ_I4SEG(n,xdum,fid,irec)

      implicit none
      integer n,fid
      integer*4 xdum(n)
      integer irec

      read(fid,rec=irec) xdum

      END

C================================================
C This subroutine writes an integer*4 1d record of length n
C================================================
      SUBROUTINE WRITE_I4SEG(n,xdum,fid,irec)

      implicit none
      integer n,fid
      integer*4 xdum(n)
      integer irec

      write(fid,rec=irec) xdum

      END

C================================================
C This subroutine reads a real*4 1d record of length n
C================================================
      SUBROUTINE READ_R4SEG(n,xdum,fid,irec)

      implicit none
      integer n,fid
      real*4 xdum(n)
      integer irec

      read(fid,rec=irec) xdum

      END

C================================================
C This subroutine writes a real*4 1d record of length n
C================================================
      SUBROUTINE WRITE_R4SEG(n,xdum,fid,irec)

      implicit none
      integer n,fid
      real*4 xdum(n)
      integer irec

      write(fid,rec=irec) xdum

      END

C     ================================================
C     This subroutine reads a real*8 1d record of length n
C     ================================================
      SUBROUTINE READ_R8SEG(n,xdum,fid,irec)

      implicit none
      integer n,fid
      real*8 xdum(n)
      integer irec

      read(fid,rec=irec) xdum

      END
      
C================================================
C This subroutine writes a real*8 1d record of length n
C================================================
      SUBROUTINE WRITE_R8SEG(n,xdum,fid,irec)

      implicit none
      integer n,fid
      real*8 xdum(n)
      integer irec

      write(fid,rec=irec) xdum

      END

C==================================================
C This subroutine reads an integer*4 array from file
C==================================================
      SUBROUTINE READ_I4_FIELD(nx,ny,iSlice,fld,filename)

      implicit none
      integer nx,ny,iSlice
      integer*4 fld(nx,ny)
      character*(*) filename

      integer i,j,irec,ioUnit

      call findunit(ioUnit)

      open(unit=ioUnit,file=filename,status='old',
     &     access='direct',recl=nx*4)
     
      do j=1,ny
        irec=j+(iSlice-1)*ny
        call read_i4seg(nx,fld(1,j),ioUnit,irec)
      enddo

      close(ioUnit)
      
      END

C==================================================
C This subroutine reads a real*4 array from file
C==================================================
      SUBROUTINE READ_R4_FIELD(nx,ny,iSlice,fld,filename)

      implicit none
      integer nx,ny,iSlice
      real*4 fld(nx,ny)
      character*(*) filename

      integer i,j,irec,ioUnit

      call findunit(ioUnit)

      open(unit=ioUnit,file=filename,status='old',
     &     access='direct',recl=nx*4)
     
      do j=1,ny
        irec=j+(iSlice-1)*ny
        call read_r4seg(nx,fld(1,j),ioUnit,irec)
      enddo

      close(ioUnit)
      
      END

C==================================================
C This subroutine reads a real*8 array from file
C==================================================
      SUBROUTINE READ_R8_FIELD(nx,ny,iSlice,fld,filename)

      implicit none
      integer nx,ny,iSlice
      real*8 fld(nx,ny)
      character*(*) filename

      integer i,j,irec,ioUnit

      call findunit(ioUnit)

      open(unit=ioUnit,file=filename,status='old',
     &     access='direct',recl=nx*8)
     
      do j=1,ny
        irec=j+(iSlice-1)*ny
        call read_r8seg(nx,fld(1,j),ioUnit,irec)
      enddo

      close(ioUnit)
      
      END

C==================================================
C This subroutine writes out an integer*4 array to file
C==================================================
      SUBROUTINE WRITE_I4_FIELD(nx,ny,iSlice,fld,filename)

      implicit none
      integer nx,ny,iSlice
      integer*4 fld(nx,ny)
      character*(*) filename

      integer i,j,irec,ioUnit

      call findunit(ioUnit)

      if (iSlice.eq.1) then  !  open file here for first time
        open(unit=ioUnit,file=filename,status='unknown',
     &       access='direct',recl=nx*4)
      else  !  open existing file for appending
        open(unit=ioUnit,file=filename,status='old',
     &       access='direct',recl=nx*4)
      endif

      do j=1,ny
        irec=j+(iSlice-1)*ny
        call write_i4seg(nx,fld(1,j),ioUnit,irec)
      enddo

      close(ioUnit)
      
      END
            
C==================================================
C This subroutine writes out a real*4 array to file
C==================================================
      SUBROUTINE WRITE_R4_FIELD(nx,ny,iSlice,fld,filename)

      implicit none
      integer nx,ny,iSlice
      real*4 fld(nx,ny)
      character*(*) filename

      integer i,j,irec,ioUnit

      call findunit(ioUnit)

      if (iSlice.eq.1) then  !  open file here for first time
        open(unit=ioUnit,file=filename,status='unknown',
     &       access='direct',recl=nx*4)
      else  !  open existing file for appending
        open(unit=ioUnit,file=filename,status='old',
     &       access='direct',recl=nx*4)
      endif

      do j=1,ny
        irec=j+(iSlice-1)*ny
        call write_r4seg(nx,fld(1,j),ioUnit,irec)
      enddo

      close(ioUnit)
      
      END

C==================================================
C This subroutine writes out a real*8 array to file
C==================================================
      SUBROUTINE WRITE_R8_FIELD(nx,ny,iSlice,fld,filename)

      implicit none
      integer nx,ny,iSlice
      real*8 fld(nx,ny)
      character*(*) filename

      integer i,j,irec,ioUnit

      call findunit(ioUnit)

      if (iSlice.eq.1) then  !  open file here for first time
        open(unit=ioUnit,file=filename,status='unknown',
     &       access='direct',recl=nx*8)
      else  !  open existing file for appending
        open(unit=ioUnit,file=filename,status='old',
     &       access='direct',recl=nx*8)
      endif

      do j=1,ny
        irec=j+(iSlice-1)*ny
        call write_r8seg(nx,fld(1,j),ioUnit,irec)
      enddo

      close(ioUnit)
      
      END

C==================================================
C This subroutine reads a sparse matrix stored in binary COO format
C Inputs: 
C  nz: number of nonzero elements
C  filename: name of binary file
C Outputs:
C  ia(nz),ja(nz): integer arrays
C  a(nz): real*8 array
C Note: The size of the binary file is 16*nz bytes. 
C The first 4*nz bytes are ia, the next 4*nz bytes are ja, 
C and the final 8*nz bytes.
C==================================================
      SUBROUTINE READ_BCOO(nz,ia,ja,a,filename)

      implicit none
      integer nz,ia(nz),ja(nz)
      real*8 a(nz)
      character*(*) filename

      call read_i4_field(nz,1,1,ia,filename)  !  read ia

      call read_i4_field(nz,1,2,ja,filename)  !  read ja

      call read_r8_field(nz,1,2,a,filename)  !  read a

      END
      
C=============================================================
C Find a free unit for I/O
C=============================================================
      SUBROUTINE FINDUNIT( iounit )
C OUT:
C     iounit   integer - unit number
C
C MDSFINDUNIT returns a valid, unused unit number for f77 I/O
C The routine stops the program is an error occurs in the process
C of searching the I/O channels.
C
C Created: 03/20/99 adcroft@mit.edu

      implicit none

C Arguments
      integer iounit
C Local
      integer ii
      logical op
      integer ios
C     ------------------------------------------------------------------

C Sweep through a valid range of unit numbers
      iounit=-1
      do ii=9,99
        if (iounit.eq.-1) then
          inquire(unit=ii,iostat=ios,opened=op)
          if (ios.ne.0) then
            stop 'ABNORMAL END: S/R MDSFINDUNIT'
          endif
          if (.NOT. op) then
            iounit=ii
          endif
        endif
      enddo

C Was there an available unit number
      if (iounit.eq.-1) then
        stop ' MDSFINDUNIT: could not find an available unit number!'
      endif

C     ------------------------------------------------------------------
      return
      END

