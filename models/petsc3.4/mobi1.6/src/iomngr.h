!====================== include file "iomngr.h" ========================

!   arrays:
!     inuse        = does iomngr know unit number is currently in use
!     hide_file    = file which because of system quirks should not
!                    be closed by relunit
!     scratch_file = file should be deleted when released
!     unix_name    = file known to unix by name [exceptions: cray_ymp
!                    "word" and "sdsalloc" files]
!     ifile        = index into fname array corresponding to unit number
!     fname        = list of file names used since ioinit
!     iunit        = unit number corresponding to file name
!   scalars:
!     iohist       = unit number of io_history file
!     nfiles       = number of files names used so far

      integer maxunit, maxfilenames
      parameter (maxunit = 1000, maxfilenames = 1000)

      character(120) :: fname
      common /iomngr_c/ fname(0:maxfilenames)

      integer ifile, iunit, iohist, nfiles
      common /iomngr_i/ ifile(1:maxunit)
      common /iomngr_i/ iunit(0:maxfilenames), iohist, nfiles

      logical inuse, hide_file, scratch_file, unix_name
      common /iomngr_l/ inuse(1:maxunit), hide_file(0:maxfilenames)
      common /iomngr_l/ scratch_file(1:maxunit), unix_name(1:maxunit)
