      subroutine ioinit

!=======================================================================
!     Ioinit initializes logical array "inuse" to show which i/o units
!     are available to the subroutine getunitnumber.  It also
!     initializes the array fname of file names associated with each
!     i/o unit, and opens a file "iohistory", if enabled.
!=======================================================================

      implicit none

      integer last_reserved_unit, i, last_reserved
      parameter (last_reserved_unit=0)

      include "stdunits.h"
      include "iomngr.h"

      do i=0,maxfilenames
        hide_file(i) = .false.
        iunit(i) = 0
        fname(i) = ' '
      enddo

      fname(0) = ' '
      fname(1) = 'reserved'
      nfiles = 1

!     reserve units 1..last_reserved_unit (making them not available)
!             units last_reserved_unit+1..maxunit=99 are available

!     unnecessary assignment to avoid WARNING on sgi f77 compiler
      last_reserved = last_reserved_unit

      do i=1,last_reserved
        inuse(i) = .true.
        ifile(i) = 1
        scratch_file(i) = .false.
        unix_name(i) = .true.
      enddo
      do i=last_reserved_unit+1,maxunit
        inuse(i) = .false.
        ifile(i) = 0
        scratch_file(i) = .false.
        unix_name(i) = .true.
      enddo

!     reserve units 5 and 6, whether or not they are stdin/stdout

      inuse(5) = .true.
      inuse(6) = .true.
      ifile(5) = 1
      ifile(6) = 1

!     reserve standard i/o files:  stdin, stdout, stderr

      call link_unit (stdin, 'stdin')
      call link_unit (stdout, 'stdout')
      if (stderr .ne. stdout) then
        call link_unit (stderr, 'stderr')
      else
        fname(ifile(stdout)) = 'stdout/stderr'
      endif

      return
      end

      subroutine iomngr_resume

!=======================================================================
!     iomngr_resume resumes normal operation of the iomngr after a
!     call release_all has closed all files cleanly to permit executing
!     a section of code that does not use the iomngr to allocate
!     units.  It reopens the file "iohistory", if enabled.
!=======================================================================

      implicit none

      include "stdunits.h"
      include "iomngr.h"

      return
      end

      subroutine getunitnumber (iounit)

!=======================================================================
!     getunitnumber assigns a Fortran i/o unit number that does not
!     conflict with any i/o unit currently in use.

!     comments:

!     see relunit which releases a Fortran i/o unit number.
!=======================================================================

      implicit none

      character(120) :: name

      integer i, iostat, iounit

      logical opened

      include "stdunits.h"
      include "iomngr.h"

      do i = 1,maxunit
        if (.NOT. inuse(i)) then
!         check if opened without getunit.
          inquire(i,OPENED=opened)
          if (.not. opened) goto 101
        endif
      enddo

!     "Normal" termination of the do loop indicates all units are in use

      write (stdout, '(a4)')
     &    'ERROR:  No i/o units available -- Execution terminated'
      call showunits
      stop '==>Getunit'

  101 continue
      iounit = i
      inuse(i) = .true.

      return
      end

      subroutine getunit (iounit, oldfilename, optionlist)

!=======================================================================
!     getunit gets an unused Fortran unit number, assigns it to the
!     dummy argument iounit, and opens it with FILE=filename and options
!     given in optionlist.  Enable cpp directive "print_io_description"
!     to see the list of options.
!=======================================================================

      implicit none

      character(*) :: oldfilename, optionlist
      character(120) :: filename, new_file_name
      character(160) :: options
      character(2) :: aiounit
      character(15) :: ssequential, sdirect
      character(15) :: sformatted, sunformatted
      character(15) :: srewind, sappend
      character(15) :: sread, swrite, sreadwrite, sunknown, sold, snew
      character(15) :: sstatus, saccess, sform, sposition, saction
      character(15) :: sieee, srec_len, sdummy
      character(15) :: sblocking, psblocking, snumeric
      character(15) :: ssave, sovfl
      character(15) :: sffio, snamelist

      logical scan_option, scan_size, exists, named
      logical skip_open

      integer errors, irec_len, iblock, iword, ibyte, inum
      integer ibuffer, maxsize, layer, i, iounit

      include "stdunits.h"
      include "iomngr.h"

      data ssequential  /'SEQUENTIAL'/
      data sdirect      /'DIRECT'/
      data sformatted   /'FORMATTED'/
      data sunformatted /'UNFORMATTED'/
      data srewind      /'REWIND'/
      data sappend      /'APPEND'/
      data sread        /'READ'/
      data swrite       /'WRITE'/
      data sreadwrite   /'READWRITE'/
      data sunknown     /'UNKNOWN'/
      data sold         /'OLD'/
      data snew         /'NEW'/

      filename = new_file_name (oldfilename)

!     set default file attributes

      saccess     = ssequential
      sform       = ' '
      sposition   = srewind
      saction     = sreadwrite
      sstatus     = sunknown
      sblocking   = ' '
      psblocking  = ' '
      sieee       = ' '
      snumeric    = ' '
      ssave       = 'save'
      sovfl       = 'novfl'
      irec_len    = 0
      iblock      = 0
      iword       = 0
      ibyte       = 0
      inum        = 0
      ibuffer     = 0
      maxsize     = 0
      layer       = 0
      sffio       = ' '
      snamelist   = ' '

!     convert options list to lowercase

      options = optionlist
      call tolower (options)

!     decode options:

      i = 1
      errors = 0
100   continue
      if (i .gt. len_trim(options)) go to 101

!       scan options by full names

        if (scan_option(options, i,
     &          'sequential', saccess, ssequential)) then
        elseif (scan_option(options, i,
     &          'direct', saccess, sdirect)) then
          call setreclen (irec_len,srec_len,sform,iblock,iword,ibyte)
          if (irec_len .eq. 0) then
            print '(a,a,/,a,a)',
     &         '=>ERROR in getunit options list:',
     &         ' ''direct'' requires a record length',
     &         ' options=',options
            errors = errors + 1
          endif
        elseif (scan_option(options, i, 'recl', sdummy, sdummy)) then
          call setreclen (irec_len,srec_len,sform,iblock,iword,ibyte)
          if (irec_len .eq. 0) then
            print '(a,a,/,a,a)',
     &         '=>ERROR in getunit options list:',
     &         ' ''recl='' requires a record length',
     &         ' options=',options
            errors = errors + 1
          endif
        elseif (scan_option(options, i,
     &          'formatted', sform, sformatted)) then
        elseif (scan_option(options, i,
     &          'unformatted', sform, sunformatted)) then
        elseif (scan_option(options, i,
     &          'readwrite', saction, sreadwrite)) then
        elseif (scan_option(options, i,
     &          'read', saction, sread)) then
          print '(a,a)', 'WARNING:  ACTION=''read'' not supported.  ',
     &                            'ACTION=''READWRITE'' is used.'
        elseif (scan_option(options, i,
     &          'write', saction, swrite)) then
          print '(a,a)', 'WARNING:  ACTION=''write'' not supported.  ',
     &                            'ACTION=''READWRITE'' is used.'
        elseif (scan_option(options, i,
     &          'rewind', sposition, srewind)) then
        elseif (scan_option(options, i,
     &          'append', sposition, sappend)) then
        elseif (scan_option(options, i, 'ieee', sieee, 'ieee')) then
        elseif (scan_option(options, i, 'f77namelist',
     &          snamelist, 'f77')) then

!       scan for blocking options

        elseif (scan_option(options, i, 'fortran',
     &          sblocking, 'f77')) then
        elseif (scan_option(options, i, 'f77',
     &          sblocking, 'f77')) then
        elseif (scan_option(options, i, 'stream',
     &          sblocking, 'null')) then
        elseif (scan_option(options, i, 'unblocked',
     &          sblocking, 'null')) then
        elseif (scan_option(options, i, 'null',
     &          sblocking, 'null')) then

!       the 'buffer' option is parsed, but ignored on many systems
!       f77 blocking is assumed

        elseif (scan_option(options, i, 'buffer', sdummy, sdummy)) then

!         force new_layer to recognize f77 blocking

          sblocking = 'f77'
          if (scan_size(options, i, ibuffer, iblock, iword, ibyte)) then
          endif
!         why assume that buffers are in bytes?   f77 buffer units=bytes.
          call set2byte (ibuffer, iblock, iword, ibyte)
            inum = ibuffer
          if (ibuffer .eq. 0) then
            print '(a,a,/,a,a)',
     &         '=>ERROR in getunit options list:',
     &         ' ''buffer='' requires a size',
     &         ' options=',options
            errors = errors + 1
          endif
!       scan options by abbreviated names

        elseif (scan_option(options, i, 's', saccess, ssequential)) then
        elseif (scan_option(options, i, 'd', saccess, sdirect)) then
          call setreclen (irec_len,srec_len,sform,iblock,iword,ibyte)
          if (irec_len .eq. 0) then
            print '(a,a,/,a,a)',
     &         '=>ERROR in getunit options list:',
     &         ' ''direct'' requires a record length',
     &         ' options=',options
            errors = errors + 1
          endif
        elseif (scan_option(options, i, 'f', sform, sformatted)) then
        elseif (scan_option(options, i, 'u', sform, sunformatted)) then
        elseif (scan_option(options, i, 'rw', saction, sreadwrite)) then
        elseif (scan_option(options, i, 'ro', saction, sread)) then
        elseif (scan_option(options, i, 'wo', saction, swrite)) then
        elseif (scan_option(options, i, 'r', sposition, srewind)) then
        elseif (scan_option(options, i, 'a', sposition, sappend)) then

!       scan for size options with no introductory keyword
!       the introductory keyword must be several options back

        elseif (scan_size(options, i, inum, iblock, iword, ibyte)) then
        else
          print '(4a)', '=>ERROR:  Getunit unable to parse options=',
     &                    options(1:i-1), '|',
     &                    options(i:len_trim(options))
          errors = errors + 1
          i = i + 1
        endif
      goto 100

!     end of scan loop
!     a do loop is not be used because i is changed within the loop

101   continue
      if (errors .gt. 0) then
        print '(/a/a,i3,a)'
     &,   'STOP in getunit'
     &,   'There were',errors,' errors scanning the list of I/O options'
        stop
      endif

!     in case neither "formatted" or "unformatted" is specified, choose
!     the proper default for sequential and direct access files.

      if (sform .eq. ' ') then
        if (saccess .eq. ssequential) then
          sform = sformatted
        elseif (saccess .eq. sdirect) then
          sform = sunformatted
        endif
      endif
!-----------------------------------------------------------------------
!     See if the file has previously been opened in this run
!-----------------------------------------------------------------------

      do i=1,nfiles
        if (fname(i) .eq. filename) then
          sstatus = sold
          if (iunit(i) .ne. 0) then

!           This happens when a call relunit was executed for this file,
!           but due to system bugs, it was necessary to merely hide the
!           file while keeping it connected to a unit to permit later
!           appending.

            iounit = iunit(i)
            ifile(iounit) = i
            inuse(iounit) = .true.
            if (sposition .eq. sappend) then
              skip_open = .true.
            else

!             A previously hidden file that is now being opened with
!             POSITION=REWIND.  It should first be closed to permit
!             reading from it.

              call close_file (iounit, filename)
              skip_open = .false.

            endif
          else

!           This file is now closed, but had previously been opened.
!           Get a new unit number and reconnect it.

            call getunitnumber (iounit)
            ifile(iounit) = i
            iunit(i) = iounit

!           STATUS=UNKNOWN may cause rewind on some computers.
!           This is incompatible with append

            if (sposition .eq. sappend) then
              sstatus = sold
            endif
            skip_open = .false.
          endif
          go to 201
        endif
      enddo

!     This file has not previously been opened in this run.

!     Get a free unit number and bind it to the file
!     unless the request is 'sdsalloc' which does not
!     use a true i/o unit number.

      if (saccess .ne. 'sdsalloc') then
        call getunitnumber (iounit)
        call link_unit (iounit, filename)
        if (sposition .eq. sappend) then

!         If file pre-exists, append to it, but if no such file
!         exists, open a new file and rewind it.

          inquire (FILE=fname(nfiles), EXIST=exists, NAMED=named)
          if (.not. exists .or. .not. named) then
            sstatus = snew
            sposition = srewind
          else
            sstatus = sold
          endif

        endif
      endif
      skip_open = .false.

201   continue

      write (aiounit,'(i2)') iounit

!     make 'stream' mean: sequential access, unformatted

      if (sblocking .eq. 'null') then
        saccess = ssequential
        sform = 'UNFORMATTED'
      endif

!-----------------------------------------------------------------------
!     select open statement compatible with options and computer platform
!-----------------------------------------------------------------------

      if (skip_open) then
        goto 301
      endif
      if (saccess .eq. ssequential .or. saccess .eq. sappend) then

!       some compilers may allow a record length specifier for
!       sequential access files.  iomngr supports this possibility

        if (irec_len .eq. 0) then
        open (iounit,
     &        FILE=filename,
     &        FORM=sform,
     &        ACCESS=saccess,
     &        STATUS=sstatus
     &,       POSITION=sposition
     &        )
        else
        open (iounit,
     &        FILE=filename,
     &        FORM=sform,
     &        ACCESS=saccess,
     &        STATUS=sstatus
     &,       POSITION=sposition
     &,       RECL=irec_len
     &        )
        endif
        unix_name(iounit) = .true.
      elseif (saccess .eq. sdirect) then
        open (iounit,
     &        FILE=filename,
     &        FORM=sform,
     &        ACCESS=saccess,
     &        RECL=irec_len,
     &        STATUS=sstatus
     &        )
!        unix_name(iounit) = .true.
        unix_name(iounit) = .false.
      endif

301   continue

      return
      end

      subroutine relunit (iounit)

!=======================================================================
!     Relunit releases a Fortran i/o unit number currently in use and
!     closes the file.
!     see Getunit which gets a Fortran i/o unit number.
!=======================================================================

      implicit none

      integer ifn, iounit

      include "stdunits.h"
      include "iomngr.h"

      ifn = ifile(iounit)
      if (iounit .ge. 1 .and. iounit .le. maxunit) then
        inuse(iounit) = .false.
        if (.not. hide_file(ifn) .and. unix_name(iounit)) then
          call close_file (iounit, fname(ifn))
          iunit(ifn) = 0
        endif
        ifile(iounit) = 0
      else
        print *,'=>ERROR: iounit is bad. iounit=',iounit
        stop 'relunit'
      endif

      return
      end

      subroutine release_all

!=======================================================================

!     Releases all user i/o units currently in use and their file names.
!     All files except stdin, stdout, and stderr are closed.
!     Hidden files are closed at this time.

!     After a call to release_all, it is safe to execute a section of
!     code that does not use the iomngr to allocate units.

!=======================================================================

      implicit none

      integer ifn

      include "stdunits.h"
      include "iomngr.h"

      do ifn=1,nfiles
        if (iunit(ifn) .ne. 0) then
          if ((inuse(iunit(ifn)) .or. hide_file(ifn)) .and.
     &        unix_name(iunit(ifn)) .and.
     &        fname(ifn) .ne. 'reserved' .and.
     &        fname(ifn) .ne. 'stdin' .and.
     &        fname(ifn) .ne. 'stdout' .and.
     &        fname(ifn) .ne. 'stderr' .and.
     &        fname(ifn) .ne. 'stdout/stderr') then
            call close_file (iunit(ifn), fname(ifn))

!           mark file as released

            inuse(iunit(ifn)) = .false.
            iunit(ifn) = 0
          endif
        endif
      enddo

      return
      end

      subroutine close_file (iounit, filename)
!=======================================================================

!     close Fortran unit iounit

!=======================================================================

      implicit none

      integer iounit
      character(120) :: filename

      close (iounit)

      return
      end

      subroutine showunits

!=======================================================================

!     lists all i/o units currently in use and their file names.

!=======================================================================

      implicit none

      integer i

      include "stdunits.h"
      include "iomngr.h"

      write(stdout, '(/,a,/)') 'i/o units in use'
      write(stdout, '(a5,tr1,a,tr1,a)') ' unit', 'hidden', 'file name'
      do i=1,maxunit
        if (inuse(i)) then
          write(stdout, '(i4,tr3,l2,tr4,a)')
     &              i, hide_file(ifile(i)),
     &              fname(ifile(i))(1:len_trim(fname(ifile(i))))
        endif
      enddo

      return
      end
!=======================================================================

      subroutine showfiles

!=======================================================================

!     lists all file names that have been used [in chronological order]

!=======================================================================

      implicit none

      integer i

      include "stdunits.h"
      include "iomngr.h"

      write(stdout, '(/,a,/)') '  Files Used in This Run'
      write(stdout, '(a5,tr1,a,tr1,a)') '     ', 'hidden', 'file name'
      do i=1,nfiles
        write(stdout, '(i4,tr3,l2,tr4,a)')
     &            i, hide_file(i),
     &              fname(i)(1:len_trim(fname(i)))
      enddo

      return
      end
!=======================================================================

      subroutine link_unit (iounit, filename)

!=======================================================================
!     links an i/o unit number to a new file name
!-----------------------------------------------------------------------

      implicit none

      character(*) :: filename

      integer iounit

      include "stdunits.h"
      include "iomngr.h"

      nfiles = nfiles + 1
      inuse(iounit) = .true.
      fname(nfiles) = filename
      ifile(iounit) = nfiles
      iunit(nfiles) = iounit

      return
      end

!=======================================================================

      subroutine setreclen (irec_len,srec_len,sform,iblock,iword,ibyte)

!-----------------------------------------------------------------------
!     calculates record length for direct access files
!-----------------------------------------------------------------------

      implicit none

      character(*) :: srec_len, sform
      character(15) :: sformatted
      data sformatted   /'FORMATTED'/

      integer irec_len, iblock, iword, ibyte

      if (sform .eq. sformatted) then
        call set2byte (irec_len, iblock, iword, ibyte)
      else
        call set2word (irec_len, iblock, iword, ibyte)
      endif

      call i2a(irec_len, srec_len)

      return
      end

!=======================================================================

      subroutine set2byte (inum, iblock, iword, ibyte)

      implicit none

      integer inum, ibyte, iword, iblock

!-----------------------------------------------------------------------
!     if first argument is zero, sets first argument to value in bytes
!-----------------------------------------------------------------------
      if (inum .eq. 0) then
        if (ibyte .ne. 0) then
          inum = ibyte
        elseif (iword .ne. 0) then
          inum = 4 * iword
        endif
      endif

      return
      end

!=======================================================================

      subroutine set2word (inum, iblock, iword, ibyte)

      implicit none

      integer inum, ibyte, iword, iblock

!-----------------------------------------------------------------------
!     if first argument is zero, sets first argument to value in bytes
!-----------------------------------------------------------------------
      if (inum .eq. 0) then
        if (ibyte .ne. 0) then
          inum = 1 + (ibyte - 1) / 4
        elseif (iword .ne. 0) then
          inum = iword
        endif
      endif

      return
      end

!=======================================================================
      subroutine tolower (string)

!-----------------------------------------------------------------------
!     converts all alphabetic characters in string to lowercase
!-----------------------------------------------------------------------

      implicit none

      character(*) :: string
      character(1) :: c

      integer lcshift, i
      parameter (lcshift=32)

      do i=1,len(string)
        c = string(i:i)
        if ('A' .le. c .and. c .le. 'Z') then
          string(i:i) = char(ichar(c) + lcshift)
        endif
      enddo

      return
      end

!=======================================================================

      subroutine i2a (i, a)

      implicit none

      character(*) :: a
      character(15) :: string

      integer i, j

      write (string, '(i15)') i
      do j=1,15
        if (string(j:j) .ne. ' ') go to 101
      enddo
101   continue
      a = string(j:)

      return
      end

!=======================================================================

      function scan_option (options, i, spattern, svar, soption)

!-----------------------------------------------------------------------
!     scans the options string for spattern starting at position i.
!     if found, svar is set to soption.
!-----------------------------------------------------------------------

      implicit none

      logical scan_option

      character(*) :: options, spattern, svar, soption
      character(1) :: c

      integer i, j

!     skip delimiters

100   continue
      if (i .le. len(options)) then
        c = options(i:i)
        if (c.eq.'=' .or. c.eq.' ' .or. c.eq.',' .or. c.eq.':') then
          i = i + 1
          goto 100
          endif
        endif

!     scan for option

      scan_option = .false.
      j = i + len_trim(spattern) - 1
      if (j .le. len(options)) then
        if (options(i:j) .eq. spattern) then
          scan_option = .true.
          svar = soption
          i = i + len_trim(spattern)
        endif
      endif

      return
      end

!=======================================================================

      function scan_size(options, i, inum, iblock, iword, ibyte)

!-----------------------------------------------------------------------
!     scans the options string for a size specification starting at
!     position i.
!     if found, one of inum, iblock, iword, ibyte is nonzero.
!-----------------------------------------------------------------------

      implicit none

      logical scan_size, scan_number, scan_option

      character(*) :: options
      character(15) :: snum, sblock, sword, sbyte, sdummy

      integer i, inum, iblock, iword, ibyte

!     scan for size options

      scan_size = .false.

      if (scan_number(options, i, inum, snum)) then
        iblock = 0
        iword = 0
        ibyte = 0
        scan_size = .true.
      elseif (scan_option(options, i, 'words', sdummy, sdummy)) then
        if (scan_number(options, i, iword, sword)) then
          inum = 0
          iblock = 0
          ibyte = 0
          scan_size = .true.
        else
          print '(a,a,/,a,a)',
     &       '=>ERROR in getunit options list:',
     &       ' ''words='' requires a size in words',
     &       ' options=',options
        endif
      elseif (scan_option(options, i, 'bytes', sdummy, sdummy)) then
        if (scan_number(options, i, ibyte, sbyte)) then
          inum = 0
          iblock = 0
          iword = 0
          scan_size = .true.
        else
          print '(a,a,/,a,a)',
     &       '=>ERROR in getunit options list:',
     &       ' ''bytes='' requires a size in bytes',
     &       ' options=',options
        endif
      endif

      return
      end

!=======================================================================

      function scan_number (options, i, n, sn)

!-----------------------------------------------------------------------
!     scans the options string for a number starting at position i.
!     if found, n is the number as an integer and sn as a string.
!-----------------------------------------------------------------------

      implicit none

      logical scan_number

      character(*) :: options, sn
      character(1) :: c

      integer i, n

!     skip delimiters

100   continue
      if (i .le. len(options)) then
        c = options(i:i)
        if (c.eq.'=' .or. c.eq.' ' .or. c.eq.',' .or. c.eq.':') then
          i = i + 1
          goto 100
        endif
      endif

!     scan for digits

      if (i .le. len(options)) then
        c = options(i:i)
        if ('0' .le. c .and. c .le. '9') then
          sn = c
          n = ichar(c) - ichar('0')
          scan_number = .true.
200       continue
          i = i + 1
          if (i .le. len(options)) then
            c = options(i:i)
            if ('0' .le. c .and. c .le. '9') then
              sn = sn(1:len_trim(sn)) // c
              n = 10 * n + ichar(c) - ichar('0')
              goto 200
            endif
          endif
        else
          scan_number = .false.
        endif
      endif

      return
      end
