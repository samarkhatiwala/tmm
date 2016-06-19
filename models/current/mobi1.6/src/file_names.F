      subroutine file_names
!=======================================================================
!     edit this routine to change file names or add directories to file
!     names. names may also be modified by editing a filenames file.
!     the filenames file will only be written to, if it already exists.
!     if you want an example of the format required, touch the file and
!     change one file below (with put names) to write the change.
!     changes from file filenames take precedence over changes below.

!     an old file name (first argument to put_names) will be replaced by
!     a new file name (second argument to put_names) when the iomngr or
!     mpp_io routines try to open a file which starts with the old file
!     name. if the start of the file matches more than one old name, the
!     first match (in the order of calls to put_names) will be used. if
!     the file name is longer than the old file name (but matches at the
!     beginning), only the matched section is replaced by the new file
!     name. if the new or old name is blank, the change is ignored.

!     caution: if an old file name is found at the start of more than
!     one file name, more than the intended file name may be changed.
!     for example the single call to put_names:
!       call put_names ('particles', 'x')
!     would change all files starting with "particles..." to "x..."
!     including changing "particles_initial..." to "x_initial..."
!     if you do not want to change "particles_initial...", you can do
!     the following (the order is important):
!       call put_names ('particles_initial', 'particles_initial')
!       call put_names ('particles', 'x')
!     you may want to use this feature to redirect many files at once.
!     for example:
!       call put_names ('t', 't_files/t')
!     would cause all files beginning with t to be read or written from
!     the directory t_files (unless a particular 't...' file was
!     renamed with put_names before this call).
!=======================================================================

      implicit none

      integer max_num_files
      parameter (max_num_files=200)

      character(120) :: old_file_names, new_file_names
      common /files_c/ old_file_names(max_num_files)
      common /files_c/ new_file_names(max_num_files)

      integer num_files
      common /files_i/ num_files

      integer iou, ios, n
      character(len=120) old_file_name, new_file_name, text
      logical exists

      num_files = 0
      text = 'old file name followed by new (see file_names.F)'

      inquire (file='filenames', exist=exists)
      if (exists) then
!       read file name changes from filenames
        call getunit (iou, 'filenames','formatted sequential rewind')
        do n=1,max_num_files
          read (iou, '(a)', IOSTAT=ios) old_file_name
          if (old_file_name /= text .and. ios == 0) then
            write(*,'(a)') '==> Warning: error reading file names'
            exit
          endif
          read (iou, '(a/a)', IOSTAT=ios) old_file_name, new_file_name
          if (ios > 0) exit
          call put_names (old_file_name, new_file_name)
        enddo
        call relunit (iou)
      endif

      print*, 'The following files will be renamed:'

      call put_names ('A_agg', 'data/A_agg')
      call put_names ('A_calb', 'data/A_calb')
      call put_names ('A_cfc', 'data/A_cfc')
      call put_names ('A_co2', 'data/A_co2')
      call put_names ('A_dc14', 'data/A_dc14')
      call put_names ('A_diff', 'data/A_diff')
      call put_names ('A_sat', 'data/A_sat')
      call put_names ('A_slat', 'data/A_slat')
      call put_names ('A_wind', 'data/A_wind')
      call put_names ('A_solar', 'data/A_solar')
      call put_names ('A_sulph', 'data/A_sulph')
      call put_names ('A_volc', 'data/A_volc')

      call put_names ('F_co2', 'data/F_co2')
      call put_names ('F_heat', 'data/F_heat')
      call put_names ('F_salt', 'data/F_salt')

      call put_names ('G_grid', 'data/G_grid')
      call put_names ('G_kmt', 'data/G_kmt')
      call put_names ('G_mskt', 'data/G_mskt')
      call put_names ('G_mskhreg', 'data/G_mskhreg')

      call put_names ('L_cropfra', 'data/L_cropfra')
      call put_names ('L_agricfra', 'data/L_agricfra')
      call put_names ('L_diurtemp', 'data/L_diurtemp')
      call put_names ('L_elev', 'data/L_elev')
      call put_names ('L_ice', 'data/L_ice')
      call put_names ('L_potveg', 'data/L_potveg')
      call put_names ('L_rivers', 'data/L_rivers')

      call put_names ('O_alk', 'data/O_alk')
      call put_names ('O_dc14', 'data/O_dc14')
      call put_names ('O_diffac', 'data/O_diffac')
      call put_names ('O_no3', 'data/O_no3')
      call put_names ('O_o2', 'data/O_o2')
      call put_names ('O_po4', 'data/O_po4')
      call put_names ('O_sal', 'data/O_sal')
      call put_names ('O_sealev', 'data/O_sealev')
      call put_names ('O_slhref', 'data/O_slhref')
      call put_names ('O_tau', 'data/O_tau')
      call put_names ('O_temp', 'data/O_temp')
      call put_names ('O_tidenrg_egbert', 'data/O_tidenrg_egbert')
      call put_names ('O_totcarb', 'data/O_totcarb')
      call put_names ('O_fe_dissolved', 'data/O_fe_dissolved')

      call put_names ('rest_embm', 'rest')
      call put_names ('restart_embm', 'data/restart')
      call put_names ('restart_2_embm', 'restart_2')
      call put_names ('tsi_embm', 'tsi')
      call put_names ('tavg_embm', 'tavg')

      call put_names ('rest_ism', 'rest')
      call put_names ('restart_ism', 'data/restart')
      call put_names ('restart_2_ism', 'restart_2')
      call put_names ('tsi_ism', 'tsi')
      call put_names ('tavg_ism', 'tavg')

      call put_names ('rest_mtlm', 'rest')
      call put_names ('restart_mtlm', 'data/restart')
      call put_names ('restart_2_mtlm', 'restart_2')
      call put_names ('tsi_mtlm', 'tsi')
      call put_names ('tavg_mtlm', 'tavg')

      call put_names ('rest_mom', 'rest')
      call put_names ('restart_mom', 'data/restart')
      call put_names ('restart_2_mom', 'restart_2')
      call put_names ('tsi_mom', 'tsi')
      call put_names ('tavg_mom', 'tavg')

      call put_names ('rest_sed', 'rest')
      call put_names ('restart_sed', 'data/restart')
      call put_names ('restart_2_sed', 'restart_2')
      call put_names ('tsi_sed', 'tsi')
      call put_names ('tavg_sed', 'tavg')

      call put_names ('G_subgrid_bathy', 'data/G_subgrid_bathy')
      call put_names ('G_subgrid_ocean_mask'
     &,               'data/G_subgrid_ocean_mask') 

      call put_names ('tsi_gsums', 'tsi')

!     write file name changes to filenames
      if (exists) then
        call getunit (iou, 'filenames','formatted sequential rewind')
        do n=1,num_files
          write (iou, '(a)') text
          write (iou, '(a/a)') old_file_names(n), new_file_names(n)
        enddo
        call relunit (iou)
      endif
      print*, ' '

      end subroutine file_names

      subroutine put_names (old_file_name, new_file_name)
!=======================================================================
!     put file names in the appropriate array location

!     input:
!     old_file_name = old file name
!     new_file_name = new file name

!     output:
!     old_file_names = updated array of old file names
!     new_file_names = updated array of new file names
!     num_files      = number of renamed files
!=======================================================================

      implicit none

      integer max_num_files
      parameter (max_num_files=200)

      character(120) :: old_file_names, new_file_names
      common /files_c/ old_file_names(max_num_files)
      common /files_c/ new_file_names(max_num_files)

      integer num_files
      common /files_i/ num_files

      integer n
      character (len=*) old_file_name, new_file_name

      do n=1,num_files
        if (trim(old_file_name) == trim(old_file_names(n))) then
          new_file_name = ''
        endif
      enddo

      if (trim(old_file_name) == '') new_file_name = ''

      if (new_file_name /= '') then
        num_files = num_files + 1
        old_file_names(num_files) = trim(old_file_name)
        new_file_names(num_files) = trim(new_file_name)
        print*, '  ',trim(old_file_name),' renamed ',trim(new_file_name)
      endif

      end subroutine put_names

      function new_file_name (old_file_name)

!=======================================================================
!     provides a new file name if the beginning of the old file name is
!     matched with a name from the old_file_names array. old and new
!     file names are set in the set_file_names routine. if a match is
!     not found, then the old file name is used.

!     input:
!     old_file_name = old file name

!     output:
!     new_file_name = new file name
!=======================================================================

      implicit none

      integer max_num_files
      parameter (max_num_files=200)

      character(120) :: old_file_names, new_file_names, new_file_name
      common /files_c/ old_file_names(max_num_files)
      common /files_c/ new_file_names(max_num_files)

      integer num_files
      common /files_i/ num_files

      integer i, j, k, m, n
      character (len=*) old_file_name

      i = 0
      new_file_name = trim(old_file_name)
      do n=1,num_files
        i = index(old_file_name, trim(old_file_names(n)))
        if (i == 1) then
          j = len_trim(old_file_name)
          k = len_trim(old_file_names(n))
          m = len_trim(new_file_names(n))
          new_file_name = trim(new_file_names(n))
          if (j > k) new_file_name(m+1:m+1+j-k) = old_file_name(k+1:j)
          exit
        endif
      enddo

      end function new_file_name
