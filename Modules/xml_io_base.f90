!
! Copyright (C) 2005-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE xml_io_base
  !----------------------------------------------------------------------------
  !
  ! ... this module contains some common subroutines used to read and write
  ! ... in XML format the data produced by Quantum ESPRESSO package
  !
  ! ... written by Carlo Sbraccia (2005)
  ! ... modified by Andrea Ferretti (2006-08)
  !
  USE iotk_module
  !
  USE kinds,     ONLY : DP
  USE io_files,  ONLY : tmp_dir, prefix, iunpun, xmlpun, &
                        current_fmt_version => qexml_version
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  USE parser,    ONLY : version_compare
  !
  IMPLICIT NONE
  PRIVATE
  !
  CHARACTER(5),  PARAMETER :: fmt_name = "QEXML"
  CHARACTER(5),  PARAMETER :: fmt_version = "1.4.0"
  !
  LOGICAL,       PARAMETER :: rho_binary = .TRUE.
  !
  CHARACTER(iotk_attlenx)  :: attr
  !
  !
  PUBLIC :: fmt_name, fmt_version
  PUBLIC :: current_fmt_version
  !
  PUBLIC :: rho_binary
  PUBLIC :: attr
  !
  PUBLIC :: create_directory, change_directory,                          &
            kpoint_dir, wfc_filename, copy_file,       &
            restart_dir, check_restartfile, check_file_exst,             &
            pp_check_file, save_history, save_print_counter,             &
            read_print_counter, set_kpoints_vars,                        &
            write_header, write_control, write_moving_cell,              &
            write_cell, write_ions, write_symmetry, write_planewaves,    &
            write_efield, write_spin, write_magnetization, write_xc,     &
            write_exx, write_occ, write_bz, write_para,                  &
            write_phonon, write_rho_xml, write_wfc, write_eig,           &
            read_wfc, read_rho_xml
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE create_directory( dirname )
      !------------------------------------------------------------------------
      !
      USE wrappers,  ONLY : f_mkdir
      USE mp,        ONLY : mp_barrier
      USE mp_global, ONLY : me_image, intra_image_comm
      USE io_files,  ONLY : check_writable
      !
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      !
      INTEGER                    :: ierr

      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      !
      IF ( ionode ) ierr = f_mkdir( TRIM( dirname ) )
      CALL mp_bcast ( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'create_directory', &
                   'unable to create directory ' // TRIM( dirname ), ierr )
      !
      ! ... syncronize all jobs (not sure it is really useful)
      !
      CALL mp_barrier( intra_image_comm )
      !
      ! ... check whether the scratch directory is writable
      !
      IF ( ionode ) ierr = check_writable ( dirname, me_image )
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'create_directory:', &
                   TRIM( dirname ) // ' non existent or non writable', ierr )
      !
      RETURN
      !
    END SUBROUTINE create_directory
    !
    !------------------------------------------------------------------------
    SUBROUTINE change_directory( dirname )
      !------------------------------------------------------------------------
      !
      USE wrappers,  ONLY : f_chdir
      USE mp,        ONLY : mp_barrier
      USE mp_global, ONLY : me_image, intra_image_comm
      !
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      !
      INTEGER                    :: ierr

      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      !
      ierr = f_chdir( TRIM( dirname ) )
      CALL mp_bcast ( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'change_directory', &
                   'unable to change to directory ' // TRIM( dirname ), ierr )
      !
      ! ... syncronize all jobs (not sure it is really useful)
      !
      CALL mp_barrier( intra_image_comm )
      !
      !
      RETURN
      !
    END SUBROUTINE change_directory
    !
    !------------------------------------------------------------------------
    FUNCTION kpoint_dir( basedir, ik )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=256)           :: kpoint_dir
      CHARACTER(LEN=*), INTENT(IN) :: basedir
      INTEGER,          INTENT(IN) :: ik
      !
      CHARACTER(LEN=256) :: kdirname
      CHARACTER(LEN=5)   :: kindex
      CHARACTER(LEN=6)   :: kindex1
      !
      IF (ik<99999) THEN
         WRITE( kindex, FMT = '( I5.5 )' ) ik     
         kdirname = TRIM( basedir ) // '/K' // kindex
      ELSEIF (ik<999999) THEN
         WRITE( kindex1, FMT = '( I6.6 )' ) ik     
         kdirname = TRIM( basedir ) // '/K' // kindex1
      ELSE
         call errore('kpoint_dir','ik too large, increase format',1)
      ENDIF
      !
      kpoint_dir = TRIM( kdirname )
      !
      RETURN
      !
    END FUNCTION kpoint_dir
    !
    !------------------------------------------------------------------------
    FUNCTION wfc_filename( basedir, name, ik, ipol, tag, extension, dir )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=256)                 :: wfc_filename
      CHARACTER(LEN=*),       INTENT(IN) :: basedir
      CHARACTER(LEN=*),       INTENT(IN) :: name
      INTEGER,                INTENT(IN) :: ik
      INTEGER,      OPTIONAL, INTENT(IN) :: ipol
      CHARACTER(*), OPTIONAL, INTENT(IN) :: tag
      CHARACTER(*), OPTIONAL, INTENT(IN) :: extension
      LOGICAL,      OPTIONAL, INTENT(IN) :: dir
      !    
      CHARACTER(LEN=256) :: filename, tag_, ext_
      LOGICAL :: dir_true
      !
      !
      filename = ' '
      tag_     = ' '
      ext_     = '.dat'
      dir_true = .true.
      !
      IF ( PRESENT( tag ) )         tag_ = '_'//TRIM(tag)
      IF ( PRESENT( extension ) )   ext_ = '.'//TRIM(extension)
      !
      IF ( PRESENT( ipol ) ) THEN
         !      
         WRITE( filename, FMT = '( I1 )' ) ipol
         !
      END IF
      IF ( PRESENT( dir )) dir_true=dir
      !
      IF (dir_true) THEN
         filename = TRIM( kpoint_dir( basedir, ik ) ) // '/' // &
                 & TRIM( name ) // TRIM( filename ) // TRIM( tag_ ) // TRIM( ext_)
      ELSE
         filename = TRIM( kpoint_dir( basedir, ik ) ) // '_' // &
                 & TRIM( name ) // TRIM( filename ) // TRIM( tag_ ) // TRIM( ext_)
      ENDIF
      !
      wfc_filename = TRIM( filename )
      !
      RETURN
      !
    END FUNCTION
    !
    !------------------------------------------------------------------------
    SUBROUTINE copy_file( file_in, file_out )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=*), INTENT(IN) :: file_in, file_out
      !
      CHARACTER(LEN=256) :: string
      INTEGER            :: iun_in, iun_out, ierr
      !
      !
      IF ( .NOT. ionode ) RETURN
      !
      CALL iotk_free_unit( iun_in,  ierr )
      CALL iotk_free_unit( iun_out, ierr )
      !
      CALL errore( 'copy_file', 'no free units available', ierr )
      !
      OPEN( UNIT = iun_in,  FILE = file_in,  STATUS = "OLD" )
      OPEN( UNIT = iun_out, FILE = file_out, STATUS = "UNKNOWN" )         
      !
      copy_loop: DO
         !
         READ( UNIT = iun_in, FMT = '(A256)', IOSTAT = ierr ) string
         !
         IF ( ierr < 0 ) EXIT copy_loop
         !
         WRITE( UNIT = iun_out, FMT = '(A)' ) TRIM( string )
         !
      END DO copy_loop
      !
      CLOSE( UNIT = iun_in )
      CLOSE( UNIT = iun_out )
      !
      RETURN
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    FUNCTION restart_dir( outdir, runit )
      !------------------------------------------------------------------------
      !
      ! KNK_nimage
      ! USE mp_global, ONLY:  my_image_id
      CHARACTER(LEN=256)           :: restart_dir
      CHARACTER(LEN=*), INTENT(IN) :: outdir
      INTEGER,          INTENT(IN) :: runit
      !
      CHARACTER(LEN=256)         :: dirname
      INTEGER                    :: strlen
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      !
      ! ... main restart directory
      !
      ! ... keep the line below ( this is the old style RESTARTXX ) !!!
      !
      ! dirname = 'RESTART' // int_to_char( runit )
      ! the next line is to have separate RESTART for each image
      ! KNK_nimage
      ! if (my_image_id > 0) dirname = trim(dirname) // '_' // trim(int_to_char( my_image_id ))
      !
      dirname = TRIM( prefix ) // '_' // TRIM( int_to_char( runit ) )// '.save'
      !
      IF ( LEN( outdir ) > 1 ) THEN
         !
         strlen = INDEX( outdir, ' ' ) - 1
         !
         dirname = outdir(1:strlen) // '/' // dirname
         !
      END IF
      !
      restart_dir = TRIM( dirname )
      !
      RETURN
      !
    END FUNCTION restart_dir
    !
    !------------------------------------------------------------------------
    FUNCTION check_restartfile( outdir, ndr )
      !------------------------------------------------------------------------
      !
      USE io_global, ONLY : ionode, ionode_id
      USE mp_global, ONLY : intra_image_comm
      !
      IMPLICIT NONE
      !
      LOGICAL                      :: check_restartfile
      INTEGER,          INTENT(IN) :: ndr
      CHARACTER(LEN=*), INTENT(IN) :: outdir
      CHARACTER(LEN=256)           :: filename
      LOGICAL                      :: lval
      !
      !
      filename = restart_dir( outdir, ndr )
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( filename ) // '/' // TRIM( xmlpun )
         !
         INQUIRE( FILE = TRIM( filename ), EXIST = lval )
         !
      END IF
      !
      CALL mp_bcast( lval, ionode_id, intra_image_comm )
      !
      check_restartfile = lval
      !
      RETURN
      !
    END FUNCTION check_restartfile
    !
    !------------------------------------------------------------------------
    FUNCTION check_file_exst( filename )
      !------------------------------------------------------------------------
      !
      USE io_global, ONLY : ionode, ionode_id
      USE mp_global, ONLY : intra_image_comm
      !
      IMPLICIT NONE
      !
      LOGICAL          :: check_file_exst
      CHARACTER(LEN=*) :: filename
      !
      LOGICAL :: lexists
      !
      IF ( ionode ) THEN 
         !
         INQUIRE( FILE = TRIM( filename ), EXIST = lexists )
         !
      ENDIF
      !
      CALL mp_bcast ( lexists, ionode_id, intra_image_comm )
      !
      check_file_exst = lexists
      RETURN
      !
    END FUNCTION check_file_exst
    !
    !------------------------------------------------------------------------
    FUNCTION pp_check_file()
      !------------------------------------------------------------------------
      !
      USE io_global,         ONLY : ionode, ionode_id
      USE mp_global,         ONLY : intra_image_comm
      USE control_flags,     ONLY : lkpoint_dir, tqr
      !
      IMPLICIT NONE
      !
      LOGICAL            :: pp_check_file
      CHARACTER(LEN=256) :: dirname, filename
      INTEGER            :: ierr
      LOGICAL            :: lval, found, back_compat
      !
      !
      dirname  = TRIM( tmp_dir ) // TRIM( prefix ) // '.save'
      filename = TRIM( dirname ) // '/' // TRIM( xmlpun )
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = filename, IERR = ierr )
      !
      CALL mp_bcast ( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'pp_check_file', 'file ' // &
                 & TRIM( dirname ) // ' not found', ierr )

      !
      ! set a flag for back compatibility (before fmt v1.4.0)
      !
      back_compat = .FALSE.
      !
      IF ( TRIM( version_compare( current_fmt_version, "1.4.0" )) == "older") &
         back_compat = .TRUE.
      !
      IF ( ionode ) THEN
         !
         IF ( .NOT. back_compat ) THEN
             !
             CALL iotk_scan_begin( iunpun, "CONTROL" ) 
             !
         ENDIF
         !
         CALL iotk_scan_dat( iunpun, "PP_CHECK_FLAG", lval, FOUND = found)
         !
         IF ( .NOT. found ) lval = .FALSE. 
         !
         CALL iotk_scan_dat( iunpun, "LKPOINT_DIR", lkpoint_dir, FOUND = found)
         !
         IF ( .NOT. found ) lkpoint_dir = .TRUE. 
         !
         CALL iotk_scan_dat( iunpun, "Q_REAL_SPACE", tqr, FOUND = found)
         !
         IF ( .NOT. found ) tqr = .FALSE. 
         !
         !
         IF ( .NOT. back_compat ) THEN
             !
             CALL iotk_scan_end( iunpun, "CONTROL" ) 
             !
         ENDIF
         !
         CALL iotk_close_read( iunpun )
         !
      END IF
      !
      CALL mp_bcast( lval, ionode_id, intra_image_comm )
      !
      CALL mp_bcast( lkpoint_dir, ionode_id, intra_image_comm )
      !
      CALL mp_bcast( tqr, ionode_id, intra_image_comm )
      !
      pp_check_file = lval
      !
      RETURN
      !
    END FUNCTION pp_check_file
    !
    !------------------------------------------------------------------------
    SUBROUTINE save_history( dirname, iter )
      !------------------------------------------------------------------------
      !
      ! ... a copy of the xml descriptor (data-file.xml) is saved in the 
      ! ... history subdir
      !
      USE io_files, ONLY : xmlpun_base
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      INTEGER,          INTENT(IN) :: iter
      !
#if defined (__VERBOSE_SAVE)
      !
      CHARACTER(LEN=256) :: filename
      CHARACTER(LEN=6)   :: hindex
      !
      CALL create_directory( TRIM( dirname ) // '/history' )
      !
      WRITE( hindex, FMT = '(I6.6)' ) iter
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( dirname ) // '/history/' // &
                  & TRIM( xmlpun_base ) // hindex // '.xml'
         !
         CALL copy_file( TRIM( dirname ) // "/" // TRIM( xmlpun ), &
                         TRIM( filename ) )
         !
      END IF
      !
#endif
      !
      RETURN
      !
    END SUBROUTINE save_history
    !
    !------------------------------------------------------------------------
    SUBROUTINE save_print_counter( iter, outdir, wunit )
      !------------------------------------------------------------------------
      !
      ! ... a counter indicating the last successful printout iteration is saved
      !
      USE io_global, ONLY : ionode, ionode_id
      USE mp_global, ONLY : intra_image_comm
      USE mp,        ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(IN) :: iter
      CHARACTER(LEN=*), INTENT(IN) :: outdir
      INTEGER,          INTENT(IN) :: wunit
      !
      INTEGER            :: ierr
      CHARACTER(LEN=256) :: filename, dirname
      !
      !
      dirname = restart_dir( outdir, wunit )
      !
      CALL create_directory( TRIM( dirname ) )
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( dirname ) // '/print_counter.xml'
         !
         CALL iotk_open_write( iunpun, FILE = filename, &
                             & ROOT = "PRINT_COUNTER",  IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'save_print_counter', &
                   'cannot open restart file for writing', ierr )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_write_begin( iunpun, "LAST_SUCCESSFUL_PRINTOUT" )
         CALL iotk_write_dat(   iunpun, "STEP", iter )
         CALL iotk_write_end(   iunpun, "LAST_SUCCESSFUL_PRINTOUT" )
         !
         CALL iotk_close_write( iunpun )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE save_print_counter
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_print_counter( nprint_nfi, outdir, runit )
      !------------------------------------------------------------------------
      !
      ! ... the counter indicating the last successful printout iteration 
      ! ... is read here
      !
      USE io_global, ONLY : ionode, ionode_id
      USE mp_global, ONLY : intra_image_comm
      USE mp,        ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: nprint_nfi
      CHARACTER(LEN=*), INTENT(IN)  :: outdir
      INTEGER,          INTENT(IN)  :: runit
      !
      INTEGER            :: ierr
      CHARACTER(LEN=256) :: filename, dirname
      !
      !
      dirname = restart_dir( outdir, runit )
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( dirname ) // '/print_counter.xml'
         !
         CALL iotk_open_read( iunpun, FILE = filename, IERR = ierr )
         !
         IF ( ierr > 0 ) THEN
            !
            nprint_nfi = -1
            !
         ELSE
            !
            CALL iotk_scan_begin( iunpun, "LAST_SUCCESSFUL_PRINTOUT" )
            CALL iotk_scan_dat(   iunpun, "STEP", nprint_nfi )
            CALL iotk_scan_end(   iunpun, "LAST_SUCCESSFUL_PRINTOUT" )
            !
            CALL iotk_close_read( iunpun )
            !
         END IF
         !
      END IF
      !
      CALL mp_bcast( nprint_nfi, ionode_id, intra_image_comm )
      !
      RETURN
      !
    END SUBROUTINE read_print_counter   
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_kpoints_vars( ik, nk, kunit, ngwl, igl, &
                                 ngroup, ikt, iks, ike, igwx, ipmask, ipsour, &
                                 ionode, root_in_group, intra_group_comm, inter_group_comm, parent_group_comm )
      !------------------------------------------------------------------------
      !
      ! ... set working variables for k-point index (ikt) and 
      ! ... k-points number (nkt)
      !
      USE mp,         ONLY : mp_sum, mp_get, mp_max, mp_rank, mp_size
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)  :: ik, nk, kunit
      INTEGER, INTENT(IN)  :: ngwl, igl(:)
      INTEGER, INTENT(OUT) :: ngroup
      INTEGER, INTENT(OUT) :: ikt, iks, ike, igwx
      INTEGER, INTENT(OUT) :: ipmask(:), ipsour
      LOGICAL, INTENT(IN)  :: ionode
      INTEGER, INTENT(IN)  :: root_in_group, intra_group_comm, inter_group_comm, parent_group_comm
      !
      INTEGER :: ierr, i
      INTEGER :: nkl, nkr, nkbl, nkt
      INTEGER :: nproc_parent, nproc_group, my_group_id, me_in_group, me_in_parent, io_in_parent
      !
      nproc_parent = mp_size( parent_group_comm )
      nproc_group  = mp_size( intra_group_comm )
      my_group_id  = mp_rank( inter_group_comm )
      me_in_group  = mp_rank( intra_group_comm )
      me_in_parent = mp_rank( parent_group_comm )
      !
      ! find the ID (io_in_parent) of the io PE ( where ionode == .true. )
      !
      io_in_parent = 0
      IF( ionode ) io_in_parent = me_in_parent
      CALL mp_sum( io_in_parent, parent_group_comm )
      !
      ikt = ik
      nkt = nk
      !
      ! ... find out the number of pools
      !
      ngroup = nproc_parent / nproc_group 
      !
      ! ... find out number of k points blocks
      !
      nkbl = nkt / kunit  
      !
      ! ... k points per pool
      !
      nkl = kunit * ( nkbl / ngroup )
      !
      ! ... find out the reminder
      !
      nkr = ( nkt - nkl * ngroup ) / kunit
      !
      ! ... Assign the reminder to the first nkr pools
      !
      IF ( my_group_id < nkr ) nkl = nkl + kunit
      !
      ! ... find out the index of the first k point in this pool
      !
      iks = nkl * my_group_id + 1
      !
      IF ( my_group_id >= nkr ) iks = iks + nkr * kunit
      !
      ! ... find out the index of the last k point in this pool
      !
      ike = iks + nkl - 1
      !
      ipmask = 0
      ipsour = io_in_parent
      !
      ! ... find out the index of the processor which collect the data 
      ! ... in the pool of ik
      !
      IF ( ngroup > 1 ) THEN
         !
         IF ( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
            !
            IF ( me_in_group == root_in_group ) ipmask( me_in_parent + 1 ) = 1
            !
         END IF
         !
         ! ... Collect the mask for all proc in the image
         !
         CALL mp_sum( ipmask, parent_group_comm )
         !
         DO i = 1, nproc_parent
            !
            IF( ipmask(i) == 1 ) ipsour = ( i - 1 )
            !
         END DO
         !
      END IF
      !
      igwx = 0
      ierr = 0
      !
      IF ( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
         !
         IF ( ngwl > SIZE( igl ) ) THEN
            !
            ierr = 1
            !
         ELSE
            !
            igwx = MAXVAL( igl(1:ngwl) )
            !
         END IF
         !
      END IF
      !
      ! ... get the maximum index within the pool
      !
      CALL mp_max( igwx, intra_group_comm )
      !
      ! ... now notify all procs if an error has been found 
      !
      CALL mp_max( ierr, parent_group_comm )
      !
      CALL errore( 'set_kpoint_vars ', 'wrong size ngl', ierr )
      !
      IF ( ipsour /= io_in_parent ) &
         CALL mp_get( igwx, igwx, me_in_parent, io_in_parent, ipsour, 1, parent_group_comm )
      !
      RETURN
      !
    END SUBROUTINE set_kpoints_vars
    !
    !
    ! ... writing subroutines
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_header( creator_name, creator_version ) 
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: creator_name, creator_version


      CALL iotk_write_begin( iunpun, "HEADER" )
      !
      CALL iotk_write_attr(attr, "NAME",TRIM(fmt_name), FIRST=.TRUE.)
      CALL iotk_write_attr(attr, "VERSION",TRIM(fmt_version) )
      CALL iotk_write_empty( iunpun, "FORMAT", ATTR=attr )
      !
      CALL iotk_write_attr(attr, "NAME",TRIM(creator_name), FIRST=.TRUE.)
      CALL iotk_write_attr(attr, "VERSION",TRIM(creator_version) )
      CALL iotk_write_empty( iunpun, "CREATOR", ATTR=attr )
      !
      CALL iotk_write_end( iunpun, "HEADER" )
      !
    END SUBROUTINE write_header
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_control( pp_check_flag, lkpoint_dir, q_real_space, beta_real_space) 
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      LOGICAL, OPTIONAL, INTENT(IN) :: pp_check_flag, lkpoint_dir, q_real_space, beta_real_space


      CALL iotk_write_begin( iunpun, "CONTROL" )
      !
      !  This flag is used to check if the file can be used for post-processing
      IF ( PRESENT( pp_check_flag ) ) &
         CALL iotk_write_dat( iunpun, "PP_CHECK_FLAG", pp_check_flag )
      !
      !  This flag says how eigenvalues are saved
      IF ( PRESENT( lkpoint_dir ) ) &
         CALL iotk_write_dat( iunpun, "LKPOINT_DIR", lkpoint_dir )
      !
      !  This flag says if Q in real space has to be used
      IF ( PRESENT( q_real_space ) ) &
         CALL iotk_write_dat( iunpun, "Q_REAL_SPACE", q_real_space )
      ! This flag says if Beta functions were treated in real space
      IF ( PRESENT( beta_real_space ) ) &
         CALL iotk_write_dat( iunpun, "BETA_REAL_SPACE", beta_real_space )
      !
      CALL iotk_write_end( iunpun, "CONTROL" )
      !
    END SUBROUTINE write_control
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_cell( ibrav, celldm, alat, a1, a2, a3, b1, b2, b3, &
                           do_mp, do_mt, do_esm )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(IN) :: ibrav
      REAL(DP),         INTENT(IN) :: celldm(6), alat
      REAL(DP),         INTENT(IN) :: a1(3), a2(3), a3(3)
      REAL(DP),         INTENT(IN) :: b1(3), b2(3), b3(3)
      LOGICAL,          INTENT(IN) :: do_mp, do_mt, do_esm
      !
      CHARACTER(LEN=256) :: bravais_lattice,es_corr
      !
      CALL iotk_write_begin( iunpun, "CELL" )
      !
      SELECT CASE ( ibrav )
        CASE(  0 )
           bravais_lattice = "free"
        CASE(  1 )
           bravais_lattice = "cubic P (sc)"
        CASE(  2 )
           bravais_lattice = "cubic F (fcc)"
        CASE(  3 )
           bravais_lattice = "cubic I (bcc)"
        CASE(  4 )
           bravais_lattice = "Hexagonal and Trigonal P"
        CASE(  5 )
           bravais_lattice = "Trigonal R"
        CASE(  6 )
           bravais_lattice = "Tetragonal P (st)"
        CASE(  7 )
           bravais_lattice = "Tetragonal I (bct)"
        CASE(  8 )
           bravais_lattice = "Orthorhombic P"
        CASE(  9 )
           bravais_lattice = "Orthorhombic base-centered(bco)"
        CASE( 10 )
           bravais_lattice = "Orthorhombic face-centered"
        CASE( 11 )
           bravais_lattice = "Orthorhombic body-centered"
        CASE( 12 )
           bravais_lattice = "Monoclinic P"
        CASE( 13 )
           bravais_lattice = "Monoclinic base-centered"
        CASE( 14 )
           bravais_lattice = "Triclinic P"
      END SELECT
      !
      IF(do_mp)THEN
        es_corr = "Makov-Payne"
      ELSE IF(do_mt) THEN
        es_corr = "Martyna-Tuckerman"
      ELSE IF(do_esm) THEN
        es_corr = "ESM"
      ELSE
        es_corr = "None"
      ENDIF
      CALL iotk_write_dat( iunpun, &
                           "NON-PERIODIC_CELL_CORRECTION", TRIM( es_corr ) )
         
      CALL iotk_write_dat( iunpun, &
                           "BRAVAIS_LATTICE", TRIM( bravais_lattice ) )
      !
      CALL iotk_write_attr( attr, "UNITS", "Bohr", FIRST = .TRUE. )
      CALL iotk_write_dat( iunpun, "LATTICE_PARAMETER", alat, ATTR = attr )
      !
      CALL iotk_write_dat( iunpun, "CELL_DIMENSIONS", celldm(1:6) )
      !
      CALL iotk_write_begin( iunpun, "DIRECT_LATTICE_VECTORS" )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_DIRECT_LATTICE_VECTORS", ATTR=attr )
      CALL iotk_write_dat(   iunpun, "a1", a1(:) * alat, COLUMNS=3 )
      CALL iotk_write_dat(   iunpun, "a2", a2(:) * alat, COLUMNS=3 )
      CALL iotk_write_dat(   iunpun, "a3", a3(:) * alat, COLUMNS=3 )
      CALL iotk_write_end(   iunpun, "DIRECT_LATTICE_VECTORS" )
      !
      CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
      CALL iotk_write_begin( iunpun, "RECIPROCAL_LATTICE_VECTORS" )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_RECIPROCAL_LATTICE_VECTORS", ATTR=attr )
      CALL iotk_write_dat(   iunpun, "b1", b1(:), COLUMNS=3 )
      CALL iotk_write_dat(   iunpun, "b2", b2(:), COLUMNS=3 )
      CALL iotk_write_dat(   iunpun, "b3", b3(:), COLUMNS=3 )
      CALL iotk_write_end(   iunpun, "RECIPROCAL_LATTICE_VECTORS" )
      !
      CALL iotk_write_end( iunpun, "CELL" )
      !
    END SUBROUTINE write_cell
    !
    SUBROUTINE write_moving_cell(lmovecell, cell_factor)

    LOGICAL, INTENT(IN) :: lmovecell
    REAL(DP), INTENT(IN) :: cell_factor

    CALL iotk_write_begin( iunpun, "MOVING_CELL" )
    CALL iotk_write_dat( iunpun, "CELL_FACTOR", cell_factor)
    CALL iotk_write_end( iunpun, "MOVING_CELL"  )

    RETURN
    END SUBROUTINE write_moving_cell
    !------------------------------------------------------------------------
    SUBROUTINE write_ions( nsp, nat, atm, ityp, psfile, &
                           pseudo_dir, amass, tau, if_pos, dirname, pos_unit )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(IN) :: nsp, nat
      INTEGER,          INTENT(IN) :: ityp(:)
      CHARACTER(LEN=*), INTENT(IN) :: atm(:)
      CHARACTER(LEN=*), INTENT(IN) :: psfile(:)
      CHARACTER(LEN=*), INTENT(IN) :: pseudo_dir
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      REAL(DP),         INTENT(IN) :: amass(:)
      REAL(DP),         INTENT(IN) :: tau(:,:)
      INTEGER,          INTENT(IN) :: if_pos(:,:)
      REAL(DP),         INTENT(IN) :: pos_unit
      !
      INTEGER            :: i, flen, flen2
      LOGICAL            :: exst
      CHARACTER(LEN=256) :: file_pseudo_in, file_pseudo_out
      !
      !
      CALL iotk_write_begin( iunpun, "IONS" )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_ATOMS", nat )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_SPECIES", nsp )
      !
      flen = LEN_TRIM( pseudo_dir )
      flen2= LEN_TRIM( dirname )
      !
      CALL iotk_write_attr ( attr, "UNITS", "a.m.u.", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_ATOMIC_MASSES", ATTR = attr )
      !
      DO i = 1, nsp
         !
         CALL iotk_write_begin( iunpun, "SPECIE"//TRIM(iotk_index(i)) )
         !
         CALL iotk_write_dat( iunpun, "ATOM_TYPE", atm(i) )
         !
         ! ... Copy PP file into the .save directory - verify that the
         ! ...  sourcefile exists and is not the same as the destination file
         !
         IF ( pseudo_dir(flen:flen) /= '/' ) THEN
            file_pseudo_in = pseudo_dir(1:flen) // '/' // TRIM(psfile(i))
         ELSE
            file_pseudo_in = pseudo_dir(1:flen) // TRIM(psfile(i))
         END IF
         !
         IF ( dirname(flen2:flen2) /= '/' ) THEN
            file_pseudo_out = dirname(1:flen2) // '/' // TRIM(psfile(i))
         ELSE
            file_pseudo_out = dirname(1:flen2) // TRIM(psfile(i))
         END IF
         !
         IF ( file_pseudo_in .ne. file_pseudo_out ) THEN
            INQUIRE ( FILE=file_pseudo_in, EXIST=exst )
            IF ( exst ) THEN
               CALL copy_file( TRIM( file_pseudo_in ), TRIM( file_pseudo_out ) )
            ELSE
               CALL infomsg( 'write_ions', &
                   'file ' // TRIM( file_pseudo_in) // ' not present' )
            END IF
         END IF
         !
         CALL iotk_write_dat( iunpun, "MASS", amass(i) )
         !
         CALL iotk_write_dat( iunpun, "PSEUDO", TRIM( psfile(i) ) )
         !
         CALL iotk_write_end( iunpun, "SPECIE"//TRIM(iotk_index(i)) )
         !
      ENDDO
      !
      ! This is the original location where PP files are read from
      !
      CALL iotk_write_dat( iunpun, "PSEUDO_DIR", TRIM( pseudo_dir) )
      !
      CALL iotk_write_attr( attr, "UNITS", "Bohr", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_ATOMIC_POSITIONS", ATTR = attr )
      !
      DO i = 1, nat
         !
         CALL iotk_write_attr( attr, "SPECIES", &
                             & atm( ityp(i) ), FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "INDEX",  ityp(i) )                     
         CALL iotk_write_attr( attr, "tau",    tau(:,i)*pos_unit )
         CALL iotk_write_attr( attr, "if_pos", if_pos(:,i) )
         CALL iotk_write_empty( iunpun, &
                              & "ATOM" // TRIM( iotk_index( i ) ), attr )
         !
      END DO
      !
      CALL iotk_write_end( iunpun, "IONS" )
      !
    END SUBROUTINE write_ions
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_symmetry( ibrav, nrot, nsym, invsym, noinv, &
                               time_reversal, no_t_rev, ft, &
                               s, sname, irt, nat, t_rev )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(IN) :: ibrav, nrot, nsym
      LOGICAL,          INTENT(IN) :: invsym, noinv, time_reversal, no_t_rev
      INTEGER,          INTENT(IN) :: s(:,:,:), irt(:,:), nat, t_rev(:)
      REAL(DP),         INTENT(IN) :: ft(:,:)
      CHARACTER(LEN=*), INTENT(IN) :: sname(:)
      !
      INTEGER  :: i
      !
      !
      CALL iotk_write_begin( iunpun, "SYMMETRIES" )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_SYMMETRIES", nsym )
      CALL iotk_write_dat( iunpun, "NUMBER_OF_BRAVAIS_SYMMETRIES", nrot )
      !
      CALL iotk_write_dat( iunpun, "INVERSION_SYMMETRY", invsym )
      !
      CALL iotk_write_dat( iunpun, "DO_NOT_USE_TIME_REVERSAL", noinv )
      !
      CALL iotk_write_dat( iunpun, "TIME_REVERSAL_FLAG", time_reversal )
      !
      CALL iotk_write_dat( iunpun, "NO_TIME_REV_OPERATIONS", no_t_rev )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_ATOMS", nat )
      !
      CALL iotk_write_attr( attr, "UNITS", "Crystal", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_SYMMETRIES", ATTR = attr )
      !
      DO i = 1, nsym
         !
         CALL iotk_write_begin( iunpun, "SYMM" // TRIM( iotk_index( i ) ) )
         !
         CALL iotk_write_attr ( attr, "NAME", TRIM( sname(i) ), FIRST=.TRUE. )
         CALL iotk_write_attr ( attr, "T_REV", t_rev(i) )
         CALL iotk_write_empty( iunpun, "INFO", ATTR = attr )
         !
         CALL iotk_write_dat( iunpun, "ROTATION", s(:,:,i), COLUMNS=3 )
         CALL iotk_write_dat( iunpun, "FRACTIONAL_TRANSLATION", ft(:,i), COLUMNS=3 )
         !
         IF ( nat > 0 ) &
            CALL iotk_write_dat( iunpun, "EQUIVALENT_IONS", irt(i,1:nat), COLUMNS=8 )
         !
         CALL iotk_write_end( iunpun, "SYMM" // TRIM( iotk_index( i ) ) )
         !
      ENDDO
      !
      ! ... the following are the symmetries of the Bravais lattice alone
      ! ... (they may be more than crystal, i.e. basis+lattice, symmetries)
      !
      DO i = nsym+1, nrot
         !
         CALL iotk_write_begin( iunpun, "SYMM" // TRIM( iotk_index( i ) ) )
         !
         CALL iotk_write_attr ( attr, "NAME", TRIM( sname(i) ), FIRST=.TRUE. )
         CALL iotk_write_empty( iunpun, "INFO", ATTR = attr )
         CALL iotk_write_dat( iunpun, "ROTATION", s(:,:,i), COLUMNS=3 )
         !
         CALL iotk_write_end( iunpun, "SYMM" // TRIM( iotk_index( i ) ) )
         !
      ENDDO
      !
      CALL iotk_write_end( iunpun, "SYMMETRIES" )
      !
    END SUBROUTINE write_symmetry
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_efield( tefield, dipfield, edir, emaxpos, eopreg, eamp )
      !------------------------------------------------------------------------
      !
      LOGICAL, INTENT(IN) :: &
           tefield,      &! if .TRUE. a finite electric field is added to the
                          ! local potential
           dipfield       ! if .TRUE. the dipole field is subtracted
      INTEGER, INTENT(IN) :: &
           edir           ! direction of the field
      REAL(DP), INTENT(IN) :: &
           emaxpos,  &! position of the maximum of the field (0<emaxpos<1)
           eopreg,   &! amplitude of the inverse region (0<eopreg<1)
           eamp       ! field amplitude (in a.u.) (1 a.u. = 51.44 10^11 V/m)
      !
      !
      CALL iotk_write_begin( iunpun, "ELECTRIC_FIELD" )
      !
      CALL iotk_write_dat( iunpun, "HAS_ELECTRIC_FIELD", tefield )
      !
      CALL iotk_write_dat( iunpun, "HAS_DIPOLE_CORRECTION", dipfield )
      !
      CALL iotk_write_dat( iunpun, "FIELD_DIRECTION", edir )
      !
      CALL iotk_write_dat( iunpun, "MAXIMUM_POSITION", emaxpos )
      !
      CALL iotk_write_dat( iunpun, "INVERSE_REGION", eopreg )
      !
      CALL iotk_write_dat( iunpun, "FIELD_AMPLITUDE", eamp )
      !
      CALL iotk_write_end( iunpun, "ELECTRIC_FIELD" )
      !
    END SUBROUTINE write_efield
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_planewaves( ecutwfc, dual, npwx, gamma_only, nr1, nr2, &
                                 nr3, ngm_g, nr1s, nr2s, nr3s, ngms_g, nr1b, &
                                 nr2b, nr3b, itmp, lgvec )
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : e2
      !
      INTEGER,  INTENT(IN) :: npwx, nr1, nr2, nr3, ngm_g, &
                              nr1s, nr2s, nr3s, ngms_g, nr1b, nr2b, nr3b
      INTEGER,  INTENT(IN) :: itmp(:,:)
      REAL(DP), INTENT(IN) :: ecutwfc, dual
      LOGICAL,  INTENT(IN) :: gamma_only, lgvec
      !
      !
      CALL iotk_write_begin( iunpun, "PLANE_WAVES" )
      !
      CALL iotk_write_attr ( attr, "UNITS", "Hartree", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_CUTOFF", ATTR = attr )
      !
      CALL iotk_write_dat( iunpun, "WFC_CUTOFF", ecutwfc / e2 )
      !
      CALL iotk_write_dat( iunpun, "RHO_CUTOFF", ecutwfc * dual / e2 )
      !
      CALL iotk_write_dat( iunpun, "MAX_NUMBER_OF_GK-VECTORS", npwx )
      !
      CALL iotk_write_dat( iunpun, "GAMMA_ONLY", gamma_only )
      !
      CALL iotk_write_attr( attr, "nr1", nr1, FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "nr2", nr2 )
      CALL iotk_write_attr( attr, "nr3", nr3 )
      CALL iotk_write_empty( iunpun, "FFT_GRID", ATTR = attr )
      !
      CALL iotk_write_dat( iunpun, "GVECT_NUMBER", ngm_g )
      !
      CALL iotk_write_attr( attr, "nr1s", nr1s, FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "nr2s", nr2s )
      CALL iotk_write_attr( attr, "nr3s", nr3s )
      CALL iotk_write_empty( iunpun, "SMOOTH_FFT_GRID", ATTR = attr )
      !
      CALL iotk_write_dat( iunpun, "SMOOTH_GVECT_NUMBER", ngms_g )
      !
      IF ( lgvec ) THEN
         !
         ! ... write the G-vectors
         !
         CALL iotk_link( iunpun, "G-VECTORS", &
                         "./gvectors.dat", CREATE = .TRUE., BINARY = .TRUE. )
         !
         CALL iotk_write_begin( iunpun, "G-VECTORS" )
         !
         CALL iotk_write_attr( attr, "nr1s", nr1s, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "nr2s", nr2s )
         CALL iotk_write_attr( attr, "nr3s", nr3s )
         CALL iotk_write_attr( attr, "gvect_number", ngm_g )
         CALL iotk_write_attr( attr, "gamma_only", gamma_only )
         CALL iotk_write_attr( attr, "units", "crystal" )
         CALL iotk_write_empty( iunpun, "INFO", ATTR = attr )
         !
         CALL iotk_write_dat  ( iunpun, "g", itmp(1:3,1:ngm_g), COLUMNS = 3 )
         CALL iotk_write_end  ( iunpun, "G-VECTORS" )
         !
      END IF
      !
      CALL iotk_write_attr( attr, "nr1b", nr1b , FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "nr2b", nr2b )
      CALL iotk_write_attr( attr, "nr3b", nr3b )
      CALL iotk_write_empty( iunpun, "SMALLBOX_FFT_GRID", ATTR = attr )
      !
      CALL iotk_write_end( iunpun, "PLANE_WAVES" )
      !
    END SUBROUTINE write_planewaves
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_spin( lsda, noncolin, npol, lspinorb, domag )
      !------------------------------------------------------------------------
      !
      LOGICAL, INTENT(IN) :: lsda, noncolin, lspinorb, domag
      INTEGER, INTENT(IN) :: npol
      !
      !
      CALL iotk_write_begin( iunpun, "SPIN" )
      !
      CALL iotk_write_dat( iunpun, "LSDA", lsda )
      !
      CALL iotk_write_dat( iunpun, "NON-COLINEAR_CALCULATION", noncolin )
      !
      IF ( noncolin ) &
         CALL iotk_write_dat( iunpun, "SPINOR_DIM", npol )
      !
      CALL iotk_write_dat( iunpun, "SPIN-ORBIT_CALCULATION", lspinorb )
      CALL iotk_write_dat( iunpun, "SPIN-ORBIT_DOMAG", domag )
      !
      CALL iotk_write_end( iunpun, "SPIN" )
      !
    END SUBROUTINE write_spin
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_magnetization(starting_magnetization, angle1, angle2, &
                                   nsp, two_fermi_energies, i_cons, mcons, bfield, &
                                   ef_up, ef_dw, nelup, neldw, lambda)
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : PI
      !
      IMPLICIT NONE
      INTEGER,  INTENT(IN) :: nsp, i_cons
      REAL(DP), INTENT(IN) :: starting_magnetization(nsp), &
                              angle1(nsp), angle2(nsp), mcons(3,nsp), &
                              bfield(3), ef_up, ef_dw, nelup, neldw, lambda
      LOGICAL,  INTENT(IN) :: two_fermi_energies
      !
      INTEGER :: i
      !
      CALL iotk_write_begin( iunpun, "MAGNETIZATION_INIT" )

      CALL iotk_write_dat(iunpun,"CONSTRAINT_MAG", i_cons)

      CALL iotk_write_dat( iunpun, "NUMBER_OF_SPECIES", nsp ) 

      DO i = 1, nsp
         !
         CALL iotk_write_begin( iunpun, "SPECIE"//TRIM(iotk_index(i)) )
         !
         CALL iotk_write_dat( iunpun, "STARTING_MAGNETIZATION",  &
                                       starting_magnetization(i) )
         CALL iotk_write_dat( iunpun, "ANGLE1", &
                                       angle1(i)*180.0_DP/PI )
         CALL iotk_write_dat( iunpun, "ANGLE2", &
                                       angle2(i)*180.0_DP/PI )
         IF (i_cons==1.OR.i_cons==2) THEN
            CALL iotk_write_dat( iunpun, "CONSTRANT_1", mcons(1,i) )
            CALL iotk_write_dat( iunpun, "CONSTRANT_2", mcons(2,i) )
            CALL iotk_write_dat( iunpun, "CONSTRANT_3", mcons(3,i) )
         ENDIF
         !
         CALL iotk_write_end( iunpun, "SPECIE"//TRIM(iotk_index(i)) )
         !
      ENDDO

      IF (i_cons==3) THEN
         !
         CALL iotk_write_dat( iunpun, "FIXED_MAGNETIZATION_1", mcons(1,1) )
         CALL iotk_write_dat( iunpun, "FIXED_MAGNETIZATION_2", mcons(2,1) )
         CALL iotk_write_dat( iunpun, "FIXED_MAGNETIZATION_3", mcons(3,1) )
         !
      ELSE IF (i_cons==4) THEN
         !
         CALL iotk_write_dat( iunpun, "MAGNETIC_FIELD_1", bfield(1) )
         CALL iotk_write_dat( iunpun, "MAGNETIC_FIELD_2", bfield(2) )
         CALL iotk_write_dat( iunpun, "MAGNETIC_FIELD_3", bfield(3) )
         !
      ENDIF
      !
      CALL iotk_write_dat(iunpun,"TWO_FERMI_ENERGIES",two_fermi_energies)
      !
      IF (two_fermi_energies) THEN
         !
         CALL iotk_write_attr ( attr, "UNITS", "Hartree", FIRST = .TRUE. )
         CALL iotk_write_empty( iunpun, "UNITS_FOR_ENERGIES", ATTR = attr )
         !
         CALL iotk_write_dat( iunpun, "FIXED_MAGNETIZATION", mcons(3,1) )
         CALL iotk_write_dat( iunpun, "ELECTRONS_UP", nelup )
         CALL iotk_write_dat( iunpun, "ELECTRONS_DOWN", neldw )
         CALL iotk_write_dat( iunpun, "FERMI_ENERGY_UP", ef_up )
         CALL iotk_write_dat( iunpun, "FERMI_ENERGY_DOWN", ef_dw )
         !
      ENDIF
      !
      IF (i_cons>0) CALL iotk_write_dat(iunpun,"LAMBDA",lambda)
      !
      CALL iotk_write_end( iunpun, "MAGNETIZATION_INIT" )
      !
    RETURN
    !
    END SUBROUTINE write_magnetization
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_xc( dft, nsp, lda_plus_u, lda_plus_u_kind, U_projection, Hubbard_lmax, &
                         Hubbard_l, Hubbard_U, Hubbard_J, Hubbard_J0, Hubbard_beta, &
                         Hubbard_alpha, inlc, vdw_table_name, pseudo_dir, dirname )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=*),   INTENT(IN) :: dft
      LOGICAL,            INTENT(IN) :: lda_plus_u
      INTEGER,  OPTIONAL, INTENT(IN) :: lda_plus_u_kind
      INTEGER,  OPTIONAL, INTENT(IN) :: nsp
      CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: U_projection
      INTEGER,  OPTIONAL, INTENT(IN) :: Hubbard_lmax
      INTEGER,  OPTIONAL, INTENT(IN) :: Hubbard_l(:)
      REAL(DP), OPTIONAL, INTENT(IN) :: Hubbard_U(:), Hubbard_J(:,:), Hubbard_alpha(:), &
                                        Hubbard_J0(:), Hubbard_beta(:)
      INTEGER,  OPTIONAL, INTENT(IN) :: inlc
      CHARACTER(LEN=*), OPTIONAL,   INTENT(IN) :: vdw_table_name, pseudo_dir, dirname
      !
      INTEGER            :: i, flen
      CHARACTER(LEN=256) :: file_table
      !
      CALL iotk_write_begin( iunpun, "EXCHANGE_CORRELATION" )
      !
      CALL iotk_write_dat( iunpun, "DFT", dft )
      !
      CALL iotk_write_dat( iunpun, "LDA_PLUS_U_CALCULATION", lda_plus_u )
      !
      IF ( lda_plus_u ) THEN
         !
         IF ( .NOT. PRESENT( Hubbard_lmax ) .OR. &
              .NOT. PRESENT( Hubbard_l )    .OR. & 
              .NOT. PRESENT( U_projection ) .OR. &
              .NOT. PRESENT( Hubbard_U )    .OR. &
              .NOT. PRESENT( Hubbard_J0 )    .OR. &
              .NOT. PRESENT( Hubbard_J )    .OR. &
              .NOT. PRESENT( nsp )          .OR. &
              .NOT. PRESENT( Hubbard_alpha )    .OR. &
              .NOT. PRESENT( Hubbard_beta ) ) &
            CALL errore( 'write_xc', &
                         ' variables for LDA+U not present', 1 )
         !
         CALL iotk_write_dat( iunpun, "NUMBER_OF_SPECIES", nsp )
         CALL iotk_write_dat( iunpun, "HUBBARD_LMAX", Hubbard_lmax )
         CALL iotk_write_dat( iunpun, "HUBBARD_L",  Hubbard_l(1:nsp) )
         CALL iotk_write_dat( iunpun, "HUBBARD_U", Hubbard_U(1:nsp) )
         CALL iotk_write_dat( iunpun, "HUBBARD_ALPHA", Hubbard_alpha(1:nsp) )
         CALL iotk_write_dat( iunpun, "U_PROJECTION_TYPE", trim(U_projection) )
         !
         CALL iotk_write_dat( iunpun, "LDA_PLUS_U_KIND", lda_plus_u_kind )
         CALL iotk_write_dat( iunpun, "HUBBARD_J", Hubbard_J(1:3,1:nsp), COLUMNS = 3 )
         CALL iotk_write_dat( iunpun, "HUBBARD_J0", Hubbard_J0(1:nsp) )
         CALL iotk_write_dat( iunpun, "HUBBARD_BETA", Hubbard_beta(1:nsp) )
         !
      END IF

      !
      ! Vdw kernel table
      !
      CALL iotk_write_dat( iunpun, "NON_LOCAL_DF", inlc )

      IF ( inlc == 1 .OR. inlc ==2 ) THEN
          
         IF ( .NOT. PRESENT( vdw_table_name ) .OR. &
              .NOT. PRESENT( pseudo_dir ) .OR. &
              .NOT. PRESENT( dirname )) &
            CALL errore( 'write_xc', &
                         ' variable vdw_table_name not present', 1 )
        
         CALL iotk_write_dat( iunpun, "VDW_KERNEL_NAME", vdw_table_name )

         !
         ! Copy the file in .save directory
         !
         flen = LEN_TRIM( pseudo_dir )
         IF ( pseudo_dir(flen:flen) /= '/' ) THEN
             !
             file_table = pseudo_dir(1:flen) // '/' // vdw_table_name
             !
         ELSE
             !
             file_table = pseudo_dir(1:flen) // vdw_table_name
             !
         END IF
         !
         CALL copy_file( TRIM( file_table ), TRIM( dirname ) // "/" // TRIM( vdw_table_name ) )
 
      END IF

      !
      CALL iotk_write_end( iunpun, "EXCHANGE_CORRELATION" )
      !
    END SUBROUTINE write_xc

    !------------------------------------------------------------------------
    SUBROUTINE write_exx( x_gamma_extrapolation, nqx1, nqx2, nqx3, &
                          exxdiv_treatment, yukawa, ecutvcut, exx_fraction, &
                          screening_parameter, exx_is_active )
      !------------------------------------------------------------------------
      !
      LOGICAL,            INTENT(IN) :: x_gamma_extrapolation, exx_is_active
      INTEGER,  OPTIONAL, INTENT(IN) :: nqx1, nqx2, nqx3
      CHARACTER(LEN=*),   INTENT(IN) :: exxdiv_treatment
      REAL(DP),           INTENT(IN) :: yukawa, ecutvcut, exx_fraction
      REAL(DP),           INTENT(IN) :: screening_parameter

      CALL iotk_write_begin(iunpun, "EXACT_EXCHANGE" )
      call iotk_write_dat(iunpun, "x_gamma_extrapolation", x_gamma_extrapolation)
      call iotk_write_dat(iunpun, "nqx1", nqx1)
      call iotk_write_dat(iunpun, "nqx2", nqx2)
      call iotk_write_dat(iunpun, "nqx3", nqx3)
      call iotk_write_dat(iunpun, "exxdiv_treatment", exxdiv_treatment)
      call iotk_write_dat(iunpun, "yukawa", yukawa)
      call iotk_write_dat(iunpun, "ecutvcut", ecutvcut)
      call iotk_write_dat(iunpun, "exx_fraction", exx_fraction)
      call iotk_write_dat(iunpun, "screening_parameter", screening_parameter)
      call iotk_write_dat(iunpun, "exx_is_active", exx_is_active)
      CALL iotk_write_end(iunpun, "EXACT_EXCHANGE" )
    END SUBROUTINE write_exx
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_occ( lgauss, ngauss, degauss, ltetra, ntetra, &
                          tetra, tfixed_occ, lsda, nstates_up, nstates_down, f_inp )
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : e2
      !
      LOGICAL,            INTENT(IN) :: lgauss, ltetra, tfixed_occ, lsda
      INTEGER,  OPTIONAL, INTENT(IN) :: ngauss, ntetra, nstates_up, nstates_down
      INTEGER,  OPTIONAL, INTENT(IN) :: tetra(:,:)
      REAL(DP), OPTIONAL, INTENT(IN) :: degauss, f_inp(:,:)      
      !
      INTEGER :: i
      !
      !
      CALL iotk_write_begin( iunpun, "OCCUPATIONS" )
      !
      CALL iotk_write_dat( iunpun, "SMEARING_METHOD", lgauss )
      !
      IF ( lgauss ) THEN
         !
         CALL iotk_write_dat( iunpun, "SMEARING_TYPE", ngauss )
         !
         CALL iotk_write_attr( attr, "UNITS", "Hartree", FIRST = .TRUE. )
         !
         CALL iotk_write_dat( iunpun, "SMEARING_PARAMETER", &
                              degauss / e2, ATTR = attr )
         !
      END IF
      !
      CALL iotk_write_dat( iunpun, "TETRAHEDRON_METHOD", ltetra )
      !
      IF ( ltetra ) THEN
         !
         CALL iotk_write_dat( iunpun, "NUMBER_OF_TETRAHEDRA", ntetra )
         !
         DO i = 1, ntetra
            !
            CALL iotk_write_dat( iunpun, "TETRAHEDRON" // &
                               & iotk_index( i ), tetra(1:4,i) )
            !
         END DO
         !
      END IF
      !
      CALL iotk_write_dat( iunpun, "FIXED_OCCUPATIONS", tfixed_occ )
      !
      IF ( tfixed_occ ) THEN
         !
         CALL iotk_write_attr( attr, "lsda" , lsda, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "nstates_up", nstates_up )
         CALL iotk_write_attr( attr, "nstates_down", nstates_down )
         !
         CALL iotk_write_empty( iunpun, 'INFO', ATTR = attr )
         !
         CALL iotk_write_dat( iunpun, "INPUT_OCC_UP", f_inp(1:nstates_up,1) )
         !
         IF ( lsda ) &
            CALL iotk_write_dat( iunpun, "INPUT_OCC_DOWN", f_inp(1:nstates_down,2) )
         !
      END IF
      !
      CALL iotk_write_end( iunpun, "OCCUPATIONS" )
      !
    END SUBROUTINE write_occ
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_bz( num_k_points, xk, wk, k1, k2, k3, nk1, nk2, nk3, &
                         qnorm, nks_start, xk_start, wk_start )
      !------------------------------------------------------------------------
      !
      INTEGER,  INTENT(IN) :: num_k_points, k1, k2, k3, nk1, nk2, nk3
      REAL(DP), INTENT(IN) :: xk(:,:), wk(:)
      REAL(DP), INTENT(IN) :: qnorm
      INTEGER,  INTENT(IN), OPTIONAL ::  nks_start
      REAL(DP), INTENT(IN), OPTIONAL :: xk_start(:,:), wk_start(:)
      !
      INTEGER :: ik, i
      !
      !
      CALL iotk_write_begin( iunpun, "BRILLOUIN_ZONE" )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_K-POINTS", num_k_points )
      !
      CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_K-POINTS", attr )
      !
      CALL iotk_write_attr( attr, "nk1", nk1, FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "nk2", nk2 )
      CALL iotk_write_attr( attr, "nk3", nk3 )
      CALL iotk_write_empty( iunpun, "MONKHORST_PACK_GRID", attr )
      CALL iotk_write_attr( attr, "k1", k1, FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "k2", k2 )
      CALL iotk_write_attr( attr, "k3", k3 )
      CALL iotk_write_empty( iunpun, "MONKHORST_PACK_OFFSET", attr )
      !
      DO ik = 1, num_k_points
         !
         CALL iotk_write_attr( attr, "XYZ", xk(:,ik), FIRST = .TRUE. )
         !            
         CALL iotk_write_attr( attr, "WEIGHT", wk(ik) )
         !
         CALL iotk_write_empty( iunpun, "K-POINT" // &
                              & TRIM( iotk_index(ik) ), attr )
         !
      END DO
      !
      ! ... these are k-points and weights in the Irreducible BZ
      !
      IF (present(nks_start).and.present(xk_start).and.present(wk_start)) THEN
         !
         CALL iotk_write_dat( iunpun, "STARTING_K-POINTS", nks_start )
         !
         DO ik = 1, nks_start
            !
            CALL iotk_write_attr( attr, "XYZ", xk_start(:,ik), FIRST = .TRUE. )
            !            
            CALL iotk_write_attr( attr, "WEIGHT", wk_start(ik) )
            !
            CALL iotk_write_empty( iunpun, "K-POINT_START" // &
                              & TRIM( iotk_index(ik) ), attr )
            !
         END DO
      ENDIF
      !
      CALL iotk_write_dat( iunpun, "NORM-OF-Q", qnorm )
      !
      CALL iotk_write_end( iunpun, "BRILLOUIN_ZONE" )
      !
    END SUBROUTINE write_bz
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_para( kunit, nproc, nproc_pool, nproc_image, &
                    ntask_groups, nproc_pot, nproc_bgrp, nproc_ortho ) 
      !------------------------------------------------------------------------
      !
      INTEGER,  INTENT(IN) :: kunit, nproc, nproc_pool, nproc_image, &
                              ntask_groups, nproc_pot, nproc_bgrp, nproc_ortho 
      !
      !
      CALL iotk_write_begin( iunpun, "PARALLELISM" )
      CALL iotk_write_dat( iunpun, &
                              "GRANULARITY_OF_K-POINTS_DISTRIBUTION", kunit )
      CALL iotk_write_dat( iunpun, "NUMBER_OF_PROCESSORS", nproc )
      CALL iotk_write_dat( iunpun, &
                              "NUMBER_OF_PROCESSORS_PER_POOL", nproc_pool )
      CALL iotk_write_dat( iunpun, &
                              "NUMBER_OF_PROCESSORS_PER_IMAGE", nproc_image )
      CALL iotk_write_dat( iunpun, "NUMBER_OF_PROCESSORS_PER_TASKGROUP", &
                                              ntask_groups )
      CALL iotk_write_dat( iunpun, "NUMBER_OF_PROCESSORS_PER_POT", &
                                              nproc_pot )
      CALL iotk_write_dat( iunpun, "NUMBER_OF_PROCESSORS_PER_BAND_GROUP", &
                                              nproc_bgrp )
      CALL iotk_write_dat( iunpun, "NUMBER_OF_PROCESSORS_PER_DIAGONALIZATION", &
                                              nproc_ortho )
      CALL iotk_write_end( iunpun, "PARALLELISM" )
      !
      !
    END SUBROUTINE write_para
    !------------------------------------------------------------------------
    SUBROUTINE write_phonon( modenum, xqq )
      !------------------------------------------------------------------------
      !
      INTEGER,  INTENT(IN) :: modenum
      REAL(DP), INTENT(IN) :: xqq(:)
      !
      !
      CALL iotk_write_begin( iunpun, "PHONON" )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_MODES", modenum )
      !
      CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_Q-POINT", attr )
      !
      CALL iotk_write_dat( iunpun, "Q-POINT", xqq(:), COLUMNS=3 )
      !
      CALL iotk_write_end( iunpun, "PHONON" )
      !
    END SUBROUTINE write_phonon
    !
    ! ... methods to write and read charge_density
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_rho_xml( rho_file_base, rho, &
                              nr1, nr2, nr3, nr1x, nr2x, ipp, npp, &
                              ionode, intra_group_comm, inter_group_comm )
      !------------------------------------------------------------------------
      !
      ! ... Writes charge density rho, one plane at a time.
      ! ... If ipp and npp are specified, planes are collected one by one from
      ! ... all processors, avoiding an overall collect of the charge density
      ! ... on a single proc.
      !
      USE io_files,  ONLY : rhounit
      USE mp,        ONLY : mp_get, mp_sum, mp_rank, mp_size
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*),  INTENT(IN) :: rho_file_base
      REAL(DP),          INTENT(IN) :: rho(:)
      INTEGER,           INTENT(IN) :: nr1, nr2, nr3
      INTEGER,           INTENT(IN) :: nr1x, nr2x
      INTEGER,           INTENT(IN) :: ipp(:)
      INTEGER,           INTENT(IN) :: npp(:)
      LOGICAL,           INTENT(IN) :: ionode
      INTEGER,           INTENT(IN) :: intra_group_comm, inter_group_comm
      !
      INTEGER               :: ierr, i, j, k, kk, ldr, ip
      CHARACTER(LEN=256)    :: rho_file
      CHARACTER(LEN=10)     :: rho_extension
      REAL(DP), ALLOCATABLE :: rho_plane(:)
      INTEGER,  ALLOCATABLE :: kowner(:)
      INTEGER               :: my_group_id, me_group, nproc_group, io_group_id, io_group
      !
      me_group    = mp_rank( intra_group_comm )
      nproc_group = mp_size( intra_group_comm )
      my_group_id = mp_rank( inter_group_comm )
      !
      rho_extension = '.dat'
      IF ( .NOT. rho_binary ) rho_extension = '.xml'
      !
      rho_file = TRIM( rho_file_base ) // TRIM( rho_extension )
      !
      IF ( ionode ) THEN 
         CALL iotk_open_write( rhounit, FILE = rho_file,  BINARY = rho_binary, IERR = ierr )
         CALL errore( 'write_rho_xml', 'cannot open' // TRIM( rho_file ) // ' file for writing', ierr )
      END IF 
      !
      IF ( ionode ) THEN
         !
         CALL iotk_write_begin( rhounit, "CHARGE-DENSITY" )
         !
         CALL iotk_write_attr( attr, "nr1", nr1, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "nr2", nr2 )
         CALL iotk_write_attr( attr, "nr3", nr3 )
         !
         CALL iotk_write_empty( rhounit, "INFO", attr )
         !
      END IF
      !
      ALLOCATE( rho_plane( nr1*nr2 ) )
      ALLOCATE( kowner( nr3 ) )
      !
      ! ... find the index of the group (pool) that will write rho
      !
      io_group_id = 0
      !
      IF ( ionode ) io_group_id = my_group_id
      !
      CALL mp_sum( io_group_id, intra_group_comm )
      CALL mp_sum( io_group_id, inter_group_comm )
      !
      ! ... find the index of the ionode within its own group (pool)
      !
      io_group = 0
      !
      IF ( ionode ) io_group = me_group
      !
      CALL mp_sum( io_group, intra_group_comm )
      !
      ! ... find out the owner of each "z" plane
      !
      DO ip = 1, nproc_group
         !
         kowner( (ipp(ip)+1):(ipp(ip)+npp(ip)) ) = ip - 1
         !
      END DO
      !
      ldr = nr1x*nr2x
      !
      DO k = 1, nr3
         !
         !  Only one subgroup write the charge density
         !
         IF( ( kowner(k) == me_group ) .AND. ( my_group_id == io_group_id ) ) THEN
            !
            kk = k - ipp( me_group + 1 )
            ! 
            DO j = 1, nr2
               !
               DO i = 1, nr1
                  !
                  rho_plane(i+(j-1)*nr1) = rho(i+(j-1)*nr1x+(kk-1)*ldr)
                  !
               END DO
               !
            END DO
            !
         END IF
         !
         IF ( kowner(k) /= io_group .AND. my_group_id == io_group_id ) &
            CALL mp_get( rho_plane, rho_plane, me_group, io_group, kowner(k), k, intra_group_comm )
         !
         IF ( ionode ) &
            CALL iotk_write_dat( rhounit, "z" // iotk_index( k ), rho_plane )
         !
      END DO
      !
      DEALLOCATE( rho_plane )
      DEALLOCATE( kowner )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_write_end( rhounit, "CHARGE-DENSITY" )
         !
         CALL iotk_close_write( rhounit )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE write_rho_xml
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_rho_xml( rho_file_base, nr1, nr2, nr3, nr1x, nr2x, &
                             ipp, npp, rho )
      !------------------------------------------------------------------------
      !
      ! ... Reads charge density rho, one plane at a time, to avoid 
      ! ... collecting the entire charge density on a single processor
      !
      USE io_files,  ONLY : rhounit
      USE io_global, ONLY : ionode, ionode_id
      USE mp_global, ONLY : intra_bgrp_comm, intra_image_comm
      USE mp,        ONLY : mp_put, mp_sum, mp_rank, mp_size
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*),  INTENT(IN)  :: rho_file_base
      INTEGER,           INTENT(IN)  :: nr1, nr2, nr3
      INTEGER,           INTENT(IN)  :: nr1x, nr2x
      REAL(DP),          INTENT(OUT) :: rho(:)
      INTEGER,           INTENT(IN)  :: ipp(:)
      INTEGER,           INTENT(IN)  :: npp(:)
      !
      INTEGER               :: ierr, i, j, k, kk, ldr, ip
      INTEGER               :: nr( 3 )
      INTEGER               :: me_group, nproc_group
      CHARACTER(LEN=256)    :: rho_file
      REAL(DP), ALLOCATABLE :: rho_plane(:)
      INTEGER,  ALLOCATABLE :: kowner(:)
      LOGICAL               :: exst
      !
      me_group     = mp_rank ( intra_bgrp_comm )
      nproc_group  = mp_size ( intra_bgrp_comm )
      !
      rho_file = TRIM( rho_file_base ) // ".dat"
      exst = check_file_exst( TRIM(rho_file) ) 
      !
      IF ( .NOT. exst ) THEN
          !
          rho_file = TRIM( rho_file_base ) // ".xml"
          exst = check_file_exst( TRIM(rho_file) ) 
          !
      ENDIF
      !
      IF ( .NOT. exst ) CALL errore('read_rho_xml', 'searching for '//TRIM(rho_file), 10)
      !
      IF ( ionode ) THEN
         CALL iotk_open_read( rhounit, FILE = rho_file, IERR = ierr )
         CALL errore( 'read_rho_xml', 'cannot open ' // TRIM( rho_file ) // ' file for reading', ierr )
      END IF
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( rhounit, "CHARGE-DENSITY" )
         !
         CALL iotk_scan_empty( rhounit, "INFO", attr )
         !
         CALL iotk_scan_attr( attr, "nr1", nr(1) )
         CALL iotk_scan_attr( attr, "nr2", nr(2) )
         CALL iotk_scan_attr( attr, "nr3", nr(3) )
         !
         IF ( nr1 /= nr(1) .OR. nr2 /= nr(2) .OR. nr3 /= nr(3) ) &
            CALL errore( 'read_rho_xml', 'dimensions do not match', 1 )
         !
      END IF
      !
      ALLOCATE( rho_plane( nr1*nr2 ) )
      ALLOCATE( kowner( nr3 ) )
      !
      DO ip = 1, nproc_group
         !
         kowner((ipp(ip)+1):(ipp(ip)+npp(ip))) = ip - 1
         !
      END DO
      !
      ldr = nr1x*nr2x
      !
      ! ... explicit initialization to zero is needed because the physical
      ! ... dimensions rho may exceed the true size of the FFT grid 
      !
      rho(:) = 0.0_DP
      !
      DO k = 1, nr3
         !
         ! ... only ionode reads the charge planes
         !
         IF ( ionode ) &
            CALL iotk_scan_dat( rhounit, "z" // iotk_index( k ), rho_plane )
         !
         ! ... planes are sent to the destination processor
         !
         CALL mp_bcast( rho_plane, ionode_id, intra_image_comm )
         !
         IF( kowner(k) == me_group ) THEN
            !
            kk = k - ipp( me_group + 1 )
            DO j = 1, nr2
               DO i = 1, nr1
                  rho(i+(j-1)*nr1x+(kk-1)*ldr) = rho_plane(i+(j-1)*nr1)
               END DO
            END DO
            !
         END IF
         !
      END DO
      !
      DEALLOCATE( rho_plane )
      DEALLOCATE( kowner )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_end( rhounit, "CHARGE-DENSITY" )
         !
         CALL iotk_close_read( rhounit )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE read_rho_xml
    !
    !------------------------------------------------------------------------
    ! ... methods to write and read wavefunctions
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_wfc( iuni, ik, nk, kunit, ispin, nspin, wf0, ngw,   &
                          gamma_only, nbnd, igl, ngwl, filename, scalef, &
                          ionode, root_in_group, intra_group_comm,       &
                          inter_group_comm, parent_group_comm )
      !------------------------------------------------------------------------
      !
      USE mp_wave,    ONLY : mergewf
      USE mp,         ONLY : mp_get, mp_size, mp_rank, mp_sum
      USE control_flags,     ONLY : lwfnscf, lwfpbe0nscf  ! Lingzhu Kong
      !
      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN) :: iuni
      INTEGER,            INTENT(IN) :: ik, nk, kunit, ispin, nspin
      COMPLEX(DP),        INTENT(IN) :: wf0(:,:)
      INTEGER,            INTENT(IN) :: ngw
      LOGICAL,            INTENT(IN) :: gamma_only
      INTEGER,            INTENT(IN) :: nbnd
      INTEGER,            INTENT(IN) :: ngwl
      INTEGER,            INTENT(IN) :: igl(:)
      CHARACTER(LEN=256), INTENT(IN) :: filename
      REAL(DP),           INTENT(IN) :: scalef    
        ! scale factor, usually 1.0 for pw and 1/SQRT( omega ) for CP
      LOGICAL,            INTENT(IN) :: ionode
      INTEGER,            INTENT(IN) :: root_in_group, intra_group_comm, inter_group_comm, parent_group_comm
      !
      INTEGER                  :: j
      INTEGER                  :: iks, ike, ikt, igwx
      INTEGER                  :: ngroup, ipsour
      INTEGER,     ALLOCATABLE :: ipmask(:)
      INTEGER                  :: me_in_group, nproc_in_group, io_in_parent, nproc_in_parent, me_in_parent
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
      !
      ngroup          = mp_size( inter_group_comm )
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      me_in_parent    = mp_rank( parent_group_comm )
      nproc_in_parent = mp_size( parent_group_comm )
      !
      ALLOCATE( ipmask( nproc_in_parent ) )
      !
      io_in_parent = 0
      IF( ionode ) io_in_parent = me_in_parent
      CALL mp_sum( io_in_parent, parent_group_comm )
      !
      CALL set_kpoints_vars( ik, nk, kunit, ngwl, igl, &
                             ngroup, ikt, iks, ike, igwx, ipmask, ipsour, &
                             ionode, root_in_group, intra_group_comm, inter_group_comm, parent_group_comm )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_open_write( iuni, FILE = TRIM( filename ), ROOT="WFC", BINARY = .TRUE. )
         !
         CALL iotk_write_attr( attr, "ngw",          ngw, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "igwx",         igwx )
         CALL iotk_write_attr( attr, "gamma_only",   gamma_only )
         CALL iotk_write_attr( attr, "nbnd",         nbnd )
         CALL iotk_write_attr( attr, "ik",           ik )
         CALL iotk_write_attr( attr, "nk",           nk )
         CALL iotk_write_attr( attr, "ispin",        ispin )
         CALL iotk_write_attr( attr, "nspin",        nspin )
         CALL iotk_write_attr( attr, "scale_factor", scalef )
         !
         CALL iotk_write_empty( iuni, "INFO", attr )
         !
      END IF
      !
      ALLOCATE( wtmp( MAX( igwx, 1 ) ) )
      !
      wtmp = 0.0_DP
      ! Next 3 lines: Lingzhu Kong
      IF ( ( index(filename,'evc0') > 0 ) .and. (lwfnscf .or. lwfpbe0nscf) )THEN
         IF ( ionode ) OPEN(60,file='cp_wf.dat',status='unknown',form='unformatted')
      ENDIF

      DO j = 1, nbnd
         !
         IF ( ngroup > 1 ) THEN
            !
            IF ( ikt >= iks .AND. ikt <= ike ) &      
               CALL mergewf( wf0(:,j), wtmp, ngwl, igl, me_in_group, &
                             nproc_in_group, root_in_group, intra_group_comm )
            !
            IF ( ipsour /= io_in_parent ) &
               CALL mp_get( wtmp, wtmp, me_in_parent, &
                            io_in_parent, ipsour, j, parent_group_comm )
            !
         ELSE
            !
            CALL mergewf( wf0(:,j), wtmp, ngwl, igl, &
                          me_in_parent, nproc_in_parent, io_in_parent, parent_group_comm )
            !
         END IF
         !
         IF ( ionode ) &
            CALL iotk_write_dat( iuni, "evc" // iotk_index( j ), wtmp(1:igwx) )
         ! Next 3 lines : Lingzhu Kong
         IF ( ( index(filename,'evc0') > 0 ) .and. (lwfnscf .or. lwfpbe0nscf) ) THEN
            IF ( ionode ) write(60)wtmp(1:igwx) 
         ENDIF
         !
      END DO
      ! Next 4 lines : Lingzhu Kong
      IF ( ( index(filename,'evc0') > 0 ) .and. (lwfnscf .or. lwfpbe0nscf) )THEN
          IF ( ionode ) close(60)   !Lingzhu Kong
          write(*,*)'done writing evc0'
      ENDIF

      IF ( ionode ) CALL iotk_close_write( iuni )
      !
      DEALLOCATE( wtmp )
      DEALLOCATE( ipmask )
      !
      RETURN
      !
    END SUBROUTINE write_wfc
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_wfc( iuni, ik, nk, kunit, ispin, &
                         nspin, wf, ngw, nbnd, igl, ngwl, filename, scalef, &
                         ionode, root_in_group, intra_group_comm, inter_group_comm, parent_group_comm, &
                         flink )
      !------------------------------------------------------------------------
      !
      USE mp_wave,   ONLY : splitwf
      USE mp,        ONLY : mp_put, mp_size, mp_rank, mp_sum
      !
      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN)    :: iuni
      COMPLEX(DP),        INTENT(OUT)   :: wf(:,:)
      INTEGER,            INTENT(IN)    :: ik, nk
      INTEGER,            INTENT(IN)    :: kunit
      INTEGER,            INTENT(INOUT) :: ngw, nbnd, ispin, nspin
      INTEGER,            INTENT(IN)    :: ngwl
      INTEGER,            INTENT(IN)    :: igl(:)
      CHARACTER(LEN=256), INTENT(IN)    :: filename
      REAL(DP),           INTENT(OUT)   :: scalef
      LOGICAL,            INTENT(IN)    :: ionode
      INTEGER,            INTENT(IN)    :: root_in_group, intra_group_comm, inter_group_comm, parent_group_comm
      LOGICAL, OPTIONAL,  INTENT(IN)    :: flink
      !
      INTEGER                  :: j
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
      INTEGER                  :: ierr
      INTEGER                  :: iks, ike, ikt
      INTEGER                  :: igwx, igwx_, ik_, nk_
      INTEGER                  :: ngroup, ipdest
      INTEGER,     ALLOCATABLE :: ipmask(:)
      LOGICAL                  :: flink_
      INTEGER                  :: me_in_group, nproc_in_group, io_in_parent, nproc_in_parent, me_in_parent
      !
      flink_ = .FALSE.
      IF( PRESENT( flink ) ) flink_ = flink
      !
      ngroup          = mp_size( inter_group_comm )
      me_in_group     = mp_rank( intra_group_comm )
      nproc_in_group  = mp_size( intra_group_comm )
      me_in_parent    = mp_rank( parent_group_comm )
      nproc_in_parent = mp_size( parent_group_comm )
      !
      ALLOCATE( ipmask( nproc_in_parent ) )
      !
      io_in_parent = 0
      IF( ionode ) io_in_parent = me_in_parent
      CALL mp_sum( io_in_parent, parent_group_comm )
      !
      CALL set_kpoints_vars( ik, nk, kunit, ngwl, igl, &
                             ngroup, ikt, iks, ike, igwx, ipmask, ipdest, &
                             ionode, root_in_group, intra_group_comm, inter_group_comm, parent_group_comm )
      !
      !  if flink = .true. we are following a link and the file is
      !  already opened for read
      !
      ierr = 0
      !
      IF ( ionode .AND. .NOT. flink_ ) &
         CALL iotk_open_read( iuni, FILE = filename, &
                              BINARY = .TRUE., IERR = ierr )
      !
      CALL mp_bcast( ierr, io_in_parent, parent_group_comm )
      !
      CALL errore( 'read_wfc ', &
                   'cannot open restart file for reading', ierr )
      !
      IF ( ionode ) THEN
          !
          CALL iotk_scan_empty( iuni, "INFO", attr )
          !
          CALL iotk_scan_attr( attr, "ngw",          ngw )
          CALL iotk_scan_attr( attr, "nbnd",         nbnd )
          CALL iotk_scan_attr( attr, "ik",           ik_ )
          CALL iotk_scan_attr( attr, "nk",           nk_ )
          CALL iotk_scan_attr( attr, "ispin",        ispin )
          CALL iotk_scan_attr( attr, "nspin",        nspin )
          CALL iotk_scan_attr( attr, "igwx",         igwx_ )
          CALL iotk_scan_attr( attr, "scale_factor", scalef )
          !
      END IF
      !
      CALL mp_bcast( ngw,    io_in_parent, parent_group_comm )
      CALL mp_bcast( nbnd,   io_in_parent, parent_group_comm )
      CALL mp_bcast( ik_,    io_in_parent, parent_group_comm )
      CALL mp_bcast( nk_,    io_in_parent, parent_group_comm )
      CALL mp_bcast( ispin,  io_in_parent, parent_group_comm )
      CALL mp_bcast( nspin,  io_in_parent, parent_group_comm )
      CALL mp_bcast( igwx_,  io_in_parent, parent_group_comm )
      CALL mp_bcast( scalef, io_in_parent, parent_group_comm )
      !
      ALLOCATE( wtmp( MAX( igwx_, igwx ) ) )
      !
      DO j = 1, nbnd
         !
         IF ( j <= SIZE( wf, 2 ) ) THEN
            !
            IF ( ionode ) THEN 
               !
               CALL iotk_scan_dat( iuni, &
                                   "evc" // iotk_index( j ), wtmp(1:igwx_) )
               !
               IF ( igwx > igwx_ ) wtmp((igwx_+1):igwx) = 0.0_DP
               ! ===========================================================
               !       Lingzhu Kong
               !IF ( j .eq. 1)write(*,'(10f12.5)')(wtmp(i),i=1,igwx_)
               ! ===========================================================
               !
            END IF
            !
            IF ( ngroup > 1 ) THEN
               !
               IF ( ipdest /= io_in_parent ) &
                  CALL mp_put( wtmp, wtmp, me_in_parent, &
                               io_in_parent, ipdest, j, parent_group_comm )
               !
               IF ( ( ikt >= iks ) .AND. ( ikt <= ike ) ) &
                  CALL splitwf( wf(:,j), wtmp, ngwl, igl, me_in_group, &
                                nproc_in_group, root_in_group, intra_group_comm )
               !
            ELSE
               !
               CALL splitwf( wf(:,j), wtmp, ngwl, igl, &
                             me_in_parent, nproc_in_parent, io_in_parent, parent_group_comm )
               !
            END IF
            !
         END IF
         !
      END DO
      !
      IF ( ionode .AND. .NOT. flink_ ) CALL iotk_close_read( iuni )
      !
      DEALLOCATE( wtmp )
      DEALLOCATE( ipmask )
      !
      RETURN
      !
    END SUBROUTINE read_wfc
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_eig( iuni, filename, nbnd, eig, energy_units, &
                          occ, ik, ispin, lkpoint_dir )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN) :: iuni
      INTEGER,            INTENT(IN) :: nbnd
      REAL(DP),           INTENT(IN) :: eig(:)
      CHARACTER(*),       INTENT(IN) :: energy_units
      REAL(DP), OPTIONAL, INTENT(IN) :: occ(:)
      INTEGER,  OPTIONAL, INTENT(IN) :: ik, ispin
      LOGICAL,  OPTIONAL, INTENT(IN) :: lkpoint_dir
      CHARACTER(LEN=256), INTENT(IN) :: filename
      LOGICAL :: lkpoint_dir0
      !
      lkpoint_dir0=.TRUE.
      IF (present(lkpoint_dir)) lkpoint_dir0=lkpoint_dir
      IF ( ionode ) THEN
         !
         if (lkpoint_dir0) CALL iotk_open_write ( iuni, &
                           FILE = TRIM( filename ), BINARY = .FALSE. )
         !
         CALL iotk_write_attr ( attr, "nbnd", nbnd, FIRST=.TRUE. )
         IF ( PRESENT( ik) )    CALL iotk_write_attr ( attr, "ik", ik )
         IF ( PRESENT( ispin) ) CALL iotk_write_attr ( attr, "ispin", ispin )
         CALL iotk_write_empty( iuni, "INFO", ATTR = attr )
         !
         CALL iotk_write_attr ( attr, "UNITS", TRIM(energy_units), FIRST = .TRUE. )
         CALL iotk_write_empty( iuni, "UNITS_FOR_ENERGIES", ATTR=attr)
         !
         CALL iotk_write_dat( iuni, "EIGENVALUES", eig(:) )
         !
         IF ( PRESENT( occ ) ) THEN
            !
            CALL iotk_write_dat( iuni, "OCCUPATIONS", occ(:) )
            !
         ENDIF
         !
         IF (lkpoint_dir0) CALL iotk_close_write ( iuni )
         !
      ENDIF
      !
    END SUBROUTINE write_eig
         
END MODULE xml_io_base
