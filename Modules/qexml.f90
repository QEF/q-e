!
! Copyright (C) 2006 WanT Group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! TB
! included monopole related variables in qexml_write_efield and qexml_read_efield
!
!----------------------------------------------------------------------------
MODULE qexml_module
  !----------------------------------------------------------------------------
  !
  ! This module contains some common subroutines used to read and write
  ! in XML format the data produced by Quantum ESPRESSO package.
  !
  ! Written by Andrea Ferretti (2006).
  ! Modified by Simone Ziraldo (2013).
  !
  ! Important parts of the implementation are taken from xml_io_base.f90
  ! (written by Carlo Sbraccia) in the Quantum ESPRESSO distribution,
  ! under the GNU-GPL licensing:
  !
  ! Copyright (C) 2005 Quantum ESPRESSO group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
#if defined(__OLDXML)
  !
  USE iotk_module
  USE kinds, ONLY : DP
  IMPLICIT NONE
  !
  PRIVATE
  SAVE
  !
  ! definitions for the fmt
  !
  CHARACTER(5), PARAMETER :: fmt_name = "QEXML"
  CHARACTER(5), PARAMETER :: fmt_version = "1.4.0"
  !
  ! some default for kinds
  !
  !INTEGER,   PARAMETER :: dbl = selected_real_kind( 14, 200 )
  !
  ! internal data to be set
  !
  CHARACTER(256)   :: datadir_in, datadir_out
  INTEGER          :: iunit, ounit
  !
  ! vars to manage back compatibility
  !
  CHARACTER(10)    :: qexml_current_version = " "
  CHARACTER(10)    :: qexml_default_version = trim( fmt_version  )
  LOGICAL          :: qexml_current_version_init = .false.
  LOGICAL          :: qexml_version_before_1_4_0 = .false.
  !
  CHARACTER(iotk_attlenx) :: attr
  !
  !
  ! end of declarations
  !
  PUBLIC :: qexml_current_version, qexml_default_version
  PUBLIC :: qexml_current_version_init
  !
  PUBLIC :: qexml_init,  qexml_openfile, qexml_closefile
  !
  PUBLIC :: qexml_write_header, qexml_write_control, qexml_write_status_cp, qexml_write_cell,  &
            qexml_write_moving_cell, qexml_write_ions, qexml_write_symmetry, qexml_write_efield, &
            qexml_write_planewaves, qexml_write_spin, qexml_write_magnetization, &
            qexml_write_xc, qexml_write_exx, qexml_write_occ, qexml_write_bz, qexml_write_para, &
            qexml_write_bands_pw,qexml_write_bands_cp, qexml_write_bands_info, qexml_write_eig, &
            qexml_write_gk, qexml_write_wfc, qexml_write_rho, qexml_write_esm
  !
  PUBLIC :: qexml_read_header, qexml_read_status_cp, qexml_read_cell, qexml_read_moving_cell, qexml_read_ions,      &
            qexml_read_symmetry, qexml_read_efield,                   &
            qexml_read_planewaves, qexml_read_spin, qexml_read_xc,    &
            qexml_read_occ, qexml_read_bz, qexml_read_phonon,         &
            qexml_read_bands_pw, qexml_read_bands_cp, qexml_read_bands_info,                  &
            qexml_read_gk, qexml_read_wfc, qexml_read_rho, qexml_read_magnetization, &
            qexml_read_exx, qexml_read_para, qexml_read_esm
  
  !
  PUBLIC :: qexml_wfc_filename, qexml_create_directory, &
            qexml_kpoint_dirname, qexml_restart_dirname
  !
CONTAINS
!
!-------------------------------------------
! ... basic (public) subroutines
!-------------------------------------------
!
    !------------------------------------------------------------------------
    SUBROUTINE qexml_init( unit_in, unit_out, dir, dir_in, dir_out, &
                           datafile, datafile_in, datafile_out )
      !------------------------------------------------------------------------
      !
      ! just init module data
      !
      IMPLICIT NONE
      INTEGER,                INTENT(in) :: unit_in
      INTEGER,      OPTIONAL, INTENT(in) :: unit_out
      CHARACTER(*), OPTIONAL, INTENT(in) :: dir
      CHARACTER(*), OPTIONAL, INTENT(in) :: dir_in, dir_out
      CHARACTER(*), OPTIONAL, INTENT(in) :: datafile
      CHARACTER(*), OPTIONAL, INTENT(in) :: datafile_in, datafile_out
      !
      iunit       = unit_in
      ounit       = unit_in
      IF ( present( unit_out ) ) ounit  = unit_out
      !
      !
      datadir_in  = "./"
      datadir_out = "./"
      !
      ! first check whether datafile is given
      !
      IF ( present( datafile ) ) THEN
          !
          datadir_in = datafile
          CALL qexml_basename ( datadir_in,  "data-file.xml")
          !
          datadir_out = datadir_in
          !
      ENDIF
      !
      IF ( present( datafile_in ) ) THEN
          !
          datadir_in = datafile_in
          CALL qexml_basename ( datadir_in,  "data-file.xml")
          !
      ENDIF
      !
      IF ( present( datafile_out ) ) THEN
          !
          datadir_out = datafile_out
          CALL qexml_basename ( datadir_out,  "data-file.xml")
          !
      ENDIF
      !
      ! the presence of directories overwrites any info
      ! about datafiles
      !
      IF ( present( dir ) ) THEN
          datadir_in  = trim(dir)
          datadir_out = trim(dir)
      ENDIF
      !
      IF ( present( dir_in ) ) THEN
          datadir_in  = trim(dir_in)
      ENDIF
      !
      IF ( present( dir_out ) ) THEN
          datadir_out  = trim(dir_out)
      ENDIF
      !
    END SUBROUTINE qexml_init
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_openfile( filename, action, binary, ierr)
      !------------------------------------------------------------------------
      !
      ! open data file
      !
      IMPLICIT NONE
      !
      CHARACTER(*),       INTENT(in)  :: filename
      CHARACTER(*),       INTENT(in)  :: action      ! ("read"|"write")
      LOGICAL, OPTIONAL,  INTENT(in)  :: binary
      INTEGER,            INTENT(out) :: ierr
      !
      LOGICAL :: binary_

      ierr = 0
      binary_ = .false.
      IF ( present(binary) ) binary_ = binary
      !
      SELECT CASE ( trim(action) )
      !
      CASE ( "read", "READ" )
          !
          CALL iotk_open_read ( iunit, FILE = trim(filename), IERR=ierr )
          IF ( ierr/=0 ) RETURN
          !
          CALL qexml_read_header( FORMAT_VERSION=qexml_current_version, IERR=ierr )
          IF ( ierr/=0 ) qexml_current_version = trim( qexml_default_version )
          !
          !
      CASE ( "write", "WRITE" )
          !
          CALL iotk_open_write( iunit, FILE = trim(filename), BINARY=binary_, IERR=ierr )
          IF ( ierr/=0 ) RETURN
          !
          qexml_current_version = trim( qexml_default_version )
          !
      CASE DEFAULT
          ierr = 1
      END SELECT
      !
      ! init logical variables for versioning
      !
      qexml_version_before_1_4_0 = .false.
      !
      IF ( trim( qexml_version_compare( qexml_current_version, "1.4.0" )) == "older" ) &
         qexml_version_before_1_4_0 = .true.
      !
      qexml_current_version_init = .true.
      !
      !
    END SUBROUTINE qexml_openfile
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_closefile( action, ierr)
      !------------------------------------------------------------------------
      !
      ! close data file
      !
      IMPLICIT NONE
      !
      CHARACTER(*),  INTENT(in)  :: action      ! ("read"|"write")
      INTEGER,       INTENT(out) :: ierr
      !
      ierr = 0
      !
      SELECT CASE ( trim(action) )
      CASE ( "read", "READ" )
          !
          CALL iotk_close_read( iunit, IERR=ierr )
          !
      CASE ( "write", "WRITE" )
          !
          CALL iotk_close_write( iunit, IERR=ierr )
          !
      CASE DEFAULT
          ierr = 2
      END SELECT
      !
    END SUBROUTINE qexml_closefile
!
!-------------------------------------------
! ... basic (private) subroutines
!-------------------------------------------
!
    !------------------------------------------------------------------------
    FUNCTION int_to_char( i )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) :: i
      CHARACTER (len=6)   :: int_to_char
      !
      !
      IF ( i < 10 ) THEN
         !
         WRITE( UNIT = int_to_char , FMT = "(I1)" ) i
         !
      ELSEIF ( i < 100 ) THEN
         !
         WRITE( UNIT = int_to_char , FMT = "(I2)" ) i
         !
       ELSEIF ( i < 1000 ) THEN
         !
         WRITE( UNIT = int_to_char , FMT = "(I3)" ) i
         !
       ELSEIF ( i < 10000 ) THEN
         !
         WRITE( UNIT = int_to_char , FMT = "(I4)" ) i
         !
       ELSE
         !
       WRITE( UNIT = int_to_char , FMT = "(I5)" ) i
       !
      ENDIF
      !
    END FUNCTION int_to_char
    !
    !
    !--------------------------------------------------------------------------
    SUBROUTINE qexml_basename( str, extension )
      !--------------------------------------------------------------------------
      !
      ! perform the basename operation on the string str, eliminating
      ! any ending (rightmost) occurrence of extension
      !
      CHARACTER(*),  INTENT(inout) :: str
      CHARACTER(*),  INTENT(in)    :: extension
      !
      INTEGER :: ind, strlen, extlen, i
      !
      IF( len_trim(extension) == 0  .or. len_trim(str) == 0 ) RETURN
      !
      strlen = len_trim( str )
      extlen = len_trim( extension )
      ind    = index( str, trim(extension), BACK=.true. )
      !
      IF ( ind <= 0 .or. ind > strlen ) RETURN
      !
      ! we want to cut only the last part of the name
      ! any intermediate matching is rejected
      !
      IF ( strlen -ind +1 /= extlen ) RETURN
      !
      DO i = ind, strlen
         str(i:i) = ' '
      ENDDO
      !
    END SUBROUTINE qexml_basename
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_version_parse(str, major, minor, patch, ierr)
      !--------------------------------------------------------------------------
      !
      ! Determine the major, minor and patch numbers from
      ! a version string with the fmt "i.j.k"
      !
      ! The ierr variable assumes the following values
      !
      ! ierr < 0     emtpy string
      ! ierr = 0     no problem
      ! ierr > 0     fatal error
      !
      IMPLICIT NONE
      CHARACTER(*),     INTENT(in)    :: str
      INTEGER,          INTENT(out)   :: major, minor, patch, ierr
      !
      INTEGER       :: i1, i2, length
      INTEGER       :: ierrtot
      CHARACTER(10) :: num(3)

      !
      major = 0
      minor = 0
      patch = 0

      length = len_trim( str )
      !
      IF ( length == 0 ) THEN
         !
         ierr = -1
         RETURN
         !
      ENDIF

      i1 = scan( str, ".")
      i2 = scan( str, ".", BACK=.true.)
      !
      IF ( i1 == 0 .or. i2 == 0 .or. i1 == i2 ) THEN
         !
         ierr = 1
         RETURN
         !
      ENDIF
      !
      num(1) = str(    1 : i1-1 )
      num(2) = str( i1+1 : i2-1 )
      num(3) = str( i2+1 : )
      !
      ierrtot = 0
      !
      READ( num(1), *, IOSTAT=ierr ) major
      IF (ierr/=0) RETURN
      !
      READ( num(2), *, IOSTAT=ierr ) minor
      IF (ierr/=0) RETURN
      !
      READ( num(3), *, IOSTAT=ierr ) patch
      IF (ierr/=0) RETURN
      !
    END SUBROUTINE qexml_version_parse
    !
    !--------------------------------------------------------------------------
    FUNCTION qexml_version_compare(str1, str2)
      !--------------------------------------------------------------------------
      !
      ! Compare two version strings; the result is
      !
      ! "newer":   str1 is newer that str2
      ! "equal":   str1 is equal   to str2
      ! "older":   str1 is older than str2
      ! " ":       str1 or str2 has a wrong format
      !
      IMPLICIT NONE
      CHARACTER(*)  :: str1, str2
      CHARACTER(10) :: qexml_version_compare
      !
      INTEGER   :: version1(3), version2(3)
      INTEGER   :: basis, icheck1, icheck2
      INTEGER   :: ierr
      !

      qexml_version_compare = " "
      !
      CALL qexml_version_parse( str1, version1(1), version1(2), version1(3), ierr)
      IF ( ierr/=0 ) RETURN
      !
      CALL qexml_version_parse( str2, version2(1), version2(2), version2(3), ierr)
      IF ( ierr/=0 ) RETURN
      !
      !
      basis = 1000
      !
      icheck1 = version1(1) * basis**2 + version1(2)* basis + version1(3)
      icheck2 = version2(1) * basis**2 + version2(2)* basis + version2(3)
      !
      IF ( icheck1 > icheck2 ) THEN
         !
         qexml_version_compare = 'newer'
         !
      ELSEIF( icheck1 == icheck2 ) THEN
         !
         qexml_version_compare = 'equal'
         !
      ELSE
         !
         qexml_version_compare = 'older'
         !
      ENDIF
      !
    END FUNCTION qexml_version_compare
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_create_directory( dirname, ierr )
      !------------------------------------------------------------------------
      !
      USE wrappers,      ONLY : f_mkdir_safe
      !
      CHARACTER(len=*), INTENT(in)  :: dirname
      INTEGER,          INTENT(out) :: ierr
      !
      INTEGER  :: iunaux
      !
      ierr = 0
      CALL iotk_free_unit( iunaux )
      !
      ierr = f_mkdir_safe( TRIM(dirname) )
      !
      IF ( ierr/=0 ) RETURN
      !
      ! ... check whether the scratch directory is writable
      !
      OPEN( iunaux , FILE = trim( dirname ) // '/test', IOSTAT = ierr )
      IF ( ierr/=0 ) RETURN
      !
      CLOSE( iunaux , STATUS = 'DELETE' )
      !
      RETURN
      !
    END SUBROUTINE qexml_create_directory
    !
    !
    !------------------------------------------------------------------------
    FUNCTION qexml_kpoint_dirname( basedir, ik )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=256)           :: qexml_kpoint_dirname
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
         call errore('qexml_kpoint_dirname','ik too large, increase format',1)
      ENDIF
      !
      qexml_kpoint_dirname = TRIM( kdirname )
      !
      RETURN
      !
    END FUNCTION qexml_kpoint_dirname
    !
    !
    !------------------------------------------------------------------------
    FUNCTION qexml_wfc_filename( basedir, name, ik, ipol, tag, extension, dir )
      !------------------------------------------------------------------------
      !
      CHARACTER(len=256)                 :: qexml_wfc_filename
      CHARACTER(len=*),       INTENT(in) :: basedir
      CHARACTER(len=*),       INTENT(in) :: name
      INTEGER,                INTENT(in) :: ik
      INTEGER,      OPTIONAL, INTENT(in) :: ipol
      CHARACTER(*), OPTIONAL, INTENT(in) :: tag
      CHARACTER(*), OPTIONAL, INTENT(in) :: extension
      LOGICAL,      OPTIONAL, INTENT(in) :: dir
      !
      CHARACTER(len=256) :: filename, tag_, ext_
      LOGICAL :: dir_true
      !
      !
      filename = ' '
      tag_     = ' '
      ext_     = '.dat'
      dir_true = .true.
      !
      IF ( present( tag ) )         tag_ = '_'//trim(tag)
      IF ( present( extension ) )   ext_ = '.'//trim(extension)
      !
      IF ( present( ipol ) ) THEN
         !
         WRITE( filename, FMT = '( I1 )' ) ipol
         !
      ENDIF
      IF (PRESENT(dir)) dir_true=dir
      !
      IF (dir_true) THEN
         filename = TRIM( qexml_kpoint_dirname( basedir, ik ) ) // '/' // &
                 & TRIM( name ) // TRIM( filename ) // TRIM( tag_ ) // TRIM( ext_)
      ELSE
         filename = TRIM( qexml_kpoint_dirname( basedir, ik ) ) // '_' // &
                 & TRIM( name ) // TRIM( filename ) // TRIM( tag_ ) // TRIM( ext_)
      ENDIF
      !
      !
      qexml_wfc_filename = trim( filename )
      !
      RETURN
      !
    END FUNCTION qexml_wfc_filename
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_copy_file( file_in, file_out, ierr )
      !------------------------------------------------------------------------
      !
      CHARACTER(len=*),  INTENT(in) :: file_in, file_out
      INTEGER,           INTENT(out):: ierr
      !
      CHARACTER(len=256) :: string
      INTEGER            :: iun_in, iun_out, ios
      !
      !
      ierr = 0
      !
      CALL iotk_free_unit( iun_in,  ierr )
      IF ( ierr /= 0) RETURN
      CALL iotk_free_unit( iun_out, ierr )
      IF ( ierr /= 0) RETURN
      !
      OPEN( UNIT = iun_in,  FILE = file_in,  STATUS = "OLD", IOSTAT=ierr )
      IF ( ierr /= 0) RETURN
      OPEN( UNIT = iun_out, FILE = file_out, STATUS = "UNKNOWN", IOSTAT=ierr )
      IF ( ierr /= 0) RETURN
      !
      copy_loop: DO
         !
         READ( UNIT = iun_in, FMT = '(A256)', IOSTAT = ios ) string
         !
         IF ( ios < 0 ) exit copy_loop
         !
         WRITE( UNIT = iun_out, FMT = '(A)' ) trim( string )
         !
      ENDDO copy_loop
      !
      CLOSE( UNIT = iun_in )
      CLOSE( UNIT = iun_out )
      !
      RETURN
      !
    END SUBROUTINE qexml_copy_file
    !
    !
    !------------------------------------------------------------------------
    FUNCTION check_file_exst( filename )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      LOGICAL          :: check_file_exst
      CHARACTER(len=*) :: filename
      !
      LOGICAL :: lexists
      !
      INQUIRE( FILE = trim( filename ), EXIST = lexists )
      !
      check_file_exst = lexists
      RETURN
      !
    END FUNCTION check_file_exst
    !
    !
    !------------------------------------------------------------------------
    FUNCTION qexml_restart_dirname( outdir, prefix, runit )
      !------------------------------------------------------------------------
      !
      CHARACTER(len=256)           :: qexml_restart_dirname
      CHARACTER(len=*), INTENT(in) :: outdir, prefix
      INTEGER,          INTENT(IN) :: runit
      !
      CHARACTER(len=256)         :: dirname
      INTEGER                    :: strlen
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      !
      ! ... main restart directory
      !
      dirname = TRIM( prefix ) // '_' // TRIM( int_to_char( runit ) )// '.save/'
      !
      IF ( len( outdir ) > 1 ) THEN
         !
         strlen = len_trim( outdir )
         IF ( outdir(strlen:strlen) == '/' ) strlen = strlen -1
         !
         dirname = outdir(1:strlen) // '/' // dirname
         !
      ENDIF
      !
      qexml_restart_dirname = trim( dirname )
      !
      RETURN
      !
    END FUNCTION qexml_restart_dirname
!
!-------------------------------------------
! ... write subroutines
!-------------------------------------------
!
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_header( creator_name, creator_version )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in) :: creator_name, creator_version
      CHARACTER(iotk_attlenx)  :: attr
      !
      CALL iotk_write_begin( ounit, "HEADER" )
      !
      CALL iotk_write_attr(attr, "NAME",trim(fmt_name), FIRST=.true.)
      CALL iotk_write_attr(attr, "VERSION",trim(fmt_version) )
      CALL iotk_write_empty( ounit, "FORMAT", ATTR=attr )
      !
      CALL iotk_write_attr(attr, "NAME",trim(creator_name), FIRST=.true.)
      CALL iotk_write_attr(attr, "VERSION",trim(creator_version) )
      CALL iotk_write_empty( ounit, "CREATOR", ATTR=attr )
      !
      CALL iotk_write_end( ounit, "HEADER" )
      !
    END SUBROUTINE qexml_write_header
    !
    !
    !------------------------------------------------------------------------  
    SUBROUTINE qexml_write_control( pp_check_flag, lkpoint_dir, q_real_space, tq_smoothing, tbeta_smoothing, beta_real_space)
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      LOGICAL, OPTIONAL, INTENT(IN) :: pp_check_flag, lkpoint_dir, q_real_space, tq_smoothing, tbeta_smoothing, beta_real_space


      CALL iotk_write_begin( ounit, "CONTROL" )
      !
      !  This flag is used to check if the file can be used for post-processing
      IF ( PRESENT( pp_check_flag ) ) &
         CALL iotk_write_dat( ounit, "PP_CHECK_FLAG", pp_check_flag )
      !
      !  This flag says how eigenvalues are saved
      IF ( PRESENT( lkpoint_dir ) ) &
         CALL iotk_write_dat( ounit, "LKPOINT_DIR", lkpoint_dir )
      !
      !  This flag says if Q in real space has to be used
      IF ( PRESENT( q_real_space ) ) &
         CALL iotk_write_dat( ounit, "Q_REAL_SPACE", q_real_space )
      ! This flag says if Beta functions were treated in real space
      IF ( PRESENT( beta_real_space ) ) &
         CALL iotk_write_dat( ounit, "BETA_REAL_SPACE", beta_real_space )
      ! This flag says if the Q are being smoothed 
      IF ( PRESENT( tq_smoothing ) ) &
         CALL iotk_write_dat( ounit, "TQ_SMOOTHING", tq_smoothing )
      ! This flag says if the beta are being smoothed 
      IF ( PRESENT( tbeta_smoothing ) ) &
         CALL iotk_write_dat( ounit, "TBETA_SMOOTHING", tbeta_smoothing )
      !
      CALL iotk_write_end( ounit, "CONTROL" )
      !
    END SUBROUTINE qexml_write_control
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_status_cp( nfi,simtime,time_units,title, &
                                  ekin, eht, esr, eself, epseu, enl, exc, vave, enthal, &
                                  energy_units)
      !------------------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: nfi
      REAL(DP), INTENT(in) :: simtime, ekin,eht,esr,eself,epseu,enl,exc,vave,enthal
      CHARACTER(len=*), INTENT(in) :: time_units, title, energy_units
      
      CALL iotk_write_begin( ounit, "STATUS" )
      !
      CALL iotk_write_attr( attr, "ITERATION", nfi, FIRST = .TRUE. )
      CALL iotk_write_empty( ounit, "STEP", attr )
      !
      CALL iotk_write_attr( attr, "UNITS", time_units, FIRST = .TRUE. )
      CALL iotk_write_dat( ounit, "TIME", simtime, ATTR = attr )
      !
      CALL iotk_write_dat( ounit, "TITLE", title )
      !
      CALL iotk_write_attr( attr, "UNITS", energy_units, FIRST = .TRUE. )
      CALL iotk_write_dat( ounit, "KINETIC_ENERGY", ekin,   ATTR = attr )
      CALL iotk_write_dat( ounit, "HARTREE_ENERGY", eht,    ATTR = attr )
      CALL iotk_write_dat( ounit, "EWALD_TERM",     esr,    ATTR = attr )
      CALL iotk_write_dat( ounit, "GAUSS_SELFINT",  eself,  ATTR = attr )
      CALL iotk_write_dat( ounit, "LPSP_ENERGY",    epseu,  ATTR = attr )
      CALL iotk_write_dat( ounit, "NLPSP_ENERGY",   enl,    ATTR = attr )
      CALL iotk_write_dat( ounit, "EXC_ENERGY",     exc,    ATTR = attr )
      CALL iotk_write_dat( ounit, "AVERAGE_POT",    vave,   ATTR = attr )
      CALL iotk_write_dat( ounit, "ENTHALPY",       enthal, ATTR = attr )
      !
      CALL iotk_write_end( ounit, "STATUS" )
      !
    END SUBROUTINE qexml_write_status_cp
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_cell( ibravais_latt, celldm, alat, &
                                 a1, a2, a3, b1, b2, b3, alat_units, a_units, b_units, &
                                 do_mp, do_mt, do_esm)
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(in) :: ibravais_latt
      REAL(DP),        INTENT(in) :: celldm(6), alat
      REAL(DP),        INTENT(in) :: a1(3), a2(3), a3(3)
      REAL(DP),        INTENT(in) :: b1(3), b2(3), b3(3)
      CHARACTER(len=*), INTENT(in) :: alat_units, a_units, b_units
      LOGICAL,          INTENT(in) :: do_mp, do_mt, do_esm
      !
      CHARACTER(len=256) :: bravais_lattice, es_corr
      !
      CALL iotk_write_begin( ounit, "CELL" )
      !
      SELECT CASE ( ibravais_latt )
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
      !
      CALL iotk_write_dat( ounit, &
                           "NON-PERIODIC_CELL_CORRECTION", TRIM( es_corr ) )
      !
      CALL iotk_write_dat( ounit, &
                           "BRAVAIS_LATTICE", trim( bravais_lattice ) )
      !
      !
      CALL iotk_write_attr( attr, "UNITS", trim(alat_units), FIRST = .true. )
      CALL iotk_write_dat( ounit, "LATTICE_PARAMETER", alat, ATTR = attr )
      !
      CALL iotk_write_dat( ounit, "CELL_DIMENSIONS", celldm(1:6) )
      !
      CALL iotk_write_attr ( attr,   "UNITS", trim(a_units), FIRST = .true. )
      CALL iotk_write_begin( ounit, "DIRECT_LATTICE_VECTORS" )
      CALL iotk_write_empty( ounit, "UNITS_FOR_DIRECT_LATTICE_VECTORS", &
                                     ATTR=attr )
      CALL iotk_write_dat(   ounit, "a1", a1(:) * alat, COLUMNS=3 )
      CALL iotk_write_dat(   ounit, "a2", a2(:) * alat, COLUMNS=3 )
      CALL iotk_write_dat(   ounit, "a3", a3(:) * alat, COLUMNS=3 )
      CALL iotk_write_end(   ounit, "DIRECT_LATTICE_VECTORS" )
      !
      CALL iotk_write_attr ( attr,   "UNITS", trim(b_units), FIRST = .true. )
      CALL iotk_write_begin( ounit, "RECIPROCAL_LATTICE_VECTORS" )
      CALL iotk_write_empty( ounit, "UNITS_FOR_RECIPROCAL_LATTICE_VECTORS", &
                                     ATTR=attr )
      CALL iotk_write_dat(   ounit, "b1", b1(:), COLUMNS=3 )
      CALL iotk_write_dat(   ounit, "b2", b2(:), COLUMNS=3 )
      CALL iotk_write_dat(   ounit, "b3", b3(:), COLUMNS=3 )
      CALL iotk_write_end(   ounit, "RECIPROCAL_LATTICE_VECTORS" )
      !
      CALL iotk_write_end( ounit, "CELL" )
      !
    END SUBROUTINE qexml_write_cell
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_moving_cell(lmovecell, cell_factor)
      !------------------------------------------------------------------------
      !
      LOGICAL, INTENT(IN) :: lmovecell
      REAL(DP), INTENT(IN) :: cell_factor
      !
      CALL iotk_write_begin( ounit, "MOVING_CELL" )
      CALL iotk_write_dat( ounit, "CELL_FACTOR", cell_factor)
      CALL iotk_write_end( ounit, "MOVING_CELL"  )
      !
    END SUBROUTINE qexml_write_moving_cell
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_ions( nsp, nat, atm, ityp, psfile, pseudo_dir,  &
                                 amass, amass_units, tau, tau_units, &
                                 if_pos, dirname, pos_unit )
      !------------------------------------------------------------------------
      !
      USE wrappers, ONLY: f_copy
      !
      INTEGER,          INTENT(in) :: nsp, nat
      INTEGER,          INTENT(in) :: ityp(:)
      CHARACTER(len=*), INTENT(in) :: atm(:)
      CHARACTER(len=*), INTENT(in) :: psfile(:)
      CHARACTER(len=*), INTENT(in) :: pseudo_dir
      CHARACTER(len=*), INTENT(in) :: dirname
      REAL(DP),        INTENT(in) :: amass(:)
      CHARACTER(len=*), INTENT(in) :: amass_units
      REAL(DP),        INTENT(in) :: tau(:,:)
      CHARACTER(len=*), INTENT(in) :: tau_units
      INTEGER,          INTENT(in) :: if_pos(:,:)
      REAL(DP),         INTENT(in) :: pos_unit
      !
      INTEGER            :: i, flen, flen2, ierrl
      CHARACTER(len=256) :: file_pseudo_in, file_pseudo_out
      LOGICAL            :: pseudo_exists
      !
      !
      CALL iotk_write_begin( ounit, "IONS" )
      !
      CALL iotk_write_dat( ounit, "NUMBER_OF_ATOMS", nat )
      !
      CALL iotk_write_dat( ounit, "NUMBER_OF_SPECIES", nsp )
      !
      flen = len_trim( pseudo_dir )
      flen2 = len_trim( dirname )
      !
      CALL iotk_write_attr ( attr, "UNITS", trim(amass_units), FIRST = .true. )
      CALL iotk_write_empty( ounit, "UNITS_FOR_ATOMIC_MASSES", ATTR = attr )
      !
      DO i = 1, nsp
         !
         CALL iotk_write_begin( ounit, "SPECIE"//trim(iotk_index(i)) )
         !
         CALL iotk_write_dat( ounit, "ATOM_TYPE", atm(i) )
         !
         CALL iotk_write_dat( ounit, "MASS", amass(i) )
         !
         CALL iotk_write_dat( ounit, "PSEUDO", trim( psfile(i) ) )
         !
         CALL iotk_write_end( ounit, "SPECIE"//trim(iotk_index(i)) )
         !
         ! copy pseudopotential file into data directory
         !
         IF ( pseudo_dir(flen:flen) /= '/' ) THEN
            file_pseudo_in = pseudo_dir(1:flen) // '/' // TRIM(psfile(i))
         ELSE
            file_pseudo_in = pseudo_dir(1:flen) // TRIM(psfile(i))
         ENDIF
         !
         IF ( dirname(flen2:flen2) /= '/' ) THEN
            file_pseudo_out = dirname(1:flen2) // '/' // TRIM(psfile(i))
         ELSE
            file_pseudo_out = dirname(1:flen2) // TRIM(psfile(i))
         END IF
         !
         IF ( file_pseudo_in .ne. file_pseudo_out ) THEN
            !
            INQUIRE ( FILE=file_pseudo_in, EXIST = pseudo_exists )
            IF ( pseudo_exists ) THEN
               ierrl = f_copy( file_pseudo_in, file_pseudo_out )
            ELSE
               CALL infomsg( 'write_ions', &
                   'file ' // TRIM( file_pseudo_in) // ' not present' )
            END IF
            !
         END IF
         !
      ENDDO
      !
      !
      CALL iotk_write_dat( ounit, "PSEUDO_DIR", trim( pseudo_dir) )
      !
      CALL iotk_write_attr( attr, "UNITS", trim(tau_units), FIRST = .true. )
      CALL iotk_write_empty( ounit, "UNITS_FOR_ATOMIC_POSITIONS", ATTR = attr )
      !
      DO i = 1, nat
         !
         CALL iotk_write_attr( attr, "SPECIES", atm( ityp(i) ), FIRST = .true. )
         CALL iotk_write_attr( attr, "INDEX",  ityp(i) )
         CALL iotk_write_attr( attr, "tau",    tau(:,i)*pos_unit )
         CALL iotk_write_attr( attr, "if_pos", if_pos(:,i) )
         CALL iotk_write_empty( ounit, "ATOM" // trim( iotk_index( i ) ), attr )
         !
      ENDDO
      !
      CALL iotk_write_end( ounit, "IONS" )
      !
    END SUBROUTINE qexml_write_ions
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_symmetry( ibrav, nrot, nsym, invsym, noinv, &
                               time_reversal, no_t_rev, ft, &
                               s, sname, s_units, irt, nat, t_rev )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(in) :: ibrav, nrot, nsym
      LOGICAL,          INTENT(in) :: invsym, noinv, time_reversal, no_t_rev
      INTEGER,          INTENT(in) :: s(:,:,:), irt(:,:), nat, t_rev(:)
      REAL(DP),         INTENT(in) :: ft(:,:)
      CHARACTER(LEN=*), INTENT(in) :: sname(:),s_units
      !
      INTEGER  :: i
      !
      !
      CALL iotk_write_begin( ounit, "SYMMETRIES" )
      !
      CALL iotk_write_dat( ounit, "NUMBER_OF_SYMMETRIES", nsym )
      CALL iotk_write_dat( ounit, "NUMBER_OF_BRAVAIS_SYMMETRIES", nrot )
      !
      CALL iotk_write_dat( ounit, "INVERSION_SYMMETRY", invsym )
      !
      CALL iotk_write_dat( ounit, "DO_NOT_USE_TIME_REVERSAL", noinv )
      !
      CALL iotk_write_dat( ounit, "TIME_REVERSAL_FLAG", time_reversal )
      !
      CALL iotk_write_dat( ounit, "NO_TIME_REV_OPERATIONS", no_t_rev )
      !
      CALL iotk_write_dat( ounit, "NUMBER_OF_ATOMS", nat )
      !
      CALL iotk_write_attr( attr, "UNITS", TRIM(s_units), FIRST = .TRUE. )
      CALL iotk_write_empty( ounit, "UNITS_FOR_SYMMETRIES", ATTR = attr )
      !
      DO i = 1, nsym
         !
         CALL iotk_write_begin( ounit, "SYMM" // TRIM( iotk_index( i ) ) )
         !
         CALL iotk_write_attr ( attr, "NAME", TRIM( sname(i) ), FIRST=.TRUE. )
         CALL iotk_write_attr ( attr, "T_REV", t_rev(i) )
         CALL iotk_write_empty( ounit, "INFO", ATTR = attr )
         !
         CALL iotk_write_dat( ounit, "ROTATION", s(:,:,i), COLUMNS=3 )
         CALL iotk_write_dat( ounit, "FRACTIONAL_TRANSLATION", ft(:,i), COLUMNS=3 )
         !
         IF ( nat > 0 ) &
            CALL iotk_write_dat( ounit, "EQUIVALENT_IONS", irt(i,1:nat), COLUMNS=8 )
         !
         CALL iotk_write_end( ounit, "SYMM" // TRIM( iotk_index( i ) ) )
         !
      ENDDO
      !
      ! ... the following are the symmetries of the Bravais lattice alone
      ! ... (they may be more than crystal, i.e. basis+lattice, symmetries)
      !
      DO i = nsym+1, nrot
         !
         CALL iotk_write_begin( ounit, "SYMM" // TRIM( iotk_index( i ) ) )
         !
         CALL iotk_write_attr ( attr, "NAME", TRIM( sname(i) ), FIRST=.TRUE. )
         CALL iotk_write_empty( ounit, "INFO", ATTR = attr )
         CALL iotk_write_dat( ounit, "ROTATION", s(:,:,i), COLUMNS=3 )
         !
         CALL iotk_write_end( ounit, "SYMM" // TRIM( iotk_index( i ) ) )
         !
      ENDDO
      !
      CALL iotk_write_end( ounit, "SYMMETRIES" )
      !
    END SUBROUTINE qexml_write_symmetry
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_efield( tefield, dipfield, edir, emaxpos, eopreg, eamp, &
                                   monopole, zmon, relaxz, block, block_1, block_2,&
                                   block_height)
     !------------------------------------------------------------------------
     !
      LOGICAL, INTENT(in)   :: tefield        ! if .TRUE. a finite electric field
                                              ! is added to the local potential
      LOGICAL, INTENT(in)   :: dipfield       ! if .TRUE. the dipole field is subtracted
      LOGICAL, INTENT(in)   :: monopole       ! if .TRUE. counter charge is represented by monopole (gate)
      LOGICAL, INTENT(in)   :: block          ! add potential barrier
      LOGICAL, INTENT(in)   :: relaxz         ! relax in z direction  
      INTEGER, INTENT(in)   :: edir           ! direction of the field
      REAL(DP), INTENT(in) :: emaxpos        ! position of the maximum of the field (0<emaxpos<1)
      REAL(DP), INTENT(in) :: eopreg         ! amplitude of the inverse region (0<eopreg<1)
      REAL(DP), INTENT(in) :: eamp           ! field amplitude (in a.u.) (1 a.u. = 51.44 10^11 V/m)
      REAL(DP), INTENT(in) :: zmon           ! position of monopole plane in units of cell vector in z direction
      REAL(DP), INTENT(in) :: block_1        ! potential barrier
      REAL(DP), INTENT(in) :: block_2
      REAL(DP), INTENT(in) :: block_height
      !
      !
      CALL iotk_write_begin( ounit, "ELECTRIC_FIELD" )
      !
      CALL iotk_write_dat( ounit, "HAS_ELECTRIC_FIELD", tefield )
      !
      CALL iotk_write_dat( ounit, "HAS_DIPOLE_CORRECTION", dipfield )
      !
      CALL iotk_write_dat( ounit, "FIELD_DIRECTION", edir )
      !
      CALL iotk_write_dat( ounit, "MAXIMUM_POSITION", emaxpos )
      !
      CALL iotk_write_dat( ounit, "INVERSE_REGION", eopreg )
      !
      CALL iotk_write_dat( ounit, "FIELD_AMPLITUDE", eamp )
      !
      CALL iotk_write_dat( ounit, "MONOPOLE_PLANE", monopole )
      !
      CALL iotk_write_dat( ounit, "MONOPOLE_POS", zmon )
      !
      CALL iotk_write_dat( ounit, "RELAX_Z", relaxz )
      !
      CALL iotk_write_dat( ounit, "BLOCK", block )
      !
      CALL iotk_write_dat( ounit, "BLOCK_1", block_1 )
      !
      CALL iotk_write_dat( ounit, "BLOCK_2", block_2 )
      !
      CALL iotk_write_dat( ounit, "BLOCK_HEIGHT", block_height )
      !
      CALL iotk_write_end( ounit, "ELECTRIC_FIELD" )
      !
    END SUBROUTINE qexml_write_efield
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_planewaves( ecutwfc, ecutrho, npwx, gamma_only, &
                                       nr1, nr2, nr3,  ngm,  nr1s, nr2s, nr3s, ngms, &
                                       nr1b, nr2b, nr3b, igv, lgvec, cutoff_units )
      !------------------------------------------------------------------------
#if defined __HDF5
      USE hdf5_qe
      USE mp_pools,  ONLY : inter_pool_comm
      USE io_files,  ONLY : tmp_dir
      USE mp_world,  ONLY : mpime
#endif

      !
      INTEGER,       INTENT(in) :: npwx, nr1, nr2, nr3, ngm, &
                                   nr1s, nr2s, nr3s, ngms, nr1b, nr2b, nr3b
      INTEGER,       INTENT(in) :: igv(:,:)
      REAL(DP),     INTENT(in) :: ecutwfc, ecutrho
      LOGICAL,       INTENT(in) :: gamma_only, lgvec
      CHARACTER(*),  INTENT(in) :: cutoff_units
#if defined __HDF5
      CHARACTER(LEN=256) :: filename_hdf5
      integer           :: ierr
#endif

      !
      !
#if defined __HDF5
      filename_hdf5=trim(tmp_dir) //"g.hdf5"
      CALL prepare_for_writing_final(g_hdf5_write,inter_pool_comm,filename_hdf5)
      CALL add_attributes_hdf5(g_hdf5_write,ecutwfc,"WFC_CUTOFF")
      CALL add_attributes_hdf5(g_hdf5_write,ecutrho,"RHO_CUTOFF")
      CALL add_attributes_hdf5(g_hdf5_write,npwx,"MAX_NUMBER_OF_GK-VECTORS")
      CALL add_attributes_hdf5(g_hdf5_write,gamma_only,"GAMMA_ONLY")
      CALL add_attributes_hdf5(g_hdf5_write,trim(cutoff_units),"UNITS_FOR_CUTOFF")
      CALL add_attributes_hdf5(g_hdf5_write,nr1,"nr1")
      CALL add_attributes_hdf5(g_hdf5_write,nr2,"nr2")
      CALL add_attributes_hdf5(g_hdf5_write,nr3,"nr3")
      CALL add_attributes_hdf5(g_hdf5_write,ngm,"GVECT_NUMBER")
      CALL add_attributes_hdf5(g_hdf5_write,nr1s,"nr1s")
      CALL add_attributes_hdf5(g_hdf5_write,nr2s,"nr2s")
      CALL add_attributes_hdf5(g_hdf5_write,nr3s,"nr3s")
      CALL add_attributes_hdf5(g_hdf5_write,ngms,"SMOOTH_GVECT_NUMBER")
      CALL add_attributes_hdf5(g_hdf5_write,nr1s,"nr1b")
      CALL add_attributes_hdf5(g_hdf5_write,nr2s,"nr2b")
      CALL add_attributes_hdf5(g_hdf5_write,nr3s,"nr3b")
#endif

      CALL iotk_write_begin( ounit, "PLANE_WAVES" )
      !
      CALL iotk_write_attr ( attr, "UNITS", trim(cutoff_units), FIRST = .true. )
      CALL iotk_write_empty( ounit, "UNITS_FOR_CUTOFF", ATTR = attr )
      !
      CALL iotk_write_dat( ounit, "WFC_CUTOFF", ecutwfc )
      !
      CALL iotk_write_dat( ounit, "RHO_CUTOFF", ecutrho )
      !
      CALL iotk_write_dat( ounit, "MAX_NUMBER_OF_GK-VECTORS", npwx )
      !
      CALL iotk_write_dat( ounit, "GAMMA_ONLY", gamma_only )
      !
      CALL iotk_write_attr( attr, "nr1", nr1, FIRST = .true. )
      CALL iotk_write_attr( attr, "nr2", nr2 )
      CALL iotk_write_attr( attr, "nr3", nr3 )
      CALL iotk_write_empty( ounit, "FFT_GRID", ATTR = attr )
      !
      CALL iotk_write_dat( ounit, "GVECT_NUMBER", ngm )
      !
      CALL iotk_write_attr( attr, "nr1s", nr1s, FIRST = .true. )
      CALL iotk_write_attr( attr, "nr2s", nr2s )
      CALL iotk_write_attr( attr, "nr3s", nr3s )
      CALL iotk_write_empty( ounit, "SMOOTH_FFT_GRID", ATTR = attr )
      !
      CALL iotk_write_dat( ounit, "SMOOTH_GVECT_NUMBER", ngms )
      !
      IF ( lgvec ) THEN
         !
         ! ... write the G-vectors
         !
#if defined __HDF5
         CALL write_g(g_hdf5_write,igv(1:3,1:ngm))
#else
         CALL iotk_link( ounit, "G-VECTORS", "./gvectors.dat", &
                         CREATE = .true., BINARY = .true. )
         !
         CALL iotk_write_begin( ounit, "G-VECTORS" )
         !
         CALL iotk_write_attr( attr, "nr1s", nr1s, FIRST = .true. )
         CALL iotk_write_attr( attr, "nr2s", nr2s )
         CALL iotk_write_attr( attr, "nr3s", nr3s )
         CALL iotk_write_attr( attr, "gvect_number", ngm )
         CALL iotk_write_attr( attr, "gamma_only", gamma_only )
         CALL iotk_write_attr( attr, "units", "crystal" )
         CALL iotk_write_empty( ounit, "INFO", ATTR = attr )
         !

         CALL iotk_write_dat  ( ounit, "g", igv(1:3,1:ngm), COLUMNS = 3 )
         CALL iotk_write_end  ( ounit, "G-VECTORS" )
#endif
      ENDIF
      !
      CALL iotk_write_attr( attr, "nr1b", nr1b , FIRST = .true. )
      CALL iotk_write_attr( attr, "nr2b", nr2b )
      CALL iotk_write_attr( attr, "nr3b", nr3b )
      CALL iotk_write_empty( ounit, "SMALLBOX_FFT_GRID", ATTR = attr )
      !
      CALL iotk_write_end( ounit, "PLANE_WAVES" )
#if defined __HDF5
      CALL h5fclose_f(g_hdf5_write%file_id,ierr)
#endif
      !
    END SUBROUTINE qexml_write_planewaves
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_gk( ik, npwk, npwkx, gamma_only, xk, k_units, index, igk )
      !------------------------------------------------------------------------
      !
      INTEGER,      INTENT(in) :: ik
      INTEGER,      INTENT(in) :: npwk, npwkx
      LOGICAL,      INTENT(in) :: gamma_only
      REAL(DP),    INTENT(in) :: xk(3)
      CHARACTER(*), INTENT(in) :: k_units
      LOGICAL,      INTENT(in) :: index(:), igk(:,:)
      !
      INTEGER        :: iunaux
      CHARACTER(256) :: filename

      CALL iotk_free_unit( iunaux )
      filename = qexml_wfc_filename( datadir_out, 'gkvectors', ik )
      !
      CALL iotk_open_write( iunaux, FILE = trim( filename ), &
                            ROOT="GK-VECTORS", BINARY = .true. )
      !
      CALL iotk_write_dat( iunaux, "NUMBER_OF_GK-VECTORS", npwk )
      CALL iotk_write_dat( iunaux, "MAX_NUMBER_OF_GK-VECTORS", npwkx )
      CALL iotk_write_dat( iunaux, "GAMMA_ONLY", gamma_only )
      !
      CALL iotk_write_attr ( attr, "UNITS", trim(k_units), FIRST = .true. )
      CALL iotk_write_dat( iunaux, "K-POINT_COORDS", xk, ATTR = attr )
      !
      CALL iotk_write_dat( iunaux, "INDEX", index(1:npwk) )
      CALL iotk_write_dat( iunaux, "GRID", igk(1:npwk,ik), COLUMNS = 3 )
      !
      CALL iotk_close_write( iunaux )
      !
    END SUBROUTINE qexml_write_gk
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_spin( lsda, noncolin, npol, lspinorb, domag )
      !------------------------------------------------------------------------
      !
      LOGICAL, INTENT(in) :: lsda, noncolin, lspinorb, domag
      INTEGER, INTENT(in) :: npol
      !
      !
      CALL iotk_write_begin( ounit, "SPIN" )
      !
      CALL iotk_write_dat( ounit, "LSDA", lsda )
      !
      CALL iotk_write_dat( ounit, "NON-COLINEAR_CALCULATION", noncolin )
      !
      IF ( noncolin ) &
         CALL iotk_write_dat( ounit, "SPINOR_DIM", npol )
      !
      CALL iotk_write_dat( ounit, "SPIN-ORBIT_CALCULATION", lspinorb )
      CALL iotk_write_dat( ounit, "SPIN-ORBIT_DOMAG", domag )
      !
      CALL iotk_write_end( ounit, "SPIN" )
      !
    END SUBROUTINE qexml_write_spin
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_magnetization(starting_magnetization, angle1, angle2, &
                                   nsp, two_fermi_energies, i_cons, mcons, bfield, &
                                   ef_up, ef_dw, nelup, neldw, lambda, energy_units)
      !------------------------------------------------------------------------
      !
      !
      IMPLICIT NONE
      INTEGER,  INTENT(IN) :: nsp, i_cons
      REAL(DP), INTENT(IN) :: starting_magnetization(nsp), &
                              angle1(nsp), angle2(nsp), mcons(3,nsp), &
                              bfield(3), ef_up, ef_dw, nelup, neldw, lambda
      LOGICAL,  INTENT(IN) :: two_fermi_energies
      CHARACTER(*),       INTENT(IN) :: energy_units
      !
      INTEGER :: i
      !
      CALL iotk_write_begin( ounit, "MAGNETIZATION_INIT" )

      CALL iotk_write_dat( ounit,"CONSTRAINT_MAG", i_cons)

      CALL iotk_write_dat( ounit, "NUMBER_OF_SPECIES", nsp ) 

      DO i = 1, nsp
         !
         CALL iotk_write_begin( ounit, "SPECIE"//TRIM(iotk_index(i)) )
         !
         CALL iotk_write_dat( ounit, "STARTING_MAGNETIZATION",  &
                                       starting_magnetization(i) )
         CALL iotk_write_dat( ounit, "ANGLE1", &
                                       angle1(i) )
         CALL iotk_write_dat( ounit, "ANGLE2", &
                                       angle2(i) )
         IF (i_cons==1.OR.i_cons==2) THEN
            CALL iotk_write_dat( ounit, "CONSTRANT_1", mcons(1,i) )
            CALL iotk_write_dat( ounit, "CONSTRANT_2", mcons(2,i) )
            CALL iotk_write_dat( ounit, "CONSTRANT_3", mcons(3,i) )
         ENDIF
         !
         CALL iotk_write_end( ounit, "SPECIE"//TRIM(iotk_index(i)) )
         !
      ENDDO
      !
      IF (i_cons==3) THEN
         !
         CALL iotk_write_dat( ounit, "FIXED_MAGNETIZATION_1", mcons(1,1) )
         CALL iotk_write_dat( ounit, "FIXED_MAGNETIZATION_2", mcons(2,1) )
         CALL iotk_write_dat( ounit, "FIXED_MAGNETIZATION_3", mcons(3,1) )
         !
      ELSE IF (i_cons==4) THEN
         !
         CALL iotk_write_dat( ounit, "MAGNETIC_FIELD_1", bfield(1) )
         CALL iotk_write_dat( ounit, "MAGNETIC_FIELD_2", bfield(2) )
         CALL iotk_write_dat( ounit, "MAGNETIC_FIELD_3", bfield(3) )
         !
      ENDIF
      !
      CALL iotk_write_dat(ounit,"TWO_FERMI_ENERGIES",two_fermi_energies)
      !
      IF (two_fermi_energies) THEN
         !
         CALL iotk_write_attr ( attr, "UNITS", trim(energy_units), FIRST = .TRUE. )
         CALL iotk_write_empty( ounit, "UNITS_FOR_ENERGIES", ATTR = attr )
         !
         CALL iotk_write_dat( ounit, "FIXED_MAGNETIZATION", mcons(3,1) )
         CALL iotk_write_dat( ounit, "ELECTRONS_UP", nelup )
         CALL iotk_write_dat( ounit, "ELECTRONS_DOWN", neldw )
         CALL iotk_write_dat( ounit, "FERMI_ENERGY_UP", ef_up )
         CALL iotk_write_dat( ounit, "FERMI_ENERGY_DOWN", ef_dw )
         !
      ENDIF
      !
      IF (i_cons>0) CALL iotk_write_dat(ounit,"LAMBDA",lambda)
      !
      CALL iotk_write_end( ounit, "MAGNETIZATION_INIT" )
      !
    RETURN
    !
    END SUBROUTINE qexml_write_magnetization
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_xc( dft, nsp, lda_plus_u, lda_plus_u_kind, U_projection, &
                         Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_J, Hubbard_J0, &
                         Hubbard_beta, Hubbard_alpha,                               &
                         inlc, vdw_table_name, pseudo_dir, acfdt_in_pw, dirname, & 
                         llondon, london_s6, london_rcut, lxdm, ts_vdw, vdw_isolated )
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
      LOGICAL, OPTIONAL,  INTENT(IN) :: acfdt_in_pw
      !
      LOGICAL,  OPTIONAL, INTENT(IN) :: llondon, lxdm, ts_vdw, vdw_isolated
      REAL(DP), OPTIONAL, INTENT(IN) :: london_s6, london_rcut

      INTEGER            :: i, flen, ierrl
      CHARACTER(LEN=256) :: file_table
      !
      CALL iotk_write_begin( ounit, "EXCHANGE_CORRELATION" )
      !
      CALL iotk_write_dat( ounit, "DFT", dft )
      !
      IF ( lda_plus_u ) THEN
         !
         IF ( .NOT. PRESENT( Hubbard_lmax ) .OR. &
              .NOT. PRESENT( Hubbard_l )    .OR. & 
              .NOT. PRESENT( Hubbard_U )    .OR. &
              .NOT. PRESENT( nsp )              )&
            CALL errore( 'write_xc', &
                         ' variables for LDA+U not present', 1 )
         !
         CALL iotk_write_dat( ounit, "LDA_PLUS_U_CALCULATION", lda_plus_u )
         CALL iotk_write_dat( ounit, "NUMBER_OF_SPECIES", nsp )
         CALL iotk_write_dat( ounit, "HUBBARD_LMAX", Hubbard_lmax )
         CALL iotk_write_dat( ounit, "HUBBARD_L", Hubbard_l(1:nsp) )
         CALL iotk_write_dat( ounit, "HUBBARD_U", Hubbard_U(1:nsp) )
         !
         IF ( PRESENT( lda_plus_u_kind ) ) THEN
            CALL iotk_write_dat( ounit, "LDA_PLUS_U_KIND", lda_plus_u_kind )
            CALL iotk_write_dat( ounit, "U_PROJECTION_TYPE", trim(U_projection) )
         END IF
         !
         IF ( PRESENT( Hubbard_J ) ) &
              CALL iotk_write_dat( ounit, "HUBBARD_J", Hubbard_J(1:3,1:nsp), COLUMNS = 3)
         !
         IF ( PRESENT( Hubbard_J0 ) ) &
         CALL iotk_write_dat( ounit, "HUBBARD_J0", Hubbard_J0(1:nsp) )
         !
         IF ( PRESENT( Hubbard_alpha ) ) &
         CALL iotk_write_dat( ounit, "HUBBARD_ALPHA", Hubbard_alpha(1:nsp) )
         !
         IF ( PRESENT( Hubbard_beta ) ) &
         CALL iotk_write_dat( ounit, "HUBBARD_BETA", Hubbard_beta(1:nsp) )
         !
      END IF
      !
      ! Vdw kernel table
      !
      IF ( present(inlc) ) THEN
         IF ( inlc > 0 ) THEN
            IF ( .NOT. PRESENT( vdw_table_name ) .OR. &
                 .NOT. PRESENT( pseudo_dir )     .OR. &
                 .NOT. PRESENT( dirname ))            &
                 CALL errore( 'write_xc', ' variable vdw_table_name not present', 1 )
        
            CALL iotk_write_dat( ounit, "NON_LOCAL_DF", inlc )
            CALL iotk_write_dat( ounit, "VDW_KERNEL_NAME", TRIM(vdw_table_name))
            !
            ! Copy the file in .save directory
            !
            flen = LEN_TRIM( pseudo_dir )
            IF ( pseudo_dir(flen:flen) /= '/' ) THEN
               file_table = pseudo_dir(1:flen) // '/' // vdw_table_name
            ELSE
               file_table = pseudo_dir(1:flen) // vdw_table_name
            END IF
            !
            CALL qexml_copy_file( TRIM( file_table ), TRIM( dirname ) // "/" // TRIM( vdw_table_name ),ierrl )
            !
         ENDIF
      ENDIF
      !
      IF ( PRESENT (llondon) ) THEN
         IF ( llondon ) THEN
            IF ( .NOT. PRESENT( london_s6 )  .OR. &
                 .NOT. PRESENT( london_rcut ) ) & 
               CALL errore( 'write_xc', &
                            ' variables for DFT+D not present', 1 )
            CALL iotk_write_begin( ounit, "DFT_D2" )
            CALL iotk_write_dat( ounit, "SCALING_FACTOR", london_s6 )
            CALL iotk_write_dat( ounit, "CUTOFF_RADIUS",  london_rcut )
            CALL iotk_write_end  ( ounit, "DFT_D2" )
         ENDIF
      ENDIF
      !
      IF ( PRESENT (lxdm) ) THEN
         IF ( lxdm) CALL iotk_write_dat( ounit, "XDM", lxdm )
      ENDIF
      !
      IF ( PRESENT (ts_vdw) ) THEN
         IF ( ts_vdw) THEN
            IF ( .NOT. PRESENT (vdw_isolated) ) &
               CALL errore( 'write_xc', &
                            ' variables for TS not present', 1 )
            CALL iotk_write_begin( ounit, "TKATCHENKO-SCHEFFLER" )
            CALL iotk_write_dat( ounit, "ISOLATED_SYSTEM", vdw_isolated )
            CALL iotk_write_end( ounit, "TKATCHENKO-SCHEFFLER" )
         END IF
      END IF
      !
      IF ( PRESENT (acfdt_in_pw) ) THEN
         CALL iotk_write_dat( ounit, "ACFDT_IN_PW", acfdt_in_pw )
      ENDIF
      !
      IF ( PRESENT (acfdt_in_pw) ) THEN
         CALL iotk_write_dat( ounit, "ACFDT_IN_PW", acfdt_in_pw )
      ENDIF
      !
      CALL iotk_write_end( ounit, "EXCHANGE_CORRELATION" )
      !
 END SUBROUTINE qexml_write_xc
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_exx( x_gamma_extrapolation, nqx1, nqx2, nqx3, &
                          exxdiv_treatment, yukawa, ecutvcut, exx_fraction, &
                          gau_parameter, screening_parameter, exx_is_active, ecutfock )
      !------------------------------------------------------------------------
      !
      LOGICAL,            INTENT(IN) :: x_gamma_extrapolation, exx_is_active
      INTEGER,            INTENT(IN) :: nqx1, nqx2, nqx3
      CHARACTER(LEN=*),   INTENT(IN) :: exxdiv_treatment
      REAL(DP),           INTENT(IN) :: yukawa, ecutvcut, exx_fraction, ecutfock
      REAL(DP),           INTENT(IN) :: screening_parameter
      REAL(DP),           INTENT(IN) :: gau_parameter
      !
      CALL iotk_write_begin(ounit, "EXACT_EXCHANGE" )
      call iotk_write_dat(ounit, "x_gamma_extrapolation", x_gamma_extrapolation)
      call iotk_write_dat(ounit, "nqx1", nqx1)
      call iotk_write_dat(ounit, "nqx2", nqx2)
      call iotk_write_dat(ounit, "nqx3", nqx3)
      call iotk_write_dat(ounit, "exxdiv_treatment", exxdiv_treatment)
      call iotk_write_dat(ounit, "yukawa", yukawa)
      call iotk_write_dat(ounit, "ecutvcut", ecutvcut)
      call iotk_write_dat(ounit, "exx_fraction", exx_fraction)
      call iotk_write_dat(ounit, "screening_parameter", screening_parameter)
      call iotk_write_dat(ounit, "gau_parameter", gau_parameter)
      call iotk_write_dat(ounit, "exx_is_active", exx_is_active)
      call iotk_write_dat(ounit, "ecutfock", ecutfock)
      CALL iotk_write_end(ounit, "EXACT_EXCHANGE" )
      !
    END SUBROUTINE qexml_write_exx
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_esm( esm_nfit, esm_efield, esm_w, esm_a, esm_bc )
      !------------------------------------------------------------------------
      !
      INTEGER,            INTENT(IN) :: esm_nfit
      REAL(DP),           INTENT(IN) :: esm_efield, esm_w, esm_a
      CHARACTER(LEN=*),   INTENT(IN) :: esm_bc
      !
      CALL iotk_write_begin(ounit, "ESM" )
      call iotk_write_dat(ounit, "esm_nfit", esm_nfit)
      call iotk_write_dat(ounit, "esm_efield", esm_efield)
      call iotk_write_dat(ounit, "esm_w", esm_w)
      call iotk_write_dat(ounit, "esm_a", esm_a)
      call iotk_write_dat(ounit, "esm_bc", esm_bc)
      CALL iotk_write_end(ounit, "ESM" )
      !
    END SUBROUTINE qexml_write_esm
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_occ( lgauss, ngauss, degauss, degauss_units, &
         ltetra, tetra_type, ntetra, tetra, tfixed_occ, &
         lsda, nstates_up, nstates_dw, input_occ )
      !------------------------------------------------------------------------
      !
      LOGICAL,                INTENT(in) :: lgauss, ltetra, tfixed_occ, lsda
      INTEGER,      OPTIONAL, INTENT(in) :: ngauss, ntetra, tetra_type, &
           nstates_up, nstates_dw
      INTEGER,      OPTIONAL, INTENT(in) :: tetra(:,:)
      REAL(DP),    OPTIONAL, INTENT(in) :: degauss, input_occ(:,:)
      CHARACTER(*), OPTIONAL, INTENT(in) :: degauss_units
      !
      INTEGER :: i
      !
      !
      CALL iotk_write_begin( ounit, "OCCUPATIONS" )
      !
      CALL iotk_write_dat( ounit, "SMEARING_METHOD", lgauss )
      !
      IF ( lgauss ) THEN
         !
         CALL iotk_write_dat( ounit, "SMEARING_TYPE", ngauss )
         !
         CALL iotk_write_attr( attr, "UNITS", trim(degauss_units), FIRST = .true. )
         !
         CALL iotk_write_dat( ounit, "SMEARING_PARAMETER", degauss , ATTR = attr )
         !
      ENDIF
      !
      CALL iotk_write_dat( ounit, "TETRAHEDRON_METHOD", ltetra )
      !
      IF ( ltetra ) THEN
         !
         CALL iotk_write_dat( ounit, "NUMBER_OF_TETRAHEDRA", ntetra )
         !
         CALL iotk_write_dat( ounit, "TETRAHEDRON_TYPE", tetra_type)
         !
         IF(tetra_type == 0) then
            !
            DO i = 1, ntetra
               !
               CALL iotk_write_dat( ounit, "TETRAHEDRON" // &
                                  & iotk_index( i ), tetra(1:4,i) )
               !
            ENDDO
            !
         END IF
         !
      ENDIF
      !
      CALL iotk_write_dat( ounit, "FIXED_OCCUPATIONS", tfixed_occ )
      !
      IF ( tfixed_occ ) THEN
         !
         CALL iotk_write_attr( attr, "lsda" , lsda, FIRST = .true. )
         CALL iotk_write_attr( attr, "nstates_up", nstates_up )
         CALL iotk_write_attr( attr, "nstates_down", nstates_dw )
         !
         CALL iotk_write_empty( ounit, 'INFO', ATTR = attr )
         !
         CALL iotk_write_dat( ounit, "INPUT_OCC_UP", input_occ(1:nstates_up,1) )
         !
         IF ( lsda ) &
            CALL iotk_write_dat( ounit, "INPUT_OCC_DOWN", input_occ(1:nstates_dw,2) )
         !
      ENDIF
      !
      CALL iotk_write_end( ounit, "OCCUPATIONS" )
      !
    END SUBROUTINE qexml_write_occ
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_bz( num_k_points, xk, wk, k1, k2, k3, &
                               nk1, nk2, nk3, k_units, qnorm, &
                               nks_start, xk_start, wk_start )
      !------------------------------------------------------------------------
      !
      INTEGER,      INTENT(in) :: num_k_points, k1, k2, k3, nk1, nk2, nk3
      REAL(DP),    INTENT(in) :: xk(:,:), wk(:)
      CHARACTER(*), INTENT(in) :: k_units
      REAL(DP), INTENT(IN) :: qnorm
      INTEGER,  INTENT(IN), OPTIONAL ::  nks_start
      REAL(DP), INTENT(IN), OPTIONAL :: xk_start(:,:), wk_start(:)
      !
      INTEGER :: ik
      !
      !
      CALL iotk_write_begin( ounit, "BRILLOUIN_ZONE" )
      !
      CALL iotk_write_dat( ounit, "NUMBER_OF_K-POINTS", num_k_points )
      !
      CALL iotk_write_attr( attr, "UNITS", trim(k_units), FIRST = .true. )
      CALL iotk_write_empty( ounit, "UNITS_FOR_K-POINTS", attr )
      !
      CALL iotk_write_attr( attr, "nk1", nk1, FIRST = .true. )
      CALL iotk_write_attr( attr, "nk2", nk2 )
      CALL iotk_write_attr( attr, "nk3", nk3 )
      CALL iotk_write_empty( ounit, "MONKHORST_PACK_GRID", attr )
      CALL iotk_write_attr( attr, "k1", k1, FIRST = .true. )
      CALL iotk_write_attr( attr, "k2", k2 )
      CALL iotk_write_attr( attr, "k3", k3 )
      CALL iotk_write_empty( ounit, "MONKHORST_PACK_OFFSET", attr )
      !
      DO ik = 1, num_k_points
         !
         CALL iotk_write_attr( attr, "XYZ", xk(:,ik), FIRST = .true. )
         !
         CALL iotk_write_attr( attr, "WEIGHT", wk(ik) )
         !
         CALL iotk_write_empty( ounit, "K-POINT" // &
                              & trim( iotk_index(ik) ), attr )
         !
      ENDDO
      !
      ! ... these are k-points and weights in the Irreducible BZ
      !
      IF (present(nks_start).and.present(xk_start).and.present(wk_start)) THEN
         !
         CALL iotk_write_dat( ounit, "STARTING_K-POINTS", nks_start )
         !
         DO ik = 1, nks_start
            !
            CALL iotk_write_attr( attr, "XYZ", xk_start(:,ik), FIRST = .TRUE. )
            !
            CALL iotk_write_attr( attr, "WEIGHT", wk_start(ik) )
            !
            CALL iotk_write_empty( ounit, "K-POINT_START" // &
                              & TRIM( iotk_index(ik) ), attr )
            !
         END DO
      ENDIF
      !
      CALL iotk_write_dat( ounit, "NORM-OF-Q", qnorm )
      !
      CALL iotk_write_end( ounit, "BRILLOUIN_ZONE" )
      !
    END SUBROUTINE qexml_write_bz
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_para( kunit, nproc, nproc_pool, nproc_image, &
                                 ntask_groups, nproc_bgrp, nproc_ortho ) 
      !------------------------------------------------------------------------
      !
      INTEGER,  INTENT(IN) :: kunit, nproc, nproc_pool, nproc_image, &
                              ntask_groups, nproc_bgrp, nproc_ortho 
      !
      !
      CALL iotk_write_begin( ounit, "PARALLELISM" )
      CALL iotk_write_dat( ounit, &
                              "GRANULARITY_OF_K-POINTS_DISTRIBUTION", kunit )
      CALL iotk_write_dat( ounit, "NUMBER_OF_PROCESSORS", nproc )
      CALL iotk_write_dat( ounit, &
                              "NUMBER_OF_PROCESSORS_PER_POOL", nproc_pool )
      CALL iotk_write_dat( ounit, &
                              "NUMBER_OF_PROCESSORS_PER_IMAGE", nproc_image )
      CALL iotk_write_dat( ounit, "NUMBER_OF_PROCESSORS_PER_TASKGROUP", &
                                              ntask_groups )
      CALL iotk_write_dat( ounit, "NUMBER_OF_PROCESSORS_PER_BAND_GROUP", &
                                              nproc_bgrp )
      CALL iotk_write_dat( ounit, "NUMBER_OF_PROCESSORS_PER_DIAGONALIZATION", &
                                              nproc_ortho )
      CALL iotk_write_end( ounit, "PARALLELISM" )
      !
      !
    END SUBROUTINE qexml_write_para
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_bands_info( num_k_points, natomwfc, &
                                       nbnd, nbnd_up, nbnd_down, &
                                       nspin, nelec, nel_up, nel_down, &
                                       energy_units, k_units, &
                                       ef, two_fermi_energies ,&
                                       ef_up, ef_down, noncolin )
      !------------------------------------------------------------------------
      !
      INTEGER,       INTENT(in) ::  num_k_points, natomwfc, nbnd, nbnd_up, nbnd_down, &
                                    nspin, nel_up, nel_down
      REAL(DP),     INTENT(in) ::   nelec
      CHARACTER(*),  INTENT(in) :: energy_units, k_units
      LOGICAL,       INTENT(in), OPTIONAL :: noncolin,two_fermi_energies
      REAL(DP),     INTENT(in), OPTIONAL :: ef,ef_up,ef_down
      !
      !
      CALL iotk_write_begin( ounit, "BAND_STRUCTURE_INFO" )
      !
      CALL iotk_write_dat  ( ounit, "NUMBER_OF_K-POINTS", num_k_points )
      !
      CALL iotk_write_dat  ( ounit, "NUMBER_OF_SPIN_COMPONENTS", nspin )
      !
      IF (present(noncolin)) CALL iotk_write_dat  ( ounit, "NON-COLINEAR_CALCULATION", noncolin )
      !
      CALL iotk_write_dat  ( ounit, "NUMBER_OF_ATOMIC_WFC", natomwfc )
      !
      IF ( nspin == 2 ) THEN
         !
         CALL iotk_write_attr( attr, "UP", nbnd_up, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "DW", nbnd_down )
         CALL iotk_write_dat( ounit, &
                              "NUMBER_OF_BANDS", nbnd, ATTR = attr )
         CALL iotk_write_attr( attr, "UP", nel_up, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "DW", nel_down )
         CALL iotk_write_dat( ounit, &
                              "NUMBER_OF_ELECTRONS", nelec, ATTR = attr )
      ELSE
         !
         CALL iotk_write_dat  ( ounit, "NUMBER_OF_BANDS", nbnd )
         CALL iotk_write_dat  ( ounit, "NUMBER_OF_ELECTRONS", nelec )
         !
      END IF
      !
      CALL iotk_write_attr ( attr, "UNITS", trim(k_units), FIRST = .TRUE. )
      CALL iotk_write_empty( ounit, "UNITS_FOR_K-POINTS", ATTR = attr )
      !
      CALL iotk_write_attr ( attr, "UNITS", trim(energy_units), FIRST = .TRUE. )
      CALL iotk_write_empty( ounit, "UNITS_FOR_ENERGIES", ATTR = attr )
      !
      !
      !
      IF (present(two_fermi_energies) ) THEN
         IF (two_fermi_energies) THEN
            !
            CALL iotk_write_dat( ounit,"TWO_FERMI_ENERGIES",two_fermi_energies)
            CALL iotk_write_dat( ounit, "ELECTRONS_UP", nel_up )
            CALL iotk_write_dat( ounit, "ELECTRONS_DOWN", nel_down )
            CALL iotk_write_dat( ounit, "FERMI_ENERGY_UP", ef_up )
            CALL iotk_write_dat( ounit, "FERMI_ENERGY_DOWN", ef_down )
            !
         ELSE
            !
            IF (present(ef)) CALL iotk_write_dat( ounit, "FERMI_ENERGY", ef )
            !
         ENDIF
      ELSE
         !
         IF (present(ef)) CALL iotk_write_dat( ounit, "FERMI_ENERGY", ef )
         !
      ENDIF
      !
      CALL iotk_write_end  ( ounit, "BAND_STRUCTURE_INFO" )
      !
      !
    END SUBROUTINE qexml_write_bands_info
    !
    !
    !------------------------------------------------------------------------  
    SUBROUTINE qexml_write_bands_pw( nbnd, num_k_points, nspin, xk, wk, wg , et, energy_units,  lkpoint_dir ,auxunit, dirname )
      !------------------------------------------------------------------------
      !
      INTEGER, INTENT(in) :: nbnd,num_k_points,nspin,auxunit
      REAL(DP), INTENT(in) :: xk(:,:),wk(:),wg(:,:),et(:,:)
      CHARACTER(len=*), INTENT(IN) :: energy_units
      LOGICAL, INTENT(in) :: lkpoint_dir
      CHARACTER(len=*), INTENT(in) :: dirname
      
      !
      REAL(DP), ALLOCATABLE :: raux(:)
      INTEGER :: ik,ispin,ik_eff
      CHARACTER(LEN=256)    :: filename
      !
      !
      CALL iotk_write_begin( ounit, "EIGENVALUES" )
      !
      ALLOCATE( raux( nbnd) )
      !
      DO ik = 1, num_k_points
         !
         !
         CALL iotk_write_begin( ounit, "K-POINT" // TRIM( iotk_index( ik ) ) )
         !
         CALL iotk_write_dat( ounit, "K-POINT_COORDS", xk(:,ik), COLUMNS=3 )
         !
         CALL iotk_write_dat( ounit, "WEIGHT", wk(ik) )
         !
         !
         IF ( nspin == 2 ) THEN
            !
            ispin = 1
            !
            IF (lkpoint_dir) THEN
               !
               filename = qexml_wfc_filename(".",'eigenval1', ik, EXTENSION='xml',&
                                     DIR=lkpoint_dir )
               !
               CALL iotk_link( ounit, "DATAFILE.1", &
                               filename, CREATE = .FALSE., BINARY = .FALSE. )
            ELSE
               CALL iotk_write_begin( auxunit, &
                             "DATA_EIG"//TRIM( iotk_index( ik ) )//"_SPIN_UP" )
            ENDIF
            !
            IF ( wk(ik) == 0.D0 ) THEN
               !
               raux = wg(:,ik)
               !
            ELSE
               !
               raux = wg(:,ik) / wk(ik)
               !
            END IF
            !
            !
            IF (lkpoint_dir) THEN
               filename = qexml_wfc_filename( dirname, 'eigenval1', ik, &
                    EXTENSION='xml',  DIR=lkpoint_dir )
               !
               CALL qexml_write_eig( auxunit, filename, nbnd, et(:, ik), &
                     trim(energy_units), OCC = raux(:), IK=ik, ISPIN=ispin )
            ELSE
               filename=' '
               CALL qexml_write_eig( auxunit, filename, nbnd, et(:, ik), &
                     trim(energy_units), OCC = raux(:), IK=ik, ISPIN=ispin,  &
                                LKPOINT_DIR=.FALSE. )
            ENDIF
            !
            ispin = 2
            !
            ik_eff = ik + num_k_points
            !
            IF (lkpoint_dir) THEN
               filename = qexml_wfc_filename( ".", 'eigenval2', ik, &
                          EXTENSION='xml',  DIR=lkpoint_dir )
               !
               CALL iotk_link( ounit, "DATAFILE.2", &
                    filename, CREATE = .FALSE., BINARY = .FALSE. )
            ELSE
               CALL iotk_write_end( auxunit, &
                    "DATA_EIG"//TRIM( iotk_index( ik ) )//"_SPIN_UP" )
               CALL iotk_write_begin( auxunit, &
                    "DATA_EIG"//TRIM( iotk_index( ik ) )//"_SPIN_DW" )
            ENDIF
            !
            IF ( wk(ik_eff) == 0.D0 ) THEN
               !
               raux = wg(:,ik_eff)
               !
            ELSE
               !
               raux = wg(:,ik_eff) / wk(ik_eff)
               !
            END IF
            !
            IF (lkpoint_dir) THEN
               filename = qexml_wfc_filename( dirname, 'eigenval2', ik, &
                    EXTENSION = 'xml',  DIR=lkpoint_dir )
               !
               CALL qexml_write_eig( auxunit, filename, nbnd, et(:, ik_eff) , &
                    trim(energy_units), OCC = raux(:), IK = ik, ISPIN = ispin)
            ELSE
               filename=' '
               CALL qexml_write_eig( auxunit, filename, nbnd, et(:, ik_eff) , &
                    trim(energy_units), OCC = raux(:), IK = ik, &
                    ISPIN = ispin, LKPOINT_DIR=.false.)
               CALL iotk_write_end( auxunit, &
                    "DATA_EIG"//TRIM( iotk_index( ik ) )//"_SPIN_DW" )
            ENDIF
               !
         ELSE
            !
            IF (lkpoint_dir) THEN
               filename = qexml_wfc_filename( ".", 'eigenval', ik, &
                    EXTENSION='xml',  DIR=lkpoint_dir )
               !
               CALL iotk_link( ounit, "DATAFILE", &
                    filename, CREATE = .FALSE., BINARY = .FALSE. )
            ELSE
               CALL iotk_write_begin( auxunit, &
                    "DATA_EIG"//TRIM( iotk_index( ik ) ) )
            ENDIF
            !
            IF ( wk(ik) == 0.D0 ) THEN
               !
               raux(:) = wg(:,ik)
               !
            ELSE
               !
               raux(:) = wg(:,ik) / wk(ik)
               !
            END IF
            !
            IF (lkpoint_dir) THEN
               filename = qexml_wfc_filename( dirname, 'eigenval', ik, &
                    EXTENSION='xml',  DIR=lkpoint_dir )
               !
               CALL qexml_write_eig( auxunit, filename, nbnd, et(:, ik), &
                               trim(energy_units), OCC = raux(:), IK = ik )
            ELSE
               filename=' '
               CALL qexml_write_eig( auxunit, filename, nbnd, et(:, ik), &
                    trim(energy_units), OCC = raux(:), IK = ik, &
                    LKPOINT_DIR=.false. )
               CALL iotk_write_end( auxunit, &
                    "DATA_EIG"//TRIM( iotk_index( ik ) ) )
            ENDIF
            !
         END IF
         !
         CALL iotk_write_end( ounit, "K-POINT" // TRIM( iotk_index( ik ) ) )
         !
      ENDDO
      !
      !
      DEALLOCATE ( raux )
      !
      !
      CALL iotk_write_end( ounit, "EIGENVALUES" )
      !
    END SUBROUTINE qexml_write_bands_pw
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_bands_cp( nbnd, num_k_points, nspin, iupdwn, nupdwn, xk, wk, et, tksw, &
         occ0, occm, energy_units, k_units, auxunit ,dirname )
      !------------------------------------------------------------------------
      !
      !
      INTEGER, INTENT(in) :: nbnd,num_k_points,nspin, iupdwn(2),nupdwn(2),auxunit
      REAL(DP), INTENT(in) :: xk(:,:),wk(:),et(:,:)
      CHARACTER(len=*), INTENT(in) :: dirname,k_units,energy_units
      LOGICAL, INTENT(in) :: tksw
      REAL(DP), INTENT(in) :: occ0(:)
      REAL(DP), INTENT(in) :: occm(:)
      !
      !
      REAL(DP), ALLOCATABLE :: dtmp(:)
      INTEGER :: iss, ik
      CHARACTER(LEN=4)     :: cspin
      CHARACTER(LEN=256)    :: filename
      !
      !
      CALL iotk_write_begin( ounit, "EIGENVALUES" )
      !
      DO ik = 1, num_k_points
         !
         CALL iotk_write_begin( ounit, "K-POINT" // TRIM( iotk_index(ik) ) )
         !
         CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
         CALL iotk_write_dat( ounit, "K-POINT_COORDS", xk(:,ik), ATTR = attr )
         !
         CALL iotk_write_dat( ounit, "WEIGHT", wk(ik) )
         !
         ALLOCATE( dtmp ( nbnd ) )
         !
         DO iss = 1, nspin
            !
            cspin = iotk_index( iss )
            !
            dtmp = 0.0d0
            !
            IF( tksw ) THEN
               !
               !
               IF( nspin == 2 ) THEN
                  IF( iss == 1 ) filename = qexml_wfc_filename( ".", 'eigenval1', ik, EXTENSION='xml' )
                  IF( iss == 2 ) filename = qexml_wfc_filename( ".", 'eigenval2', ik, EXTENSION='xml' )
                  !
                  IF( iss == 1 ) CALL iotk_link( ounit, "DATAFILE.1", &
                       filename, CREATE = .FALSE., BINARY = .FALSE. )
                  IF( iss == 2 ) CALL iotk_link( ounit, "DATAFILE.2", &
                       filename, CREATE = .FALSE., BINARY = .FALSE. )
                  
                  IF( iss == 1 ) filename = qexml_wfc_filename( dirname, 'eigenval1', ik, EXTENSION='xml' )
                  IF( iss == 2 ) filename = qexml_wfc_filename( dirname, 'eigenval2', ik, EXTENSION='xml' )
               ELSE
                  filename = qexml_wfc_filename( ".", 'eigenval', ik, EXTENSION='xml' )
                  CALL iotk_link( ounit, "DATAFILE", filename, CREATE = .FALSE., BINARY = .FALSE. )
                  filename = qexml_wfc_filename( dirname, 'eigenval', ik, EXTENSION='xml' )
               END IF
               
               dtmp ( 1:nupdwn( iss ) ) = occ0( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) / wk(ik)
               !
               CALL qexml_write_eig( auxunit, filename, nbnd, et( 1:nbnd, iss) , energy_units, &
                    OCC = dtmp(:), IK=ik, ISPIN=iss )
            END IF
               !
            CALL iotk_write_dat( ounit, "OCC0"  // TRIM( cspin ), &
                                    occ0( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) )
            !
            CALL iotk_write_dat( ounit, "OCCM" // TRIM( cspin ), &
                 occm( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) )
            !
         END DO
         !
         DEALLOCATE( dtmp )
         !
         CALL iotk_write_end( ounit, "K-POINT" // TRIM( iotk_index(ik) ) )
         !
      END DO
      !
      CALL iotk_write_end( ounit, "EIGENVALUES" )
      !
      !
    END SUBROUTINE qexml_write_bands_cp
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_eig( iuni, filename, nbnd, eig, energy_units, &
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
      !
      !
      IF (lkpoint_dir0) CALL iotk_open_write ( iuni, &
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
      !
    END SUBROUTINE qexml_write_eig
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_wfc( nbnd, nkpts, nspin, ik, ispin, ipol, igk, ngw, igwx, &
                                gamma_only, wf, wf_kindip, scale_factor )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER,                INTENT(in) :: nbnd, nkpts, nspin
      INTEGER,                INTENT(in) :: ik
      INTEGER,      OPTIONAL, INTENT(in) :: ispin, ipol
      INTEGER,                INTENT(in) :: ngw, igwx
      LOGICAL,                INTENT(in) :: gamma_only
      INTEGER,      OPTIONAL, INTENT(in) :: igk(:)
      COMPLEX(DP), OPTIONAL, INTENT(in) :: wf(:,:)
      COMPLEX(DP), OPTIONAL, INTENT(in) :: wf_kindip(:,:)
      REAL(DP),    OPTIONAL, INTENT(in) :: scale_factor
      !
      INTEGER         :: iunaux, ierr
      INTEGER         :: ig, ib
      CHARACTER(256)  :: filename
      COMPLEX(DP),  ALLOCATABLE :: wtmp(:)

      ierr = 0
      !
      IF ( present( ispin ) .and. present( ipol )  ) THEN
         !
         ierr = 1
         RETURN
         !
      ENDIF
      !
      ! open the file to write
      !
      CALL iotk_free_unit( iunaux )
      !
      IF ( present( ispin ) ) THEN
         !
         filename = trim( qexml_wfc_filename( datadir_out, 'evc', ik, ispin ) )
         !
      ELSEIF ( present( ipol )  ) THEN
         !
         filename = trim( qexml_wfc_filename( datadir_out, 'evc', ik, ipol ) )
         !
      ELSE
         !
         filename = trim( qexml_wfc_filename( datadir_out, 'evc', ik ) )
         !
      ENDIF
      !
      CALL iotk_open_write ( iunaux, FILE = trim(filename), ROOT="WFC", BINARY=.true., IERR=ierr )
      IF (ierr/=0)  RETURN
      !
      !
      CALL iotk_write_attr( attr, "ngw",          ngw, FIRST = .true. )
      CALL iotk_write_attr( attr, "igwx",         igwx )
      CALL iotk_write_attr( attr, "gamma_only",   gamma_only )
      CALL iotk_write_attr( attr, "nbnd",         nbnd )
      CALL iotk_write_attr( attr, "ik",           ik )
      CALL iotk_write_attr( attr, "nk",           nkpts )
      CALL iotk_write_attr( attr, "ispin",        ispin )
      CALL iotk_write_attr( attr, "nspin",        nspin )
      IF ( present( scale_factor) ) CALL iotk_write_attr( attr, "scale_factor", scale_factor )
      !
      CALL iotk_write_empty( iunaux, "INFO", attr )
      !
      !
      IF ( present( wf ) ) THEN
         !
         ! write wfcs without any G-reordering
         !
         DO ib = 1, nbnd
            !
            CALL iotk_write_dat( iunaux, "evc" // trim(iotk_index( ib )), wf( 1: ngw, ib) )
            !
         ENDDO
         !
      ENDIF
      !
      !
      IF ( present( wf_kindip ) ) THEN
         !
         ! we need to reorder wfcs in terms of G-vectors
         ! we need the igk map
         !
         IF ( .not. present( igk ) ) THEN
            ierr = 71
            RETURN
         ENDIF
         !
         ALLOCATE( wtmp( ngw ) )
         !
         DO ib = 1, nbnd
            !
            DO ig = 1, ngw
               !
               wtmp( ig ) = wf_kindip( igk(ig), ib)
               !
            ENDDO
            !
            CALL iotk_write_dat( iunaux, "evc" // trim(iotk_index( ib )), wtmp( 1: ngw) )
            !
         ENDDO
         !
         DEALLOCATE( wtmp )
         !
      ENDIF
      !
      !
      CALL iotk_close_write( iunaux )
      !
    END SUBROUTINE qexml_write_wfc
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_rho( nr1, nr2, nr3, rho, nr1x, nr2x, rhov, binary )
      !------------------------------------------------------------------------
      !
      ! Writes charge density rho, one plane at a time.
      !
      IMPLICIT NONE
      !
      INTEGER,             INTENT(in) :: nr1, nr2, nr3
      INTEGER,   OPTIONAL, INTENT(in) :: nr1x, nr2x
      REAL(DP), OPTIONAL, INTENT(in) :: rho(:,:,:), rhov(:)
      LOGICAL,   OPTIONAL, INTENT(in) :: binary
      !
      INTEGER        :: iunaux, nr1x_, nr2x_, ip, i1, i2, i
      LOGICAL        :: binary_
      CHARACTER(256) :: filename
      REAL(DP), ALLOCATABLE :: plane(:,:)
      !
      !
      CALL iotk_free_unit( iunaux )
      !
      binary_ = .true.
      IF ( present (binary) ) binary_ = binary
      !
      IF ( binary_ ) THEN
         !
         filename = trim( datadir_out ) // '/' //'charge-density.dat'
         !
      ELSE
         !
         filename = trim( datadir_out ) // '/' //'charge-density.xml'
         !
      ENDIF
      !
      CALL iotk_open_write( iunaux, FILE = trim(filename), BINARY=binary_ )
      !
      !
      CALL iotk_write_begin( iunaux, "CHARGE-DENSITY" )
      !
      CALL iotk_write_attr( attr, "nr1", nr1, FIRST = .true. )
      CALL iotk_write_attr( attr, "nr2", nr2 )
      CALL iotk_write_attr( attr, "nr3", nr3 )
      !
      CALL iotk_write_empty( iunaux, "INFO", attr )
      !
      !
      IF ( present( rho ) ) THEN
         !
         DO ip = 1, nr3
            !
            CALL iotk_write_dat( iunaux, "z"//trim(iotk_index(ip)), rho(1:nr1,1:nr2,ip) )
            !
         ENDDO
         !
      ELSEIF ( present( rhov ) ) THEN
         !
         nr1x_ = nr1
         IF ( present( nr1x )) nr1x_ = nr1x
         nr2x_ = nr2
         IF ( present( nr2x )) nr2x_ = nr2x
         !
         IF ( nr1x_ /= nr1 .or. nr2x_ /= nr2 ) THEN
            !
            ! we need to separately reconstruct the rho-plane
            !
            ALLOCATE( plane(nr1, nr2 ) )
            !
            DO ip = 1, nr3
               !
               DO i2 = 1, nr2
               DO i1 = 1, nr1
                   !
                   i = (nr1x_ * nr2x_) * ( ip -1 ) + nr1x_ * ( i2 -1 ) + i1
                   !
                   plane( i1, i2) = rhov( i )
                   !
               ENDDO
               ENDDO
               !
               CALL iotk_write_dat( iunaux, "z"//trim(iotk_index(ip)), plane )
               !
            ENDDO
            !
            DEALLOCATE( plane )
            !
         ELSE
            !
            DO ip = 1, nr3
               !
               i1  = ( nr1 * nr2 ) * ( ip -1 ) + 1
               i2  = ( nr1 * nr2 ) * ip
               !
               CALL iotk_write_dat( iunaux, "z"//trim(iotk_index(ip)), rhov(i1:i2) )
               !
            ENDDO
            !
         ENDIF
         !
      ENDIF
      !
      !
      CALL iotk_write_end( iunaux, "CHARGE-DENSITY" )
      !
      CALL iotk_close_write( iunaux )
      !
      !
    END SUBROUTINE qexml_write_rho
!
!-------------------------------------------
! ... read subroutines
!-------------------------------------------
!
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_header( creator_name, creator_version, &
                                  format_name, format_version, ierr )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      CHARACTER(len=*),  OPTIONAL, INTENT(out) :: creator_name, creator_version
      CHARACTER(len=*),  OPTIONAL, INTENT(out) :: format_name, format_version
      INTEGER,           INTENT(out) :: ierr
      !
      CHARACTER(256) :: creator_name_, creator_version_
      CHARACTER(256) :: format_name_,     format_version_
      !
      ierr = 0
      !
      !
      CALL iotk_scan_begin( iunit, "HEADER", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_empty( iunit, "FORMAT", ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_attr(attr, "NAME", format_name_, IERR=ierr)
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr(attr, "VERSION", format_version_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_empty( iunit, "CREATOR", ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_attr(attr, "NAME", creator_name_, IERR=ierr)
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr(attr, "VERSION", creator_version_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_end( iunit, "HEADER", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( present(creator_name) )     creator_name    = trim(creator_name_)
      IF ( present(creator_version) )  creator_version = trim(creator_version_)
      IF ( present(format_name) )      format_name     = trim(format_name_)
      IF ( present(format_version) )   format_version  = trim(format_version_)
      !
    END SUBROUTINE qexml_read_header
    !
    !
    !------------------------------------------------------------------------  
    SUBROUTINE qexml_read_status_cp( nfi,simtime,time_units,title, &
                                  ekin, eht, esr, eself, epseu, enl, exc, vave, enthal, &
                                  energy_units, found, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER, OPTIONAL, INTENT(OUT) :: nfi
      REAL(DP), OPTIONAL, INTENT(OUT) :: simtime, ekin,eht,esr,eself,epseu,enl,exc,vave,enthal
      CHARACTER(len=*), OPTIONAL, INTENT(OUT) :: time_units, title, energy_units
      LOGICAL,INTENT(OUT) :: found
      INTEGER, INTENT(OUT) :: ierr
      !
      INTEGER :: nfi_
      REAL(DP) :: simtime_, ekin_,eht_,esr_,eself_,epseu_,enl_,exc_,vave_,enthal_
      CHARACTER(len=256) :: time_units_, title_, energy_units_
      !
      CALL iotk_scan_begin( iunit, "STATUS", ATTR=attr, FOUND = found )
      IF ( .NOT.found ) RETURN
      !
      CALL iotk_scan_empty( iunit, "STEP", ATTR = attr, IERR = ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_attr( attr, "ITERATION", nfi_, IERR = ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "TIME", simtime_, ATTR = attr, IERR = ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_attr( attr, "UNITS", time_units_, IERR = ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "TITLE", title_, IERR = ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "KINETIC_ENERGY", ekin_,   ATTR = attr, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_attr( attr, "UNITS", energy_units_, IERR = ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "HARTREE_ENERGY", eht_,   IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "EWALD_TERM",     esr_,   IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "GAUSS_SELFINT",  eself_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "LPSP_ENERGY",    epseu_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "NLPSP_ENERGY",   enl_,   IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "EXC_ENERGY",     exc_,   IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "AVERAGE_POT",    vave_,  IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "ENTHALPY",       enthal_,IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_end( iunit, "STATUS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF (present(nfi)) nfi = nfi_
      IF (present(simtime)) simtime = simtime_
      IF (present(time_units)) time_units = time_units_
      IF (present(title)) title = title_
      IF (present(ekin)) ekin = ekin_
      IF (present(eht)) eht = eht_
      IF (present(esr)) esr = esr_
      IF (present(eself)) eself = eself_
      IF (present(epseu)) epseu = epseu_
      IF (present(enl)) enl = enl_
      IF (present(exc)) exc = exc_
      IF (present(vave)) vave = vave_
      IF (present(enthal)) enthal = enthal_
      IF (present(energy_units)) energy_units = energy_units_
      !
    END SUBROUTINE qexml_read_status_cp
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_cell( bravais_lattice, celldm, alat, &
                                a1, a2, a3, b1, b2, b3, alat_units, a_units, b_units, es_corr, ierr )
      !------------------------------------------------------------------------
      !
      CHARACTER(len=*),  OPTIONAL, INTENT(out) :: bravais_lattice
      REAL(DP),         OPTIONAL, INTENT(out) :: celldm(6), alat
      REAL(DP),         OPTIONAL, INTENT(out) :: a1(3), a2(3), a3(3)
      REAL(DP),         OPTIONAL, INTENT(out) :: b1(3), b2(3), b3(3)
      CHARACTER(len=*),  OPTIONAL, INTENT(out) :: alat_units, a_units, b_units
      CHARACTER(len=*),  OPTIONAL, INTENT(out) :: es_corr
      INTEGER,                     INTENT(out) :: ierr
      !
      CHARACTER(256)     :: bravais_lattice_
      CHARACTER(256)     :: alat_units_, a_units_, b_units_,es_corr_
      REAL(DP)          :: celldm_(6), alat_
      REAL(DP)          :: a1_(3), a2_(3), a3_(3)
      REAL(DP)          :: b1_(3), b2_(3), b3_(3)
      !

      ierr=0
      !
      !
      CALL iotk_scan_begin( iunit, "CELL" )
      !
      CALL iotk_scan_dat( iunit, "BRAVAIS_LATTICE", bravais_lattice_, IERR=ierr )
      IF ( ierr /= 0 ) RETURN
      !
      !
      CALL iotk_scan_dat( iunit, "NON-PERIODIC_CELL_CORRECTION", es_corr_, IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      !
      CALL iotk_scan_dat( iunit, "LATTICE_PARAMETER", alat_, ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      CALL iotk_scan_attr( attr, "UNITS", alat_units_, IERR=ierr )
      IF ( ierr /= 0 ) RETURN
      !
      !
      CALL iotk_scan_dat( iunit, "CELL_DIMENSIONS", celldm_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      CALL iotk_scan_begin( iunit, "DIRECT_LATTICE_VECTORS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      CALL iotk_scan_empty( iunit, "UNITS_FOR_DIRECT_LATTICE_VECTORS", &
                            ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "UNITS", a_units_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_dat(   iunit, "a1", a1_(:), ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_dat(   iunit, "a2", a2_(:), IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_dat(   iunit, "a3", a3_(:), IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_end(   iunit, "DIRECT_LATTICE_VECTORS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_begin( iunit, "RECIPROCAL_LATTICE_VECTORS", IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_empty( iunit, "UNITS_FOR_RECIPROCAL_LATTICE_VECTORS", &
                            ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "UNITS", b_units_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_dat(   iunit, "b1", b1_(:), ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_dat(   iunit, "b2", b2_(:), IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_dat(   iunit, "b3", b3_(:), IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_end(   iunit, "RECIPROCAL_LATTICE_VECTORS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_end( iunit, "CELL", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( present(bravais_lattice) )  bravais_lattice = bravais_lattice_
      IF ( present(celldm) )        celldm       = celldm_
      IF ( present(alat) )          alat         = alat_
      IF ( present(a1) )            a1           = a1_
      IF ( present(a2) )            a2           = a2_
      IF ( present(a3) )            a3           = a3_
      IF ( present(b1) )            b1           = b1_
      IF ( present(b2) )            b2           = b2_
      IF ( present(b3) )            b3           = b3_
      IF ( present(alat_units) )    alat_units   = trim(alat_units_)
      IF ( present(a_units) )       a_units      = trim(a_units_)
      IF ( present(b_units) )       b_units      = trim(b_units_)
      IF ( present(es_corr) )       es_corr      = trim(es_corr_)
      !
    END SUBROUTINE qexml_read_cell
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_moving_cell(lmovecell, cell_factor, ierr)
      !------------------------------------------------------------------------
      !
      LOGICAL,  INTENT(OUT) :: lmovecell
      REAL(DP), INTENT(OUT) :: cell_factor
      INTEGER,  INTENT(OUT) :: ierr
      !
      LOGICAL :: found
      !
      CALL iotk_scan_begin( iunit, "MOVING_CELL", found=lmovecell, IERR=ierr )
      !
      IF (lmovecell) THEN
         CALL iotk_scan_dat( iunit, "CELL_FACTOR", cell_factor, IERR=ierr)
         IF (ierr/=0) RETURN
         !
         CALL iotk_scan_end( iunit, "MOVING_CELL", IERR=ierr )
         IF (ierr/=0) RETURN
      END IF
      !
    END SUBROUTINE qexml_read_moving_cell
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_ions( nsp, nat, atm, ityp, psfile, amass, amass_units, &
                                tau, tau_units, if_pos, pseudo_dir, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,          OPTIONAL, INTENT(out) :: nsp, nat
      INTEGER,          OPTIONAL, INTENT(out) :: ityp(:)
      CHARACTER(len=*), OPTIONAL, INTENT(out) :: atm(:)
      CHARACTER(len=*), OPTIONAL, INTENT(out) :: psfile(:)
      REAL(DP),         OPTIONAL, INTENT(out) :: amass(:)
      CHARACTER(len=*), OPTIONAL, INTENT(out) :: amass_units
      REAL(DP),         OPTIONAL, INTENT(out) :: tau(:,:)
      INTEGER,          OPTIONAL, INTENT(out) :: if_pos(:,:)
      CHARACTER(len=*), OPTIONAL, INTENT(out) :: tau_units
      CHARACTER(len=*), OPTIONAL, INTENT(out) :: pseudo_dir
      INTEGER,                    INTENT(out) :: ierr
      !
      INTEGER                     :: nat_, nsp_
      CHARACTER(256)              :: tau_units_, amass_units_
      INTEGER,        ALLOCATABLE :: ityp_(:)
      CHARACTER(3),   ALLOCATABLE :: atm_(:)
      CHARACTER(256), ALLOCATABLE :: psfile_(:)
      CHARACTER(256)              :: pseudo_dir_
      REAL(DP),       ALLOCATABLE :: amass_(:)
      REAL(DP),       ALLOCATABLE :: tau_(:,:)
      INTEGER,        ALLOCATABLE :: if_pos_(:,:)
      !
      INTEGER :: i
      !
      ierr=0
      !
      CALL iotk_scan_begin( iunit, "IONS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "NUMBER_OF_ATOMS", nat_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "NUMBER_OF_SPECIES", nsp_ )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_empty( iunit, "UNITS_FOR_ATOMIC_MASSES", ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "UNITS", amass_units_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      IF ( present(nat) )   nat = nat_
      IF ( present(nsp) )   nsp = nsp_
      !
      ALLOCATE( atm_(nsp_) )
      ALLOCATE( amass_(nsp_) )
      ALLOCATE( psfile_(nsp_) )
      !
      DO i = 1, nsp_
         !
         IF ( qexml_version_before_1_4_0 ) THEN
            !
            CALL iotk_scan_dat( iunit, "ATOM_TYPE", atm_(i), IERR=ierr )
            IF (ierr/=0) RETURN
            CALL iotk_scan_dat( iunit, trim( atm_(i) ) // "_MASS", amass_(i), IERR=ierr )
            IF (ierr/=0) RETURN
            CALL iotk_scan_dat( iunit, "PSEUDO_FOR_" // trim( atm_(i) ), &
                                psfile_(i), IERR=ierr )
            IF (ierr/=0) RETURN
            !
         ELSE
            !
            ! current version
            !
            CALL iotk_scan_begin( iunit, "SPECIE"//trim(iotk_index(i)), IERR=ierr )
            IF (ierr/=0) RETURN
            !
            CALL iotk_scan_dat( iunit, "ATOM_TYPE", atm_(i), IERR=ierr )
            IF (ierr/=0) RETURN
            CALL iotk_scan_dat( iunit, "MASS", amass_(i), IERR=ierr )
            IF (ierr/=0) RETURN
            CALL iotk_scan_dat( iunit, "PSEUDO", psfile_(i), IERR=ierr )
            IF (ierr/=0) RETURN
            !
            CALL iotk_scan_end( iunit, "SPECIE"//trim(iotk_index(i)), IERR=ierr )
            IF (ierr/=0) RETURN
            !
         ENDIF
         !
      ENDDO
      !
      CALL iotk_scan_dat( iunit, "PSEUDO_DIR", pseudo_dir_ )
      !
      CALL iotk_scan_empty( iunit, "UNITS_FOR_ATOMIC_POSITIONS", ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "UNITS", tau_units_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      ALLOCATE( ityp_(nat_) )
      ALLOCATE( tau_(3,nat_) )
      ALLOCATE( if_pos_(3,nat_) )
      !
      DO i = 1, nat_
         !
         CALL iotk_scan_empty( iunit, &
                               "ATOM" // trim( iotk_index(i) ), ATTR=attr, IERR=ierr )
         IF (ierr/=0) RETURN
         !
         CALL iotk_scan_attr( attr, "INDEX",  ityp_(i), IERR=ierr )
         IF (ierr/=0) RETURN
         CALL iotk_scan_attr( attr, "tau",    tau_(:,i), IERR=ierr )
         IF (ierr/=0) RETURN
         CALL iotk_scan_attr( attr, "if_pos", if_pos_(:,i), IERR=ierr )
         IF (ierr/=0) RETURN
         !
      ENDDO
      !
      CALL iotk_scan_end( iunit, "IONS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( present(nsp) )         nsp    = nsp_
      IF ( present(nat) )         nat    = nat_
      IF ( present(atm) )         atm(1:nsp_)    = atm_
      IF ( present(amass) )       amass(1:nsp_)  = amass_
      IF ( present(amass_units) ) amass_units    = trim(amass_units_)
      IF ( present(psfile) )      psfile(1:nsp_) = psfile_(1:nsp_)
      IF ( present(ityp) )        ityp(1:nat_)   = ityp_
      IF ( present(tau_units) )   tau_units      = trim(tau_units_)
      IF ( present(tau) )         tau(1:3, 1:nat_)    = tau_
      IF ( present(if_pos) )      if_pos(1:3, 1:nat_) = if_pos_
      IF ( present(pseudo_dir) )  pseudo_dir = pseudo_dir_
      !
      DEALLOCATE( atm_ )
      DEALLOCATE( amass_ )
      DEALLOCATE( psfile_ )
      DEALLOCATE( ityp_ )
      DEALLOCATE( tau_ )
      DEALLOCATE( if_pos_ )
      !
    END SUBROUTINE qexml_read_ions
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_magnetization(starting_magnetization, angle1, angle2, &
                                   nsp, two_fermi_energies, i_cons, mcons, bfield, &
                                   ef_up, ef_dw, nelup, neldw, lambda, energy_units, found, ierr)
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      INTEGER,      OPTIONAL, INTENT(OUT) :: nsp, i_cons
      REAL(DP),     OPTIONAL, INTENT(OUT) :: starting_magnetization(:), &
                                             angle1(:), angle2(:), mcons(:,:), &
                                             bfield(:), ef_up, ef_dw, nelup, neldw, lambda
      LOGICAL,      OPTIONAL, INTENT(OUT) :: two_fermi_energies
      LOGICAL,      OPTIONAL, INTENT(OUT) :: found
      CHARACTER(*), OPTIONAL, INTENT(OUT) :: energy_units
      INTEGER,                INTENT(OUT) :: ierr
      !
      INTEGER  :: i
      INTEGER  :: nsp_, i_cons_
      LOGICAL  :: two_fermi_energies_,found_,found2
      REAL(DP) :: ef_up_, ef_dw_, nelup_, neldw_, lambda_, bfield_(3)
      REAL(DP), ALLOCATABLE :: angle1_(:), angle2_(:)
      REAL(DP), ALLOCATABLE :: mcons_(:,:), starting_magnetization_(:) 
      CHARACTER(256) :: energy_units_
      !
      !
      CALL iotk_scan_begin( iunit, "MAGNETIZATION_INIT", IERR=ierr, FOUND = found_ )
      !
      IF (found_) THEN
         !
         CALL iotk_scan_dat(iunit,"CONSTRAINT_MAG", i_cons_)
         !
         CALL iotk_scan_dat( iunit, "NUMBER_OF_SPECIES", nsp_ )
         !
         ALLOCATE( starting_magnetization_(nsp_) )
         ALLOCATE( angle1_(nsp_) )
         ALLOCATE( angle2_(nsp_) )

         IF ( i_cons_ ==1 .OR. i_cons_ ==2 ) ALLOCATE( mcons_(3,nsp_) )
         !
         DO i=1,nsp_
            !
            CALL iotk_scan_begin( iunit, "SPECIE"//TRIM(iotk_index(i)) )
            !
            CALL iotk_scan_dat( iunit, "STARTING_MAGNETIZATION", &
                 starting_magnetization_(i) )
            CALL iotk_scan_dat( iunit, "ANGLE1", angle1_(i) )
            CALL iotk_scan_dat( iunit, "ANGLE2", angle2_(i) )
            !
            !
            IF (i_cons_==1.OR.i_cons_==2) THEN
               !
               CALL iotk_scan_dat( iunit, "CONSTRANT_1", mcons_(1,i) )
               CALL iotk_scan_dat( iunit, "CONSTRANT_2", mcons_(2,i) )
               CALL iotk_scan_dat( iunit, "CONSTRANT_3", mcons_(3,i) )
               !
            ENDIF
            !
            CALL iotk_scan_end( iunit, "SPECIE"//TRIM(iotk_index(i)) )
            !
         ENDDO
         !
         IF ( i_cons_ ==1 .OR. i_cons_ ==2 ) THEN
            !
            mcons_(1:3,1:nsp_) = mcons_
            !
            DEALLOCATE( mcons_ )
            !
         ENDIF
         !
         IF (i_cons_==3) THEN
            !
            ALLOCATE( mcons_(3,1) )
            !
            CALL iotk_scan_dat( iunit, "FIXED_MAGNETIZATION_1", mcons_(1,1) )
            CALL iotk_scan_dat( iunit, "FIXED_MAGNETIZATION_2", mcons_(2,1) )
            CALL iotk_scan_dat( iunit, "FIXED_MAGNETIZATION_3", mcons_(3,1) )
            !
            IF (present(mcons) ) mcons(1:3,1:1) = mcons_
            !
            DEALLOCATE( mcons_)
            !
         ELSE IF (i_cons_==4) THEN
            ! 
            CALL iotk_scan_dat( iunit, "MAGNETIC_FIELD_1", bfield_(1) )
            CALL iotk_scan_dat( iunit, "MAGNETIC_FIELD_2", bfield_(2) )
            CALL iotk_scan_dat( iunit, "MAGNETIC_FIELD_3", bfield_(3) )
            !
            IF (present(bfield)) bfield(1:3) = bfield_(1:3)
            !
         ENDIF
         !
         CALL iotk_scan_dat(iunit,"TWO_FERMI_ENERGIES", &
              two_fermi_energies_,FOUND=found2 )
         IF ( .not. found2 ) two_fermi_energies_=.FALSE.
         !
         IF (two_fermi_energies_) THEN
            !
            CALL iotk_scan_empty( iunit, "UNITS_FOR_ENERGIES", ATTR=attr, IERR=ierr )
            IF (ierr/=0) RETURN
            CALL iotk_scan_attr( attr, "UNITS", energy_units_, IERR=ierr )
            !
            ALLOCATE( mcons_(3,1) )
            !
            CALL iotk_scan_dat( iunit, "FIXED_MAGNETIZATION", mcons_(3,1) )
            CALL iotk_scan_dat( iunit, "ELECTRONS_UP", nelup_ )
            CALL iotk_scan_dat( iunit, "ELECTRONS_DOWN", neldw_ )
            CALL iotk_scan_dat( iunit, "FERMI_ENERGY_UP", ef_up_ )
            CALL iotk_scan_dat( iunit, "FERMI_ENERGY_DOWN", ef_dw_ )
            !
            IF (present(mcons) ) mcons(3,1) = mcons_(3,1)
            IF (present(ef_up) ) ef_up = ef_up_
            IF (present(ef_dw) ) ef_dw = ef_dw_
            IF (present(nelup) ) nelup = nelup_
            IF (present(neldw) ) neldw = neldw_
            IF (present(energy_units) ) energy_units = trim(energy_units_)
            !
            DEALLOCATE( mcons_)
            !
         ENDIF
         !
         lambda_ = 0.0d0
         IF (i_cons_ > 0) CALL iotk_scan_dat(iunit,"LAMBDA",lambda_)
         !
         CALL iotk_scan_end( iunit, "MAGNETIZATION_INIT" )
         !
         IF (present(nsp)) nsp = nsp_
         IF (present(two_fermi_energies)) two_fermi_energies = two_fermi_energies_
         IF (present(i_cons)) i_cons = i_cons_
         !
         IF (present(lambda) ) lambda = lambda_
         IF (present(starting_magnetization) ) starting_magnetization(1:nsp_) = starting_magnetization_
         IF (present(angle1) ) angle1(1:nsp_) = angle1_(1:nsp_)
         IF (present(angle2) ) angle2(1:nsp_) = angle2_(1:nsp_)
         !
      END IF
      !
      IF (present(found)) found = found_
      IF ( (.NOT. present(found)) .AND. ( .NOT. found_) ) ierr = 1
      !
      !
    END SUBROUTINE qexml_read_magnetization
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_symmetry( nsym, nrot, invsym, noinv, time_reversal, no_t_rev, &
                                    trasl, s, sname, s_units, t_rev, &
                                    irt, nat, found, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,          OPTIONAL, INTENT(out) :: nsym, nrot
      LOGICAL,          OPTIONAL, INTENT(out) :: invsym, noinv, time_reversal, no_t_rev
      INTEGER,          OPTIONAL, INTENT(out) :: s(:,:,:)
      REAL(DP),         OPTIONAL, INTENT(out) :: trasl(:,:)
      CHARACTER(len=*), OPTIONAL, INTENT(out) :: sname(:)
      CHARACTER(len=*), OPTIONAL, INTENT(out) :: s_units
      INTEGER,          OPTIONAL, INTENT(out) :: t_rev(:)
      INTEGER,          OPTIONAL, INTENT(out) :: irt(:,:), nat
      LOGICAL,                    INTENT(out) :: found
      INTEGER,                    INTENT(out) :: ierr
      !
      INTEGER              :: nsym_
      INTEGER              :: nrot_
      CHARACTER(256)       :: sname_(48), s_units_
      LOGICAL              :: invsym_, noinv_, time_reversal_, no_t_rev_
      INTEGER              :: s_(3,3,48)
      REAL(DP)             :: trasl_(3,48)
      INTEGER              :: t_rev_(48)
      INTEGER              :: nat_
      INTEGER, ALLOCATABLE :: irt_(:,:)
      INTEGER              :: i
      LOGICAL              :: found_tmp
      !
      !
      ierr=0
      !
      !
      CALL iotk_scan_begin( iunit, "SYMMETRIES", FOUND=found ,IERR=ierr )
      IF ((ierr/=0).OR.(.NOT.found)) RETURN
      !
      CALL iotk_scan_dat( iunit, "NUMBER_OF_SYMMETRIES", nsym_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "NUMBER_OF_BRAVAIS_SYMMETRIES", &
                                 nrot_, FOUND=found_tmp, IERR=ierr )
      IF (ierr/=0) RETURN
      IF (.NOT. found_tmp) nrot_ = nsym_
      !
      CALL iotk_scan_dat( iunit, "INVERSION_SYMMETRY", invsym_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "DO_NOT_USE_TIME_REVERSAL", &
                                  noinv_, FOUND = found_tmp, IERR=ierr )
      IF (ierr/=0) RETURN
      IF (.NOT. found_tmp) noinv_ = .FALSE.
      !
      CALL iotk_scan_dat( iunit, "TIME_REVERSAL_FLAG", &
                                  time_reversal_, FOUND = found_tmp, IERR=ierr )
      IF (ierr/=0) RETURN
      IF (.NOT. found_tmp) time_reversal_ = .TRUE.
      !
      CALL iotk_scan_dat( iunit, "NO_TIME_REV_OPERATIONS", &
                                   no_t_rev_, FOUND = found_tmp, IERR=ierr )
      IF (ierr/=0) RETURN
      IF (.NOT. found_tmp) no_t_rev_ = .FALSE.
      !
      CALL iotk_scan_dat( iunit, "NUMBER_OF_ATOMS", nat_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      ALLOCATE( irt_(48, nat_) )
      !
      CALL iotk_scan_empty( iunit, "UNITS_FOR_SYMMETRIES", ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "UNITS", s_units_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      DO i = 1, nsym_
          !
          CALL iotk_scan_begin( iunit, "SYMM"//trim( iotk_index( i ) ), IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_empty( iunit, "INFO", ATTR=attr, IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_attr( attr, "NAME", sname_(i), IERR=ierr )
          IF (ierr/=0) RETURN
          CALL iotk_scan_attr( attr, "T_REV", t_rev_(i), IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_dat( iunit, "ROTATION", s_(1:3,1:3,i), IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_dat( iunit, "FRACTIONAL_TRANSLATION", trasl_(1:3,i), IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_dat( iunit, "EQUIVALENT_IONS", irt_(i,1:nat_), IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_end( iunit, "SYMM"//trim( iotk_index( i ) ), IERR=ierr )
          IF (ierr/=0) RETURN
          !
      ENDDO
      !
      DO i = nsym_+1, nrot_
         !    
         CALL iotk_scan_begin( iunit, "SYMM" // TRIM( iotk_index( i ) ), IERR=ierr )
         IF (ierr/=0) RETURN
         !
         CALL iotk_scan_empty( iunit, "INFO", ATTR = attr, IERR=ierr )
         IF (ierr/=0) RETURN
         !
         CALL iotk_scan_attr( attr, "NAME",  sname_(i), IERR=ierr )
         IF (ierr/=0) RETURN
         !
         CALL iotk_scan_dat( iunit, "ROTATION", s_(1:3,1:3,i), IERR=ierr )
         IF (ierr/=0) RETURN
         !
         CALL iotk_scan_end( iunit, "SYMM" // TRIM( iotk_index( i ) ), IERR=ierr )
         IF (ierr/=0) RETURN
         !
      END DO
      !
      CALL iotk_scan_end( iunit, "SYMMETRIES", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( present(nsym) )        nsym          = nsym_
      IF ( present(nrot) )        nrot          = nrot_
      IF ( present(invsym) )      invsym        = invsym_
      IF ( present(noinv) )       noinv         = noinv_
      IF ( present(time_reversal) ) time_reversal = time_reversal_
      IF ( present(no_t_rev) )    no_t_rev      = no_t_rev_
      IF ( present(nat) )         nat           = nat_
      IF ( present(trasl) )       trasl(1:3, 1:nsym_)   = trasl_(1:3, 1:nsym_)
      IF ( present(s) )           s(1:3, 1:3, 1:nrot_)  = s_(1:3, 1:3, 1:nrot_)
      IF ( present(irt) )         irt(1:nsym_, 1:nat_)  = irt_(1:nsym_, 1:nat_)
      IF ( present(sname) )  THEN
          DO i = 1, nrot_
                                  sname( i )            = trim( sname_( i ) )
          ENDDO
      ENDIF
      IF ( present(s_units) )     s_units               = trim( s_units_ )
      IF ( present(t_rev) )       t_rev( 1:nsym_ )      = t_rev_( 1:nsym_ )
      !
      DEALLOCATE( irt_ )
      !
    END SUBROUTINE qexml_read_symmetry
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_efield( tefield, dipfield, edir, emaxpos, eopreg, eamp, &
                                  monopole, zmon, relaxz, block, block_1, block_2,&
                                  block_height, found, ierr )
      !----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      LOGICAL,   OPTIONAL, INTENT(out) :: tefield, dipfield, monopole, relaxz, block
      INTEGER,   OPTIONAL, INTENT(out) :: edir
      REAL(DP),  OPTIONAL, INTENT(out) :: emaxpos, eopreg, eamp, zmon, block_1, block_2, block_height
      LOGICAL,             INTENT(out) :: found
      INTEGER,             INTENT(out) :: ierr
      !
      LOGICAL   :: tefield_, dipfield_, monopole_, block_, relaxz_
      INTEGER   :: edir_
      REAL(DP)  :: emaxpos_, eopreg_, eamp_, zmon_, block_1_, block_2_, block_height_
      !
      ierr = 0
      !
      CALL iotk_scan_begin( iunit, "ELECTRIC_FIELD", FOUND=found, IERR=ierr )
      IF ( ( .NOT. found ).OR.( ierr /= 0 ) ) RETURN
      !
      !
      CALL iotk_scan_dat( iunit, "HAS_ELECTRIC_FIELD", tefield_, IERR=ierr )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "HAS_DIPOLE_CORRECTION", dipfield_, IERR=ierr )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "FIELD_DIRECTION", edir_, IERR=ierr )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "MAXIMUM_POSITION", emaxpos_, IERR=ierr )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "INVERSE_REGION", eopreg_, IERR=ierr )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "FIELD_AMPLITUDE", eamp_, IERR=ierr )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "MONOPOLE_PLANE", monopole_ )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "MONOPOLE_POS", zmon_ )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "RELAX_Z", relaxz_ )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "BLOCK", block_ )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "BLOCK_1", block_1_ )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "BLOCK_2", block_2_ )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "BLOCK_HEIGHT", block_height_ )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_end( iunit, "ELECTRIC_FIELD", IERR=ierr )
      IF ( ierr /= 0 ) RETURN
      !
      !
      IF ( present(tefield) )        tefield      = tefield_
      IF ( present(dipfield) )       dipfield     = dipfield_
      IF ( present(edir) )           edir         = edir_
      IF ( present(emaxpos) )        emaxpos      = emaxpos_
      IF ( present(eopreg) )         eopreg       = eopreg_
      IF ( present(eamp) )           eamp         = eamp_
      IF ( present(monopole) )       monopole     = monopole_
      IF ( present(zmon) )           zmon         = zmon_
      IF ( present(relaxz) )         relaxz       = relaxz_
      IF ( present(block) )          block        = block_
      IF ( present(block_1) )        block_1      = block_1_
      IF ( present(block_2) )        block_2      = block_2_
      IF ( present(block_height) )   block_height = block_height_
      !
    END SUBROUTINE qexml_read_efield
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_exx( x_gamma_extrapolation, nqx1, nqx2, nqx3, &
                          exxdiv_treatment, yukawa, ecutvcut, exx_fraction, &
                          screening_parameter, gau_parameter, exx_is_active, ecutfock, &
                          found, ierr )
      !----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      LOGICAL,          OPTIONAL, INTENT(OUT) :: x_gamma_extrapolation, exx_is_active
      INTEGER,          OPTIONAL, INTENT(OUT) :: nqx1, nqx2, nqx3
      CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: exxdiv_treatment
      REAL(DP),         OPTIONAL, INTENT(OUT) :: yukawa, ecutvcut, exx_fraction
      REAL(DP),         OPTIONAL, INTENT(OUT) :: screening_parameter, ecutfock
      REAL(DP),         OPTIONAL, INTENT(OUT) :: gau_parameter
      LOGICAL,                    INTENT(out) :: found
      INTEGER,                    INTENT(out) :: ierr
      !
      LOGICAL  :: x_gamma_extrapolation_, exx_is_active_
      INTEGER  :: nqx1_, nqx2_, nqx3_
      REAL(DP) :: yukawa_, ecutvcut_, exx_fraction_
      REAL(DP) :: screening_parameter_, ecutfock_
      REAL(DP) :: gau_parameter_
      CHARACTER(LEN=80) :: exxdiv_treatment_
      !
      !
      ierr = 0
      !
      CALL iotk_scan_begin( iunit, "EXACT_EXCHANGE", FOUND=found, IERR=ierr )
      IF ( ( .NOT. found ).OR.( ierr /= 0 ) ) RETURN
      !
      call iotk_scan_dat(iunit, "x_gamma_extrapolation", x_gamma_extrapolation_, IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      call iotk_scan_dat(iunit, "nqx1", nqx1_, IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      call iotk_scan_dat(iunit, "nqx2", nqx2_, IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      call iotk_scan_dat(iunit, "nqx3", nqx3_, IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      call iotk_scan_dat(iunit, "exxdiv_treatment", exxdiv_treatment_, IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      call iotk_scan_dat(iunit, "yukawa", yukawa_, IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      call iotk_scan_dat(iunit, "ecutvcut", ecutvcut_, IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      call iotk_scan_dat(iunit, "exx_fraction", exx_fraction_, IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      call iotk_scan_dat(iunit, "screening_parameter", screening_parameter_, IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      ! Check if existing, for back-compatibility
      call iotk_scan_dat(iunit, "gau_parameter", gau_parameter_, FOUND=found, IERR=ierr)
      IF ( .NOT. found )  gau_parameter_=0.0_dp
      IF ( ierr /= 0 ) RETURN
      !
      call iotk_scan_dat(iunit, "exx_is_active", exx_is_active_, IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      ! Check if existing, for back-compatibility
      call iotk_scan_dat(iunit, "ecutfock", ecutfock_, IERR=ierr, FOUND=found)
      IF ( .NOT. found )  ecutfock_=-1.0_dp
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_end(iunit, "EXACT_EXCHANGE", IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      !
      IF ( present(x_gamma_extrapolation) ) x_gamma_extrapolation = x_gamma_extrapolation_
      IF ( present(nqx1) )                                   nqx1 = nqx1_
      IF ( present(nqx2) )                                   nqx2 = nqx2_
      IF ( present(nqx3) )                                   nqx3 = nqx3_
      IF ( present(exxdiv_treatment) )           exxdiv_treatment = exxdiv_treatment_
      IF ( present(yukawa) )                               yukawa = yukawa_
      IF ( present(ecutvcut) )                           ecutvcut = ecutvcut_
      IF ( present(exx_fraction) )                   exx_fraction = exx_fraction_
      IF ( present(screening_parameter) )     screening_parameter = screening_parameter_
      IF ( present(ecutfock) )                           ecutfock = ecutfock_
      IF ( present(gau_parameter) )                 gau_parameter = gau_parameter_
      IF ( present(exx_is_active) )                 exx_is_active = exx_is_active_
      !
      found = .TRUE.
      !
    END SUBROUTINE qexml_read_exx
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_esm( esm_nfit, esm_efield, esm_w, esm_a, esm_bc, ierr )
      !----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER,          OPTIONAL, INTENT(OUT) :: esm_nfit
      REAL(DP),         OPTIONAL, INTENT(OUT) :: esm_efield, esm_w, esm_a
      CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: esm_bc
      INTEGER,                    INTENT(out) :: ierr
      !
      INTEGER  :: esm_nfit_
      REAL(DP) :: esm_efield_, esm_w_, esm_a_
      CHARACTER(LEN=3) :: esm_bc_
      !
      !
      ierr = 0
      !
      CALL iotk_scan_begin( iunit, "ESM", IERR=ierr )
      IF ( ierr /= 0 ) RETURN
      !
      call iotk_scan_dat(iunit, "esm_nfit", esm_nfit_, IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      call iotk_scan_dat(iunit, "esm_efield", esm_efield_, IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      call iotk_scan_dat(iunit, "esm_w", esm_w_, IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      call iotk_scan_dat(iunit, "esm_a", esm_a_, IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      call iotk_scan_dat(iunit, "esm_bc", esm_bc_, IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_end(iunit, "ESM", IERR=ierr)
      IF ( ierr /= 0 ) RETURN
      !
      !
      IF ( present(esm_nfit) )    esm_nfit    = esm_nfit_
      IF ( present(esm_efield) )  esm_efield  = esm_efield_
      IF ( present(esm_w) )       esm_w       = esm_w_
      IF ( present(esm_a) )       esm_a       = esm_a_
      IF ( present(esm_bc) )      esm_bc      = esm_bc_
      !
    END SUBROUTINE qexml_read_esm
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_planewaves( ecutwfc, ecutrho, npwx, gamma_only, &
                                      nr1, nr2, nr3,  ngm,  nr1s, nr2s, nr3s, ngms, &
                                      nr1b, nr2b, nr3b,  igv, cutoff_units, ierr )
      !------------------------------------------------------------------------
      !
      !
      INTEGER,      OPTIONAL, INTENT(out) :: npwx, nr1, nr2, nr3, ngm, &
                                             nr1s, nr2s, nr3s, ngms, nr1b, nr2b, nr3b
      INTEGER,      OPTIONAL, INTENT(out) :: igv(:,:)
      REAL(DP),     OPTIONAL, INTENT(out) :: ecutwfc, ecutrho
      LOGICAL,      OPTIONAL, INTENT(out) :: gamma_only
      CHARACTER(*), OPTIONAL, INTENT(out) :: cutoff_units
      INTEGER,                INTENT(out) :: ierr
      !
      INTEGER        :: npwx_, nr1_, nr2_, nr3_, ngm_, &
                        nr1s_, nr2s_, nr3s_, ngms_, nr1b_, nr2b_, nr3b_
      REAL(DP)       :: ecutwfc_, ecutrho_
      CHARACTER(256) :: cutoff_units_
      LOGICAL        :: gamma_only_
      !
      ierr = 0
      !
      CALL iotk_scan_begin( iunit, "PLANE_WAVES", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_empty( iunit, "UNITS_FOR_CUTOFF", ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "UNITS", cutoff_units_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "WFC_CUTOFF", ecutwfc_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "RHO_CUTOFF", ecutrho_ , IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "MAX_NUMBER_OF_GK-VECTORS", npwx_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "GAMMA_ONLY", gamma_only_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_empty( iunit, "FFT_GRID", ATTR = attr, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_attr( attr, "nr1", nr1_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "nr2", nr2_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "nr3", nr3_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "GVECT_NUMBER", ngm_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_empty( iunit, "SMOOTH_FFT_GRID", ATTR = attr, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_attr( attr, "nr1s", nr1s_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "nr2s", nr2s_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "nr3s", nr3s_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "SMOOTH_GVECT_NUMBER", ngms_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( present( igv ) ) THEN
          !
          CALL iotk_scan_begin( iunit, "G-VECTORS", IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_dat( iunit, "g", igv(1:3,1:ngm_), IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_end( iunit, "G-VECTORS", IERR=ierr )
          IF (ierr/=0) RETURN
          !
      ENDIF
      !
      !
      CALL iotk_scan_empty( iunit, "SMALLBOX_FFT_GRID", ATTR = attr, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_attr( attr, "nr1b", nr1b_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "nr2b", nr2b_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "nr3b", nr3b_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_end( iunit, "PLANE_WAVES", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( present( ecutwfc ) )           ecutwfc      = ecutwfc_
      IF ( present( ecutrho ) )           ecutrho      = ecutrho_
      IF ( present( npwx ) )              npwx         = npwx_
      IF ( present( gamma_only ) )        gamma_only   = gamma_only_
      IF ( present( nr1 ) )               nr1          = nr1_
      IF ( present( nr2 ) )               nr2          = nr2_
      IF ( present( nr3 ) )               nr3          = nr3_
      IF ( present( ngm ) )               ngm          = ngm_
      IF ( present( nr1s ) )              nr1s         = nr1s_
      IF ( present( nr2s ) )              nr2s         = nr2s_
      IF ( present( nr3s ) )              nr3s         = nr3s_
      IF ( present( ngms ) )              ngms         = ngms_
      IF ( present( nr1b ) )              nr1b         = nr1b_
      IF ( present( nr2b ) )              nr2b         = nr2b_
      IF ( present( nr3b ) )              nr3b         = nr3b_
      IF ( present( cutoff_units ) )      cutoff_units = trim( cutoff_units_ )
      !
    END SUBROUTINE qexml_read_planewaves
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_gk( ik, npwk, npwkx, gamma_only, xk, k_units, index, igk, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,                INTENT(in)  :: ik
      INTEGER,      OPTIONAL, INTENT(out) :: npwk, npwkx
      LOGICAL,      OPTIONAL, INTENT(out) :: gamma_only
      REAL(DP),     OPTIONAL, INTENT(out) :: xk(3)
      CHARACTER(*), OPTIONAL, INTENT(out) :: k_units
      INTEGER,      OPTIONAL, INTENT(out) :: igk(:,:), index(:)
      INTEGER,                INTENT(out) :: ierr
      !
      CHARACTER(256) :: filename, k_units_
      INTEGER   :: npwk_, npwkx_
      LOGICAL   :: gamma_only_
      REAL(DP) :: xk_(3)
      INTEGER   :: iunaux
      !

      ierr = 0
      !
      CALL iotk_free_unit( iunaux )
      filename = qexml_wfc_filename( datadir_in, 'gkvectors', ik )
      !
      CALL iotk_open_read ( iunaux, FILE = trim(filename), IERR=ierr )
      IF (ierr/=0)  RETURN
      !
      CALL iotk_scan_dat( iunaux, 'NUMBER_OF_GK-VECTORS', npwk_, IERR=ierr)
      IF (ierr/=0)  RETURN
      !
      CALL iotk_scan_dat( iunaux, 'MAX_NUMBER_OF_GK-VECTORS', npwkx_, IERR=ierr)
      IF (ierr/=0)  RETURN
      !
      IF ( qexml_version_before_1_4_0 ) THEN
         !
         IF ( present( gamma_only ) ) THEN
             !
             CALL qexml_read_planewaves( GAMMA_ONLY=gamma_only_, IERR=ierr)
             IF (ierr/=0)  RETURN
             !
         ENDIF
         !
      ELSE
         !
         CALL iotk_scan_dat( iunaux, 'GAMMA_ONLY', gamma_only_, IERR=ierr)
         IF (ierr/=0)  RETURN
         !
      ENDIF
      !
      CALL iotk_scan_dat( iunaux, 'K-POINT_COORDS', xk_, ATTR=attr, IERR=ierr)
      IF (ierr/=0)  RETURN
      CALL iotk_scan_attr( attr, 'UNITS', k_units_, IERR=ierr)
      IF (ierr/=0)  RETURN
      !
      IF ( present( index ) ) THEN
          !
          CALL iotk_scan_dat( iunaux, 'INDEX', index(1:npwk_), IERR=ierr)
          IF (ierr/=0)  RETURN
          !
      ENDIF
      !
      IF ( present( igk ) ) THEN
          !
          CALL iotk_scan_dat( iunaux, 'GRID', igk(1:3, 1:npwk_), IERR=ierr)
          IF (ierr/=0)  RETURN
          !
      ENDIF
      !
      CALL iotk_close_read ( iunaux, IERR=ierr )
      IF (ierr/=0)  RETURN
      !
      !
      IF ( present( npwk ) )       npwk          = npwk_
      IF ( present( npwkx ) )      npwkx         = npwkx_
      IF ( present( gamma_only ) ) gamma_only    = gamma_only_
      IF ( present( xk ) )         xk(1:3)       = xk_(1:3)
      IF ( present( k_units ) )    k_units       = trim(k_units_)
      !
    END SUBROUTINE qexml_read_gk
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_spin( lsda, noncolin, npol, lspinorb, domag, ierr )
      !------------------------------------------------------------------------
      !
      LOGICAL, OPTIONAL, INTENT(out) :: lsda, noncolin, lspinorb, domag
      INTEGER, OPTIONAL, INTENT(out) :: npol
      INTEGER,           INTENT(out) :: ierr
      !
      LOGICAL   :: lsda_, noncolin_, lspinorb_, domag_,found
      INTEGER   :: npol_
      !

      ierr = 0
      !
      CALL iotk_scan_begin( iunit, "SPIN", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "LSDA", lsda_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "NON-COLINEAR_CALCULATION", noncolin_, IERR=ierr, FOUND=found )
      IF (ierr/=0) RETURN
      IF ( .not. found ) noncolin_ = .FALSE.
      !
      npol_ = 1
      !
      IF ( noncolin_ ) THEN
          !
          CALL iotk_scan_dat( iunit, "SPINOR_DIM", npol_, IERR=ierr )
          IF (ierr/=0) RETURN
          !
      ENDIF
      !
      CALL iotk_scan_dat( iunit, "SPIN-ORBIT_CALCULATION", lspinorb_, IERR=ierr, FOUND=found )
      IF (ierr/=0) RETURN
      IF ( .NOT. found ) lspinorb_ = .FALSE.
      !
      CALL iotk_scan_dat( iunit, "SPIN-ORBIT_DOMAG", domag_, IERR=ierr, FOUND=found )
      IF ( .NOT. found ) domag_ = .FALSE.
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_end( iunit, "SPIN", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( present( lsda ) )       lsda      = lsda_
      IF ( present( noncolin ) )   noncolin  = noncolin_
      IF ( present( npol ) )       npol      = npol_
      IF ( present( lspinorb ) )   lspinorb  = lspinorb_
      IF ( present( domag ) )      domag     = domag_
      !

    END SUBROUTINE qexml_read_spin
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_xc( dft, lda_plus_u, lda_plus_u_kind, U_projection, &
                              Hubbard_lmax, Hubbard_l, nsp, Hubbard_U, Hubbard_J,&
                              Hubbard_J0, Hubbard_alpha, Hubbard_beta, &
                              inlc, vdw_table_name, acfdt_in_pw, llondon, london_s6, &
                              london_rcut, lxdm, ts_vdw, vdw_isolated, ierr )
      !----------------------------------------------------------------------
      !
      CHARACTER(len=*), OPTIONAL, INTENT(out) :: dft
      LOGICAL,          OPTIONAL, INTENT(out) :: lda_plus_u
      INTEGER,          OPTIONAL, INTENT(out) :: lda_plus_u_kind
      !
      INTEGER,          OPTIONAL, INTENT(out) :: Hubbard_lmax
      INTEGER,          OPTIONAL, INTENT(out) :: Hubbard_l(:)
      INTEGER,          OPTIONAL, INTENT(out) :: nsp
      REAL(DP),         OPTIONAL, INTENT(out) :: Hubbard_U(:), Hubbard_J(:,:),&
                                                 Hubbard_alpha(:), &
                                                 Hubbard_J0(:), Hubbard_beta(:)
      INTEGER,          OPTIONAL, INTENT(out) :: inlc
      CHARACTER(LEN=*), OPTIONAL, INTENT(out) :: U_projection
      CHARACTER(LEN=*), OPTIONAL, INTENT(out) :: vdw_table_name
      LOGICAL,          OPTIONAL, INTENT(out) :: acfdt_in_pw
      LOGICAL,  OPTIONAL, INTENT(out) :: llondon, lxdm, ts_vdw, vdw_isolated
      REAL(DP), OPTIONAL, INTENT(out) :: london_s6, london_rcut
      !
      INTEGER,                    INTENT(out) :: ierr
      !
      CHARACTER(LEN=256)      :: dft_, vdw_table_name_, U_projection_
      LOGICAL                 :: lda_plus_u_, found
      LOGICAL                 :: acfdt_in_pw_
      INTEGER                 :: Hubbard_lmax_, nsp_,lda_plus_u_kind_, inlc_
      INTEGER,    ALLOCATABLE :: Hubbard_l_(:)
      REAL(DP),   ALLOCATABLE :: Hubbard_U_(:), Hubbard_J_(:,:)
      REAL(DP),   ALLOCATABLE :: Hubbard_alpha_(:), Hubbard_J0_(:), Hubbard_beta_(:)
      LOGICAL                 :: llondon_, lxdm_, ts_vdw_, vdw_isolated_
      REAL(DP)                :: london_s6_=0._dp, london_rcut_=0._dp
      !
      ierr = 0
      !
      !
      CALL iotk_scan_begin( iunit, "EXCHANGE_CORRELATION", IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "DFT", dft_, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "LDA_PLUS_U_CALCULATION", lda_plus_u_, FOUND=found, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      IF ( .NOT. found ) lda_plus_u_ = .FALSE.
      !
      IF ( lda_plus_u_ ) THEN
         !
         CALL iotk_scan_dat( iunit, "NUMBER_OF_SPECIES", nsp_, IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
         CALL iotk_scan_dat( iunit, "HUBBARD_LMAX", Hubbard_lmax_, IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
         ALLOCATE( Hubbard_l_(1:nsp_) )
         ALLOCATE( Hubbard_U_(nsp_) )
         !
         CALL iotk_scan_dat( iunit, "HUBBARD_L", Hubbard_l_, IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
         CALL iotk_scan_dat( iunit, "HUBBARD_U", Hubbard_U_, IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
         IF ( PRESENT (lda_plus_u_kind) ) THEN
            CALL iotk_scan_dat( iunit, "LDA_PLUS_U_KIND", lda_plus_u_kind_, FOUND = found, IERR=ierr )
            IF ( ierr/=0 ) RETURN
         END IF
         !
         IF ( PRESENT (U_projection) ) THEN
            CALL iotk_scan_dat( iunit, "U_PROJECTION_TYPE", U_projection_, FOUND = found, IERR=ierr )
            IF ( ierr/=0 ) RETURN
            IF ( .NOT. found ) U_projection_='atomic' ! for compatibility
         END IF
         !
         IF ( PRESENT (Hubbard_J) ) THEN
            ALLOCATE( Hubbard_J_(3,nsp_) )
            CALL iotk_scan_dat( iunit, "HUBBARD_J", Hubbard_J_, FOUND = found, IERR=ierr )
            IF ( ierr/=0 ) RETURN
         END IF
         !
         IF ( PRESENT (Hubbard_J0) ) THEN
            ALLOCATE( Hubbard_J0_(nsp_) )
            CALL iotk_scan_dat( iunit, "HUBBARD_J0", Hubbard_J0_, FOUND = found, IERR=ierr )
            IF ( ierr/=0 ) RETURN
         END IF
         !
         IF ( PRESENT (Hubbard_alpha) ) THEN
            ALLOCATE( Hubbard_alpha_(nsp_) )
            CALL iotk_scan_dat( iunit, "HUBBARD_ALPHA", Hubbard_alpha_, FOUND = found, IERR=ierr )
            IF ( ierr/=0 ) RETURN
         END IF
         !
         IF ( PRESENT (Hubbard_beta) ) THEN
            ALLOCATE( Hubbard_beta_(nsp_) )
            CALL iotk_scan_dat( iunit, "HUBBARD_BETA", Hubbard_beta_, FOUND = found, IERR=ierr )
            IF ( ierr/=0 ) RETURN
         END IF
         !
      ENDIF
      !
      CALL iotk_scan_dat( iunit, "NON_LOCAL_DF", inlc_, FOUND = found )
      IF ( found ) THEN
         !
         IF ( inlc_ > 0 ) CALL iotk_scan_dat( iunit, "VDW_KERNEL_NAME", vdw_table_name_ )
         !
      ELSE
         !
         inlc_ = 0
         vdw_table_name_ = ' '
         !
      ENDIF
      !
      CALL iotk_scan_dat( iunit, "ACFDT_IN_PW", acfdt_in_pw_, FOUND=found, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      IF ( .NOT. found ) acfdt_in_pw_ = .FALSE.
!      IF (acfdt_in_pw) dft_name = 'NOX NOC NOGX NOGC'
      !
      CALL iotk_scan_dat( iunit, "ACFDT_IN_PW", acfdt_in_pw_, FOUND=found, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      IF ( .NOT. found ) acfdt_in_pw_ = .FALSE.
!      IF (acfdt_in_pw) dft_name = 'NOX NOC NOGX NOGC'
      !
      CALL iotk_scan_begin( iunit, "DFT_D2", FOUND=found, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      llondon_ = found
      IF ( llondon_ ) THEN
         CALL iotk_scan_dat( iunit, "SCALING_FACTOR", london_s6_ )
         CALL iotk_scan_dat( iunit, "CUTOFF_RADIUS",  london_rcut_)
         CALL iotk_scan_end( iunit, "DFT_D2" )
      ENDIF
      !
      CALL iotk_scan_dat( iunit, "XDM", lxdm_, FOUND = found, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      IF (.NOT. found) lxdm_ = .FALSE.
      !
      CALL iotk_scan_begin( iunit,"TKATCHENKO-SCHEFFLER", FOUND=found, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      ts_vdw_ = found
      IF ( ts_vdw_ ) THEN
         CALL iotk_scan_dat( iunit, "ISOLATED_SYSTEM", vdw_isolated_ )
         CALL iotk_scan_end( iunit, "TKATCHENKO-SCHEFFLER" )
      END IF
      !
      CALL iotk_scan_end( iunit, "EXCHANGE_CORRELATION", IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      IF ( present( dft ) )           dft           = dft_
      IF ( present( lda_plus_u ) )    lda_plus_u    = lda_plus_u_
      !
      IF ( lda_plus_u_ )  THEN
         !
         IF ( present( nsp ) )             nsp                   = nsp_
         IF ( present( lda_plus_u_kind ) ) lda_plus_u_kind       = lda_plus_u_kind_
         IF ( present( U_projection ) )    U_projection          = U_projection_
         IF ( present( Hubbard_lmax ) )    Hubbard_lmax          = Hubbard_lmax_
         IF ( present( Hubbard_l ) )       Hubbard_l(1:nsp_)     = Hubbard_l_(:)
         IF ( present( Hubbard_U ) )       Hubbard_U(1:nsp_)     = Hubbard_U_(1:nsp_)
         IF ( present( Hubbard_J ) )       Hubbard_J(1:3,1:nsp_) = Hubbard_J_(1:3,1:nsp_)
         IF ( present( Hubbard_J0 ) )      Hubbard_J0(1:nsp_)    = Hubbard_J0_(1:nsp_)
         IF ( present( Hubbard_alpha ) )   Hubbard_alpha(1:nsp_) = Hubbard_alpha_(1:nsp_)
         IF ( present( Hubbard_beta )  )   Hubbard_beta(1:nsp_)  = Hubbard_beta_(1:nsp_)
         !
         DEALLOCATE( Hubbard_l_ )
         DEALLOCATE( Hubbard_U_ )
         DEALLOCATE( Hubbard_J_ )
         DEALLOCATE( Hubbard_J0_ )
         DEALLOCATE( Hubbard_alpha_ )
         DEALLOCATE( Hubbard_beta_ )
         !
      ENDIF
      !
      IF (present(inlc) ) inlc = inlc_
      IF (present( vdw_table_name) )  vdw_table_name =  vdw_table_name_
      !
      IF ( present( acfdt_in_pw ) )    acfdt_in_pw    = acfdt_in_pw_
      IF (present(llondon) ) THEN
         llondon = llondon_
         IF (present(london_s6) )   london_s6   = london_s6_
         IF (present(london_rcut) ) london_rcut = london_rcut_
      ELSE IF (present(lxdm) ) THEN
         lxdm = lxdm_
      ELSE IF (present(ts_vdw) ) THEN
         ts_vdw = ts_vdw_
         IF (present(vdw_isolated) ) vdw_isolated=vdw_isolated_
      END IF
      !
    END SUBROUTINE qexml_read_xc
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_occ( lgauss, ngauss, degauss, degauss_units, &
                               ltetra, tetra_type, ntetra, tetra, tfixed_occ, &
                               nstates_up, nstates_dw, input_occ, ierr )
      !------------------------------------------------------------------------
      !
      LOGICAL,      OPTIONAL, INTENT(out) :: lgauss, ltetra, tfixed_occ
      INTEGER,      OPTIONAL, INTENT(out) :: ngauss, ntetra, tetra_type
      INTEGER,      OPTIONAL, INTENT(out) :: tetra(:,:)
      INTEGER,      OPTIONAL, INTENT(out) :: nstates_up, nstates_dw
      REAL(DP),     OPTIONAL, INTENT(out) :: degauss, input_occ(:,:)
      CHARACTER(*), OPTIONAL, INTENT(out) :: degauss_units
      INTEGER,                INTENT(out) :: ierr
      !
      LOGICAL        :: lgauss_, ltetra_, tfixed_occ_
      INTEGER        :: ngauss_, ntetra_, nstates_up_, nstates_dw_, tetra_type_
      LOGICAL        :: lsda_
      REAL(DP)      :: degauss_
      CHARACTER(256) :: degauss_units_
      INTEGER,  ALLOCATABLE :: tetra_(:,:)
      INTEGER :: i
      LOGICAL :: lfound,found
      !
      ierr = 0
      !
      CALL iotk_scan_begin( iunit, "OCCUPATIONS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "SMEARING_METHOD", lgauss_, FOUND=found, IERR=ierr )
      IF (ierr/=0) RETURN
      IF ( .NOT. found ) lgauss_ = .FALSE.
      !
      !
      ngauss_=0
      degauss_=-1.0d0 
      degauss_units_="none"
      !
      IF ( lgauss_ ) THEN
         !
         CALL iotk_scan_dat( iunit, "SMEARING_TYPE", ngauss_, IERR=ierr )
         IF (ierr/=0) RETURN
         !
         CALL iotk_scan_dat( iunit, "SMEARING_PARAMETER", degauss_ , &
                                     ATTR=attr, IERR=ierr )
         IF (ierr/=0) RETURN
         !
         CALL iotk_scan_attr( ATTR, "UNITS", degauss_units_, IERR=ierr )
         IF (ierr/=0) RETURN
         !
      ENDIF
      !
      CALL iotk_scan_dat( iunit, "TETRAHEDRON_METHOD", ltetra_, FOUND=found, IERR=ierr )
      IF (ierr/=0) RETURN
      IF ( .NOT. found ) ltetra_ = .FALSE.
      !
      !
      ntetra_ = 0
      IF ( ltetra_ ) THEN
         !
         CALL iotk_scan_dat( iunit, "NUMBER_OF_TETRAHEDRA", ntetra_, IERR=ierr )
         IF (ierr/=0) RETURN
         !
         CALL iotk_scan_dat( iunit, "TETRAHEDRON_TYPE", tetra_type_, IERR=ierr )
         !
         IF(tetra_type_ == 0) then
            !
            ALLOCATE( tetra_(4, ntetra_) )
            !
            DO i = 1, ntetra_
               !
               CALL iotk_scan_dat( iunit, "TETRAHEDRON"//iotk_index(i), &
                                               tetra_(1:4,i), IERR=ierr )
               IF (ierr/=0) RETURN
               !
            ENDDO
            !
         END IF
         !
      ENDIF
      !
      CALL iotk_scan_dat( iunit, "FIXED_OCCUPATIONS", tfixed_occ_, FOUND=found, IERR=ierr )
      IF (ierr/=0) RETURN
      IF ( .NOT. found ) tfixed_occ_ = .FALSE.
      !
      nstates_up_=0.0d0
      nstates_dw_=0.0d0
      !
      IF ( tfixed_occ_  .and. ( present(input_occ)   .or. &
                                present(nstates_up)  .or. present(nstates_dw) ) ) THEN
         !
         CALL iotk_scan_empty( iunit, "INFO", ATTR=attr, IERR=ierr)
         IF (ierr /=0 ) RETURN
         !
         CALL iotk_scan_attr( attr, "lsda", lsda_, IERR=ierr )
         IF (ierr /=0 ) RETURN
         !
         IF ( qexml_version_before_1_4_0 ) THEN
            !
            CALL iotk_scan_attr( attr, "nelup", nstates_up_, IERR=ierr )
            IF (ierr /=0 ) RETURN
            CALL iotk_scan_attr( attr, "neldw", nstates_dw_, IERR=ierr )
            IF (ierr /=0 ) RETURN
            !
         ELSE
            !
            ! current version
            !
            CALL iotk_scan_attr( attr, "nstates_up", nstates_up_, IERR=ierr )
            IF (ierr /=0 ) RETURN
            CALL iotk_scan_attr( attr, "nstates_down", nstates_dw_, IERR=ierr )
            IF (ierr /=0 ) RETURN
            !
         ENDIF
         !
         IF ( present( input_occ ) ) THEN
            !
            !
            CALL iotk_scan_dat( iunit, "INPUT_OCC_UP", input_occ(1:nstates_up_,1), IERR=ierr )
            IF (ierr/=0) RETURN
            !
            IF ( lsda_  ) THEN
               !
               CALL iotk_scan_dat( iunit, "INPUT_OCC_DOWN", input_occ(1:nstates_dw_,2), IERR=ierr )
               IF (ierr/=0) RETURN
               !
            ENDIF
            !
         ENDIF
         !
      ENDIF
      !
      CALL iotk_scan_end( iunit, "OCCUPATIONS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( present( lgauss ))           lgauss      = lgauss_
      IF ( present( ltetra ))           ltetra      = ltetra_
      IF ( present( tfixed_occ ))       tfixed_occ  = tfixed_occ_
      IF ( present( ngauss ))           ngauss      = ngauss_
      IF ( present( ntetra ))           ntetra      = ntetra_
      IF ( present( tetra_type ))       tetra_type  = tetra_type_
      IF ( present( degauss ))          degauss     = degauss_
      IF ( present( degauss_units ))    degauss_units  = trim(degauss_units_)
      IF ( present( nstates_up ))       nstates_up  = nstates_up_
      IF ( present( nstates_dw ))       nstates_dw  = nstates_dw_
      !
      IF ( ltetra_ .and. (tetra_type_ == 0)) THEN
         !
         IF ( present( tetra ) )         tetra(1:4, 1:ntetra_)  = tetra_
         !
         DEALLOCATE( tetra_ )
         !
      ENDIF

    END SUBROUTINE qexml_read_occ
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_bz( num_k_points, xk, wk, k1, k2, k3, nk1, nk2, nk3, &
                              nks_start, xk_start, wk_start, qnorm, &
                              k_units, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,       OPTIONAL, INTENT(out) :: num_k_points, k1, k2, k3, nk1, nk2, nk3, &
                                              nks_start
      REAL(DP),     OPTIONAL, INTENT(out) :: xk(:,:), wk(:), qnorm
      REAL(DP),     OPTIONAL, ALLOCATABLE, INTENT(out) :: xk_start(:,:), wk_start(:)
      CHARACTER(*),  OPTIONAL, INTENT(out) :: k_units
      INTEGER,                 INTENT(out) :: ierr
      !
      INTEGER                :: num_k_points_, k1_, k2_, k3_, nk1_, nk2_, nk3_,nks_start_
      CHARACTER(256)         :: k_units_
      REAL(DP)               :: qnorm_
      REAL(DP), ALLOCATABLE  :: xk_(:,:), wk_(:)
      REAL(DP), ALLOCATABLE  :: xk_start_(:,:), wk_start_(:)
      !
      INTEGER :: ik, i
      LOGICAL :: found
      !
      ierr = 0
      !
      CALL iotk_scan_begin( iunit, "BRILLOUIN_ZONE", IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "NUMBER_OF_K-POINTS", num_k_points_, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      !
      CALL iotk_scan_empty( iunit, "UNITS_FOR_K-POINTS", ATTR=attr, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      CALL iotk_scan_attr( attr, "UNITS", k_units_, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      CALL iotk_scan_empty( iunit, "MONKHORST_PACK_GRID", ATTR=attr, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      CALL iotk_scan_attr( attr, "nk1", nk1_, IERR=ierr  )
      IF ( ierr/=0 ) RETURN
      CALL iotk_scan_attr( attr, "nk2", nk2_, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      CALL iotk_scan_attr( attr, "nk3", nk3_, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      !
      CALL iotk_scan_empty( iunit, "MONKHORST_PACK_OFFSET", ATTR=attr, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      CALL iotk_scan_attr( attr, "k1", k1_, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      CALL iotk_scan_attr( attr, "k2", k2_, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      CALL iotk_scan_attr( attr, "k3", k3_, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      !
      ALLOCATE( xk_( 3, num_k_points_ ) )
      ALLOCATE( wk_(    num_k_points_ ) )
      !
      DO ik = 1, num_k_points_
         !
         CALL iotk_scan_empty( iunit, "K-POINT" // trim( iotk_index(ik) ), &
                               ATTR=attr, IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
         CALL iotk_scan_attr( attr, "XYZ", xk_(:,ik), IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
         CALL iotk_scan_attr( attr, "WEIGHT", wk_(ik), IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
      ENDDO
      !
      nks_start_=0
      !
      IF ( present(nks_start) .or. present(xk_start) .or. present(wk_start) ) THEN
         !
         CALL iotk_scan_dat( iunit, "STARTING_K-POINTS", nks_start_, &
              FOUND = found )
         IF (.NOT. found) nks_start_=0
         !
         IF (nks_start_ > 0 ) THEN
            !
            ALLOCATE( xk_start_(3,nks_start_) )
            ALLOCATE( wk_start_(nks_start_) )
            !
         END IF
         !
         DO ik = 1, nks_start_
            !
            CALL iotk_scan_empty( iunit, "K-POINT_START" // &
                 & TRIM( iotk_index( ik ) ), ATTR=attr )
            !
            CALL iotk_scan_attr( attr, "XYZ", xk_start_(:,ik) )
            !
            CALL iotk_scan_attr( attr, "WEIGHT", wk_start_(ik) )
            !
         END DO
         !
      END IF
      !   
      CALL iotk_scan_dat( iunit, "NORM-OF-Q", qnorm_, FOUND = found )
      IF (.not. found) qnorm_=0.0_DP
      !
      !
      CALL iotk_scan_end( iunit, "BRILLOUIN_ZONE", IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      !
      IF ( present( num_k_points ) )       num_k_points  = num_k_points_
      IF ( present( nk1 ) )                nk1           = nk1_
      IF ( present( nk2 ) )                nk2           = nk2_
      IF ( present( nk3 ) )                nk3           = nk3_
      IF ( present( k1 ) )                 k1            =  k1_
      IF ( present( k2 ) )                 k2            =  k2_
      IF ( present( k3 ) )                 k3            =  k3_
      IF ( present( k_units ) )            k_units       =  trim(k_units_)
      IF ( present( xk ) )                 xk(1:3,1:num_k_points_) = xk_(:,:)
      IF ( present( wk ) )                 wk(1:num_k_points_)     = wk_(:)
      IF ( present( nks_start ) )          nks_start     =  nks_start_
      !
      IF ( nks_start_>0 .AND. present( xk_start ) ) THEN
         IF (.NOT.ALLOCATED(xk_start)) ALLOCATE(xk_start(3,nks_start_))
         xk_start(1:3,1:nks_start_) =  xk_start_(:,:)
      ENDIF
      IF ( nks_start_>0 .AND. present( wk_start ) ) THEN
         IF (.NOT.ALLOCATED(wk_start)) ALLOCATE(wk_start(nks_start_))
         wk_start(1:nks_start_) =  wk_start_(:)
      ENDIF
      !
      IF ( present( qnorm ) )              qnorm         = qnorm_
      !
      DEALLOCATE( xk_ )
      DEALLOCATE( wk_ )
      IF (ALLOCATED(xk_start_)) DEALLOCATE(xk_start_)
      IF (ALLOCATED(wk_start_)) DEALLOCATE(wk_start_)
      !
    END SUBROUTINE qexml_read_bz
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_para( kunit, nproc, nproc_pool, nproc_image, &
                    ntask_groups, nproc_bgrp, nproc_ortho, found, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER, OPTIONAL, INTENT(OUT) :: kunit, nproc, nproc_pool, nproc_image, &
           ntask_groups, nproc_bgrp, nproc_ortho
      LOGICAL, INTENT(OUT) :: found
      INTEGER, INTENT(OUT) :: ierr
      !
      INTEGER :: kunit_, nproc_, nproc_pool_, nproc_image_, ntask_groups_, &
                 nproc_bgrp_, nproc_ortho_
      !
      LOGICAL :: found2
      !
      !
      CALL iotk_scan_begin( iunit, "PARALLELISM", FOUND=found,IERR=ierr )
      IF ((.NOT. found ) .OR. (ierr /= 0 ) ) RETURN
      !
      CALL iotk_scan_dat( iunit, &
                              "GRANULARITY_OF_K-POINTS_DISTRIBUTION", kunit_ )
      !
      CALL iotk_scan_dat( iunit, "NUMBER_OF_PROCESSORS", nproc_, FOUND=found2 )
      IF ( .NOT. found2) nproc_=1  !compatibility
      !
      CALL iotk_scan_dat( iunit, &
                              "NUMBER_OF_PROCESSORS_PER_POOL", nproc_pool_, FOUND=found2 )
      IF ( .NOT. found2) nproc_pool_=1 ! compatibility
      !
      CALL iotk_scan_dat( iunit, &
                              "NUMBER_OF_PROCESSORS_PER_IMAGE", nproc_image_, FOUND=found2 )
      IF ( .NOT. found2) nproc_image_=1 ! compatibility
      !
      CALL iotk_scan_dat( iunit, "NUMBER_OF_PROCESSORS_PER_TASKGROUP", &
                                              ntask_groups_, FOUND=found2 )
      IF ( .NOT. found2) ntask_groups_=1 ! compatibility
      !
      CALL iotk_scan_dat( iunit, "NUMBER_OF_PROCESSORS_PER_BAND_GROUP", &
                                              nproc_bgrp_, FOUND=found2 )
      IF ( .NOT. found2) nproc_bgrp_=1 ! compatibility
      !
      CALL iotk_scan_dat( iunit, "NUMBER_OF_PROCESSORS_PER_DIAGONALIZATION", &
                                              nproc_ortho_, FOUND=found2 )
      IF ( .NOT. found2) nproc_ortho_=1 ! compatibility
      !
      CALL iotk_scan_end( iunit, "PARALLELISM" )
      !
      !
      IF (present(kunit)) kunit = kunit_
      IF (present(nproc)) nproc = nproc_
      IF (present(nproc_pool)) nproc_pool = nproc_pool_
      IF (present(nproc_image)) nproc_image = nproc_image_
      IF (present(ntask_groups)) ntask_groups = ntask_groups_
      IF (present(nproc_bgrp)) nproc_bgrp = nproc_bgrp_
      IF (present(nproc_ortho)) nproc_ortho = nproc_ortho_
      !
    END SUBROUTINE qexml_read_para
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_phonon( modenum, xqq, q_units, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,       OPTIONAL, INTENT(out) :: modenum
      REAL(DP),     OPTIONAL, INTENT(out) :: xqq(:)
      CHARACTER(*),  OPTIONAL, INTENT(out) :: q_units
      INTEGER,                 INTENT(out) :: ierr
      !
      INTEGER         :: modenum_
      CHARACTER(256)  :: q_units_
      !

      ierr = 0
      !
      CALL iotk_scan_begin( iunit, "PHONON", IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "NUMBER_OF_MODES", modenum_, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      CALL iotk_scan_empty( iunit, "UNITS_FOR_Q-POINT", attr, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      CALL iotk_scan_attr( attr, "UNITS", q_units_, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      IF ( present (xqq) ) THEN
         !
         CALL iotk_scan_dat( iunit, "Q-POINT", xqq(:), IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
      ENDIF
      !
      CALL iotk_scan_end( iunit, "PHONON", IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      !
      IF ( present (modenum) )      modenum = modenum_
      IF ( present (q_units) )      q_units = trim(q_units_)
      !
    END SUBROUTINE qexml_read_phonon
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_bands_info( num_k_points, natomwfc, &
                                      nbnd, nbnd_up, nbnd_down, &
                                      nspin, nelec, nel_up, nel_down, &
                                      ef, two_fermi_energies, &
                                      ef_up, ef_dw ,energy_units, k_units, &
                                      noncolin, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,      OPTIONAL, INTENT(out) :: nbnd, nbnd_up, nbnd_down, num_k_points, nspin, natomwfc, nel_up,nel_down
      LOGICAL,      OPTIONAL, INTENT(out) :: noncolin, two_fermi_energies
      REAL(DP),     OPTIONAL, INTENT(out) :: ef, nelec, ef_up, ef_dw
      CHARACTER(*), OPTIONAL, INTENT(out) :: energy_units, k_units
      INTEGER,                INTENT(out) :: ierr
      !
      INTEGER        :: nbnd_, nbnd_up_, nbnd_down_, num_k_points_, nspin_, natomwfc_, nel_up_, nel_down_
      LOGICAL        :: noncolin_, two_fermi_energies_
      REAL(DP)       :: ef_, nelec_, ef_up_, ef_dw_
      CHARACTER(256) :: energy_units_, k_units_
      !
      LOGICAL :: found
      !
      ierr = 0
      !
      !
      CALL iotk_scan_begin( iunit, "BAND_STRUCTURE_INFO", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      CALL iotk_scan_dat  ( iunit, "NUMBER_OF_K-POINTS", num_k_points_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat  ( iunit, "NUMBER_OF_SPIN_COMPONENTS", nspin_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat  ( iunit, "NON-COLINEAR_CALCULATION", noncolin_, FOUND=found, IERR=ierr )
      IF (ierr/=0) RETURN
      IF (.NOT. found) noncolin_ = .FALSE.
      !
      CALL iotk_scan_dat  ( iunit, "NUMBER_OF_ATOMIC_WFC", natomwfc_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      nbnd_up_   = 0
      nbnd_down_ = 0
      nel_up_    = 0.0d0
      nel_down_  = 0.0d0
      !
      IF ( nspin_ == 2 ) THEN
         !
         CALL iotk_scan_dat( iunit, &
                              "NUMBER_OF_BANDS", nbnd_, ATTR = attr, IERR = ierr )
         IF (ierr/=0) RETURN
         !
         CALL iotk_scan_attr( attr, "UP", nbnd_up_, IERR = ierr )
         IF (ierr/=0) RETURN
         !
         CALL iotk_scan_attr( attr, "DW", nbnd_down_, IERR = ierr )
         IF (ierr/=0) RETURN
         !
         CALL iotk_scan_dat( iunit, &
                              "NUMBER_OF_ELECTRONS", nelec_, ATTR = attr, IERR = ierr )
         IF (ierr/=0) RETURN
         !
         CALL iotk_scan_attr( attr, "UP", nel_up_, IERR = ierr )
         IF (ierr/=0) RETURN
         !
         CALL iotk_scan_attr( attr, "DW", nel_down_, IERR = ierr )
         IF (ierr/=0) RETURN
         !
      ELSE
         !
         CALL iotk_scan_dat( iunit, "NUMBER_OF_BANDS", nbnd_, IERR = ierr )
         IF (ierr/=0) RETURN
         !
         CALL iotk_scan_dat( iunit, "NUMBER_OF_ELECTRONS", nelec_, IERR = ierr )
         IF (ierr/=0) RETURN
         !
      END IF
      !
      CALL iotk_scan_empty( iunit, "UNITS_FOR_K-POINTS", ATTR = attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr ( attr,   "UNITS", k_units_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_empty( iunit, "UNITS_FOR_ENERGIES", ATTR = attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr ( attr,   "UNITS", energy_units_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "TWO_FERMI_ENERGIES", two_fermi_energies_, FOUND = found)
      IF ( .not. found ) two_fermi_energies_=.FALSE.
      !
      ef_up_  =0.0d0
      ef_dw_  =0.0d0
      ef_     =0.0d0
      !
      IF ( two_fermi_energies_ ) THEN
         !
         CALL iotk_scan_dat( iunit, "FERMI_ENERGY_UP", ef_up_, IERR=ierr )
         IF (ierr/=0) RETURN
         CALL iotk_scan_dat( iunit, "FERMI_ENERGY_DOWN", ef_dw_, IERR=ierr )
         IF (ierr/=0) RETURN
         CALL iotk_scan_dat( iunit, "ELECTRONS_UP", nel_up_, IERR=ierr )
         IF (ierr/=0) RETURN
         CALL iotk_scan_dat( iunit, "ELECTRONS_DOWN", nel_down_, IERR=ierr )
         IF (ierr/=0) RETURN
         !
      ELSE
         !
         CALL iotk_scan_dat  ( iunit, "FERMI_ENERGY", ef_ , FOUND=found )
         IF (ierr/=0) RETURN
         !
         !
      END IF
      !
      CALL iotk_scan_end( iunit, "BAND_STRUCTURE_INFO", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( present( nbnd ) )             nbnd           = nbnd_
      IF ( present( nbnd_up ) )          nbnd_up        = nbnd_up_
      IF ( present( nbnd_down ) )        nbnd_down      = nbnd_down_
      IF ( present( num_k_points ) )     num_k_points   = num_k_points_
      IF ( present( nspin ) )            nspin          = nspin_
      IF ( present( noncolin ) )         noncolin       = noncolin_
      IF ( present( natomwfc ) )         natomwfc       = natomwfc_
      IF ( present( nelec ) )            nelec          = nelec_
      IF ( present( nel_up ) )           nel_up         = nel_up_
      IF ( present( nel_down ) )         nel_down       = nel_down_
      IF ( present( ef ) )               ef             = ef_
      IF ( present( two_fermi_energies ) ) two_fermi_energies = two_fermi_energies_
      IF ( present( ef_up ) )            ef_up          = ef_up_
      IF ( present( ef_dw ) )            ef_dw          = ef_dw_
      IF ( present( energy_units ) )     energy_units   = trim( energy_units_ )
      IF ( present( k_units ) )          k_units        = trim( k_units_ )
      !
      !
    END SUBROUTINE qexml_read_bands_info
    !
    !
    !--------------------------------------------------------------------------
    SUBROUTINE qexml_read_bands_pw( num_k_points, nbnd, nkstot, lsda, lkpoint_dir, &
                                    filename, isk, et, wg , ierr )
      !------------------------------------------------------------------------ 
      !
      INTEGER,   INTENT(in) :: num_k_points, nbnd, nkstot
      LOGICAL,   INTENT(in) :: lsda, lkpoint_dir
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER,   INTENT(out), OPTIONAL :: isk(:)
      REAL(DP),  INTENT(out), OPTIONAL :: et(:,:), wg(:,:)
      INTEGER,   INTENT(out):: ierr
      !
      INTEGER  :: ik, ik_eff, iunaux
      INTEGER  :: isk_(nkstot)
      REAL(DP) :: et_(nbnd, nkstot),wg_(nbnd, nkstot)
      LOGICAL  :: found
      !
      !
      IF ( .NOT. lkpoint_dir) THEN
         !
         CALL iotk_free_unit( iunaux )
         !
         CALL iotk_open_read ( iunaux, FILE = trim(filename), IERR=ierr )
         IF (ierr/=0)  RETURN
         !
      END IF
      !
      !
      CALL iotk_scan_begin( iunit, "EIGENVALUES", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      k_points_loop: DO ik = 1, num_k_points
         !
         CALL iotk_scan_begin( iunit, "K-POINT" // TRIM( iotk_index( ik ) ) )
         !
         IF ( lsda ) THEN
            !
            isk_(ik) = 1
            !
            IF (lkpoint_dir) THEN
               CALL iotk_scan_begin(iunit, "DATAFILE"//TRIM(iotk_index(1)) , FOUND = found)
               IF (.NOT. found ) GO TO 10 ! workaround: PW-CP compatibility
               CALL iotk_scan_dat  ( iunit, "EIGENVALUES", et_(:,ik)  )
               CALL iotk_scan_dat  ( iunit, "OCCUPATIONS", wg_(:,ik) )
               CALL iotk_scan_end  ( iunit, "DATAFILE"//TRIM(iotk_index(1)) )
            ELSE
               CALL iotk_scan_begin( iunaux, &
                    "DATA_EIG"//TRIM( iotk_index(ik) )//"_SPIN_UP", FOUND=found )
               IF (.NOT. found ) GO TO 10 ! workaround: PW-CP compatibility
               CALL iotk_scan_dat  ( iunaux, "EIGENVALUES", et_(:,ik)  )
               CALL iotk_scan_dat  ( iunaux, "OCCUPATIONS", wg_(:,ik) )
               CALL iotk_scan_end( iunaux, &
                    "DATA_EIG"//TRIM( iotk_index( ik ) )//"_SPIN_UP")
            ENDIF
            !
10          CONTINUE
            !
            ik_eff = ik + num_k_points
            isk_(ik_eff) = 2
            !
            IF (lkpoint_dir) THEN
               CALL iotk_scan_begin(iunit,"DATAFILE"//TRIM(iotk_index(2)) &
                    , FOUND = found)
               IF (.NOT. found ) GO TO 20 ! workaround: PW-CP compatibility
               CALL iotk_scan_dat  ( iunit, "EIGENVALUES", et_(:,ik_eff) )
               CALL iotk_scan_dat  ( iunit, "OCCUPATIONS", wg_(:,ik_eff) )
               CALL iotk_scan_end( iunit, "DATAFILE"//TRIM(iotk_index(2)) )
            ELSE
               CALL iotk_scan_begin( iunaux, &
               "DATA_EIG"//TRIM( iotk_index(ik) )//"_SPIN_DW", FOUND=found )
               IF (.NOT. found ) GO TO 20 ! workaround: PW-CP compatibility
               CALL iotk_scan_dat  ( iunaux, "EIGENVALUES", et_(:,ik_eff) )
               CALL iotk_scan_dat  ( iunaux, "OCCUPATIONS", wg_(:,ik_eff) )
               CALL iotk_scan_end( iunaux, &
                    "DATA_EIG"//TRIM( iotk_index( ik ) )//"_SPIN_DW")
            ENDIF
            !
20          CONTINUE
            !
         ELSE
            !
            isk_(ik) = 1
            !
            IF (lkpoint_dir) THEN
               CALL iotk_scan_begin( iunit, "DATAFILE" , FOUND = found)
               IF (.NOT. found ) GO TO 15 ! workaround: PW-CP compatibility
               CALL iotk_scan_dat  ( iunit, "EIGENVALUES", et_(:,ik) )
               CALL iotk_scan_dat  ( iunit, "OCCUPATIONS", wg_(:,ik) )
               CALL iotk_scan_end  ( iunit, "DATAFILE" )
            ELSE
               CALL iotk_scan_begin( iunaux, &
                    "DATA_EIG"//TRIM( iotk_index(ik) ), FOUND = found )
               IF (.NOT. found ) GO TO 15 ! workaround: PW-CP compatibility
               CALL iotk_scan_dat  ( iunaux, "EIGENVALUES", et_(:,ik) )
               CALL iotk_scan_dat  ( iunaux, "OCCUPATIONS", wg_(:,ik) )
               CALL iotk_scan_end( iunaux, &
                    "DATA_EIG"//TRIM( iotk_index( ik ) ))
            ENDIF
15          CONTINUE
            !
         END IF
         !
         CALL iotk_scan_end( iunit, "K-POINT" // TRIM( iotk_index( ik ) ) )
         !
      END DO k_points_loop
      !
      CALL iotk_scan_end( iunit, "EIGENVALUES", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      IF (.NOT.lkpoint_dir) THEN
         CALL iotk_close_read ( iunaux, IERR=ierr )
         IF (ierr/=0)  RETURN
      END IF
      !
      IF ( present( isk ) )              isk( 1:nkstot )         = isk_(:)
      IF ( present( et ) )                et( 1:nbnd, 1:nkstot ) = et_(:,:)
      IF ( present( wg ) )                wg( 1:nbnd, 1:nkstot ) = wg_(:,:)
      !
    END SUBROUTINE qexml_read_bands_pw
    !
    !
    !-----------------------------------------------------------------------------
    SUBROUTINE qexml_read_bands_cp( num_k_points, nbnd_tot, nudx , nspin, iupdwn, &
      nupdwn, occ0, occm, ierr )
    !-----------------------------------------------------------------------------
      !
      INTEGER,  INTENT(OUT) :: ierr
      REAL(DP), INTENT(OUT) :: occ0(:)
      REAL(DP), INTENT(OUT) :: occm(:)
      !
      INTEGER, INTENT(in)  :: num_k_points, nspin, nbnd_tot, nudx
      INTEGER,               INTENT(IN) :: iupdwn(:)
      INTEGER,               INTENT(IN) :: nupdwn(:)
      !
      INTEGER :: ik, iss, ik_eff
      CHARACTER(LEN=4)      :: cspin
      REAL(DP), ALLOCATABLE :: occ_(:)
      REAL(DP)              :: wk_
      LOGICAL :: found
      !
#if !defined __HDF5
      CALL iotk_scan_begin( iunit, "EIGENVALUES", IERR=ierr )
#endif
      IF (ierr /= 0) RETURN
      !
      !
      k_points_loop1: DO ik = 1, num_k_points
         !
#if !defined __HDF5
         CALL iotk_scan_begin( iunit, "K-POINT" // TRIM( iotk_index(ik) ) )
#endif
         !
         CALL iotk_scan_dat( iunit, "WEIGHT", wk_ )
         !
         !
         DO iss = 1, nspin
            !
            cspin = iotk_index( iss )
            !
            ik_eff = ik + ( iss - 1 ) * num_k_points
            !
            ALLOCATE( occ_ ( MAX( nudx , nbnd_tot ) ) )
            !
            occ_ = 0.0d0
            !
            CALL iotk_scan_dat( iunit, "OCC0" // TRIM( cspin ), occ_ ( 1 :  nupdwn( iss ) ), FOUND = found )
            !
            IF( .NOT. found ) THEN
               !
               IF( nspin == 1 ) THEN
                  CALL iotk_scan_begin( iunit, "DATAFILE", FOUND = found )
               ELSE
                  CALL iotk_scan_begin( iunit, "DATAFILE"//TRIM(cspin), FOUND = found )
               END IF
               !
               CALL iotk_scan_dat  ( iunit, "OCCUPATIONS", occ_( 1:nbnd_tot ) )
               IF( nspin == 1 ) THEN
                  CALL iotk_scan_end( iunit, "DATAFILE" )
               ELSE
                  CALL iotk_scan_end( iunit, "DATAFILE"//TRIM(cspin) )
               END IF
               !
               IF( found ) THEN
                  occ0( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) = occ_ ( 1:nupdwn( iss ) ) * wk_
                  occm( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) = occ_ ( 1:nupdwn( iss ) ) * wk_
               END IF
               !
            ELSE
               !
               occ0( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) = occ_ ( 1:nupdwn( iss ) )
               !
               CALL iotk_scan_dat( iunit, "OCCM" // TRIM( cspin ), occ_ ( 1 :  nupdwn( iss ) ), FOUND = found )
               !
               IF( found ) THEN
                  occm( iupdwn( iss ) : iupdwn( iss ) + nupdwn( iss ) - 1 ) = occ_ ( 1:nupdwn( iss ) )
               END IF
               !
            END IF
            !
            DEALLOCATE ( occ_ )
            !
            IF( .NOT. found ) THEN
               ierr = 1
               RETURN
            END IF
            !
         END DO
         !
#if !defined __HDF5
         CALL iotk_scan_end( iunit, "K-POINT" // TRIM( iotk_index(ik) ), IERR = ierr )
#endif
         IF (ierr /= 0) RETURN
         !
      END DO k_points_loop1
      !
#if !defined __HDF5
      CALL iotk_scan_end  ( iunit, "EIGENVALUES", IERR = ierr )
#endif
      IF (ierr /= 0) RETURN
      !
      !
    END SUBROUTINE qexml_read_bands_cp
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_wfc( ibnds, ibnde, ik, ispin, ipol, igk, ngw, igwx, &
                               gamma_only, wf, wf_kindip, ierr )
      !------------------------------------------------------------------------
      !
      ! read wfc from IBNDS to IBNDE, for kpt IK and spin ISPIN
      ! WF is the wfc on its proper k+g grid, while WF_KINDIP is the same wfc
      ! but on a truncated rho grid (k-point indipendent)
      !
      INTEGER,                 INTENT(in)  :: ibnds, ibnde, ik
      INTEGER,       OPTIONAL, INTENT(in)  :: ispin, ipol
      INTEGER,       OPTIONAL, INTENT(in)  :: igk(:)
      INTEGER,       OPTIONAL, INTENT(out) :: ngw, igwx
      LOGICAL,       OPTIONAL, INTENT(out) :: gamma_only
      COMPLEX(DP),  OPTIONAL, INTENT(out) :: wf(:,:), wf_kindip(:,:)
      INTEGER,                 INTENT(out) :: ierr
      !
      INTEGER :: iunaux
      INTEGER :: ngw_, igwx_, ig, ib, lindex
      LOGICAL :: gamma_only_
      COMPLEX(DP),  ALLOCATABLE :: wf_(:)
      CHARACTER(256)             :: filename

      ierr = 0
      !
      !
      ! few check
      !
      IF ( present( ispin ) .and. present( ipol )  ) THEN
         !
         ierr = 1
         RETURN
         !
      ENDIF
      !
      !
      ! read the main data
      !
      CALL iotk_free_unit( iunaux )
      !
      IF ( present( ispin ) ) THEN
         !
         filename = trim( qexml_wfc_filename( datadir_in, 'evc', ik, ispin ) )
         !
      ELSEIF ( present( ipol )  ) THEN
         !
         filename = trim( qexml_wfc_filename( datadir_in, 'evc', ik, ipol ) )
         !
      ELSE
         !
         filename = trim( qexml_wfc_filename( datadir_in, 'evc', ik ) )
         !
      ENDIF
      !
      CALL iotk_open_read ( iunaux, FILE = trim(filename), IERR=ierr )
      IF (ierr/=0)  RETURN
      !
      !
      CALL iotk_scan_empty( iunaux, "INFO", ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_attr( attr, "ngw",  ngw_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "igwx", igwx_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( qexml_version_before_1_4_0 ) THEN
         !
         IF ( present( gamma_only ) ) THEN
             !
             CALL qexml_read_planewaves( GAMMA_ONLY=gamma_only_, IERR=ierr)
             IF (ierr/=0)  RETURN
             !
         ENDIF
         !
      ELSE
         !
         CALL iotk_scan_attr( attr, 'gamma_only', gamma_only_, IERR=ierr)
         IF (ierr/=0)  RETURN
         !
      ENDIF
      !
      !
      IF ( present( wf )  )  THEN
          !
          lindex = 0
          !
          DO ib = ibnds, ibnde
              !
              lindex = lindex + 1
              !
              CALL iotk_scan_dat( iunaux, "evc"//trim(iotk_index(ib)), &
                                  wf( 1:igwx_, lindex ), IERR=ierr )
              IF (ierr/=0) RETURN
              !
          ENDDO
          !
      ENDIF
      !
      IF ( present( wf_kindip )  )  THEN
          !
          ALLOCATE( wf_(igwx_ ), STAT=ierr )
          IF (ierr/=0) RETURN
          !
          IF ( .not. present( igk ) ) THEN
              ierr = 3
              RETURN
          ENDIF
          !
          IF ( maxval( igk( 1: igwx_ ) ) > size( wf_kindip, 1)  ) THEN
              ierr = 4
              RETURN
          ENDIF
          !
          !
          lindex = 0
          !
          DO ib = ibnds, ibnde
              !
              lindex = lindex + 1
              !
              CALL iotk_scan_dat( iunaux, "evc"//trim(iotk_index( ib ) ), &
                                           wf_(1:igwx_), IERR=ierr )
              IF (ierr/=0) RETURN
              !
              ! use the igk map to do the transformation
              !
              wf_kindip(:, lindex) = 0.0_DP
              !
              DO ig = 1, igwx_
                  !
                  wf_kindip( igk( ig ), lindex ) = wf_( ig )
                  !
              ENDDO
              !
          ENDDO
          !
          DEALLOCATE( wf_, STAT=ierr )
          IF (ierr/=0) RETURN
          !
      ENDIF
      !
      CALL iotk_close_read ( iunaux, IERR=ierr )
      IF (ierr/=0)  RETURN
      !
      !
      IF ( present( ngw ) )                 ngw   = ngw_
      IF ( present( igwx ) )               igwx   = igwx_
      IF ( present( gamma_only ) )   gamma_only   = gamma_only_
      !
    END SUBROUTINE qexml_read_wfc
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_rho( nr1, nr2, nr3, rho, ip, rhoz, ierr  )
      !------------------------------------------------------------------------
      !
      ! Reads charge density rho, as a whole or one plane at a time.
      ! if RHO is specified, the whole charge density is read;
      ! if RHOZ is specified only the IP-th plane is read
      !
      IMPLICIT NONE
      !
      INTEGER,   OPTIONAL, INTENT(out) :: nr1, nr2, nr3
      INTEGER,   OPTIONAL, INTENT(in)  :: ip
      REAL(DP), OPTIONAL, INTENT(out) :: rho(:,:,:), rhoz(:)
      INTEGER,             INTENT(out) :: ierr
      !
      INTEGER        :: nr1_, nr2_, nr3_, ip_
      INTEGER        :: iunaux
      LOGICAL        :: lexists
      CHARACTER(256) :: filename

      ierr = 0
      !
      !
      CALL iotk_free_unit( iunaux )
      !
      filename = trim( datadir_in ) // '/' // 'charge-density.dat'
      lexists  = check_file_exst( trim(filename) )
      !
      IF ( .not. lexists ) THEN
          !
          filename = trim( datadir_in ) // '/' // 'charge-density.xml'
          lexists  = check_file_exst( trim(filename) )
          !
      ENDIF
      !
      IF ( .not. lexists ) THEN
          !
          ierr = -1
          RETURN
          !
      ENDIF
      !
      CALL iotk_open_read( iunaux, FILE = filename, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      !
      CALL iotk_scan_begin( iunaux, "CHARGE-DENSITY", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_empty( iunaux, "INFO", ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_attr( attr, "nr1", nr1_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "nr2", nr2_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "nr3", nr3_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( present( rhoz ) ) THEN
         !
         IF ( .not. present( ip ) ) THEN
            ierr = 71
            RETURN
         ENDIF
         !
         CALL iotk_scan_dat( iunaux, "z"//trim(iotk_index(ip)), rhoz, IERR=ierr)
         IF (ierr/=0) RETURN
         !
      ENDIF
      !
      !
      IF ( present( rho ) ) THEN
         !
         DO ip_ = 1, nr3_
            !
            CALL iotk_scan_dat( iunaux, "z"//trim(iotk_index(ip_)), rho(1:nr1_,1:nr2_,ip_), &
                                IERR=ierr)
            IF (ierr/=0) RETURN
            !
         ENDDO
         !
      ENDIF
      !
      CALL iotk_scan_end( iunaux, "CHARGE-DENSITY", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_close_read( iunaux, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( present( nr1 ) )     nr1 = nr1_
      IF ( present( nr2 ) )     nr2 = nr2_
      IF ( present( nr3 ) )     nr3 = nr3_
      !
    END SUBROUTINE qexml_read_rho
    !
#endif
    !
END MODULE qexml_module
