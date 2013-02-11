!
! Copyright (C) 2006 WanT Group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE qexml_module
  !----------------------------------------------------------------------------
  !
  ! This module contains some common subroutines used to read and write
  ! in XML format the data produced by Quantum ESPRESSO package.
  !
  ! Written by Andrea Ferretti (2006).
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
  !
  USE iotk_module
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
  INTEGER,   PARAMETER :: dbl = selected_real_kind( 14, 200 )
  REAL(dbl), PARAMETER :: e2 = 2.0_dbl
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
  PUBLIC :: qexml_write_header, qexml_write_cell, qexml_write_ions,   &
            qexml_write_symmetry, qexml_write_efield,                 &
            qexml_write_planewaves, qexml_write_spin, qexml_write_xc, &
            qexml_write_occ, qexml_write_bz, qexml_write_phonon,      &
            qexml_write_bands, qexml_write_bands_info,                &
            qexml_write_gk, qexml_write_wfc, qexml_write_rho
  !
  PUBLIC :: qexml_read_header, qexml_read_cell, qexml_read_ions,      &
            qexml_read_symmetry, qexml_read_efield,                   &
            qexml_read_planewaves, qexml_read_spin, qexml_read_xc,    &
            qexml_read_occ, qexml_read_bz, qexml_read_phonon,         &
            qexml_read_bands, qexml_read_bands_info,                  &
            qexml_read_gk, qexml_read_wfc, qexml_read_rho

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
      ! the presence of directories overwirtes any info
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
    SUBROUTINE create_directory( dirname, ierr )
      !------------------------------------------------------------------------
      !
      CHARACTER(len=*), INTENT(in)  :: dirname
      INTEGER,          INTENT(out) :: ierr
      !
      INTEGER, EXTERNAL    :: c_mkdir_int
      INTEGER  :: iunaux, i, ilen
      INTEGER, ALLOCATABLE :: istr(:)
      !
      ierr = 0
      CALL iotk_free_unit( iunaux )
      !
      ilen = len_trim( dirname )
      ALLOCATE( istr( ilen ) )
      DO i = 1, ilen
         istr(i) = ichar( dirname(i:i) )
      ENDDO
      !
      ierr = c_mkdir_int( istr, ilen )
      !
      DEALLOCATE( istr )
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
    END SUBROUTINE create_directory
    !
    !
    !------------------------------------------------------------------------
    FUNCTION kpoint_dirname( basedir, ik )
      !------------------------------------------------------------------------
      !
      CHARACTER(len=256)           :: kpoint_dirname
      CHARACTER(len=*), INTENT(in) :: basedir
      INTEGER,          INTENT(in) :: ik
      !
      CHARACTER(len=256) :: kdirname
      CHARACTER(len=5)   :: kindex
      !
      WRITE( kindex, FMT = '( I5.5 )' ) ik
      !
      kdirname = trim( basedir ) // '/K' // kindex
      !
      kpoint_dirname = trim( kdirname )
      !
      RETURN
      !
    END FUNCTION kpoint_dirname
    !
    !
    !------------------------------------------------------------------------
    FUNCTION wfc_filename( basedir, name, ik, ipol, tag, extension )
      !------------------------------------------------------------------------
      !
      CHARACTER(len=256)                 :: wfc_filename
      CHARACTER(len=*),       INTENT(in) :: basedir
      CHARACTER(len=*),       INTENT(in) :: name
      INTEGER,                INTENT(in) :: ik
      INTEGER,      OPTIONAL, INTENT(in) :: ipol
      CHARACTER(*), OPTIONAL, INTENT(in) :: tag
      CHARACTER(*), OPTIONAL, INTENT(in) :: extension
      !
      CHARACTER(len=256) :: filename, tag_, ext_
      !
      !
      filename = ''
      tag_     = ''
      ext_     = '.dat'
      !
      IF ( present( tag ) )         tag_ = '_'//trim(tag)
      IF ( present( extension ) )   ext_ = '.'//trim(extension)
      !
      IF ( present( ipol ) ) THEN
         !
         WRITE( filename, FMT = '( I1 )' ) ipol
         !
      ENDIF
      !
      filename = trim( kpoint_dirname( basedir, ik ) ) // '/' // &
                 & trim( name ) // trim( filename ) // trim( tag_ ) // trim( ext_)
      !
      wfc_filename = trim( filename )
      !
      RETURN
      !
    END FUNCTION wfc_filename
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE copy_file( file_in, file_out, ierr )
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
    END SUBROUTINE copy_file
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
    FUNCTION restart_dirname( outdir, prefix )
      !------------------------------------------------------------------------
      !
      CHARACTER(len=256)           :: restart_dirname
      CHARACTER(len=*), INTENT(in) :: outdir, prefix
      !
      CHARACTER(len=256)         :: dirname
      INTEGER                    :: strlen
      !
      ! ... main restart directory
      !
      dirname = trim( prefix ) // '.save'
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
      restart_dirname = trim( dirname )
      !
      RETURN
      !
    END FUNCTION restart_dirname
    !
    !
!
!-------------------------------------------
! ... write subroutines
!-------------------------------------------
!
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_header( creator_name, creator_version )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in) :: creator_name, creator_version
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
    SUBROUTINE qexml_write_cell( ibravais_latt, symm_type, celldm, alat, &
                                 a1, a2, a3, b1, b2, b3, alat_units, a_units, b_units )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(in) :: ibravais_latt
      CHARACTER(len=*), INTENT(in) :: symm_type
      REAL(dbl),        INTENT(in) :: celldm(6), alat
      REAL(dbl),        INTENT(in) :: a1(3), a2(3), a3(3)
      REAL(dbl),        INTENT(in) :: b1(3), b2(3), b3(3)
      CHARACTER(len=*), INTENT(in) :: alat_units, a_units, b_units
      !
      CHARACTER(len=256) :: bravais_lattice
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
        CASE( -5 )
           bravais_lattice = "Trigonal R (3-fold axis <111>)"
        CASE( -9 )
           bravais_lattice = "Orthorhombic base-centered(bco), alt. axes"
        CASE( -12 )
           bravais_lattice = "Monoclinic P, alt. axis (unique axis b)"
        CASE DEFAULT
           CALL infomsg('qexml_write_cell',&
                'Unexpected value for ibrav, setting ibrav=0 in the XML')
           bravais_lattice = "free"
      END SELECT
      !
      CALL iotk_write_dat( ounit, &
                           "BRAVAIS_LATTICE", trim( bravais_lattice ) )
      !
      IF ( LEN_TRIM( symm_type) > 0 ) &
          CALL iotk_write_dat( ounit, "CELL_SYMMETRY", symm_type )
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
    SUBROUTINE qexml_write_ions( nsp, nat, atm, ityp, psfile, pseudo_dir,  &
                                 amass, amass_units, tau, tau_units, &
                                 if_pos, dirname )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(in) :: nsp, nat
      INTEGER,          INTENT(in) :: ityp(:)
      CHARACTER(len=*), INTENT(in) :: atm(:)
      CHARACTER(len=*), INTENT(in) :: psfile(:)
      CHARACTER(len=*), INTENT(in) :: pseudo_dir
      CHARACTER(len=*), INTENT(in) :: dirname
      REAL(dbl),        INTENT(in) :: amass(:)
      CHARACTER(len=*), INTENT(in) :: amass_units
      REAL(dbl),        INTENT(in) :: tau(:,:)
      CHARACTER(len=*), INTENT(in) :: tau_units
      INTEGER,          INTENT(in) :: if_pos(:,:)
      !
      INTEGER            :: i, flen, ierrl
      CHARACTER(len=256) :: file_pseudo
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
         IF ( pseudo_dir(flen:flen) /= '/' ) THEN
            !
            file_pseudo = pseudo_dir(1:flen) // '/' // psfile(i)
            !
         ELSE
            !
            file_pseudo = pseudo_dir(1:flen) // psfile(i)
            !
         ENDIF
         !
         INQUIRE( FILE = trim( dirname ) // "/" &
                       & // trim( psfile(i) ), EXIST = pseudo_exists )
         !
         IF ( .not. pseudo_exists ) THEN
               CALL copy_file( trim( file_pseudo ), &
                               trim( dirname ) // "/" // trim( psfile(i) ), ierrl )
         ENDIF
         !
         CALL iotk_write_dat( ounit, "MASS", amass(i) )
         !
         CALL iotk_write_dat( ounit, "PSEUDO", trim( psfile(i) ) )
         !
         !
         CALL iotk_write_end( ounit, "SPECIE"//trim(iotk_index(i)) )
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
         CALL iotk_write_attr( attr, "tau",    tau(:,i) )
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
    SUBROUTINE qexml_write_symmetry( nsym, invsym, trasl, s, sname, s_units, &
                                     t_rev, irt, nat )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(in) :: nsym
      LOGICAL,          INTENT(in) :: invsym
      INTEGER,          INTENT(in) :: s(:,:,:)
      REAL(dbl),        INTENT(in) :: trasl(:,:)
      CHARACTER(len=*), INTENT(in) :: sname(:)
      CHARACTER(len=*), INTENT(in) :: s_units
      INTEGER,          INTENT(in) :: t_rev(:)
      INTEGER,          INTENT(in) :: irt(:,:), nat
      !
      INTEGER   :: i
      REAL(dbl) :: tmp(3)
      !
      !
      CALL iotk_write_begin( ounit, "SYMMETRIES" )
      !
      CALL iotk_write_dat( ounit, "NUMBER_OF_SYMMETRIES", nsym )
      !
      CALL iotk_write_dat( ounit, "INVERSION_SYMMETRY", invsym )
      !
      CALL iotk_write_dat( ounit, "NUMBER_OF_ATOMS", nat )
      !
      CALL iotk_write_attr( attr, "UNITS", trim(s_units), FIRST = .true. )
      CALL iotk_write_empty( ounit, "UNITS_FOR_SYMMETRIES", ATTR = attr )
      !
      DO i = 1, nsym
         !
         CALL iotk_write_begin( ounit, "SYMM" // trim( iotk_index( i ) ) )
         !
         CALL iotk_write_attr ( attr, "NAME", trim( sname(i) ), FIRST=.true. )
         CALL iotk_write_attr ( attr, "T_REV", t_rev(i) )
         CALL iotk_write_empty( ounit, "INFO", ATTR = attr )
         !
         tmp(1) = trasl(1,i)
         tmp(2) = trasl(2,i)
         tmp(3) = trasl(3,i)
         !
         CALL iotk_write_dat( ounit, "ROTATION", s(:,:,i), COLUMNS=3 )
         CALL iotk_write_dat( ounit, "FRACTIONAL_TRANSLATION", tmp(1:3), COLUMNS=3 )
         CALL iotk_write_dat( ounit, "EQUIVALENT_IONS", irt(i,1:nat), COLUMNS=8 )
         !
         CALL iotk_write_end( ounit, "SYMM" // trim( iotk_index( i ) ) )
         !
      ENDDO
      !
      CALL iotk_write_end( ounit, "SYMMETRIES" )
      !
    END SUBROUTINE qexml_write_symmetry
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_efield( tefield, dipfield, edir, emaxpos, eopreg, eamp )
      !------------------------------------------------------------------------
      !
      LOGICAL, INTENT(in)   :: tefield        ! if .TRUE. a finite electric field
                                              ! is added to the local potential
      LOGICAL, INTENT(in)   :: dipfield       ! if .TRUE. the dipole field is subtracted
      INTEGER, INTENT(in)   :: edir           ! direction of the field
      REAL(dbl), INTENT(in) :: emaxpos        ! position of the maximum of the field (0<emaxpos<1)
      REAL(dbl), INTENT(in) :: eopreg         ! amplitude of the inverse region (0<eopreg<1)
      REAL(dbl), INTENT(in) :: eamp           ! field amplitude (in a.u.) (1 a.u. = 51.44 10^11 V/m)
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
      !
      INTEGER,       INTENT(in) :: npwx, nr1, nr2, nr3, ngm, &
                                   nr1s, nr2s, nr3s, ngms, nr1b, nr2b, nr3b
      INTEGER,       INTENT(in) :: igv(:,:)
      REAL(dbl),     INTENT(in) :: ecutwfc, ecutrho
      LOGICAL,       INTENT(in) :: gamma_only, lgvec
      CHARACTER(*),  INTENT(in) :: cutoff_units
      !
      !
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
         CALL iotk_link( ounit, "G-VECTORS", "./gvectors.dat", &
                         CREATE = .true., BINARY = .true. )
         !
         CALL iotk_write_begin( ounit, "G-VECTORS" )
         !
         CALL iotk_write_attr( attr, "nr1s", nr1s, FIRST = .true. )
         CALL iotk_write_attr( attr, "nr2s", nr2s )
         CALL iotk_write_attr( attr, "nr3s", nr3s )
         CALL iotk_write_attr( attr, "gamma_only", gamma_only )
         CALL iotk_write_attr( attr, "units", "crystal" )
         CALL iotk_write_empty( ounit, "INFO", ATTR = attr )
         !
         CALL iotk_write_dat  ( ounit, "g", igv(1:3,1:ngm), COLUMNS = 3 )
         CALL iotk_write_end  ( ounit, "G-VECTORS" )
         !
      ENDIF
      !
      CALL iotk_write_attr( attr, "nr1b", nr1b , FIRST = .true. )
      CALL iotk_write_attr( attr, "nr2b", nr2b )
      CALL iotk_write_attr( attr, "nr3b", nr3b )
      CALL iotk_write_empty( ounit, "SMALLBOX_FFT_GRID", ATTR = attr )
      !
      CALL iotk_write_end( ounit, "PLANE_WAVES" )
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
      REAL(dbl),    INTENT(in) :: xk(3)
      CHARACTER(*), INTENT(in) :: k_units
      LOGICAL,      INTENT(in) :: index(:), igk(:,:)
      !
      INTEGER        :: iunaux
      CHARACTER(256) :: filename

      CALL iotk_free_unit( iunaux )
      filename = wfc_filename( datadir_out, 'gkvectors', ik )
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
    SUBROUTINE qexml_write_xc( dft, lda_plus_u, Hubbard_lmax, Hubbard_l, &
                               nsp, Hubbard_U, Hubbard_alpha )
      !------------------------------------------------------------------------
      !
      CHARACTER(len=*),    INTENT(in) :: dft
      LOGICAL,             INTENT(in) :: lda_plus_u
      INTEGER,   OPTIONAL, INTENT(in) :: nsp
      INTEGER,   OPTIONAL, INTENT(in) :: Hubbard_lmax
      INTEGER,   OPTIONAL, INTENT(in) :: Hubbard_l(:)
      REAL(dbl), OPTIONAL, INTENT(in) :: Hubbard_U(:), Hubbard_alpha(:)
      !
      !
      CALL iotk_write_begin( ounit, "EXCHANGE_CORRELATION" )
      !
      CALL iotk_write_dat( ounit, "DFT", dft )
      !
      CALL iotk_write_dat( ounit, "LDA_PLUS_U_CALCULATION", lda_plus_u )
      !
      IF ( lda_plus_u ) THEN
         !
         IF ( .not. present( Hubbard_lmax ) .or. &
              .not. present( Hubbard_l )    .or. &
              .not. present( Hubbard_U )    .or. &
              .not. present( nsp )          .or. &
              .not. present( Hubbard_alpha ) ) &
            CALL errore( 'write_exchange_correlation', &
                         ' variables for LDA+U not present', 1 )
         !
         CALL iotk_write_dat( ounit, "NUMBER_OF_SPECIES", nsp )
         !
         CALL iotk_write_dat( ounit, "HUBBARD_LMAX", Hubbard_lmax )
         !
         CALL iotk_write_dat( ounit, "HUBBARD_L", &
                              Hubbard_l(1:nsp) )
         !
         CALL iotk_write_dat( ounit, "HUBBARD_U", Hubbard_U(1:nsp) )
         !
         CALL iotk_write_dat( ounit, "HUBBARD_ALPHA", Hubbard_alpha(1:nsp) )
         !
      ENDIF
      !
      CALL iotk_write_end( ounit, "EXCHANGE_CORRELATION" )
      !
    END SUBROUTINE qexml_write_xc
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_occ( lgauss, ngauss, degauss, degauss_units, ltetra, ntetra, tetra, &
                                tfixed_occ, lsda, nstates_up, nstates_dw, input_occ )
      !------------------------------------------------------------------------
      !
      LOGICAL,                INTENT(in) :: lgauss, ltetra, tfixed_occ, lsda
      INTEGER,      OPTIONAL, INTENT(in) :: ngauss, ntetra, nstates_up, nstates_dw
      INTEGER,      OPTIONAL, INTENT(in) :: tetra(:,:)
      REAL(dbl),    OPTIONAL, INTENT(in) :: degauss, input_occ(:,:)
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
         DO i = 1, ntetra
            !
            CALL iotk_write_dat( ounit, "TETRAHEDRON" // &
                               & iotk_index( i ), tetra(1:4,i) )
            !
         ENDDO
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
                               nk1, nk2, nk3, k_units )
      !------------------------------------------------------------------------
      !
      INTEGER,      INTENT(in) :: num_k_points, k1, k2, k3, nk1, nk2, nk3
      REAL(dbl),    INTENT(in) :: xk(:,:), wk(:)
      CHARACTER(*), INTENT(in) :: k_units
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
      CALL iotk_write_end( ounit, "BRILLOUIN_ZONE" )
      !
    END SUBROUTINE qexml_write_bz
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_phonon( modenum, xqq, q_units )
      !------------------------------------------------------------------------
      !
      INTEGER,      INTENT(in) :: modenum
      REAL(dbl),    INTENT(in) :: xqq(:)
      CHARACTER(*), INTENT(in) :: q_units
      !
      !
      CALL iotk_write_begin( ounit, "PHONON" )
      !
      CALL iotk_write_dat( ounit, "NUMBER_OF_MODES", modenum )
      !
      CALL iotk_write_attr( attr, "UNITS", trim(q_units), FIRST = .true. )
      CALL iotk_write_empty( ounit, "UNITS_FOR_Q-POINT", attr )
      !
      CALL iotk_write_dat( ounit, "Q-POINT", xqq(:), COLUMNS=3 )
      !
      CALL iotk_write_end( ounit, "PHONON" )
      !
    END SUBROUTINE qexml_write_phonon
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_bands_info( nbnd, num_k_points, nspin, noncolin, natomwfc, &
                                       nelec, ef, energy_units, k_units )
      !------------------------------------------------------------------------
      !
      INTEGER,       INTENT(in) :: nbnd, num_k_points, nspin, natomwfc
      LOGICAL,       INTENT(in) :: noncolin
      REAL(dbl),     INTENT(in) :: ef, nelec
      CHARACTER(*),  INTENT(in) :: energy_units, k_units
      !
      !
      CALL iotk_write_begin( ounit, "BAND_STRUCTURE_INFO" )
      !
      CALL iotk_write_dat  ( ounit, "NUMBER_OF_BANDS", nbnd )
      !
      CALL iotk_write_dat  ( ounit, "NUMBER_OF_K-POINTS", num_k_points )
      !
      CALL iotk_write_dat  ( ounit, "NUMBER_OF_SPIN_COMPONENTS", nspin )
      !
      CALL iotk_write_dat  ( ounit, "NON-COLINEAR_CALCULATION", noncolin )
      !
      CALL iotk_write_dat  ( ounit, "NUMBER_OF_ATOMIC_WFC", natomwfc )
      !
      CALL iotk_write_dat  ( ounit, "NUMBER_OF_ELECTRONS", nelec )
      !
      CALL iotk_write_attr ( attr, "UNITS", trim(k_units), FIRST = .true. )
      CALL iotk_write_empty( ounit, "UNITS_FOR_K-POINTS", ATTR = attr )
      !
      CALL iotk_write_attr ( attr, "UNITS", trim(energy_units), FIRST = .true. )
      CALL iotk_write_empty( ounit, "UNITS_FOR_ENERGIES", ATTR = attr )
      !
      CALL iotk_write_dat  ( ounit, "FERMI_ENERGY", ef )
      !
      CALL iotk_write_end  ( ounit, "BAND_STRUCTURE_INFO" )
      !
      RETURN
      !
    END SUBROUTINE qexml_write_bands_info
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_bands( ik, ispin, nbnd, eig, energy_units, occ, ef )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER,             INTENT(in) :: ik, nbnd
      INTEGER, OPTIONAL,   INTENT(in) :: ispin
      REAL(dbl),           INTENT(in) :: eig(:)
      CHARACTER(*),        INTENT(in) :: energy_units
      REAL(dbl), OPTIONAL, INTENT(in) :: occ(:), ef
      !
      INTEGER :: iunaux
      CHARACTER(len=256) :: filename

      !
      IF ( present( ispin) ) THEN
         !
         filename= trim( wfc_filename( datadir_out, 'eigenval', &
                                       ik, ispin, EXTENSION="xml") )
         !
      ELSE
         !
         filename= trim( wfc_filename( datadir_out, 'eigenval', &
                                       ik, EXTENSION="xml") )
         !
      ENDIF
      !
      CALL iotk_free_unit( iunaux )
      CALL iotk_open_write ( iunaux, FILE = trim( filename ), BINARY = .false. )
      !
      CALL iotk_write_attr ( attr, "nbnd", nbnd, FIRST=.true. )
      CALL iotk_write_attr ( attr, "ik", ik )
      !
      IF ( present( ispin) ) CALL iotk_write_attr ( attr, "ispin", ispin )
      !
      CALL iotk_write_empty( iunaux, "INFO", ATTR = attr )
      !
      CALL iotk_write_attr ( attr, "UNITS", trim(energy_units), FIRST = .true. )
      CALL iotk_write_empty( iunaux, "UNITS_FOR_ENERGIES", ATTR=attr)
      !
      IF ( present( ef ) ) THEN
         !
         CALL iotk_write_dat( iunaux, "FERMI_ENERGY", ef)
         !
      ENDIF
      !
      CALL iotk_write_dat( iunaux, "EIGENVALUES", eig(:) )
      !
      IF ( present( occ ) ) THEN
         !
         CALL iotk_write_dat( iunaux, "OCCUPATIONS", occ(:) )
         !
      ENDIF
      !
      CALL iotk_close_write ( iunaux )
      !
    END SUBROUTINE qexml_write_bands
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
      COMPLEX(dbl), OPTIONAL, INTENT(in) :: wf(:,:)
      COMPLEX(dbl), OPTIONAL, INTENT(in) :: wf_kindip(:,:)
      REAL(dbl),    OPTIONAL, INTENT(in) :: scale_factor
      !
      INTEGER         :: iunaux, ierr
      INTEGER         :: ig, ib
      CHARACTER(256)  :: filename
      COMPLEX(dbl),  ALLOCATABLE :: wtmp(:)

      ierr = 0
      !
      IF ( present( ispin ) .and. present( ipol )  ) THEN
         !
         ierr = 1
         RETURN
         !
      ENDIF
      !
      !
      ! open the file to write
      !
      CALL iotk_free_unit( iunaux )
      !
      IF ( present( ispin ) ) THEN
         !
         filename = trim( wfc_filename( datadir_out, 'evc', ik, ispin ) )
         !
      ELSEIF ( present( ipol )  ) THEN
         !
         filename = trim( wfc_filename( datadir_out, 'evc', ik, ipol ) )
         !
      ELSE
         !
         filename = trim( wfc_filename( datadir_out, 'evc', ik ) )
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
      REAL(dbl), OPTIONAL, INTENT(in) :: rho(:,:,:), rhov(:)
      LOGICAL,   OPTIONAL, INTENT(in) :: binary
      !
      INTEGER        :: iunaux, nr1x_, nr2x_, ip, i1, i2, i
      LOGICAL        :: binary_
      CHARACTER(256) :: filename
      REAL(dbl), ALLOCATABLE :: plane(:,:)
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

!
!-------------------------------------------
! ... read subroutines
!-------------------------------------------
!
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_header( creator_name, creator_version, &
                                  format_name, format_version, ierr )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      CHARACTER(len=*),  OPTIONAL, INTENT(out) :: creator_name, creator_version
      CHARACTER(len=*),  OPTIONAL, INTENT(out) :: format_name, format_version
      INTEGER,           OPTIONAL, INTENT(out) :: ierr

      CHARACTER(256) :: creator_name_, creator_version_
      CHARACTER(256) :: format_name_,     format_version_

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
    SUBROUTINE qexml_read_cell( bravais_latt, symm_type, celldm, alat, &
                                a1, a2, a3, b1, b2, b3, alat_units, a_units, b_units, ierr )
      !------------------------------------------------------------------------
      !
      CHARACTER(len=*),  OPTIONAL, INTENT(out) :: bravais_latt
      CHARACTER(len=*),  OPTIONAL, INTENT(out) :: symm_type
      REAL(dbl),         OPTIONAL, INTENT(out) :: celldm(6), alat
      REAL(dbl),         OPTIONAL, INTENT(out) :: a1(3), a2(3), a3(3)
      REAL(dbl),         OPTIONAL, INTENT(out) :: b1(3), b2(3), b3(3)
      CHARACTER(len=*),  OPTIONAL, INTENT(out) :: alat_units, a_units, b_units
      INTEGER,                     INTENT(out) :: ierr
      !
      CHARACTER(256)     :: bravais_latt_,
      CHARACTER(256)     :: alat_units_, a_units_, b_units_
      REAL(dbl)          :: celldm_(6), alat_
      REAL(dbl)          :: a1_(3), a2_(3), a3_(3)
      REAL(dbl)          :: b1_(3), b2_(3), b3_(3)
      !

      ierr=0
      !
      !
      CALL iotk_scan_begin( iunit, "CELL" )
      !
      CALL iotk_scan_dat( iunit, "BRAVAIS_LATTICE", bravais_latt_, IERR=ierr )
      IF ( ierr /= 0 ) RETURN
      !
      IF ( PRESENT( symm_type ) ) THEN
          CALL iotk_scan_dat( iunit, "CELL_SYMMETRY", symm_type, IERR=ierr )
          IF ( ierr /= 0 ) RETURN
      ENDIF
      !
      CALL iotk_scan_dat( iunit, "LATTICE_PARAMETER", alat_, ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "UNITS", alat_units_, IERR=ierr )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "CELL_DIMENSIONS", celldm_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_begin( iunit, "DIRECT_LATTICE_VECTORS", IERR=ierr )
      IF (ierr/=0) RETURN
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
      IF ( present(bravais_latt) )  bravais_latt = bravais_latt_
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

    END SUBROUTINE qexml_read_cell
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_ions( nsp, nat, atm, ityp, psfile, amass, amass_units, &
                                tau, tau_units, if_pos, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,          OPTIONAL, INTENT(out) :: nsp, nat
      INTEGER,          OPTIONAL, INTENT(out) :: ityp(:)
      CHARACTER(len=*), OPTIONAL, INTENT(out) :: atm(:)
      CHARACTER(len=*), OPTIONAL, INTENT(out) :: psfile(:)
      REAL(dbl),        OPTIONAL, INTENT(out) :: amass(:)
      CHARACTER(len=*), OPTIONAL, INTENT(out) :: amass_units
      REAL(dbl),        OPTIONAL, INTENT(out) :: tau(:,:)
      INTEGER,          OPTIONAL, INTENT(out) :: if_pos(:,:)
      CHARACTER(len=*), OPTIONAL, INTENT(out) :: tau_units
      INTEGER,                    INTENT(out) :: ierr
      !
      INTEGER                     :: nat_, nsp_
      CHARACTER(256)              :: tau_units_, amass_units_
      INTEGER,        ALLOCATABLE :: ityp_(:)
      CHARACTER(3),   ALLOCATABLE :: atm_(:)
      CHARACTER(256), ALLOCATABLE :: psfile_(:)
      REAL(dbl),      ALLOCATABLE :: amass_(:)
      REAL(dbl),      ALLOCATABLE :: tau_(:,:)
      INTEGER,        ALLOCATABLE :: if_pos_(:,:)
      !
      INTEGER            :: i

      !
      ierr=0
      !
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
    SUBROUTINE qexml_read_symmetry( nsym, invsym, trasl, s, sname, s_units, t_rev, &
                                    irt, nat, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,          OPTIONAL, INTENT(out) :: nsym
      LOGICAL,          OPTIONAL, INTENT(out) :: invsym
      INTEGER,          OPTIONAL, INTENT(out) :: s(:,:,:)
      REAL(dbl),        OPTIONAL, INTENT(out) :: trasl(:,:)
      CHARACTER(len=*), OPTIONAL, INTENT(out) :: sname(:)
      CHARACTER(len=*), OPTIONAL, INTENT(out) :: s_units
      INTEGER,          OPTIONAL, INTENT(out) :: t_rev(:)
      INTEGER,          OPTIONAL, INTENT(out) :: irt(:,:), nat
      INTEGER,                    INTENT(out) :: ierr
      !
      INTEGER              :: nsym_
      CHARACTER(256)       :: sname_(48), s_units_
      LOGICAL              :: invsym_
      INTEGER              :: s_(3,3,48)
      REAL(dbl)            :: trasl_(3,48)
      INTEGER              :: t_rev_(48)
      INTEGER              :: nat_
      INTEGER, ALLOCATABLE :: irt_(:,:)
      !
      INTEGER             :: i

      !
      ierr=0
      !
      !
      CALL iotk_scan_begin( iunit, "SYMMETRIES", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "NUMBER_OF_SYMMETRIES", nsym_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "INVERSION_SYMMETRY", invsym_, IERR=ierr )
      IF (ierr/=0) RETURN
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
      CALL iotk_scan_end( iunit, "SYMMETRIES", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( present(nsym) )        nsym          = nsym_
      IF ( present(invsym) )      invsym        = invsym_
      IF ( present(nat) )         nat           = nat_
      IF ( present(trasl) )       trasl(1:3, 1:nsym_)   = trasl_(1:3, 1:nsym_)
      IF ( present(s) )           s(1:3, 1:3, 1:nsym_)  = s_(1:3, 1:3, 1:nsym_)
      IF ( present(irt) )         irt(1:nsym_, 1:nat_)  = irt_(1:nsym_, 1:nat_)
      IF ( present(sname) )  THEN
          DO i = 1, nsym_
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
    SUBROUTINE qexml_read_efield( tefield, dipfield, edir, emaxpos, eopreg, eamp, ierr )
      !----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      LOGICAL,   OPTIONAL, INTENT(out) :: tefield, dipfield
      INTEGER,   OPTIONAL, INTENT(out) :: edir
      REAL(dbl), OPTIONAL, INTENT(out) :: emaxpos, eopreg, eamp
      INTEGER,             INTENT(out) :: ierr
      !
      LOGICAL   :: tefield_, dipfield_
      INTEGER   :: edir_
      REAL(dbl) :: emaxpos_, eopreg_, eamp_
      !

      ierr = 0
      !
      CALL iotk_scan_begin( iunit, "ELECTRIC_FIELD", IERR=ierr )
      IF ( ierr /= 0 ) RETURN
      !
      !
      CALL iotk_scan_dat( iunit, "HAS_ELECTRIC_FIELD", tefield_ )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "HAS_DIPOLE_CORRECTION", dipfield_ )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "FIELD_DIRECTION", edir_ )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "MAXIMUM_POSITION", emaxpos_ )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "INVERSE_REGION", eopreg_ )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunit, "FIELD_AMPLITUDE", eamp_ )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_end( iunit, "ELECTRIC_FIELD" )
      IF ( ierr /= 0 ) RETURN
      !
      !
      IF ( present(tefield) )        tefield      = tefield_
      IF ( present(dipfield) )       dipfield     = dipfield_
      IF ( present(edir) )           edir         = edir_
      IF ( present(emaxpos) )        emaxpos      = emaxpos_
      IF ( present(eopreg) )         eopreg       = eopreg_
      IF ( present(eamp) )           eamp         = eamp_
      !
    END SUBROUTINE qexml_read_efield
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
      REAL(dbl),    OPTIONAL, INTENT(out) :: ecutwfc, ecutrho
      LOGICAL,      OPTIONAL, INTENT(out) :: gamma_only
      CHARACTER(*), OPTIONAL, INTENT(out) :: cutoff_units
      INTEGER,                INTENT(out) :: ierr
      !
      INTEGER        :: npwx_, nr1_, nr2_, nr3_, ngm_, &
                        nr1s_, nr2s_, nr3s_, ngms_, nr1b_, nr2b_, nr3b_
      REAL(dbl)      :: ecutwfc_, ecutrho_
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
      REAL(dbl),    OPTIONAL, INTENT(out) :: xk(3)
      CHARACTER(*), OPTIONAL, INTENT(out) :: k_units
      INTEGER,      OPTIONAL, INTENT(out) :: igk(:,:), index(:)
      INTEGER,                INTENT(out) :: ierr
      !
      CHARACTER(256) :: filename, k_units_
      INTEGER   :: npwk_, npwkx_
      LOGICAL   :: gamma_only_
      REAL(dbl) :: xk_(3)
      INTEGER   :: iunaux
      !

      ierr = 0
      !
      CALL iotk_free_unit( iunaux )
      filename = wfc_filename( datadir_in, 'gkvectors', ik )
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
    SUBROUTINE qexml_read_spin( lsda, noncolin, npol, lspinorb, ierr )
      !------------------------------------------------------------------------
      !
      LOGICAL, OPTIONAL, INTENT(out) :: lsda, noncolin, lspinorb
      INTEGER, OPTIONAL, INTENT(out) :: npol
      INTEGER,           INTENT(out) :: ierr
      !
      LOGICAL   :: lsda_, noncolin_, lspinorb_
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
      CALL iotk_scan_dat( iunit, "NON-COLINEAR_CALCULATION", noncolin_, IERR=ierr )
      IF (ierr/=0) RETURN
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
      CALL iotk_scan_dat( iunit, "SPIN-ORBIT_CALCULATION", lspinorb_, IERR=ierr )
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
      !

    END SUBROUTINE qexml_read_spin
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_xc( dft, lda_plus_u, Hubbard_lmax, Hubbard_l, &
                              nsp, Hubbard_U, Hubbard_alpha, ierr )
      !------------------------------------------------------------------------
      !
      CHARACTER(len=*), OPTIONAL, INTENT(out) :: dft
      LOGICAL,          OPTIONAL, INTENT(out) :: lda_plus_u
      INTEGER,          OPTIONAL, INTENT(out) :: Hubbard_lmax
      INTEGER,          OPTIONAL, INTENT(out) :: Hubbard_l(:)
      INTEGER,          OPTIONAL, INTENT(out) :: nsp
      REAL(dbl),        OPTIONAL, INTENT(out) :: Hubbard_U(:), Hubbard_alpha(:)
      INTEGER,                    INTENT(out) :: ierr
      !
      CHARACTER(256) :: dft_
      LOGICAL        :: lda_plus_u_
      INTEGER        :: Hubbard_lmax_, nsp_
      INTEGER,    ALLOCATABLE :: Hubbard_l_(:)
      REAL(dbl),  ALLOCATABLE :: Hubbard_U_(:)
      REAL(dbl),  ALLOCATABLE :: Hubbard_alpha_(:)
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
      CALL iotk_scan_dat( iunit, "LDA_PLUS_U_CALCULATION", lda_plus_u_, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      IF ( lda_plus_u_ ) THEN
         !
         CALL iotk_scan_dat( iunit, "NUMBER_OF_SPECIES", nsp_, IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
         CALL iotk_scan_dat( iunit, "HUBBARD_LMAX", Hubbard_lmax_, IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
         ALLOCATE( Hubbard_l_(nsp_) )
         ALLOCATE( Hubbard_U_(nsp_) )
         ALLOCATE( Hubbard_alpha_(nsp_) )
         !
         CALL iotk_scan_dat( iunit, "HUBBARD_L", Hubbard_l_, IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
         CALL iotk_scan_dat( iunit, "HUBBARD_U", Hubbard_U_, IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
         CALL iotk_scan_dat( iunit, "HUBBARD_ALPHA", Hubbard_alpha_, IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
      ENDIF
      !
      CALL iotk_scan_end( iunit, "EXCHANGE_CORRELATION", IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      !
      IF ( present( dft ) )           dft           = dft_
      IF ( present( lda_plus_u ) )    lda_plus_u    = lda_plus_u_
      !
      IF ( lda_plus_u_ )  THEN
         !
         IF ( present( nsp ) )             nsp                   = nsp_
         IF ( present( Hubbard_lmax ) )    Hubbard_lmax          = Hubbard_lmax_
         IF ( present( Hubbard_l ) )       Hubbard_l(1:nsp_)     = Hubbard_l_(1:nsp_)
         IF ( present( Hubbard_U ) )       Hubbard_U(1:nsp_)     = Hubbard_U_(1:nsp_)
         IF ( present( Hubbard_alpha ) )   Hubbard_alpha(1:nsp_) = Hubbard_alpha_(1:nsp_)
         !
         DEALLOCATE( Hubbard_l_ )
         DEALLOCATE( Hubbard_U_ )
         DEALLOCATE( Hubbard_alpha_ )
         !
      ENDIF

    END SUBROUTINE qexml_read_xc
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_occ( lgauss, ngauss, degauss, degauss_units, &
                               ltetra, ntetra, tetra, tfixed_occ,      &
                               nstates_up, nstates_dw, input_occ, ierr )
      !------------------------------------------------------------------------
      !
      LOGICAL,      OPTIONAL, INTENT(out) :: lgauss, ltetra, tfixed_occ
      INTEGER,      OPTIONAL, INTENT(out) :: ngauss, ntetra
      INTEGER,      OPTIONAL, INTENT(out) :: tetra(:,:)
      INTEGER,      OPTIONAL, INTENT(out) :: nstates_up, nstates_dw
      REAL(dbl),    OPTIONAL, INTENT(out) :: degauss, input_occ(:,:)
      CHARACTER(*), OPTIONAL, INTENT(out) :: degauss_units
      INTEGER,                INTENT(out) :: ierr
      !
      LOGICAL        :: lgauss_, ltetra_, tfixed_occ_
      INTEGER        :: ngauss_, ntetra_, nstates_up_, nstates_dw_
      LOGICAL        :: lsda_
      REAL(dbl)      :: degauss_
      CHARACTER(256) :: degauss_units_
      INTEGER,  ALLOCATABLE :: tetra_(:,:)
      INTEGER :: i
      LOGICAL :: lfound
      !
      ierr = 0
      !
      CALL iotk_scan_begin( iunit, "OCCUPATIONS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunit, "SMEARING_METHOD", lgauss_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
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
         CALL iotk_scan_attr( ATTR, "UNITS", degauss_units_ , IERR=ierr )
         IF (ierr/=0) RETURN
         !
      ENDIF
      !
      CALL iotk_scan_dat( iunit, "TETRAHEDRON_METHOD", ltetra_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( ltetra_ ) THEN
         !
         CALL iotk_scan_dat( iunit, "NUMBER_OF_TETRAHEDRA", ntetra_, IERR=ierr )
         IF (ierr/=0) RETURN
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
      ENDIF
      !
      CALL iotk_scan_dat( iunit, "FIXED_OCCUPATIONS", tfixed_occ_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( tfixed_occ_  .and. ( present( input_occ ) .or. &
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
      IF ( present( degauss ))          degauss     = degauss
      IF ( present( degauss_units ))    degauss_units  = trim(degauss_units_)
      IF ( present( nstates_up ))       nstates_up  = nstates_up_
      IF ( present( nstates_dw ))       nstates_dw  = nstates_dw_
      !
      IF ( ltetra_ ) THEN
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
                              k_units, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,       OPTIONAL, INTENT(out) :: num_k_points, k1, k2, k3, nk1, nk2, nk3
      REAL(dbl),     OPTIONAL, INTENT(out) :: xk(:,:), wk(:)
      CHARACTER(*),  OPTIONAL, INTENT(out) :: k_units
      INTEGER,                 INTENT(out) :: ierr
      !
      INTEGER                :: num_k_points_, k1_, k2_, k3_, nk1_, nk2_, nk3_
      CHARACTER(256)         :: k_units_
      REAL(dbl), ALLOCATABLE :: xk_(:,:), wk_(:)
      !
      INTEGER :: ik
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
      !
      DEALLOCATE( xk_ )
      DEALLOCATE( wk_ )
      !
    END SUBROUTINE qexml_read_bz
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_phonon( modenum, xqq, q_units, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,       OPTIONAL, INTENT(out) :: modenum
      REAL(dbl),     OPTIONAL, INTENT(out) :: xqq(:)
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
    SUBROUTINE qexml_read_bands_info( nbnd, num_k_points, nspin, noncolin, natomwfc, &
                                      nelec, ef, energy_units, k_units, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,      OPTIONAL, INTENT(out) :: nbnd, num_k_points, nspin, natomwfc
      LOGICAL,      OPTIONAL, INTENT(out) :: noncolin
      REAL(dbl),    OPTIONAL, INTENT(out) :: ef, nelec
      CHARACTER(*), OPTIONAL, INTENT(out) :: energy_units, k_units
      INTEGER,                INTENT(out) :: ierr
      !
      INTEGER        :: nbnd_, num_k_points_, nspin_, natomwfc_
      LOGICAL        :: noncolin_
      REAL(dbl)      :: ef_, nelec_
      CHARACTER(256) :: energy_units_, k_units_

      ierr = 0
      !
      !
      CALL iotk_scan_begin( iunit, "BAND_STRUCTURE_INFO", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat  ( iunit, "NUMBER_OF_BANDS", nbnd_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat  ( iunit, "NUMBER_OF_K-POINTS", num_k_points_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat  ( iunit, "NUMBER_OF_SPIN_COMPONENTS", nspin_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat  ( iunit, "NON-COLINEAR_CALCULATION", noncolin_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat  ( iunit, "NUMBER_OF_ATOMIC_WFC", natomwfc_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat  ( iunit, "NUMBER_OF_ELECTRONS", nelec_, IERR=ierr )
      IF (ierr/=0) RETURN
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
      CALL iotk_scan_dat  ( iunit, "FERMI_ENERGY", ef_ , IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_end( iunit, "BAND_STRUCTURE_INFO", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( present( nbnd ) )             nbnd           = nbnd_
      IF ( present( num_k_points ) )     num_k_points   = num_k_points_
      IF ( present( nspin ) )            nspin          = nspin_
      IF ( present( noncolin ) )         noncolin       = noncolin_
      IF ( present( natomwfc ) )         natomwfc       = natomwfc_
      IF ( present( nelec ) )            nelec          = nelec_
      IF ( present( ef ) )               ef             = ef_
      IF ( present( energy_units ) )     energy_units   = trim( energy_units_ )
      IF ( present( k_units ) )          k_units        = trim( k_units_ )
      !
    END SUBROUTINE qexml_read_bands_info
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_bands( ik, ispin, nbnd, eig, energy_units, occ, ef, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,                INTENT(in)  :: ik
      INTEGER,      OPTIONAL, INTENT(in)  :: ispin
      INTEGER,      OPTIONAL, INTENT(out) :: nbnd
      REAL(dbl),    OPTIONAL, INTENT(out) :: eig(:)
      CHARACTER(*), OPTIONAL, INTENT(out) :: energy_units
      REAL(dbl),    OPTIONAL, INTENT(out) :: occ(:)
      REAL(dbl),    OPTIONAL, INTENT(out) :: ef
      INTEGER,                INTENT(out) :: ierr
      !
      INTEGER        :: iunaux
      INTEGER        :: nbnd_
      CHARACTER(256) :: energy_units_
      CHARACTER(256) :: filename
      REAL(dbl), ALLOCATABLE :: occ_(:), eig_(:)
      !

      ierr = 0
      !
      !
      ! read the main data
      !
      CALL iotk_free_unit( iunaux )
      !
      IF ( present( ispin) ) THEN
         !
         filename= trim( wfc_filename( datadir_in, 'eigenval', &
                                       ik, ispin, EXTENSION="xml") )
         !
      ELSE
         !
         filename= trim( wfc_filename( datadir_in, 'eigenval', &
                                       ik, EXTENSION="xml") )
         !
      ENDIF
      !
      !
      CALL iotk_open_read ( iunaux, FILE = trim(filename), IERR=ierr )
      IF (ierr/=0)  RETURN
      !
      CALL iotk_scan_empty( iunaux, "INFO", ATTR = attr, IERR=ierr )
      IF (ierr/=0)  RETURN
      CALL iotk_scan_attr( attr, "nbnd", nbnd_, IERR=ierr )
      IF (ierr/=0)  RETURN
      !
      CALL iotk_scan_empty( iunaux, "UNITS_FOR_ENERGIES", ATTR = attr, IERR=ierr )
      IF (ierr/=0)  RETURN
      CALL iotk_scan_attr( attr, "UNITS", energy_units_, IERR=ierr )
      IF (ierr/=0)  RETURN
      !
      IF ( present( ef )) THEN
         !
         CALL iotk_scan_dat( iunaux, "FERMI_ENERGY", ef, IERR=ierr )
         IF (ierr/=0)  RETURN
         !
      ENDIF
      !
      !
      ! Allocations
      !
      ALLOCATE(  eig_ ( nbnd_ ) )
      ALLOCATE(  occ_ ( nbnd_ ) )
      !
      CALL iotk_scan_dat( iunaux, "EIGENVALUES", eig_(:), IERR=ierr)
      IF (ierr/=0)  RETURN
      !
      CALL iotk_scan_dat( iunaux, "OCCUPATIONS", occ_(:), IERR=ierr)
      IF (ierr/=0)  RETURN
      !
      CALL iotk_close_read ( iunaux, IERR=ierr )
      IF (ierr/=0)  RETURN
      !
      !
      IF ( present( nbnd ) )             nbnd             = nbnd_
      IF ( present( energy_units ) )     energy_units     = trim( energy_units_ )
      IF ( present( occ ) )              occ  (1:nbnd_ )  = occ_(:)
      IF ( present( eig ) )              eig  (1:nbnd_ )  = eig_(:)
      !
      DEALLOCATE( occ_ )
      DEALLOCATE( eig_ )
      !
    END SUBROUTINE qexml_read_bands
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
      COMPLEX(dbl),  OPTIONAL, INTENT(out) :: wf(:,:), wf_kindip(:,:)
      INTEGER,                 INTENT(out) :: ierr
      !
      INTEGER :: iunaux
      INTEGER :: ngw_, igwx_, ig, ib, lindex
      LOGICAL :: gamma_only_
      COMPLEX(dbl),  ALLOCATABLE :: wf_(:)
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
         filename = trim( wfc_filename( datadir_in, 'evc', ik, ispin ) )
         !
      ELSEIF ( present( ipol )  ) THEN
         !
         filename = trim( wfc_filename( datadir_in, 'evc', ik, ipol ) )
         !
      ELSE
         !
         filename = trim( wfc_filename( datadir_in, 'evc', ik ) )
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
              wf_kindip(:, lindex) = 0.0_dbl
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
      REAL(dbl), OPTIONAL, INTENT(out) :: rho(:,:,:), rhoz(:)
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
    !
END MODULE qexml_module

