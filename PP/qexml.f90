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
  ! in XML format the data produced by Quantum-ESPRESSO package.
  !
  ! Written by Andrea Ferretti (2006).
  !
  ! Important parts of the implementation are taken from xml_io_base.f90
  ! (written by Carlo Sbraccia) in the Quantum-ESPRESSO distribution,
  ! under the GNU-GPL licensing:
  !
  ! Copyright (C) 2005 Quantum-ESPRESSO group
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
  INTEGER,   PARAMETER :: dbl = SELECTED_REAL_KIND( 14, 200 )
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
  CHARACTER(10)    :: qexml_default_version = TRIM( fmt_version  )
  LOGICAL          :: qexml_current_version_init = .FALSE.
  LOGICAL          :: qexml_version_before_1_4_0 = .FALSE.
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
    SUBROUTINE qexml_init( unit_in, unit_out, dir, dir_in, dir_out )
      !------------------------------------------------------------------------
      !
      ! just init module data
      !
      IMPLICIT NONE
      INTEGER,                INTENT(IN) :: unit_in
      INTEGER,      OPTIONAL, INTENT(IN) :: unit_out
      CHARACTER(*), OPTIONAL, INTENT(IN) :: dir    
      CHARACTER(*), OPTIONAL, INTENT(IN) :: dir_in, dir_out   
      !
      iunit       = unit_in
      ounit       = unit_in
      IF ( PRESENT( unit_out ) ) ounit  = unit_out
      !
      !
      datadir_in  = "./"
      datadir_out = "./"
      !
      IF ( PRESENT( dir ) ) THEN
          datadir_in  = TRIM(dir)
          datadir_out = TRIM(dir)
      ENDIF
      !
      IF ( PRESENT( dir_in ) ) THEN
          datadir_in  = TRIM(dir_in)
      ENDIF
      !
      IF ( PRESENT( dir_out ) ) THEN
          datadir_out  = TRIM(dir_out)
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
      CHARACTER(*),       INTENT(IN)  :: filename
      CHARACTER(*),       INTENT(IN)  :: action      ! ("read"|"write")
      LOGICAL, OPTIONAL,  INTENT(IN)  :: binary
      INTEGER,            INTENT(OUT) :: ierr
      !
      LOGICAL :: binary_

      ierr = 0
      binary_ = .FALSE.
      IF ( PRESENT(binary) ) binary_ = binary 
      !
      SELECT CASE ( TRIM(action) )
      CASE ( "read", "READ" )
          !
          CALL iotk_open_read ( iunit, FILE = TRIM(filename), IERR=ierr )
          IF ( ierr/=0 ) RETURN
          !
          CALL qexml_read_header( FORMAT_VERSION=qexml_current_version, IERR=ierr )
          IF ( ierr/=0 ) qexml_current_version = TRIM( qexml_default_version )
          !
          !
      CASE ( "write", "WRITE" )
          !
          CALL iotk_open_write( iunit, FILE = TRIM(filename), BINARY=binary_, IERR=ierr )
          IF ( ierr/=0 ) RETURN
          !
          qexml_current_version = TRIM( qexml_default_version )
          !
      CASE DEFAULT
          ierr = 1
      END SELECT
      
      !    
      ! init logical variables for versioning
      !    
      qexml_version_before_1_4_0 = .FALSE.
      !    
      IF ( TRIM( qexml_version_compare( qexml_current_version, "1.4.0" )) == "older" ) &
         qexml_version_before_1_4_0 = .TRUE.
      !
      qexml_current_version_init = .TRUE.
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
      CHARACTER(*),  INTENT(IN)  :: action      ! ("read"|"write")
      INTEGER,       INTENT(OUT) :: ierr
      !
      ierr = 0
      !
      SELECT CASE ( TRIM(action) )
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
      INTEGER, INTENT(IN) :: i
      CHARACTER (LEN=6)   :: int_to_char
      !
      !
      IF ( i < 10 ) THEN
         !
         WRITE( UNIT = int_to_char , FMT = "(I1)" ) i
         !
      ELSE IF ( i < 100 ) THEN
         !
         WRITE( UNIT = int_to_char , FMT = "(I2)" ) i
         !
       ELSE IF ( i < 1000 ) THEN
         !
         WRITE( UNIT = int_to_char , FMT = "(I3)" ) i
         !
       ELSE IF ( i < 10000 ) THEN
         !
         WRITE( UNIT = int_to_char , FMT = "(I4)" ) i
         !
       ELSE
         !
       WRITE( UNIT = int_to_char , FMT = "(I5)" ) i
       !
      END IF
      !
    END FUNCTION int_to_char
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

      length = LEN_TRIM( str )
      !
      IF ( length == 0 ) THEN
         !
         ierr = -1
         RETURN
         !
      ENDIF
  
      i1 = SCAN( str, ".")
      i2 = SCAN( str, ".", BACK=.TRUE.)
      !
      IF ( i1 == 0 .OR. i2 == 0 .OR. i1 == i2 ) THEN
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
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      INTEGER, EXTERNAL          :: c_mkdir
      INTEGER  :: iunaux
      !
      !
      ierr = 0
      CALL iotk_free_unit( iunaux )
      !
      ierr = c_mkdir( TRIM( dirname ), LEN_TRIM( dirname ) )
      IF ( ierr/=0 ) RETURN

      !
      ! ... check whether the scratch directory is writable
      !
      OPEN( iunaux , FILE = TRIM( dirname ) // '/test', IOSTAT = ierr )
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
      CHARACTER(LEN=256)           :: kpoint_dirname
      CHARACTER(LEN=*), INTENT(IN) :: basedir
      INTEGER,          INTENT(IN) :: ik
      !
      CHARACTER(LEN=256) :: kdirname
      CHARACTER(LEN=5)   :: kindex
      !
      WRITE( kindex, FMT = '( I5.5 )' ) ik     
      !
      kdirname = TRIM( basedir ) // '/K' // kindex
      !
      kpoint_dirname = TRIM( kdirname )
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
      CHARACTER(LEN=256)                 :: wfc_filename
      CHARACTER(LEN=*),       INTENT(IN) :: basedir
      CHARACTER(LEN=*),       INTENT(IN) :: name
      INTEGER,                INTENT(IN) :: ik
      INTEGER,      OPTIONAL, INTENT(IN) :: ipol
      CHARACTER(*), OPTIONAL, INTENT(IN) :: tag
      CHARACTER(*), OPTIONAL, INTENT(IN) :: extension
      !    
      CHARACTER(LEN=256) :: filename, tag_, ext_
      !
      !
      filename = ''
      tag_     = ''
      ext_     = '.dat'
      !
      IF ( PRESENT( tag ) )         tag_ = '_'//TRIM(tag)
      IF ( PRESENT( extension ) )   ext_ = '.'//TRIM(extension)
      !
      IF ( PRESENT( ipol ) ) THEN
         !      
         WRITE( filename, FMT = '( I1 )' ) ipol
         !
      END IF
      !
      filename = TRIM( kpoint_dirname( basedir, ik ) ) // '/' // &
                 & TRIM( name ) // TRIM( filename ) // TRIM( tag_ ) // TRIM( ext_)
      !
      wfc_filename = TRIM( filename )
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
      CHARACTER(LEN=*),  INTENT(IN) :: file_in, file_out
      INTEGER,           INTENT(OUT):: ierr
      !
      CHARACTER(LEN=256) :: string
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
         IF ( ios < 0 ) EXIT copy_loop
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
      CHARACTER(LEN=*) :: filename
      !    
      LOGICAL :: lexists
      !    
      INQUIRE( FILE = TRIM( filename ), EXIST = lexists )
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
      CHARACTER(LEN=256)           :: restart_dirname
      CHARACTER(LEN=*), INTENT(IN) :: outdir, prefix
      !
      CHARACTER(LEN=256)         :: dirname
      INTEGER                    :: strlen
      !
      ! ... main restart directory
      !
      dirname = TRIM( prefix ) // '.save'
      !
      IF ( LEN( outdir ) > 1 ) THEN
         !
         strlen = LEN_TRIM( outdir )
         IF ( outdir(strlen:strlen) == '/' ) strlen = strlen -1
         !
         dirname = outdir(1:strlen) // '/' // dirname
         !
      END IF
      !
      restart_dirname = TRIM( dirname )
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
      CHARACTER(LEN=*), INTENT(IN) :: creator_name, creator_version
      !
      CALL iotk_write_begin( ounit, "HEADER" )
      !
      CALL iotk_write_attr(attr, "NAME",TRIM(fmt_name), FIRST=.TRUE.)
      CALL iotk_write_attr(attr, "VERSION",TRIM(fmt_version) )
      CALL iotk_write_empty( ounit, "FORMAT", ATTR=attr )
      !
      CALL iotk_write_attr(attr, "NAME",TRIM(creator_name), FIRST=.TRUE.)
      CALL iotk_write_attr(attr, "VERSION",TRIM(creator_version) )
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
      INTEGER,          INTENT(IN) :: ibravais_latt
      CHARACTER(LEN=*), INTENT(IN) :: symm_type
      REAL(dbl),        INTENT(IN) :: celldm(6), alat
      REAL(dbl),        INTENT(IN) :: a1(3), a2(3), a3(3)
      REAL(dbl),        INTENT(IN) :: b1(3), b2(3), b3(3)
      CHARACTER(LEN=*), INTENT(IN) :: alat_units, a_units, b_units
      !
      CHARACTER(LEN=256) :: bravais_lattice
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
      CALL iotk_write_dat( ounit, &
                           "BRAVAIS_LATTICE", TRIM( bravais_lattice ) )
      !
      CALL iotk_write_dat( ounit, "CELL_SYMMETRY", symm_type )
      !
      CALL iotk_write_attr( attr, "UNITS", TRIM(alat_units), FIRST = .TRUE. )
      CALL iotk_write_dat( ounit, "LATTICE_PARAMETER", alat, ATTR = attr )
      !
      CALL iotk_write_dat( ounit, "CELL_DIMENSIONS", celldm(1:6) )
      !
      CALL iotk_write_attr ( attr,   "UNITS", TRIM(a_units), FIRST = .TRUE. )
      CALL iotk_write_begin( ounit, "DIRECT_LATTICE_VECTORS" )
      CALL iotk_write_empty( ounit, "UNITS_FOR_DIRECT_LATTICE_VECTORS", &
                                     ATTR=attr )
      CALL iotk_write_dat(   ounit, "a1", a1(:) * alat, COLUMNS=3 )
      CALL iotk_write_dat(   ounit, "a2", a2(:) * alat, COLUMNS=3 )
      CALL iotk_write_dat(   ounit, "a3", a3(:) * alat, COLUMNS=3 )
      CALL iotk_write_end(   ounit, "DIRECT_LATTICE_VECTORS" )
      !
      CALL iotk_write_attr ( attr,   "UNITS", TRIM(b_units), FIRST = .TRUE. )
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
      INTEGER,          INTENT(IN) :: nsp, nat
      INTEGER,          INTENT(IN) :: ityp(:)
      CHARACTER(LEN=*), INTENT(IN) :: atm(:)
      CHARACTER(LEN=*), INTENT(IN) :: psfile(:)
      CHARACTER(LEN=*), INTENT(IN) :: pseudo_dir
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      REAL(dbl),        INTENT(IN) :: amass(:)
      CHARACTER(LEN=*), INTENT(IN) :: amass_units
      REAL(dbl),        INTENT(IN) :: tau(:,:)
      CHARACTER(LEN=*), INTENT(IN) :: tau_units
      INTEGER,          INTENT(IN) :: if_pos(:,:)
      !
      INTEGER            :: i, flen, ierrl
      CHARACTER(LEN=256) :: file_pseudo
      LOGICAL            :: pseudo_exists
      !
      !
      CALL iotk_write_begin( ounit, "IONS" )
      !
      CALL iotk_write_dat( ounit, "NUMBER_OF_ATOMS", nat )
      !
      CALL iotk_write_dat( ounit, "NUMBER_OF_SPECIES", nsp )
      !
      flen = LEN_TRIM( pseudo_dir )
      !
      CALL iotk_write_attr ( attr, "UNITS", TRIM(amass_units), FIRST = .TRUE. )
      CALL iotk_write_empty( ounit, "UNITS_FOR_ATOMIC_MASSES", ATTR = attr )
      !
      DO i = 1, nsp
         !
         CALL iotk_write_begin( ounit, "SPECIE"//TRIM(iotk_index(i)) )
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
         END IF
         !
         INQUIRE( FILE = TRIM( dirname ) // "/" &
                       & // TRIM( psfile(i) ), EXIST = pseudo_exists )
         !
         IF ( .NOT. pseudo_exists ) THEN 
               CALL copy_file( TRIM( file_pseudo ), &
                               TRIM( dirname ) // "/" // TRIM( psfile(i) ), ierrl )
         ENDIF
         !
         CALL iotk_write_dat( ounit, "MASS", amass(i) )
         !
         CALL iotk_write_dat( ounit, "PSEUDO", TRIM( psfile(i) ) )
         !
         !
         CALL iotk_write_end( ounit, "SPECIE"//TRIM(iotk_index(i)) )
         !
      ENDDO
      !
      !
      CALL iotk_write_dat( ounit, "PSEUDO_DIR", TRIM( pseudo_dir) )
      !
      CALL iotk_write_attr( attr, "UNITS", TRIM(tau_units), FIRST = .TRUE. )
      CALL iotk_write_empty( ounit, "UNITS_FOR_ATOMIC_POSITIONS", ATTR = attr )
      !
      DO i = 1, nat
         !
         CALL iotk_write_attr( attr, "SPECIES", atm( ityp(i) ), FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "INDEX",  ityp(i) )                     
         CALL iotk_write_attr( attr, "tau",    tau(:,i) )
         CALL iotk_write_attr( attr, "if_pos", if_pos(:,i) )
         CALL iotk_write_empty( ounit, "ATOM" // TRIM( iotk_index( i ) ), attr )
         !
      END DO
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
      INTEGER,          INTENT(IN) :: nsym
      LOGICAL,          INTENT(IN) :: invsym
      INTEGER,          INTENT(IN) :: s(:,:,:) 
      REAL(dbl),        INTENT(IN) :: trasl(:,:)
      CHARACTER(LEN=*), INTENT(IN) :: sname(:)
      CHARACTER(LEN=*), INTENT(IN) :: s_units
      INTEGER,          INTENT(IN) :: t_rev(:)
      INTEGER,          INTENT(IN) :: irt(:,:), nat
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
         tmp(1) = trasl(1,i) 
         tmp(2) = trasl(2,i) 
         tmp(3) = trasl(3,i) 
         !
         CALL iotk_write_dat( ounit, "ROTATION", s(:,:,i), COLUMNS=3 )
         CALL iotk_write_dat( ounit, "FRACTIONAL_TRANSLATION", tmp(1:3), COLUMNS=3 )
         CALL iotk_write_dat( ounit, "EQUIVALENT_IONS", irt(i,1:nat), COLUMNS=8 )
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
    SUBROUTINE qexml_write_efield( tefield, dipfield, edir, emaxpos, eopreg, eamp )
      !------------------------------------------------------------------------
      !
      LOGICAL, INTENT(IN)   :: tefield        ! if .TRUE. a finite electric field
                                              ! is added to the local potential
      LOGICAL, INTENT(IN)   :: dipfield       ! if .TRUE. the dipole field is subtracted
      INTEGER, INTENT(IN)   :: edir           ! direction of the field
      REAL(dbl), INTENT(IN) :: emaxpos        ! position of the maximum of the field (0<emaxpos<1)
      REAL(dbl), INTENT(IN) :: eopreg         ! amplitude of the inverse region (0<eopreg<1)
      REAL(dbl), INTENT(IN) :: eamp           ! field amplitude (in a.u.) (1 a.u. = 51.44 10^11 V/m)
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
      INTEGER,       INTENT(IN) :: npwx, nr1, nr2, nr3, ngm, &
                                   nr1s, nr2s, nr3s, ngms, nr1b, nr2b, nr3b
      INTEGER,       INTENT(IN) :: igv(:,:)
      REAL(dbl),     INTENT(IN) :: ecutwfc, ecutrho
      LOGICAL,       INTENT(IN) :: gamma_only, lgvec
      CHARACTER(*),  INTENT(IN) :: cutoff_units
      !
      !
      CALL iotk_write_begin( ounit, "PLANE_WAVES" )
      !
      CALL iotk_write_attr ( attr, "UNITS", TRIM(cutoff_units), FIRST = .TRUE. )
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
      CALL iotk_write_attr( attr, "nr1", nr1, FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "nr2", nr2 )
      CALL iotk_write_attr( attr, "nr3", nr3 )
      CALL iotk_write_empty( ounit, "FFT_GRID", ATTR = attr )
      !
      CALL iotk_write_dat( ounit, "GVECT_NUMBER", ngm )
      !
      CALL iotk_write_attr( attr, "nr1s", nr1s, FIRST = .TRUE. )
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
                         CREATE = .TRUE., BINARY = .TRUE. )
         !
         CALL iotk_write_begin( ounit, "G-VECTORS" )
         CALL iotk_write_empty( ounit, "INFO", ATTR = attr )
         CALL iotk_write_dat  ( ounit, "g", igv(1:3,1:ngm), COLUMNS = 3 )
         CALL iotk_write_end  ( ounit, "G-VECTORS" )
         !
      END IF
      !
      CALL iotk_write_attr( attr, "nr1b", nr1b , FIRST = .TRUE. )
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
    SUBROUTINE qexml_write_gk( ik, npwk, npwkx, xk, k_units, index, igk )
      !------------------------------------------------------------------------
      !
      INTEGER,      INTENT(IN) :: ik
      INTEGER,      INTENT(IN) :: npwk, npwkx
      REAL(dbl),    INTENT(IN) :: xk(3)
      CHARACTER(*), INTENT(IN) :: k_units
      LOGICAL,      INTENT(IN) :: index(:), igk(:,:)
      !
      INTEGER        :: iunaux
      CHARACTER(256) :: filename

      CALL iotk_free_unit( iunaux )
      filename = wfc_filename( datadir_out, 'gkvectors', ik )
      !
      CALL iotk_open_write( iunaux, FILE = TRIM( filename ), &
                            ROOT="GK-VECTORS", BINARY = .TRUE. )
      !
      CALL iotk_write_dat( iunaux, "NUMBER_OF_GK-VECTORS", npwk )
      CALL iotk_write_dat( iunaux, "MAX_NUMBER_OF_GK-VECTORS", npwkx )
      !
      CALL iotk_write_attr ( attr, "UNITS", TRIM(k_units), FIRST = .TRUE. )
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
      LOGICAL, INTENT(IN) :: lsda, noncolin, lspinorb, domag
      INTEGER, INTENT(IN) :: npol
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
      CHARACTER(LEN=*),    INTENT(IN) :: dft
      LOGICAL,             INTENT(IN) :: lda_plus_u
      INTEGER,   OPTIONAL, INTENT(IN) :: nsp
      INTEGER,   OPTIONAL, INTENT(IN) :: Hubbard_lmax
      INTEGER,   OPTIONAL, INTENT(IN) :: Hubbard_l(:)
      REAL(dbl), OPTIONAL, INTENT(IN) :: Hubbard_U(:), Hubbard_alpha(:)
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
         IF ( .NOT. PRESENT( Hubbard_lmax ) .OR. &
              .NOT. PRESENT( Hubbard_l )    .OR. & 
              .NOT. PRESENT( Hubbard_U )    .OR. &
              .NOT. PRESENT( nsp )          .OR. &
              .NOT. PRESENT( Hubbard_alpha ) ) &
            CALL errore( 'write_exchange_correlation', &
                         ' variables for LDA+U not present', 1 )
         !
         CALL iotk_write_dat( ounit, "NUMBER_OF_SPECIES", nsp )
         !
         CALL iotk_write_dat( ounit, "HUBBARD_LMAX", Hubbard_lmax )
         !
         CALL iotk_write_dat( ounit, "HUBBARD_L", &
                              Hubbard_l(1:Hubbard_lmax) )
         !
         CALL iotk_write_dat( ounit, "HUBBARD_U", Hubbard_U(1:nsp) )
         !
         CALL iotk_write_dat( ounit, "HUBBARD_ALPHA", Hubbard_alpha(1:nsp) )
         !
      END IF
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
      LOGICAL,                INTENT(IN) :: lgauss, ltetra, tfixed_occ, lsda
      INTEGER,      OPTIONAL, INTENT(IN) :: ngauss, ntetra, nstates_up, nstates_dw
      INTEGER,      OPTIONAL, INTENT(IN) :: tetra(:,:)
      REAL(dbl),    OPTIONAL, INTENT(IN) :: degauss, input_occ(:,:)      
      CHARACTER(*), OPTIONAL, INTENT(IN) :: degauss_units
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
         CALL iotk_write_attr( attr, "UNITS", TRIM(degauss_units), FIRST = .TRUE. )
         !
         CALL iotk_write_dat( ounit, "SMEARING_PARAMETER", degauss , ATTR = attr )
         !
      END IF
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
         END DO
         !
      END IF
      !
      CALL iotk_write_dat( ounit, "FIXED_OCCUPATIONS", tfixed_occ )
      !
      IF ( tfixed_occ ) THEN
         !
         CALL iotk_write_attr( attr, "lsda" , lsda, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "nstates_up", nstates_up )
         CALL iotk_write_attr( attr, "nstates_dw", nstates_dw )
         !
         CALL iotk_write_empty( ounit, 'INFO', ATTR = attr )
         !
         CALL iotk_write_dat( ounit, "INPUT_OCC_UP", input_occ(1:nstates_up,1) )
         !
         IF ( lsda ) &
            CALL iotk_write_dat( ounit, "INPUT_OCC_DOWN", input_occ(1:nstates_dw,2) )
         !
      END IF
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
      INTEGER,      INTENT(IN) :: num_k_points, k1, k2, k3, nk1, nk2, nk3
      REAL(dbl),    INTENT(IN) :: xk(:,:), wk(:)
      CHARACTER(*), INTENT(IN) :: k_units
      !
      INTEGER :: ik
      !
      !
      CALL iotk_write_begin( ounit, "BRILLOUIN_ZONE" )
      !
      CALL iotk_write_dat( ounit, "NUMBER_OF_K-POINTS", num_k_points )
      !
      CALL iotk_write_attr( attr, "UNITS", TRIM(k_units), FIRST = .TRUE. )
      CALL iotk_write_empty( ounit, "UNITS_FOR_K-POINTS", attr )
      !
      CALL iotk_write_attr( attr, "nk1", nk1, FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "nk2", nk2 )
      CALL iotk_write_attr( attr, "nk3", nk3 )
      CALL iotk_write_empty( ounit, "MONKHORST_PACK_GRID", attr )
      CALL iotk_write_attr( attr, "k1", k1, FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "k2", k2 )
      CALL iotk_write_attr( attr, "k3", k3 )
      CALL iotk_write_empty( ounit, "MONKHORST_PACK_OFFSET", attr )
      !
      DO ik = 1, num_k_points
         !
         CALL iotk_write_attr( attr, "XYZ", xk(:,ik), FIRST = .TRUE. )
         !            
         CALL iotk_write_attr( attr, "WEIGHT", wk(ik) )
         !
         CALL iotk_write_empty( ounit, "K-POINT" // &
                              & TRIM( iotk_index(ik) ), attr )
         !
      END DO
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
      INTEGER,      INTENT(IN) :: modenum
      REAL(dbl),    INTENT(IN) :: xqq(:)
      CHARACTER(*), INTENT(IN) :: q_units
      !
      !
      CALL iotk_write_begin( ounit, "PHONON" )
      !
      CALL iotk_write_dat( ounit, "NUMBER_OF_MODES", modenum )
      !
      CALL iotk_write_attr( attr, "UNITS", TRIM(q_units), FIRST = .TRUE. )
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
      INTEGER,       INTENT(IN) :: nbnd, num_k_points, nspin, natomwfc
      LOGICAL,       INTENT(IN) :: noncolin
      REAL(dbl),     INTENT(IN) :: ef, nelec
      CHARACTER(*),  INTENT(IN) :: energy_units, k_units
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
      CALL iotk_write_attr ( attr, "UNITS", TRIM(k_units), FIRST = .TRUE. )
      CALL iotk_write_empty( ounit, "UNITS_FOR_K-POINTS", ATTR = attr )
      !
      CALL iotk_write_attr ( attr, "UNITS", TRIM(energy_units), FIRST = .TRUE. )
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
      INTEGER,             INTENT(IN) :: ik, nbnd
      INTEGER, OPTIONAL,   INTENT(IN) :: ispin
      REAL(dbl),           INTENT(IN) :: eig(:)
      CHARACTER(*),        INTENT(IN) :: energy_units
      REAL(dbl), OPTIONAL, INTENT(IN) :: occ(:), ef
      !
      INTEGER :: iunaux
      CHARACTER(LEN=256) :: filename

      !
      IF ( PRESENT( ispin) ) THEN
         !
         filename= TRIM( wfc_filename( datadir_out, 'eigenval', &
                                       ik, ispin, EXTENSION="xml") )
         !
      ELSE
         !
         filename= TRIM( wfc_filename( datadir_out, 'eigenval', &
                                       ik, EXTENSION="xml") )
         !
      ENDIF
      !
      CALL iotk_free_unit( iunaux )
      CALL iotk_open_write ( iunaux, FILE = TRIM( filename ), BINARY = .FALSE. )
      !
      CALL iotk_write_attr ( attr, "nbnd", nbnd, FIRST=.TRUE. )
      CALL iotk_write_attr ( attr, "ik", ik )
      !
      IF ( PRESENT( ispin) ) CALL iotk_write_attr ( attr, "ispin", ispin )
      !
      CALL iotk_write_empty( iunaux, "INFO", ATTR = attr )
      !
      CALL iotk_write_attr ( attr, "UNITS", TRIM(energy_units), FIRST = .TRUE. )
      CALL iotk_write_empty( iunaux, "UNITS_FOR_ENERGIES", ATTR=attr)
      !
      IF ( PRESENT( ef ) ) THEN
         !
         CALL iotk_write_dat( iunaux, "FERMI_ENERGY", ef)
         !
      ENDIF
      !
      CALL iotk_write_dat( iunaux, "EIGENVALUES", eig(:) )
      !
      IF ( PRESENT( occ ) ) THEN
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
    SUBROUTINE qexml_write_wfc( nbnd, nkpts, nspin, ik, ispin, ipol, igk, ngw, &
                                igwx, wf, wf_kindip, scale_factor )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER,                INTENT(IN) :: nbnd, nkpts, nspin
      INTEGER,                INTENT(IN) :: ik
      INTEGER,      OPTIONAL, INTENT(IN) :: ispin, ipol
      INTEGER,                INTENT(IN) :: ngw, igwx
      INTEGER,      OPTIONAL, INTENT(IN) :: igk(:)
      COMPLEX(dbl), OPTIONAL, INTENT(IN) :: wf(:,:)
      COMPLEX(dbl), OPTIONAL, INTENT(IN) :: wf_kindip(:,:)
      REAL(dbl),    OPTIONAL, INTENT(IN) :: scale_factor
      !
      INTEGER         :: iunaux, ierr
      INTEGER         :: ig, ib
      CHARACTER(256)  :: filename
      COMPLEX(dbl),  ALLOCATABLE :: wtmp(:)

      ierr = 0
      !
      IF ( PRESENT( ispin ) .AND. PRESENT( ipol )  ) THEN
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
      IF ( PRESENT( ispin ) ) THEN
         !
         filename = TRIM( wfc_filename( datadir_out, 'evc', ik, ispin ) )
         !
      ELSEIF ( PRESENT( ipol )  ) THEN
         !
         filename = TRIM( wfc_filename( datadir_out, 'evc', ik, ipol ) )
         !
      ELSE
         !
         filename = TRIM( wfc_filename( datadir_out, 'evc', ik ) )
         !
      ENDIF
      !
      CALL iotk_open_write ( iunaux, FILE = TRIM(filename), ROOT="WFC", BINARY=.TRUE., IERR=ierr )
      IF (ierr/=0)  RETURN
      !
      !
      CALL iotk_write_attr( attr, "ngw",          ngw, FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "igwx",         igwx )
      CALL iotk_write_attr( attr, "nbnd",         nbnd )
      CALL iotk_write_attr( attr, "ik",           ik )
      CALL iotk_write_attr( attr, "nk",           nkpts )
      CALL iotk_write_attr( attr, "ispin",        ispin )
      CALL iotk_write_attr( attr, "nspin",        nspin )
      IF ( PRESENT( scale_factor) ) CALL iotk_write_attr( attr, "scale_factor", scale_factor )
      !
      CALL iotk_write_empty( iunaux, "INFO", attr )
      !
      !
      IF ( PRESENT( wf ) ) THEN
         !
         ! write wfcs without any G-reordering
         !
         DO ib = 1, nbnd
            !
            CALL iotk_write_dat( iunaux, "evc" // TRIM(iotk_index( ib )), wf( 1: ngw, ib) )
            !
         ENDDO
         !
      ENDIF
      !
      !
      IF ( PRESENT( wf_kindip ) ) THEN
         !
         ! we need to reorder wfcs in terms of G-vectors
         ! we need the igk map
         !
         IF ( .NOT. PRESENT( igk ) ) THEN
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
            CALL iotk_write_dat( iunaux, "evc" // TRIM(iotk_index( ib )), wtmp( 1: ngw) )
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
      INTEGER,             INTENT(IN) :: nr1, nr2, nr3
      INTEGER,   OPTIONAL, INTENT(IN) :: nr1x, nr2x
      REAL(dbl), OPTIONAL, INTENT(IN) :: rho(:,:,:), rhov(:)
      LOGICAL,   OPTIONAL, INTENT(IN) :: binary
      !
      INTEGER        :: iunaux, nr1x_, nr2x_, ip, i1, i2, i
      LOGICAL        :: binary_
      CHARACTER(256) :: filename
      REAL(dbl), ALLOCATABLE :: plane(:,:)
      !
      !
      CALL iotk_free_unit( iunaux )
      !
      binary_ = .TRUE.
      IF ( PRESENT (binary) ) binary_ = binary
      !
      IF ( binary_ ) THEN
         !
         filename = TRIM( datadir_out ) // '/' //'charge-density.dat'
         !
      ELSE
         !
         filename = TRIM( datadir_out ) // '/' //'charge-density.xml'
         !
      ENDIF
      !
      CALL iotk_open_write( iunaux, FILE = TRIM(filename), BINARY=binary_ )
      !
      !
      CALL iotk_write_begin( iunaux, "CHARGE-DENSITY" )
      !
      CALL iotk_write_attr( attr, "nr1", nr1, FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "nr2", nr2 )
      CALL iotk_write_attr( attr, "nr3", nr3 )
      !
      CALL iotk_write_empty( iunaux, "INFO", attr )
      !
      !
      IF ( PRESENT( rho ) ) THEN
         !
         DO ip = 1, nr3 
            !
            CALL iotk_write_dat( iunaux, "z"//TRIM(iotk_index(ip)), rho(1:nr1,1:nr2,ip) )
            !
         ENDDO
         !
      ELSEIF ( PRESENT( rhov ) ) THEN
         !
         nr1x_ = nr1
         IF ( PRESENT( nr1x )) nr1x_ = nr1x
         nr2x_ = nr2
         IF ( PRESENT( nr2x )) nr2x_ = nr2x
         !
         IF ( nr1x_ /= nr1 .OR. nr2x_ /= nr2 ) THEN
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
               CALL iotk_write_dat( iunaux, "z"//TRIM(iotk_index(ip)), plane )
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
               CALL iotk_write_dat( iunaux, "z"//TRIM(iotk_index(ip)), rhov(i1:i2) )
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
      CHARACTER(LEN=*),  OPTIONAL, INTENT(OUT) :: creator_name, creator_version
      CHARACTER(LEN=*),  OPTIONAL, INTENT(OUT) :: format_name, format_version
      INTEGER,           OPTIONAL, INTENT(OUT) :: ierr

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
      IF ( PRESENT(creator_name) )     creator_name    = TRIM(creator_name_)
      IF ( PRESENT(creator_version) )  creator_version = TRIM(creator_version_)
      IF ( PRESENT(format_name) )      format_name     = TRIM(format_name_)
      IF ( PRESENT(format_version) )   format_version  = TRIM(format_version_)
      !
    END SUBROUTINE qexml_read_header
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_cell( bravais_latt, symm_type, celldm, alat, &
                                a1, a2, a3, b1, b2, b3, alat_units, a_units, b_units, ierr )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=*),  OPTIONAL, INTENT(OUT) :: bravais_latt
      CHARACTER(LEN=*),  OPTIONAL, INTENT(OUT) :: symm_type
      REAL(dbl),         OPTIONAL, INTENT(OUT) :: celldm(6), alat
      REAL(dbl),         OPTIONAL, INTENT(OUT) :: a1(3), a2(3), a3(3)
      REAL(dbl),         OPTIONAL, INTENT(OUT) :: b1(3), b2(3), b3(3)
      CHARACTER(LEN=*),  OPTIONAL, INTENT(OUT) :: alat_units, a_units, b_units
      INTEGER,                     INTENT(OUT) :: ierr
      !
      CHARACTER(256)     :: bravais_latt_, symm_type_
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
      CALL iotk_scan_dat( iunit, "CELL_SYMMETRY", symm_type_, IERR=ierr )
      IF ( ierr /= 0 ) RETURN
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
      IF ( PRESENT(bravais_latt) )  bravais_latt = bravais_latt_
      IF ( PRESENT(celldm) )        symm_type    = symm_type_
      IF ( PRESENT(symm_type) )     celldm       = celldm_
      IF ( PRESENT(alat) )          alat         = alat_
      IF ( PRESENT(a1) )            a1           = a1_
      IF ( PRESENT(a2) )            a2           = a2_
      IF ( PRESENT(a3) )            a3           = a3_
      IF ( PRESENT(b1) )            b1           = b1_
      IF ( PRESENT(b2) )            b2           = b2_
      IF ( PRESENT(b3) )            b3           = b3_
      IF ( PRESENT(alat_units) )    alat_units   = TRIM(alat_units_)
      IF ( PRESENT(a_units) )       a_units      = TRIM(a_units_)
      IF ( PRESENT(b_units) )       b_units      = TRIM(b_units_)

    END SUBROUTINE qexml_read_cell
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_ions( nsp, nat, atm, ityp, psfile, amass, amass_units, &
                                tau, tau_units, if_pos, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,          OPTIONAL, INTENT(OUT) :: nsp, nat
      INTEGER,          OPTIONAL, INTENT(OUT) :: ityp(:)
      CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: atm(:)
      CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: psfile(:)
      REAL(dbl),        OPTIONAL, INTENT(OUT) :: amass(:)
      CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: amass_units
      REAL(dbl),        OPTIONAL, INTENT(OUT) :: tau(:,:)
      INTEGER,          OPTIONAL, INTENT(OUT) :: if_pos(:,:)
      CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: tau_units
      INTEGER,                    INTENT(OUT) :: ierr
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
      IF ( PRESENT(nat) )   nat = nat_
      IF ( PRESENT(nsp) )   nsp = nsp_
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
            CALL iotk_scan_dat( iunit, TRIM( atm_(i) ) // "_MASS", amass_(i), IERR=ierr )
            IF (ierr/=0) RETURN
            CALL iotk_scan_dat( iunit, "PSEUDO_FOR_" // TRIM( atm_(i) ), &
                                psfile_(i), IERR=ierr )
            IF (ierr/=0) RETURN
            !
         ELSE
            !
            ! current version
            !
            CALL iotk_scan_begin( iunit, "SPECIE"//TRIM(iotk_index(i)), IERR=ierr )
            IF (ierr/=0) RETURN
            !
            CALL iotk_scan_dat( iunit, "ATOM_TYPE", atm_(i), IERR=ierr )
            IF (ierr/=0) RETURN
            CALL iotk_scan_dat( iunit, "MASS", amass_(i), IERR=ierr )
            IF (ierr/=0) RETURN
            CALL iotk_scan_dat( iunit, "PSEUDO", psfile_(i), IERR=ierr )
            IF (ierr/=0) RETURN
            !
            CALL iotk_scan_end( iunit, "SPECIE"//TRIM(iotk_index(i)), IERR=ierr )
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
                               "ATOM" // TRIM( iotk_index(i) ), ATTR=attr, IERR=ierr )
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
      IF ( PRESENT(nsp) )         nsp    = nsp_
      IF ( PRESENT(nat) )         nat    = nat_
      IF ( PRESENT(atm) )         atm(1:nsp_)    = atm_
      IF ( PRESENT(amass) )       amass(1:nsp_)  = amass_
      IF ( PRESENT(amass_units) ) amass_units    = TRIM(amass_units_)
      IF ( PRESENT(psfile) )      psfile(1:nsp_) = psfile_(1:nsp_)
      IF ( PRESENT(ityp) )        ityp(1:nat_)   = ityp_
      IF ( PRESENT(tau_units) )   tau_units      = TRIM(tau_units_)
      IF ( PRESENT(tau) )         tau(1:3, 1:nat_)    = tau_
      IF ( PRESENT(if_pos) )      if_pos(1:3, 1:nat_) = if_pos_
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
      INTEGER,          OPTIONAL, INTENT(OUT) :: nsym
      LOGICAL,          OPTIONAL, INTENT(OUT) :: invsym
      INTEGER,          OPTIONAL, INTENT(OUT) :: s(:,:,:)
      REAL(dbl),        OPTIONAL, INTENT(OUT) :: trasl(:,:)
      CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: sname(:)
      CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: s_units
      INTEGER,          OPTIONAL, INTENT(OUT) :: t_rev(:)
      INTEGER,          OPTIONAL, INTENT(OUT) :: irt(:,:), nat
      INTEGER,                    INTENT(OUT) :: ierr
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
          CALL iotk_scan_begin( iunit, "SYMM"//TRIM( iotk_index( i ) ), IERR=ierr )
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
          CALL iotk_scan_end( iunit, "SYMM"//TRIM( iotk_index( i ) ), IERR=ierr )
          IF (ierr/=0) RETURN
          !
      ENDDO
      !
      CALL iotk_scan_end( iunit, "SYMMETRIES", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( PRESENT(nsym) )        nsym          = nsym_
      IF ( PRESENT(invsym) )      invsym        = invsym_
      IF ( PRESENT(nat) )         nat           = nat_
      IF ( PRESENT(trasl) )       trasl(1:3, 1:nsym_)   = trasl_(1:3, 1:nsym_)
      IF ( PRESENT(s) )           s(1:3, 1:3, 1:nsym_)  = s_(1:3, 1:3, 1:nsym_)
      IF ( PRESENT(irt) )         irt(1:nsym_, 1:nat_)  = irt_(1:nsym_, 1:nat_)
      IF ( PRESENT(sname) )  THEN     
          DO i = 1, nsym_
                                  sname( i )            = TRIM( sname_( i ) )
          ENDDO
      ENDIF       
      IF ( PRESENT(s_units) )     s_units               = TRIM( s_units_ )
      IF ( PRESENT(t_rev) )       t_rev( 1:nsym_ )      = t_rev_( 1:nsym_ )
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
      LOGICAL,   OPTIONAL, INTENT(OUT) :: tefield, dipfield    
      INTEGER,   OPTIONAL, INTENT(OUT) :: edir
      REAL(dbl), OPTIONAL, INTENT(OUT) :: emaxpos, eopreg, eamp
      INTEGER,             INTENT(OUT) :: ierr
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
      IF ( PRESENT(tefield) )        tefield      = tefield_
      IF ( PRESENT(dipfield) )       dipfield     = dipfield_
      IF ( PRESENT(edir) )           edir         = edir_
      IF ( PRESENT(emaxpos) )        emaxpos      = emaxpos_
      IF ( PRESENT(eopreg) )         eopreg       = eopreg_
      IF ( PRESENT(eamp) )           eamp         = eamp_
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
      INTEGER,      OPTIONAL, INTENT(OUT) :: npwx, nr1, nr2, nr3, ngm, &
                                             nr1s, nr2s, nr3s, ngms, nr1b, nr2b, nr3b
      INTEGER,      OPTIONAL, INTENT(OUT) :: igv(:,:)
      REAL(dbl),    OPTIONAL, INTENT(OUT) :: ecutwfc, ecutrho
      LOGICAL,      OPTIONAL, INTENT(OUT) :: gamma_only
      CHARACTER(*), OPTIONAL, INTENT(OUT) :: cutoff_units
      INTEGER,                INTENT(OUT) :: ierr
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
      IF ( PRESENT( igv ) ) THEN
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
      IF ( PRESENT( ecutwfc ) )           ecutwfc      = ecutwfc_
      IF ( PRESENT( ecutrho ) )           ecutrho      = ecutrho_
      IF ( PRESENT( npwx ) )              npwx         = npwx_
      IF ( PRESENT( gamma_only ) )        gamma_only   = gamma_only_
      IF ( PRESENT( nr1 ) )               nr1          = nr1_
      IF ( PRESENT( nr2 ) )               nr2          = nr2_
      IF ( PRESENT( nr3 ) )               nr3          = nr3_
      IF ( PRESENT( ngm ) )               ngm          = ngm_
      IF ( PRESENT( nr1s ) )              nr1s         = nr1s_
      IF ( PRESENT( nr2s ) )              nr2s         = nr2s_
      IF ( PRESENT( nr3s ) )              nr3s         = nr3s_
      IF ( PRESENT( ngms ) )              ngms         = ngms_
      IF ( PRESENT( nr1b ) )              nr1b         = nr1b_
      IF ( PRESENT( nr2b ) )              nr2b         = nr2b_
      IF ( PRESENT( nr3b ) )              nr3b         = nr3b_
      IF ( PRESENT( cutoff_units ) )      cutoff_units = TRIM( cutoff_units_ )
      !
    END SUBROUTINE qexml_read_planewaves
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_gk( ik, npwk, npwkx, xk, k_units, index, igk, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,                INTENT(IN)  :: ik
      INTEGER,      OPTIONAL, INTENT(OUT) :: npwk, npwkx
      REAL(dbl),    OPTIONAL, INTENT(OUT) :: xk(3)
      CHARACTER(*), OPTIONAL, INTENT(OUT) :: k_units
      INTEGER,      OPTIONAL, INTENT(OUT) :: igk(:,:), index(:)
      INTEGER,                INTENT(OUT) :: ierr
      !
      CHARACTER(256) :: filename, k_units_
      INTEGER   :: npwk_, npwkx_
      REAL(dbl) :: xk_(3)
      INTEGER   :: iunaux
      !

      ierr = 0
      !
      CALL iotk_free_unit( iunaux )
      filename = wfc_filename( datadir_in, 'gkvectors', ik )
      !
      CALL iotk_open_read ( iunaux, FILE = TRIM(filename), IERR=ierr )
      IF (ierr/=0)  RETURN
      !
      CALL iotk_scan_dat( iunaux, 'NUMBER_OF_GK-VECTORS', npwk_, IERR=ierr)
      IF (ierr/=0)  RETURN
      !
      CALL iotk_scan_dat( iunaux, 'MAX_NUMBER_OF_GK-VECTORS', npwkx_, IERR=ierr)
      IF (ierr/=0)  RETURN
      !
      CALL iotk_scan_dat( iunaux, 'K-POINT_COORDS', xk_, ATTR=attr, IERR=ierr)
      IF (ierr/=0)  RETURN
      CALL iotk_scan_attr( attr, 'UNITS', k_units_, IERR=ierr)
      IF (ierr/=0)  RETURN
      !
      IF ( PRESENT( index ) ) THEN
          !
          CALL iotk_scan_dat( iunaux, 'INDEX', index(1:npwk_), IERR=ierr)
          IF (ierr/=0)  RETURN
          !
      ENDIF
      !
      IF ( PRESENT( igk ) ) THEN
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
      IF ( PRESENT( npwk ) )       npwk    = npwk_
      IF ( PRESENT( npwkx ) )      npwkx   = npwkx_
      IF ( PRESENT( xk ) )         xk(1:3) = xk_(1:3)
      IF ( PRESENT( k_units ) )    k_units = TRIM(k_units_)
      !
    END SUBROUTINE qexml_read_gk
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_spin( lsda, noncolin, npol, lspinorb, ierr )
      !------------------------------------------------------------------------
      !
      LOGICAL, OPTIONAL, INTENT(OUT) :: lsda, noncolin, lspinorb
      INTEGER, OPTIONAL, INTENT(OUT) :: npol
      INTEGER,           INTENT(OUT) :: ierr
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
      IF ( PRESENT( lsda ) )       lsda      = lsda_
      IF ( PRESENT( noncolin ) )   noncolin  = noncolin_
      IF ( PRESENT( npol ) )       npol      = npol_
      IF ( PRESENT( lspinorb ) )   lspinorb  = lspinorb_
      !

    END SUBROUTINE qexml_read_spin
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_xc( dft, lda_plus_u, Hubbard_lmax, Hubbard_l, &
                              nsp, Hubbard_U, Hubbard_alpha, ierr )
    !------------------------------------------------------------------------
      !
      CHARACTER(LEN=*), OPTIONAL, INTENT(OUT) :: dft
      LOGICAL,          OPTIONAL, INTENT(OUT) :: lda_plus_u
      INTEGER,          OPTIONAL, INTENT(OUT) :: Hubbard_lmax
      INTEGER,          OPTIONAL, INTENT(OUT) :: Hubbard_l(:)
      INTEGER,          OPTIONAL, INTENT(OUT) :: nsp
      REAL(dbl),        OPTIONAL, INTENT(OUT) :: Hubbard_U(:), Hubbard_alpha(:)
      INTEGER,                    INTENT(OUT) :: ierr
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
         ALLOCATE( Hubbard_l_(1:Hubbard_lmax_) )
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
      IF ( PRESENT( dft ) )           dft           = dft_
      IF ( PRESENT( lda_plus_u ) )    lda_plus_u    = lda_plus_u_
      !
      IF ( lda_plus_u_ )  THEN
         !
         IF ( PRESENT( nsp ) )             nsp                   = nsp_
         IF ( PRESENT( Hubbard_lmax ) )    Hubbard_lmax          = Hubbard_lmax_
         IF ( PRESENT( Hubbard_l ) )       Hubbard_l(1:Hubbard_lmax_)   = Hubbard_l_(:)
         IF ( PRESENT( Hubbard_U ) )       Hubbard_U(1:nsp_)     = Hubbard_U_(1:nsp_)
         IF ( PRESENT( Hubbard_alpha ) )   Hubbard_alpha(1:nsp_) = Hubbard_alpha_(1:nsp_)
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
      LOGICAL,      OPTIONAL, INTENT(OUT) :: lgauss, ltetra, tfixed_occ
      INTEGER,      OPTIONAL, INTENT(OUT) :: ngauss, ntetra
      INTEGER,      OPTIONAL, INTENT(OUT) :: tetra(:,:)
      INTEGER,      OPTIONAL, INTENT(OUT) :: nstates_up, nstates_dw
      REAL(dbl),    OPTIONAL, INTENT(OUT) :: degauss, input_occ(:,:)
      CHARACTER(*), OPTIONAL, INTENT(OUT) :: degauss_units
      INTEGER,                INTENT(OUT) :: ierr
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
      IF ( tfixed_occ_  .AND. ( PRESENT( input_occ ) .OR. &
                                PRESENT(nstates_up)  .OR. PRESENT(nstates_dw) ) ) THEN
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
            CALL iotk_scan_attr( attr, "nstates_dw", nstates_dw_, IERR=ierr )
            IF (ierr /=0 ) RETURN
            !
         ENDIF
         ! 
         IF ( PRESENT( input_occ ) ) THEN
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
      IF ( PRESENT( lgauss ))           lgauss      = lgauss_
      IF ( PRESENT( ltetra ))           ltetra      = ltetra_
      IF ( PRESENT( tfixed_occ ))       tfixed_occ  = tfixed_occ_
      IF ( PRESENT( ngauss ))           ngauss      = ngauss_
      IF ( PRESENT( ntetra ))           ntetra      = ntetra_
      IF ( PRESENT( degauss ))          degauss     = degauss
      IF ( PRESENT( degauss_units ))    degauss_units  = TRIM(degauss_units_)
      IF ( PRESENT( nstates_up ))       nstates_up  = nstates_up_
      IF ( PRESENT( nstates_dw ))       nstates_dw  = nstates_dw_
      !
      IF ( ltetra_ ) THEN
         !
         IF ( PRESENT( tetra ) )         tetra(1:4, 1:ntetra_)  = tetra_
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
      INTEGER,       OPTIONAL, INTENT(OUT) :: num_k_points, k1, k2, k3, nk1, nk2, nk3
      REAL(dbl),     OPTIONAL, INTENT(OUT) :: xk(:,:), wk(:)
      CHARACTER(*),  OPTIONAL, INTENT(OUT) :: k_units
      INTEGER,                 INTENT(OUT) :: ierr
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
         CALL iotk_scan_empty( iunit, "K-POINT" // TRIM( iotk_index(ik) ), &
                               ATTR=attr, IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
         CALL iotk_scan_attr( attr, "XYZ", xk_(:,ik), IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !            
         CALL iotk_scan_attr( attr, "WEIGHT", wk_(ik), IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
      END DO
      !
      CALL iotk_scan_end( iunit, "BRILLOUIN_ZONE", IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      !
      IF ( PRESENT( num_k_points ) )       num_k_points  = num_k_points_
      IF ( PRESENT( nk1 ) )                nk1           = nk1_
      IF ( PRESENT( nk2 ) )                nk2           = nk2_
      IF ( PRESENT( nk3 ) )                nk3           = nk3_
      IF ( PRESENT( k1 ) )                 k1            =  k1_
      IF ( PRESENT( k2 ) )                 k2            =  k2_
      IF ( PRESENT( k3 ) )                 k3            =  k3_
      IF ( PRESENT( k_units ) )            k_units       =  TRIM(k_units_)
      IF ( PRESENT( xk ) )                 xk(1:3,1:num_k_points_) = xk_(:,:)
      IF ( PRESENT( wk ) )                 wk(1:num_k_points_)     = wk_(:)
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
      INTEGER,       OPTIONAL, INTENT(OUT) :: modenum
      REAL(dbl),     OPTIONAL, INTENT(OUT) :: xqq(:)
      CHARACTER(*),  OPTIONAL, INTENT(OUT) :: q_units
      INTEGER,                 INTENT(OUT) :: ierr
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
      IF ( PRESENT (xqq) ) THEN
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
      IF ( PRESENT (modenum) )      modenum = modenum_
      IF ( PRESENT (q_units) )      q_units = TRIM(q_units_)
      !
    END SUBROUTINE qexml_read_phonon
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_bands_info( nbnd, num_k_points, nspin, noncolin, natomwfc, & 
                                      nelec, ef, energy_units, k_units, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,      OPTIONAL, INTENT(OUT) :: nbnd, num_k_points, nspin, natomwfc
      LOGICAL,      OPTIONAL, INTENT(OUT) :: noncolin
      REAL(dbl),    OPTIONAL, INTENT(OUT) :: ef, nelec
      CHARACTER(*), OPTIONAL, INTENT(OUT) :: energy_units, k_units
      INTEGER,                INTENT(OUT) :: ierr
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
      IF ( PRESENT( nbnd ) )             nbnd           = nbnd_
      IF ( PRESENT( num_k_points ) )     num_k_points   = num_k_points_
      IF ( PRESENT( nspin ) )            nspin          = nspin_
      IF ( PRESENT( noncolin ) )         noncolin       = noncolin_
      IF ( PRESENT( natomwfc ) )         natomwfc       = natomwfc_
      IF ( PRESENT( nelec ) )            nelec          = nelec_
      IF ( PRESENT( ef ) )               ef             = ef_
      IF ( PRESENT( energy_units ) )     energy_units   = TRIM( energy_units_ )
      IF ( PRESENT( k_units ) )          k_units        = TRIM( k_units_ )
      !
    END SUBROUTINE qexml_read_bands_info
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_bands( ik, ispin, nbnd, eig, energy_units, occ, ef, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,                INTENT(IN)  :: ik
      INTEGER,      OPTIONAL, INTENT(IN)  :: ispin
      INTEGER,      OPTIONAL, INTENT(OUT) :: nbnd
      REAL(dbl),    OPTIONAL, INTENT(OUT) :: eig(:)
      CHARACTER(*), OPTIONAL, INTENT(OUT) :: energy_units
      REAL(dbl),    OPTIONAL, INTENT(OUT) :: occ(:)
      REAL(dbl),    OPTIONAL, INTENT(OUT) :: ef
      INTEGER,                INTENT(OUT) :: ierr
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
      IF ( PRESENT( ispin) ) THEN
         !
         filename= TRIM( wfc_filename( datadir_in, 'eigenval', &
                                       ik, ispin, EXTENSION="xml") ) 
         !
      ELSE
         !
         filename= TRIM( wfc_filename( datadir_in, 'eigenval', &
                                       ik, EXTENSION="xml") ) 
         !
      ENDIF
      !
      !
      CALL iotk_open_read ( iunaux, FILE = TRIM(filename), IERR=ierr )
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
      IF ( PRESENT( ef )) THEN
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
      IF ( PRESENT( nbnd ) )             nbnd             = nbnd_
      IF ( PRESENT( energy_units ) )     energy_units     = TRIM( energy_units_ )
      IF ( PRESENT( occ ) )              occ  (1:nbnd_ )  = occ_(:)
      IF ( PRESENT( eig ) )              eig  (1:nbnd_ )  = eig_(:)
      !
      DEALLOCATE( occ_ )
      DEALLOCATE( eig_ )
      !
    END SUBROUTINE qexml_read_bands
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_wfc( ibnds, ibnde, ik, ispin, ipol, igk, ngw, igwx, &
                               wf, wf_kindip, ierr )
      !------------------------------------------------------------------------
      !
      ! read wfc from IBNDS to IBNDE, for kpt IK and spin ISPIN
      ! WF is the wfc on its proper k+g grid, while WF_KINDIP is the same wfc
      ! but on a truncated rho grid (k-point indipendent)
      !
      INTEGER,                 INTENT(IN)  :: ibnds, ibnde, ik
      INTEGER,       OPTIONAL, INTENT(IN)  :: ispin, ipol
      INTEGER,       OPTIONAL, INTENT(IN)  :: igk(:)
      INTEGER,       OPTIONAL, INTENT(OUT) :: ngw, igwx
      COMPLEX(dbl),  OPTIONAL, INTENT(OUT) :: wf(:,:), wf_kindip(:,:)
      INTEGER,                 INTENT(OUT) :: ierr
      !
      INTEGER :: iunaux
      INTEGER :: ngw_, igwx_, ig, ib, lindex
      COMPLEX(dbl),  ALLOCATABLE :: wf_(:)
      CHARACTER(256)             :: filename

      ierr = 0
      !
      !
      ! few check
      !
      IF ( PRESENT( ispin ) .AND. PRESENT( ipol )  ) THEN
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
      IF ( PRESENT( ispin ) ) THEN
         !
         filename = TRIM( wfc_filename( datadir_in, 'evc', ik, ispin ) ) 
         !
      ELSEIF ( PRESENT( ipol )  ) THEN
         !
         filename = TRIM( wfc_filename( datadir_in, 'evc', ik, ipol ) ) 
         !
      ELSE
         !
         filename = TRIM( wfc_filename( datadir_in, 'evc', ik ) ) 
         !
      ENDIF
      !
      CALL iotk_open_read ( iunaux, FILE = TRIM(filename), IERR=ierr )
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
      IF ( PRESENT( wf )  )  THEN
          !
          lindex = 0
          !
          DO ib = ibnds, ibnde
              !
              lindex = lindex + 1
              !
              CALL iotk_scan_dat( iunaux, "evc"//TRIM(iotk_index(ib)), &
                                  wf( 1:igwx_, lindex ), IERR=ierr )
              IF (ierr/=0) RETURN
              !
          ENDDO
          !
      ENDIF
      !
      IF ( PRESENT( wf_kindip )  )  THEN
          !
          ALLOCATE( wf_(igwx_ ), STAT=ierr )
          IF (ierr/=0) RETURN
          !
          IF ( .NOT. PRESENT( igk ) ) THEN
              ierr = 3
              RETURN
          ENDIF
          !
          IF ( MAXVAL( igk( 1: igwx_ ) ) > SIZE( wf_kindip, 1)  ) THEN
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
              CALL iotk_scan_dat( iunaux, "evc"//TRIM(iotk_index( ib ) ), &
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
      IF ( PRESENT( ngw ) )     ngw    = ngw_
      IF ( PRESENT( igwx ) )    igwx   = igwx_
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
      INTEGER,   OPTIONAL, INTENT(OUT) :: nr1, nr2, nr3
      INTEGER,   OPTIONAL, INTENT(IN)  :: ip
      REAL(dbl), OPTIONAL, INTENT(OUT) :: rho(:,:,:), rhoz(:)
      INTEGER,             INTENT(OUT) :: ierr
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
      filename = TRIM( datadir_in ) // '/' // 'charge-density.dat'
      lexists  = check_file_exst( TRIM(filename) )
      ! 
      IF ( .NOT. lexists ) THEN
          !
          filename = TRIM( datadir_in ) // '/' // 'charge-density.xml'
          lexists  = check_file_exst( TRIM(filename) )
          !
      ENDIF
      !
      IF ( .NOT. lexists ) THEN
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
      IF ( PRESENT( rhoz ) ) THEN
         !
         IF ( .NOT. PRESENT( ip ) ) THEN
            ierr = 71
            RETURN
         ENDIF
         !
         CALL iotk_scan_dat( iunaux, "z"//TRIM(iotk_index(ip)), rhoz, IERR=ierr)
         IF (ierr/=0) RETURN
         !
      ENDIF
      !
      !
      IF ( PRESENT( rho ) ) THEN
         !
         DO ip_ = 1, nr3_
            !
            CALL iotk_scan_dat( iunaux, "z"//TRIM(iotk_index(ip_)), rho(1:nr1_,1:nr2_,ip_), &
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
      IF ( PRESENT( nr1 ) )     nr1 = nr1_
      IF ( PRESENT( nr2 ) )     nr2 = nr2_
      IF ( PRESENT( nr3 ) )     nr3 = nr3_
      !
    END SUBROUTINE qexml_read_rho
    !
    !
END MODULE qexml_module

