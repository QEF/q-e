!
! Copyright (C) 2005 Quantum-ESPRESSO group
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
  ! in XML format the data produced by Quantum-ESPRESSO package
  !
  ! written by Andrea Ferretti (2006)
  ! using large part of implementation by Carlo Sbraccia (2005)
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
  CHARACTER(5), PARAMETER :: fmt_version = "1.2.0"
  !
  ! some default for kinds
  !
  INTEGER,   PARAMETER :: dbl = SELECTED_REAL_KIND( 14, 200 )
  REAL(dbl), PARAMETER :: e2 = 2.0_dbl
  !
  ! internal data to be set
  !
  CHARACTER(256)   :: datadir
  INTEGER          :: iunpun, rhounit
  LOGICAL          :: rho_binary
  !
  ! end of declarations
  !
  PUBLIC :: fmt_name, fmt_version
  PUBLIC :: iunpun, rhounit, rho_binary
  !
  PUBLIC :: qexml_init,  qexml_openfile, qexml_closefile
  !
  PUBLIC :: qexml_write_header, qexml_write_cell, qexml_write_ions,   &
            qexml_write_symmetry, qexml_write_planewaves,             &
            qexml_write_spin, qexml_write_xc,                         &
            qexml_write_occ, qexml_write_bz, qexml_write_phonon,      &
            qexml_write_rho, qexml_write_wfc
  !
  PUBLIC :: qexml_read_header, qexml_read_cell, qexml_read_ions,      &
            qexml_read_symmetry, qexml_read_planewaves,               &
            qexml_read_spin,  qexml_read_xc,                          &
            qexml_read_occ,  qexml_read_bz,  qexml_read_phonon,       &
            qexml_read_bands, qexml_read_gk, qexml_read_wfc
!            qexml_read_rho

  !
  CHARACTER(iotk_attlenx) :: attr
  !
CONTAINS
  !
!
!-------------------------------------------
! ... basic (public) subroutines
!-------------------------------------------
!
    !------------------------------------------------------------------------
    SUBROUTINE qexml_init( iunpun_, datadir_, rhounit_, rho_binary_ )
      !------------------------------------------------------------------------
      !
      ! just init module data
      !
      IMPLICIT NONE
      INTEGER,           INTENT(IN) :: iunpun_
      INTEGER,           INTENT(IN) :: rhounit_
      CHARACTER(*),      INTENT(IN) :: datadir_    
      LOGICAL, OPTIONAL, INTENT(IN) :: rho_binary_
      !
      iunpun  = iunpun_
      rhounit = rhounit_
      datadir = TRIM(datadir_)
      !
      rho_binary = .TRUE.
      IF ( PRESENT(rho_binary_) ) rho_binary = rho_binary_
      !
    END SUBROUTINE qexml_init


    !------------------------------------------------------------------------
    SUBROUTINE qexml_openfile( filename, action, binary, ierr)
      !------------------------------------------------------------------------
      !
      ! open data file
      !
      IMPLICIT NONE
      !
      CHARACTER(*),  INTENT(IN)  :: filename
      CHARACTER(*),  INTENT(IN)  :: action      ! ("read"|"write")
      LOGICAL,       INTENT(IN)  :: binary
      INTEGER,       INTENT(OUT) :: ierr
      !
      ierr = 0
      !
      SELECT CASE ( TRIM(action) )
      CASE ( "read", "READ" )
          !
          CALL iotk_open_read ( iunpun, FILE = TRIM(filename), &
                                BINARY=binary, IERR=ierr )
          !
      CASE ( "write", "WRITE" )
          !
          CALL iotk_open_write( iunpun, FILE = TRIM(filename), &
                                BINARY=binary, IERR=ierr )
          !
      CASE DEFAULT
          ierr = 1
      END SELECT
          
    END SUBROUTINE qexml_openfile
      

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
          CALL iotk_close_read( iunpun, IERR=ierr )
          !
      CASE ( "write", "WRITE" )
          !
          CALL iotk_close_write( iunpun, IERR=ierr )
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
    FUNCTION int_to_char( int )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: int
      CHARACTER (LEN=6)   :: int_to_char
      !
      !
      IF ( int < 10 ) THEN
         !
         WRITE( UNIT = int_to_char , FMT = "(I1)" ) int
         !
      ELSE IF ( int < 100 ) THEN
         !
         WRITE( UNIT = int_to_char , FMT = "(I2)" ) int
         !
       ELSE IF ( int < 1000 ) THEN
         !
         WRITE( UNIT = int_to_char , FMT = "(I3)" ) int
         !
       ELSE IF ( int < 10000 ) THEN
         !
         WRITE( UNIT = int_to_char , FMT = "(I4)" ) int
         !
       ELSE
         !
       WRITE( UNIT = int_to_char , FMT = "(I5)" ) int
       !
      END IF
      !
    END FUNCTION int_to_char

    !------------------------------------------------------------------------
    SUBROUTINE create_directory( dirname )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      !
      INTEGER                    :: ierr
      INTEGER, EXTERNAL          :: c_mkdir
      !
      ierr = c_mkdir( TRIM( dirname ), LEN_TRIM( dirname ) )
      !
      CALL errore( 'create_directory', &
                   'unable to create directory ' // TRIM( dirname ), ierr )
      !
      ! ... check whether the scratch directory is writable
      !
      OPEN( UNIT = 4, FILE = TRIM( dirname ) // '/test', &
            STATUS = 'UNKNOWN', IOSTAT = ierr )
      CLOSE( UNIT = 4, STATUS = 'DELETE' )
      !
      CALL errore( 'create_directory:', &
                   TRIM( dirname ) // ' non existent or non writable', ierr )
      !
      RETURN
      !
    END SUBROUTINE create_directory
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
      !
      WRITE( kindex, FMT = '( I5.5 )' ) ik     
      !
      kdirname = TRIM( basedir ) // '/K' // kindex
      !
      kpoint_dir = TRIM( kdirname )
      !
      RETURN
      !
    END FUNCTION kpoint_dir
    !
    !------------------------------------------------------------------------
    FUNCTION wfc_filename( basedir, name, ik, ipol, tag )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=256)              :: wfc_filename
      CHARACTER(LEN=*),    INTENT(IN) :: basedir
      CHARACTER(LEN=*),    INTENT(IN) :: name
      INTEGER,             INTENT(IN) :: ik
      INTEGER,   OPTIONAL, INTENT(IN) :: ipol
      CHARACTER, OPTIONAL, INTENT(IN) :: tag
      !    
      CHARACTER(LEN=256) :: filename
      !
      !
      filename = ''
      !
      IF ( PRESENT( ipol ) ) THEN
         !      
         WRITE( filename, FMT = '( I1 )' ) ipol
         !
      END IF
      !
      IF ( PRESENT( tag ) ) THEN
         !
         filename = TRIM( kpoint_dir( basedir, ik ) ) // '/' // &
                  & TRIM( name ) // TRIM( filename ) // '_' // tag // '.dat'
         !
      ELSE
         !
         filename = TRIM( kpoint_dir( basedir, ik ) ) // '/' // &
                  & TRIM( name ) // TRIM( filename ) // '.dat'
         !
      END IF
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
!-------------------------------------------
! ... writing subroutines
!-------------------------------------------
!
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_header( creator_name, creator_version )
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
    END SUBROUTINE qexml_write_header
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_cell( ibrav, symm_type, &
                                celldm, alat, a1, a2, a3, b1, b2, b3 )
      !------------------------------------------------------------------------
      !
      INTEGER,           OPTIONAL, INTENT(IN) :: ibrav
      CHARACTER(LEN=*),  OPTIONAL, INTENT(IN) :: symm_type
      REAL(dbl),         OPTIONAL, INTENT(IN) :: celldm(6), alat
      REAL(dbl),         OPTIONAL, INTENT(IN) :: a1(3), a2(3), a3(3)
      REAL(dbl),         OPTIONAL, INTENT(IN) :: b1(3), b2(3), b3(3)
      !
      CHARACTER(LEN=256) :: bravais_lattice
      !
      CALL iotk_write_begin( iunpun, "CELL" )
      !
      IF ( PRESENT(ibrav) ) THEN
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
          CALL iotk_write_dat( iunpun, &
                               "BRAVAIS_LATTICE", TRIM( bravais_lattice ) )
      ENDIF
      !
      IF (PRESENT(symm_type)) CALL iotk_write_dat( iunpun, "CELL_SYMMETRY", symm_type )
      !
      IF ( PRESENT(alat) ) THEN
          CALL iotk_write_attr( attr, "UNITS", "Bohr", FIRST = .TRUE. )
          CALL iotk_write_dat( iunpun, "LATTICE_PARAMETER", alat, ATTR = attr )
      ENDIF
      !
      IF (PRESENT(celldm)) CALL iotk_write_dat( iunpun, "CELL_DIMENSIONS", celldm(1:6) )
      !
      CALL iotk_write_begin( iunpun, "DIRECT_LATTICE_VECTORS" )
      IF (PRESENT(a1)) CALL iotk_write_dat(   iunpun, "a1", a1(:) * alat, ATTR = attr )
      IF (PRESENT(a2)) CALL iotk_write_dat(   iunpun, "a2", a2(:) * alat, ATTR = attr )
      IF (PRESENT(a3)) CALL iotk_write_dat(   iunpun, "a3", a3(:) * alat, ATTR = attr )
      CALL iotk_write_end(   iunpun, "DIRECT_LATTICE_VECTORS" )
      !
      CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
      CALL iotk_write_begin( iunpun, "RECIPROCAL_LATTICE_VECTORS" )
      IF (PRESENT(b1)) CALL iotk_write_dat(   iunpun, "b1", b1(:), ATTR = attr )
      IF (PRESENT(b2)) CALL iotk_write_dat(   iunpun, "b2", b2(:), ATTR = attr )
      IF (PRESENT(b3)) CALL iotk_write_dat(   iunpun, "b3", b3(:), ATTR = attr )
      CALL iotk_write_end(   iunpun, "RECIPROCAL_LATTICE_VECTORS" )
      !
      CALL iotk_write_end( iunpun, "CELL" )
      !
    END SUBROUTINE qexml_write_cell
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_ions( nsp, nat, atm, ityp, psfile, &
                                pseudo_dir, amass, tau, if_pos, dirname, pos_unit )
      !------------------------------------------------------------------------
      !
      INTEGER,           INTENT(IN) :: nsp, nat
      INTEGER,           INTENT(IN) :: ityp(:)
      CHARACTER(LEN=*),  INTENT(IN) :: atm(:)
      CHARACTER(LEN=*),  INTENT(IN) :: psfile(:)
      CHARACTER(LEN=*),  INTENT(IN) :: pseudo_dir
      CHARACTER(LEN=*),  INTENT(IN) :: dirname
      REAL(dbl),         INTENT(IN) :: amass(:)
      REAL(dbl),         INTENT(IN) :: tau(:,:)
      INTEGER,           INTENT(IN) :: if_pos(:,:)
      REAL(dbl),         INTENT(IN) :: pos_unit
      !
      INTEGER            :: i, flen
      CHARACTER(LEN=256) :: file_pseudo
      LOGICAL            :: pseudo_exists
      !
      !
      CALL iotk_write_begin( iunpun, "IONS" )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_ATOMS", nat )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_SPECIES", nsp )
      !
      DO i = 1, nsp
         !
         CALL iotk_write_dat( iunpun, "ATOM_TYPE", atm(i) )
         !
         flen = LEN_TRIM( pseudo_dir )
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
         IF ( .NOT. pseudo_exists ) &
            CALL copy_file( TRIM( file_pseudo ), &
                            TRIM( dirname ) // "/" // TRIM( psfile(i) ) )
         !
         CALL iotk_write_attr( attr, "UNITS", "a.m.u.", FIRST = .TRUE. )
         CALL iotk_write_dat( iunpun, TRIM( atm(i) ) // "_MASS", &
                              amass(i), ATTR = attr )
         !
         CALL iotk_write_dat( iunpun, "PSEUDO_FOR_" // &
                            & TRIM( atm(i) ), TRIM( psfile(i) ) )
         !
      END DO
      !
      CALL iotk_write_attr( attr, "UNITS", "Bohr", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_ATOMIC_POSITIONS", attr )
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
    END SUBROUTINE qexml_write_ions
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_symmetry( ibrav, symm_type, nsym, &
                                    invsym, nr1, nr2, nr3, ftau, s, sname, irt, t_rev )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(IN) :: ibrav, nsym,  nr1, nr2, nr3
      CHARACTER(LEN=*), INTENT(IN) :: symm_type
      LOGICAL,          INTENT(IN) :: invsym
      INTEGER,          INTENT(IN) :: s(:,:,:), ftau(:,:)
      CHARACTER(LEN=*), INTENT(IN) :: sname(:)
      INTEGER,          INTENT(IN) :: irt(:,:), t_rev(:)
      !
      INTEGER   :: i
      REAL(dbl) :: tmp(3)
      !
      !
      CALL iotk_write_begin( iunpun, "SYMMETRIES" )
      !
      IF ( ibrav == 0 ) &
         CALL iotk_write_dat( iunpun, "CELL_SYMMETRY", symm_type )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_SYMMETRIES", nsym )
      !
      CALL iotk_write_dat( iunpun, "INVERSION_SYMMETRY", invsym )
      !
      DO i = 1, nsym
         !
         CALL iotk_write_attr( attr, "UNITS", "Crystal", FIRST = .TRUE. )
         !
         tmp(1) = ftau(1,i) / DBLE( nr1 )
         tmp(2) = ftau(2,i) / DBLE( nr2 )
         tmp(3) = ftau(3,i) / DBLE( nr3 )
         !
         CALL iotk_write_attr( attr, "ROT", s(:,:,i) )
         CALL iotk_write_attr( attr, "T_REV", t_rev(i) )
         CALL iotk_write_attr( attr, "FRAC_TRANS", tmp(:) )
         CALL iotk_write_attr( attr, "NAME", TRIM( sname(i) ) )
         CALL iotk_write_attr( attr, "EQ_IONS", irt(i,:) )
         !
         CALL iotk_write_empty( iunpun, &
                                "SYMM" // TRIM( iotk_index( i ) ), attr )
         !
      END DO
      !
      CALL iotk_write_end( iunpun, "SYMMETRIES" )
      !
    END SUBROUTINE qexml_write_symmetry
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_planewaves( ecutwfc, dual, npwx, gamma_only, nr1, nr2, &
                                      nr3, ngm_g, nr1s, nr2s, nr3s, ngms_g, nr1b, &
                                      nr2b, nr3b, itmp, lgvec )
      !------------------------------------------------------------------------
      !
      INTEGER,   INTENT(IN) :: npwx, nr1, nr2, nr3, ngm_g, &
                               nr1s, nr2s, nr3s, ngms_g, nr1b, nr2b, nr3b
      INTEGER,   INTENT(IN) :: itmp(:,:)
      REAL(dbl), INTENT(IN) :: ecutwfc, dual
      LOGICAL,   INTENT(IN) :: gamma_only, lgvec
      !
      !
      CALL iotk_write_begin( iunpun, "PLANE_WAVES" )
      !
      CALL iotk_write_attr( attr, "UNITS", "Hartree", FIRST = .TRUE. )
      !
      CALL iotk_write_dat( iunpun, "WFC_CUTOFF", ecutwfc / e2, ATTR = attr )
      !
      CALL iotk_write_dat( iunpun, "RHO_CUTOFF", &
                           ecutwfc * dual / e2, ATTR = attr )
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
                         "gvectors.dat", CREATE = .TRUE., BINARY = .TRUE. )
         !
         CALL iotk_write_begin( iunpun, "G-VECTORS", ATTR = attr )
         CALL iotk_write_dat( iunpun, "g", itmp(1:3,1:ngm_g), COLUMNS = 3 )
         CALL iotk_write_end( iunpun, "G-VECTORS" )
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
    END SUBROUTINE qexml_write_planewaves
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_spin( lsda, noncolin, npol, lspinorb, domag )
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
    END SUBROUTINE qexml_write_spin
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_xc( dft, nsp, lda_plus_u, &
                              Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_alpha )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=*),    INTENT(IN) :: dft
      LOGICAL,             INTENT(IN) :: lda_plus_u
      INTEGER,             INTENT(IN) :: nsp
      INTEGER,   OPTIONAL, INTENT(IN) :: Hubbard_lmax
      INTEGER,   OPTIONAL, INTENT(IN) :: Hubbard_l(:)
      REAL(dbl), OPTIONAL, INTENT(IN) :: Hubbard_U(:), Hubbard_alpha(:)
      !
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
              .NOT. PRESENT( Hubbard_U )    .OR. &
              .NOT. PRESENT( Hubbard_alpha ) ) &
            CALL errore( 'write_exchange_correlation', &
                         ' variables for LDA+U not present', 1 )
         !
         CALL iotk_write_dat( iunpun, "HUBBARD_LMAX", Hubbard_lmax )
         !
         CALL iotk_write_dat( iunpun, "HUBBARD_L", &
                              Hubbard_l(1:Hubbard_lmax) )
         !
         CALL iotk_write_dat( iunpun, "HUBBARD_U", Hubbard_U(1:nsp) )
         !
         CALL iotk_write_dat( iunpun, "HUBBARD_ALPHA", Hubbard_alpha(1:nsp) )
         !
      END IF
      !
      CALL iotk_write_end( iunpun, "EXCHANGE_CORRELATION" )
      !
    END SUBROUTINE qexml_write_xc
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_occ( lgauss, ngauss, degauss, ltetra, ntetra, &
                               tetra, tfixed_occ, lsda, nelup, neldw, occupations )
      !------------------------------------------------------------------------
      !
      LOGICAL,             INTENT(IN) :: lgauss, ltetra, tfixed_occ, lsda
      INTEGER,   OPTIONAL, INTENT(IN) :: ngauss, ntetra, nelup, neldw
      INTEGER,   OPTIONAL, INTENT(IN) :: tetra(:,:)
      REAL(dbl), OPTIONAL, INTENT(IN) :: degauss, occupations(:,:)      
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
         CALL iotk_write_dat( iunpun, "INPUT_OCC_UP", occupations(1:nelup,1) )
         !
         IF ( lsda ) &
            CALL iotk_write_dat( iunpun, "INPUT_OCC_DOWN", occupations(1:neldw,2) )
         !
      END IF
      !
      CALL iotk_write_end( iunpun, "OCCUPATIONS" )
      !
    END SUBROUTINE qexml_write_occ
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_bz( num_k_points, xk, wk, k1, k2, k3, nk1, nk2, nk3 )
      !------------------------------------------------------------------------
      !
      INTEGER,   INTENT(IN) :: num_k_points, k1, k2, k3, nk1, nk2, nk3
      REAL(dbl), INTENT(IN) :: xk(:,:), wk(:)
      !
      INTEGER :: ik
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
      CALL iotk_write_end( iunpun, "BRILLOUIN_ZONE" )
      !
    END SUBROUTINE qexml_write_bz
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_phonon( modenum, xqq )
      !------------------------------------------------------------------------
      !
      INTEGER,   INTENT(IN) :: modenum
      REAL(dbl), INTENT(IN) :: xqq(:)
      !
      !
      CALL iotk_write_begin( iunpun, "PHONON" )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_MODES", modenum )
      !
      CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_Q-POINT", attr )
      !
      CALL iotk_write_dat( iunpun, "Q-POINT", xqq(:) )
      !
      CALL iotk_write_end( iunpun, "PHONON" )
      !
    END SUBROUTINE qexml_write_phonon
    !
    ! ... methods to write and read charge_density
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_rho( rho_file_base, rho, &
                          nr1, nr2, nr3, nr1x, nr2x, ipp, npp )
      !------------------------------------------------------------------------
      !
      ! ... Writes charge density rho, one plane at a time.
      ! ... If ipp and npp are specified, planes are collected one by one from
      ! ... all processors, avoiding an overall collect of the charge density
      ! ... on a single proc.
      !
#ifdef __HAVE_RHO_WRITE 
!      USE mp_global, ONLY : me_image, intra_image_comm, me_pool, nproc_pool, &
!                            intra_pool_comm, my_pool_id
!      USE mp,        ONLY : mp_get
      !
#endif
      IMPLICIT NONE
      !
      CHARACTER(LEN=*),   INTENT(IN) :: rho_file_base
      INTEGER,            INTENT(IN) :: nr1, nr2, nr3
      INTEGER,            INTENT(IN) :: nr1x, nr2x
      REAL(dbl),          INTENT(IN) :: rho(:)
      INTEGER, OPTIONAL,  INTENT(IN) :: ipp(:)
      INTEGER, OPTIONAL,  INTENT(IN) :: npp(:)
#ifdef __HAVE_RHO_WRITE 
      !
      INTEGER                :: ierr, i, j, k, kk, ldr, ip
      CHARACTER(LEN=256)     :: rho_file
      REAL(dbl), ALLOCATABLE :: rho_plane(:)
      INTEGER,   ALLOCATABLE :: kowner(:)
      INTEGER                :: iopool_id, ionode_pool
      !
     
      !
      rho_file = TRIM( rho_file_base ) // '.xml'
      !
      IF ( ionode ) &
         CALL iotk_open_write( rhounit, FILE = rho_file, &
                               BINARY = rho_binary, IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'write_rho', 'cannot open' // &
                 & TRIM( rho_file ) // ' file for writing', ierr )
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
      ! ... find the index of the pool that will write rho
      !
      IF ( ionode ) iopool_id = my_pool_id
      !
      CALL mp_bcast( iopool_id, ionode_id, intra_image_comm )
      !
      ! ... find the index of the ionode within its own pool
      !
      IF ( ionode ) ionode_pool = me_pool
      !
      CALL mp_bcast( ionode_pool, ionode_id, intra_image_comm )
      !
      ! ... find out the owner of each "z" plane
      !
      IF ( PRESENT( ipp ) .AND. PRESENT( npp ) ) THEN
         !
         DO ip = 1, nproc_pool
            !
            kowner( (ipp(ip)+1):(ipp(ip)+npp(ip)) ) = ip - 1
            !
         END DO
         !
      ELSE
         !
         kowner = ionode_id
         !
      END IF
      !
      ldr = nr1x*nr2x
      !
      DO k = 1, nr3
         !
         IF( kowner(k) == me_pool ) THEN
            !
            kk = k
            !
            IF ( PRESENT( ipp ) ) kk = k - ipp(me_pool+1)
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
         IF ( kowner(k) /= ionode_pool .AND. my_pool_id == iopool_id ) &
            CALL mp_get( rho_plane, rho_plane, &
                         me_pool, ionode_pool, kowner(k), k, intra_pool_comm )
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
#endif
      RETURN
      !
    END SUBROUTINE qexml_write_rho
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_rho( rho_file_base, rho, &
                              nr1, nr2, nr3, nr1x, nr2x, ipp, npp )
      !------------------------------------------------------------------------
      !
      ! ... Writes charge density rho, one plane at a time.
      ! ... If ipp and npp are specified, planes are collected one by one from
      ! ... all processors, avoiding an overall collect of the charge density
      ! ... on a single proc.
      !
#ifdef __HAVE_RHO_READ
!      USE io_global, ONLY : ionode, ionode_id
!      USE mp_global, ONLY : me_image, intra_image_comm, me_pool, nproc_pool, &
!                            intra_pool_comm, my_pool_id, npool
!      USE mp,        ONLY : mp_put
#endif
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*),   INTENT(IN)  :: rho_file_base
      INTEGER,            INTENT(IN)  :: nr1, nr2, nr3
      INTEGER,            INTENT(IN)  :: nr1x, nr2x
      REAL(dbl),          INTENT(OUT) :: rho(:)
      INTEGER, OPTIONAL,  INTENT(IN)  :: ipp(:)
      INTEGER, OPTIONAL,  INTENT(IN)  :: npp(:)
      !
      INTEGER                :: ierr, i, j, k, kk, ldr, ip
      INTEGER                :: nr( 3 )
      CHARACTER(LEN=256)     :: rho_file
      REAL(dbl), ALLOCATABLE :: rho_plane(:)
      INTEGER,   ALLOCATABLE :: kowner(:)
      INTEGER                :: iopool_id, ionode_pool
      !
      !
#ifdef __HAVE_RHO_READ

      rho_file = TRIM( rho_file_base ) // '.xml'
      !
      IF ( ionode ) &
         CALL iotk_open_read( rhounit, FILE = rho_file, &
                              BINARY = rho_binary, IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'read_rho', 'cannot open ' // &
                 & TRIM( rho_file ) // ' file for reading', ierr )
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
      END IF
      !
      CALL mp_bcast( nr, ionode_id, intra_image_comm )
      !
      IF ( nr1 /= nr(1) .OR. nr2 /= nr(2) .OR. nr3 /= nr(3) ) &
         CALL errore( 'read_rho', 'dimensions do not match', 1 )
      !
      ALLOCATE( rho_plane( nr1*nr2 ) )
      ALLOCATE( kowner( nr3 ) )
      !
      ! ... find the index of the pool that will write rho
      !
      IF ( ionode ) iopool_id = my_pool_id
      !
      CALL mp_bcast( iopool_id, ionode_id, intra_image_comm )
      !
      ! ... find the index of the ionode within its own pool
      !
      IF ( ionode ) ionode_pool = me_pool
      !
      CALL mp_bcast( ionode_pool, ionode_id, intra_image_comm )
      !
      ! ... find out the owner of each "z" plane
      !
      IF ( PRESENT( ipp ) .AND. PRESENT( npp ) ) THEN
         !
         DO ip = 1, nproc_pool
            !
            kowner((ipp(ip)+1):(ipp(ip)+npp(ip))) = ip - 1
            !
         END DO
         !
      ELSE
         !
         kowner = ionode_id
         !
      END IF
      !
      ldr = nr1x*nr2x
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
         IF( npool > 1 ) THEN
            !
            !  send to all proc/pools
            !
            CALL mp_bcast( rho_plane, ionode_id, intra_image_comm )
            !
         ELSE
            !
            !  send to the destination proc
            !
            IF ( kowner(k) /= ionode_id ) &
               CALL mp_put( rho_plane, rho_plane, me_image, &
                            ionode_id, kowner(k), k, intra_image_comm )
            !
         END IF
         !
         IF( kowner(k) == me_pool ) THEN
            !
            kk = k
            !
            IF ( PRESENT( ipp ) ) kk = k - ipp(me_pool+1)
            ! 
            DO j = 1, nr2
               !
               DO i = 1, nr1
                  !
                  rho(i+(j-1)*nr1x+(kk-1)*ldr) = rho_plane(i+(j-1)*nr1)
                  !
               END DO
               !
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
#endif
      RETURN
      !
    END SUBROUTINE qexml_read_rho
    !
    ! ... methods to write and read wavefunctions
    !
    !------------------------------------------------------------------------
    SUBROUTINE qexml_write_wfc( iuni, ik, nk, kunit, ispin, &
                          nspin, wf0, ngw, nbnd, igl, ngwl, filename, scalef )
      !------------------------------------------------------------------------
      !
#ifdef __HAVE_WRITE_WFC
!      USE mp_wave,    ONLY : mergewf
!      USE mp,         ONLY : mp_get
!      USE mp_global,  ONLY : me_pool, nproc_image, nproc_pool, &
!                             root_pool, intra_pool_comm, me_image, &
!                             intra_image_comm
#endif
      !
      IMPLICIT NONE
      !
      INTEGER,             INTENT(IN) :: iuni
      INTEGER,             INTENT(IN) :: ik, nk, kunit, ispin, nspin
      COMPLEX(dbl),        INTENT(IN) :: wf0(:,:)
      INTEGER,             INTENT(IN) :: ngw
      INTEGER,             INTENT(IN) :: nbnd
      INTEGER,             INTENT(IN) :: ngwl
      INTEGER,             INTENT(IN) :: igl(:)
      CHARACTER(LEN=256),  INTENT(IN) :: filename
      REAL(dbl),           INTENT(IN) :: scalef    
        ! scale factor, usually 1.0 for pw and 1/SQRT( omega ) for CP
#ifdef __HAVE_WRITE_WFC
      !
      INTEGER                   :: i, j, ierr
      INTEGER                   :: iks, ike, ikt, igwx
      INTEGER                   :: npool, ipmask(nproc_image), ipsour
      COMPLEX(dbl), ALLOCATABLE :: wtmp(:)
      !
      !
      CALL set_kpoints_vars( ik, nk, kunit, ngwl, igl, &
                             npool, ikt, iks, ike, igwx, ipmask, ipsour )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_open_write( iuni, FILE = TRIM( filename ), BINARY = .TRUE. )
         !
         CALL iotk_write_attr( attr, "ngw",          ngw, FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "nbnd",         nbnd )
         CALL iotk_write_attr( attr, "ik",           ik )
         CALL iotk_write_attr( attr, "nk",           nk )
         CALL iotk_write_attr( attr, "kunit",        kunit )
         CALL iotk_write_attr( attr, "ispin",        ispin )
         CALL iotk_write_attr( attr, "nspin",        nspin )
         CALL iotk_write_attr( attr, "igwx",         igwx )
         CALL iotk_write_attr( attr, "scale_factor", scalef )
         !
         CALL iotk_write_empty( iuni, "INFO", attr )
         !
      END IF
      !
      ALLOCATE( wtmp( MAX( igwx, 1 ) ) )
      !
      wtmp = 0.D0
      !
      DO j = 1, nbnd
         !
         IF ( npool > 1 ) THEN
            !
            IF ( ikt >= iks .AND. ikt <= ike ) &      
               CALL mergewf( wf0(:,j), wtmp, ngwl, igl, me_pool, &
                             nproc_pool, root_pool, intra_pool_comm )
            !
            IF ( ipsour /= ionode_id ) &
               CALL mp_get( wtmp, wtmp, me_image, ionode_id, ipsour, j, intra_image_comm )
            !
         ELSE
            !
            CALL mergewf( wf0(:,j), wtmp, ngwl, igl, me_image, nproc_image, & 
                          ionode_id, intra_image_comm )
            !
         END IF
         !
         IF ( ionode ) &
            CALL iotk_write_dat( iuni, "evc" // iotk_index( j ), wtmp(1:igwx) )
         !
      END DO
      !
      IF ( ionode ) CALL iotk_close_write( iuni )
      !
      DEALLOCATE( wtmp )
      !
#endif
      RETURN
      !
    END SUBROUTINE qexml_write_wfc
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
      CALL iotk_scan_begin( iunpun, "HEADER", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_empty( iunpun, "FORMAT", ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_attr(attr, "NAME", format_name_, IERR=ierr)
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr(attr, "VERSION", format_version_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_empty( iunpun, "CREATOR", ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_attr(attr, "NAME", creator_name_, IERR=ierr)
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr(attr, "VERSION", creator_version_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_end( iunpun, "HEADER", IERR=ierr )
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
      CALL iotk_scan_begin( iunpun, "CELL" )
      !
      CALL iotk_scan_dat( iunpun, "BRAVAIS_LATTICE", bravais_latt_, IERR=ierr )
      IF ( ierr /= 0 ) RETURN
      !
      IF ( TRIM( bravais_latt_ ) == "Trigonal R" .OR. &
           TRIM( bravais_latt_ ) == "Hexagonal and Trigonal P" ) THEN
         !
         symm_type_ = 'hexagonal'
         !
      ELSE
         !
         symm_type_ = 'cubic'
         !
      END IF
      !
      CALL iotk_scan_dat( iunpun, "LATTICE_PARAMETER", alat_, ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "UNITS", alat_units_, IERR=ierr )
      IF ( ierr /= 0 ) RETURN
      !
      CALL iotk_scan_dat( iunpun, "CELL_DIMENSIONS", celldm_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_begin( iunpun, "DIRECT_LATTICE_VECTORS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_empty( iunpun, "UNITS_FOR_DIRECT_LATTICE_VECTORS", ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "UNITS", a_units_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_dat(   iunpun, "a1", a1_(:), ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_dat(   iunpun, "a2", a2_(:), IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_dat(   iunpun, "a3", a3_(:), IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_end(   iunpun, "DIRECT_LATTICE_VECTORS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_begin( iunpun, "RECIPROCAL_LATTICE_VECTORS", IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_empty( iunpun, "UNITS_FOR_RECIPROCAL_LATTICE_VECTORS", ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "UNITS", b_units_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_dat(   iunpun, "b1", b1_(:), ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_dat(   iunpun, "b2", b2_(:), IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_dat(   iunpun, "b3", b3_(:), IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_end(   iunpun, "RECIPROCAL_LATTICE_VECTORS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_end( iunpun, "CELL", IERR=ierr )
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
      CALL iotk_scan_begin( iunpun, "IONS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "NUMBER_OF_ATOMS", nat_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "NUMBER_OF_SPECIES", nsp_ )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_empty( iunpun, "UNITS_FOR_ATOMIC_MASSES", ATTR=attr, IERR=ierr )
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
         CALL iotk_scan_dat( iunpun, "ATOM_TYPE", atm_(i), IERR=ierr )
         IF (ierr/=0) RETURN
         CALL iotk_scan_dat( iunpun, TRIM( atm_(i) ) // "_MASS", amass_(i), IERR=ierr )
         IF (ierr/=0) RETURN
         CALL iotk_scan_dat( iunpun, "PSEUDO_FOR_" // TRIM( atm_(i) ), psfile_(i), IERR=ierr )
         IF (ierr/=0) RETURN
         !
      ENDDO
      !
      CALL iotk_scan_empty( iunpun, "UNITS_FOR_ATOMIC_POSITIONS", ATTR=attr, IERR=ierr )
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
         CALL iotk_scan_empty( iunpun, &
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
      CALL iotk_scan_end( iunpun, "IONS", IERR=ierr )
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
      CALL iotk_scan_begin( iunpun, "SYMMETRIES", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "NUMBER_OF_SYMMETRIES", nsym_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "INVERSION_SYMMETRY", invsym_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "NUMBER_OF_ATOMS", nat_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      ALLOCATE( irt_(48, nat_) )
      !
      CALL iotk_scan_empty( iunpun, "UNITS_FOR_SYMMETRIES", ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "UNITS", s_units_, IERR=ierr )
      IF (ierr/=0) RETURN
      ! 
      DO i = 1, nsym_
          !
          CALL iotk_scan_begin( iunpun, "SYMM"//TRIM( iotk_index( i ) ), IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_empty( iunpun, "INFO", ATTR=attr, IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_attr( attr, "NAME", sname_(i), IERR=ierr )
          IF (ierr/=0) RETURN
          CALL iotk_scan_attr( attr, "T_REV", t_rev_(i), IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_dat( iunpun, "ROTATION", s_(:,:,i), IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_dat( iunpun, "FRACTIONAL_TRANSLATION", trasl_(:,i), IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_dat( iunpun, "EQUIVALENT_IONS", irt_(i,:), IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_end( iunpun, "SYMM"//TRIM( iotk_index( i ) ), IERR=ierr )
          IF (ierr/=0) RETURN
          !
      ENDDO
      !
      CALL iotk_scan_end( iunpun, "SYMMETRIES", IERR=ierr )
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


    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_planewaves( ecutwfc, ecutrho, npwx, gamma_only, nr1, nr2,  &
                                      nr3, ngm, nr1s, nr2s, nr3s, ngms, nr1b, &
                                      nr2b, nr3b, igv, cutoff_units, ierr )
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
      CALL iotk_scan_begin( iunpun, "PLANE_WAVES", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_empty( iunpun, "UNITS_FOR_CUTOFF", ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "UNITS", cutoff_units_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "WFC_CUTOFF", ecutwfc_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "RHO_CUTOFF", ecutrho_ , IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "MAX_NUMBER_OF_GK-VECTORS", npwx_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "GAMMA_ONLY", gamma_only_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_empty( iunpun, "FFT_GRID", ATTR = attr, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_attr( attr, "nr1", nr1_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "nr2", nr2_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "nr3", nr3_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "GVECT_NUMBER", ngm_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_empty( iunpun, "SMOOTH_FFT_GRID", ATTR = attr, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_attr( attr, "nr1s", nr1s_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "nr2s", nr2s_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "nr3s", nr3s_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "SMOOTH_GVECT_NUMBER", ngms_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( PRESENT( igv ) ) THEN
          !
          CALL iotk_scan_begin( iunpun, "G-VECTORS", IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_dat( iunpun, "g", igv(1:3,1:ngm_), IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_end( iunpun, "G-VECTORS", IERR=ierr )          
          IF (ierr/=0) RETURN
          !
      ENDIF
      !
      !
      CALL iotk_scan_empty( iunpun, "SMALLBOX_FFT_GRID", ATTR = attr, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_attr( attr, "nr1b", nr1b_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "nr2b", nr2b_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "nr3b", nr3b_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_end( iunpun, "PLANE_WAVES", IERR=ierr )
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


    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_gk( ik, npwk, npwkx, idx, igk, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,                INTENT(IN)  :: ik
      INTEGER,      OPTIONAL, INTENT(OUT) :: npwk, npwkx
      INTEGER,      OPTIONAL, INTENT(OUT) :: igk(:,:), idx(:)
      INTEGER,                INTENT(OUT) :: ierr
      !
      INTEGER :: npwk_, npwkx_

      ierr = 0

      !
      !
      CALL iotk_scan_begin( iunpun, "BAND_STRUCTURE", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_begin( iunpun, "EIGENVALUES_AND_EIGENVECTORS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "MAX_NUMBER_OF_GK-VECTORS", npwkx_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_begin( iunpun, "K-POINT" //TRIM(iotk_index(ik)), IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "NUMBER_OF_GK-VECTORS", npwk_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_begin( iunpun, "GK-VECTORS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_begin( iunpun, "K-POINT" //TRIM(iotk_index(ik)), IERR=ierr )
      IF (ierr/=0) RETURN
      !
      IF ( PRESENT( idx )  )  THEN
          !
          CALL iotk_scan_dat( iunpun, "INDEX", idx(1: npwk_) , IERR=ierr )
          IF (ierr/=0) RETURN
          !
      ENDIF
      !
      IF ( PRESENT( igk )  )  THEN
          !
          CALL iotk_scan_dat( iunpun, "GRID", igk(1:3, 1:npwk_) , IERR=ierr )
          IF (ierr/=0) RETURN
          !
      ENDIF
      !
      CALL iotk_scan_end( iunpun, "K-POINT" //TRIM(iotk_index(ik)), IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_end( iunpun, "GK-VECTORS", IERR=ierr )
      IF (ierr/=0) RETURN
      !    
      !    
      CALL iotk_scan_end( iunpun, "K-POINT" //TRIM(iotk_index(ik)), IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_end( iunpun, "EIGENVALUES_AND_EIGENVECTORS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_end( iunpun, "BAND_STRUCTURE", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( PRESENT( npwk) )     npwk  = npwk_
      IF ( PRESENT( npwkx) )    npwkx = npwkx_
      !
    END SUBROUTINE qexml_read_gk


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
      CALL iotk_scan_begin( iunpun, "SPIN", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "LSDA", lsda_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "NON-COLINEAR_CALCULATION", noncolin_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      npol_ = 1
      !
      IF ( noncolin_ ) THEN
          !
          CALL iotk_scan_dat( iunpun, "SPINOR_DIM", npol_, IERR=ierr )
          IF (ierr/=0) RETURN
          !
      ENDIF
      !
      CALL iotk_scan_dat( iunpun, "SPIN-ORBIT_CALCULATION", lspinorb_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_end( iunpun, "SPIN", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( PRESENT( lsda ) )       lsda      = lsda_
      IF ( PRESENT( noncolin ) )   noncolin  = noncolin_
      IF ( PRESENT( npol ) )       npol      = npol_
      IF ( PRESENT( lspinorb ) )   lspinorb  = lspinorb_
      !

    END SUBROUTINE qexml_read_spin


    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_xc( dft, lda_plus_u,  &
                              Hubbard_lmax, Hubbard_l, nsp, Hubbard_U, Hubbard_alpha, ierr )
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
      CALL iotk_scan_begin( iunpun, "EXCHANGE_CORRELATION", IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      CALL iotk_scan_dat( iunpun, "DFT", dft_, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      CALL iotk_scan_dat( iunpun, "LDA_PLUS_U_CALCULATION", lda_plus_u_, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      IF ( lda_plus_u_ ) THEN
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_SPECIES", nsp_, IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
         CALL iotk_scan_dat( iunpun, "HUBBARD_LMAX", Hubbard_lmax_, IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
         ALLOCATE( Hubbard_l_(1:Hubbard_lmax_) )
         ALLOCATE( Hubbard_U_(nsp_) )
         ALLOCATE( Hubbard_alpha_(nsp_) )
         !
         CALL iotk_scan_dat( iunpun, "HUBBARD_L", Hubbard_l_, IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
         CALL iotk_scan_dat( iunpun, "HUBBARD_U", Hubbard_U_, IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
         CALL iotk_scan_dat( iunpun, "HUBBARD_ALPHA", Hubbard_alpha_, IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
      ENDIF
      !
      CALL iotk_scan_end( iunpun, "EXCHANGE_CORRELATION", IERR=ierr )
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


    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_occ( lgauss, ngauss, degauss, degauss_units, ltetra, ntetra, &
                               tetra, tfixed_occ, occupations, ierr )
      !------------------------------------------------------------------------
      !
      LOGICAL,      OPTIONAL, INTENT(OUT) :: lgauss, ltetra, tfixed_occ
      INTEGER,      OPTIONAL, INTENT(OUT) :: ngauss, ntetra
      INTEGER,      OPTIONAL, INTENT(OUT) :: tetra(:,:)
      REAL(dbl),    OPTIONAL, INTENT(OUT) :: degauss, occupations(:,:)
      CHARACTER(*), OPTIONAL, INTENT(OUT) :: degauss_units
      INTEGER,                INTENT(OUT) :: ierr
      !
      LOGICAL        :: lgauss_, ltetra_, tfixed_occ_
      INTEGER        :: ngauss_, ntetra_
      REAL(dbl)      :: degauss_
      CHARACTER(256) :: degauss_units_
      INTEGER,  ALLOCATABLE :: tetra_(:,:)
      INTEGER :: i
      LOGICAL :: lfound
      !
      ierr = 0 
      !
      CALL iotk_scan_begin( iunpun, "OCCUPATIONS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "SMEARING_METHOD", lgauss_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( lgauss_ ) THEN
         !
         CALL iotk_scan_dat( iunpun, "SMEARING_TYPE", ngauss_, IERR=ierr )
         IF (ierr/=0) RETURN
         !
         CALL iotk_scan_dat( iunpun, "SMEARING_PARAMETER", degauss_ , ATTR=attr, IERR=ierr )
         IF (ierr/=0) RETURN
         !
         CALL iotk_scan_attr( ATTR, "UNITS", degauss_units_ , IERR=ierr )
         IF (ierr/=0) RETURN
         !
      ENDIF
      !
      CALL iotk_scan_dat( iunpun, "TETRAHEDRON_METHOD", ltetra_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( ltetra_ ) THEN
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_TETRAHEDRA", ntetra_, IERR=ierr )
         IF (ierr/=0) RETURN
         !
         ALLOCATE( tetra_(4, ntetra_) )
         !
         DO i = 1, ntetra_
            !
            CALL iotk_scan_dat( iunpun, "TETRAHEDRON"//iotk_index(i), tetra_(1:4,i), IERR=ierr )
            IF (ierr/=0) RETURN
            !
         ENDDO
         !
      ENDIF
      !
      CALL iotk_scan_dat( iunpun, "FIXED_OCCUPATIONS", tfixed_occ_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      IF ( tfixed_occ_  .AND. PRESENT( occupations ) ) THEN
         !
         CALL iotk_scan_dat( iunpun, "INPUT_OCC_UP", occupations(:,1), IERR=ierr )
         IF (ierr/=0) RETURN
         !
         IF ( SIZE(occupations, 2) >= 2  ) THEN
            !
            CALL iotk_scan_dat( iunpun, "INPUT_OCC_DOWN", occupations(:,2), FOUND=lfound, IERR=ierr )
            IF (ierr/=0) RETURN
            !
         ENDIF
         !
      ENDIF
      !
      CALL iotk_scan_end( iunpun, "OCCUPATIONS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( PRESENT( lgauss ))           lgauss     = lgauss_
      IF ( PRESENT( ltetra ))           ltetra     = ltetra_
      IF ( PRESENT( tfixed_occ ))       tfixed_occ = tfixed_occ_
      IF ( PRESENT( ngauss ))           ngauss     = ngauss_
      IF ( PRESENT( ntetra ))           ntetra     = ntetra_
      IF ( PRESENT( degauss ))          degauss    = degauss
      IF ( PRESENT( degauss_units ))    degauss_units  = TRIM(degauss_units_)
      !
      IF ( ltetra_ ) THEN
         !
         IF ( PRESENT( tetra ) )         tetra(1:4, 1:ntetra_)  = tetra_
         !
         DEALLOCATE( tetra_ )
         !
      ENDIF

    END SUBROUTINE qexml_read_occ


    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_bz( num_k_points, xk, wk, k1, k2, k3, nk1, nk2, nk3, k_units, ierr )
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
      CALL iotk_scan_begin( iunpun, "BRILLOUIN_ZONE", IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      CALL iotk_scan_dat( iunpun, "NUMBER_OF_K-POINTS", num_k_points_, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      !
      CALL iotk_scan_empty( iunpun, "UNITS_FOR_K-POINTS", ATTR=attr, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      CALL iotk_scan_attr( attr, "UNITS", k_units_, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      CALL iotk_scan_empty( iunpun, "MONKHORST_PACK_GRID", ATTR=attr, IERR=ierr )
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
      CALL iotk_scan_empty( iunpun, "MONKHORST_PACK_OFFSET", ATTR=attr, IERR=ierr )
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
         CALL iotk_scan_empty( iunpun, "K-POINT" // TRIM( iotk_index(ik) ), ATTR=attr, IERR=ierr )
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
      CALL iotk_scan_end( iunpun, "BRILLOUIN_ZONE", IERR=ierr )
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
      CALL iotk_scan_begin( iunpun, "PHONON", IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      CALL iotk_scan_dat( iunpun, "NUMBER_OF_MODES", modenum_, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      CALL iotk_scan_empty( iunpun, "UNITS_FOR_Q-POINT", attr, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      CALL iotk_scan_attr( attr, "UNITS", q_units_, IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      IF ( PRESENT (xqq) ) THEN
         !
         CALL iotk_scan_dat( iunpun, "Q-POINT", xqq(:), IERR=ierr )
         IF ( ierr/=0 ) RETURN
         !
      ENDIF
      !
      CALL iotk_scan_end( iunpun, "PHONON", IERR=ierr )
      IF ( ierr/=0 ) RETURN
      !
      !
      IF ( PRESENT (modenum) )      modenum = modenum_
      IF ( PRESENT (q_units) )      q_units = TRIM(q_units_)
      !
    END SUBROUTINE qexml_read_phonon


    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_bands( nbnd, num_k_points, nspin, natomwfc, ef, nelec, &
                                 xk, wk, occ, occ_s, eig, eig_s, energy_units, k_units, ierr )
      !------------------------------------------------------------------------
      !
      INTEGER,      OPTIONAL, INTENT(OUT) :: nbnd, num_k_points, nspin, natomwfc
      REAL(dbl),    OPTIONAL, INTENT(OUT) :: ef, nelec, xk(:,:), wk(:), occ(:,:), occ_s(:,:,:)
      REAL(dbl),    OPTIONAL, INTENT(OUT) :: eig(:,:), eig_s(:,:,:)
      CHARACTER(*), OPTIONAL, INTENT(OUT) :: energy_units, k_units
      INTEGER,                INTENT(OUT) :: ierr
      !
      INTEGER        :: nbnd_, num_k_points_, nspin_, natomwfc_
      REAL(dbl)      :: ef_, nelec_ 
      CHARACTER(256) :: energy_units_, k_units_
      INTEGER        :: ik, ispin
      REAL(dbl), ALLOCATABLE :: xk_(:,:), wk_(:), occ_s_(:,:,:), eig_s_(:,:,:)
      !
      
      ierr = 0
      !
      ! open the main section
      !
      CALL iotk_scan_begin( iunpun, "BAND_STRUCTURE", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "NUMBER_OF_ELECTRONS", nelec_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "NUMBER_OF_ATOMIC_WFC", natomwfc_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "NUMBER_OF_BANDS", nbnd_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "NUMBER_OF_K-POINTS", num_k_points_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "FERMI_ENERGY", ef_, ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "UNITS", energy_units_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_dat( iunpun, "NUMBER_OF_SPIN_COMPONENTS", nspin_, IERR=ierr )
      IF (ierr/=0) RETURN

      !
      !
      ! Allocations
      !
      ALLOCATE(  xk_    (3, num_k_points_ ) )
      ALLOCATE(  wk_    (   num_k_points_ ) )
      ALLOCATE(  eig_s_ ( nbnd_, num_k_points_, nspin_ ) )
      ALLOCATE(  occ_s_ ( nbnd_, num_k_points_, nspin_ ) )
      !
      CALL iotk_scan_begin( iunpun, "EIGENVALUES_AND_EIGENVECTORS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      DO ik = 1, num_k_points_
          !
          CALL iotk_scan_begin( iunpun, "K-POINT" // TRIM(iotk_index( ik )), IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_dat( iunpun, "K-POINT_COORDS", xk_(:,ik), ATTR=attr, IERR=ierr )
          IF (ierr/=0) RETURN
          CALL iotk_scan_attr( attr, "UNITS", k_units_, IERR=ierr )
          IF (ierr/=0) RETURN
          !
          CALL iotk_scan_dat( iunpun, "WEIGHT", wk_(ik), IERR=ierr )
          IF (ierr/=0) RETURN
          !
          DO ispin = 1, nspin_
              !
              CALL iotk_scan_dat( iunpun, "ET" // TRIM( iotk_index( ispin ) ), &
                                  eig_s_(:,ik,ispin), IERR=ierr )
              IF (ierr/=0) RETURN
              !
              CALL iotk_scan_dat( iunpun, "OCC" // TRIM( iotk_index( ispin ) ), &
                                  occ_s_(:,ik,ispin), IERR=ierr )
              IF (ierr/=0) RETURN
              !
          ENDDO
          !
          CALL iotk_scan_end( iunpun, "K-POINT" // TRIM(iotk_index( ik )), IERR=ierr )
          IF (ierr/=0) RETURN
          !
      ENDDO    
      !         
      CALL iotk_scan_end( iunpun, "EIGENVALUES_AND_EIGENVECTORS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_end( iunpun, "BAND_STRUCTURE", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( PRESENT( nbnd ) )             nbnd           = nbnd_
      IF ( PRESENT( num_k_points ) )     num_k_points   = num_k_points_
      IF ( PRESENT( nspin ) )            nspin          = nspin_
      IF ( PRESENT( nelec ) )            nelec          = nelec_
      IF ( PRESENT( natomwfc ) )         natomwfc       = natomwfc_
      IF ( PRESENT( ef ) )               ef             = ef_
      IF ( PRESENT( energy_units ) )     energy_units   = TRIM( energy_units_ )
      IF ( PRESENT( k_units ) )          k_units        = TRIM( k_units_ )
      IF ( PRESENT( xk ) )               xk(1:3, 1:num_k_points_)  = xk_(:,:)
      IF ( PRESENT( wk ) )               wk(     1:num_k_points_)  = wk_(:)
      IF ( PRESENT( occ ) )              occ  (1:nbnd_, 1:num_k_points_ )           = occ_s_(:,:,1)
      IF ( PRESENT( occ_s ) )            occ_s(1:nbnd_, 1:num_k_points_, 1:nspin_ ) = occ_s_(:,:,:)
      IF ( PRESENT( eig ) )              eig  (1:nbnd_, 1:num_k_points_ )           = eig_s_(:,:,1)
      IF ( PRESENT( eig_s ) )            eig_s(1:nbnd_, 1:num_k_points_, 1:nspin_ ) = eig_s_(:,:,:)
      !
      DEALLOCATE( xk_)
      DEALLOCATE( wk_)
      DEALLOCATE( occ_s_ )
      DEALLOCATE( eig_s_ )
      !
    END SUBROUTINE qexml_read_bands


    !------------------------------------------------------------------------
    SUBROUTINE qexml_read_wfc( ibnds, ibnde, ik, ispin, npwk, igk, ngw, igwx, wf, wf_kindip, ierr )
      !------------------------------------------------------------------------
      !
      ! read wfc from IBNDS to IBNDE, for kpt IK and spin ISPIN
      ! WF is the wfc on itsproper k+g grid, while WF_KINDIP is the same wfc
      ! but on a truncated rho grid (k-point indipendent)
      !
      INTEGER,                 INTENT(IN)  :: ibnds, ibnde, ik, ispin
      INTEGER,       OPTIONAL, INTENT(IN)  :: npwk, igk(:)
      INTEGER,       OPTIONAL, INTENT(OUT) :: ngw, igwx
      COMPLEX(dbl),  OPTIONAL, INTENT(OUT) :: wf(:,:), wf_kindip(:,:)
      INTEGER,                 INTENT(OUT) :: ierr
      !
      INTEGER :: ngw_, igwx_, ig, ib, lindex
      COMPLEX(dbl),  ALLOCATABLE :: wf_(:)

      ierr = 0
      !
      ! read some dimensions
      !
      CALL iotk_scan_begin( iunpun, "BAND_STRUCTURE", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_begin( iunpun, "EIGENVALUES_AND_EIGENVECTORS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_begin( iunpun, "K-POINT" //TRIM(iotk_index(ik)), IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_begin( iunpun, "WFC" // TRIM(iotk_index( ispin) ) , IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      CALL iotk_scan_empty( iunpun, "INFO", ATTR=attr, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_attr( attr, "ngw",  ngw_, IERR=ierr )
      IF (ierr/=0) RETURN
      CALL iotk_scan_attr( attr, "igwx", igwx_, IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( PRESENT( wf_kindip )  )  THEN
          !
          lindex = 0
          !
          DO ib = ibnds, ibnde
              !
              lindex = lindex + 1
              !
              CALL iotk_scan_dat( iunpun, "evc"//TRIM(iotk_index(ib)), &
                                  wf_kindip( 1:igwx_, lindex ), IERR=ierr )
              IF (ierr/=0) RETURN
              !
          ENDDO
          !
      ENDIF
      !
      IF ( PRESENT( wf )  )  THEN
          !
          ALLOCATE( wf_(igwx_ ), STAT=ierr )
          IF (ierr/=0) RETURN
          !
          IF ( .NOT. PRESENT( igk ) .OR. .NOT. PRESENT( npwk ) ) THEN
              ierr = 3
              RETURN
          ENDIF
          !
          IF ( MAXVAL( igk( 1:npwk ) ) > igwx_ ) THEN 
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
              CALL iotk_scan_dat( iunpun, "evc"//TRIM(iotk_index( ib ) ), wf_(1:igwx_), IERR=ierr )
              IF (ierr/=0) RETURN
              !
              ! use the igk map to do the transformation
              !
              DO ig = 1, npwk
                  !
                  wf( ig, lindex ) = wf_( igk( ig )  )
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
      !
      CALL iotk_scan_end( iunpun, "WFC" // TRIM( iotk_index(ispin) ), IERR=ierr )
      IF (ierr/=0) RETURN
      !    
      CALL iotk_scan_end( iunpun, "K-POINT" //TRIM(iotk_index(ik)), IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_end( iunpun, "EIGENVALUES_AND_EIGENVECTORS", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      CALL iotk_scan_end( iunpun, "BAND_STRUCTURE", IERR=ierr )
      IF (ierr/=0) RETURN
      !
      !
      IF ( PRESENT( ngw ) )     ngw    = ngw_
      IF ( PRESENT( igwx ) )    igwx   = igwx_
      !
    END SUBROUTINE qexml_read_wfc
    !
    !
END MODULE qexml_module
