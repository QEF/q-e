!
! Copyright (C) 2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE xml_io_base
  !----------------------------------------------------------------------------
  !
  ! ... this module contains subroutines to read and write in XML format 
  ! ... the data produced by Quantum-ESPRESSO package
  !
  ! ... written by Carlo Sbraccia (2005)
  !
  USE iotk_module
  !
  USE kinds,     ONLY : DP
  USE io_files,  ONLY : tmp_dir, prefix, iunpun, xmlpun
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast  
  !
  IMPLICIT NONE
  !
  CHARACTER(iotk_attlenx) :: attr
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE create_directory( dirname )
      !------------------------------------------------------------------------
      !
      USE mp,        ONLY : mp_barrier
      USE mp_global, ONLY : mpime
      USE parser,    ONLY : int_to_char
      !
      CHARACTER(LEN=256) :: dirname
      INTEGER            :: ierr, ik
      !
      INTEGER,  EXTERNAL :: c_mkdir
      !
      !
      IF ( ionode ) &
         ierr = c_mkdir( TRIM( dirname ), LEN_TRIM( dirname ) )
      !
      IF ( ierr > 0 ) &
        CALL errore( 'create_directory', &
                     'unable to create directory ' // TRIM( dirname ), ierr )
      !
      ! ... all jobs are syncronized
      !
      CALL mp_barrier()
      !
      ! ... each job checks whether the scratch directory is accessible
      ! ... or not
      !
      OPEN( UNIT = 4, FILE = TRIM( dirname ) // '/test' // &
            TRIM( int_to_char( mpime ) ), STATUS = 'UNKNOWN', IOSTAT = ierr )
      CLOSE( UNIT = 4, STATUS = 'DELETE' )
      !
      IF ( ierr /= 0 ) &
         CALL errore( ' create_directory: ', TRIM( dirname ) // &
                    & ' non existent or non writable', 1 )
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
    END FUNCTION
    !
    !------------------------------------------------------------------------
    SUBROUTINE copy_file( file_in, file_out )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=*), INTENT(IN) :: file_in, file_out
      CHARACTER(LEN=256)           :: string
      INTEGER                      :: ierr
      !
      !
      OPEN( UNIT = 776, FILE = file_in,  STATUS = "OLD" )
      OPEN( UNIT = 777, FILE = file_out, STATUS = "UNKNOWN" )         
      !
      copy_loop: DO
         !
         READ( UNIT = 776, FMT = '(A256)', IOSTAT = ierr ) string
         !
         IF ( ierr < 0 ) EXIT copy_loop
         !
         WRITE( UNIT = 777, FMT = * ) TRIM( string )
         !
      END DO copy_loop
      !
      CLOSE( UNIT = 776 )
      CLOSE( UNIT = 777 )
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    FUNCTION restart_dir( scradir, runit )
      !------------------------------------------------------------------------
      !
      USE parser, ONLY: int_to_char
      !
      CHARACTER(LEN=256)           :: restart_dir
      CHARACTER(LEN=*), INTENT(IN) :: scradir
      INTEGER,          INTENT(IN) :: runit
      !
      CHARACTER(LEN=256) :: dirname
      INTEGER            :: strlen
      !
      ! ... main restart directory
      !
      dirname = 'RESTART' // int_to_char( runit )
      !
      IF ( LEN( scradir ) > 1 ) THEN
         !
         strlen = INDEX( scradir, ' ' ) - 1
         !
         dirname = scradir(1:strlen) // '/' // dirname
         !
      END IF
      !
      restart_dir = TRIM( dirname )
      !
    END FUNCTION restart_dir
    !
    ! ... writing
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_cell( ibrav, symm_type, &
                           celldm, alat, a1, a2, a3, b1, b2, b3 )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(IN) :: ibrav
      CHARACTER(LEN=*), INTENT(IN) :: symm_type
      REAL(KIND=DP),    INTENT(IN) :: celldm(6), alat
      REAL(KIND=DP),    INTENT(IN) :: a1(3), a2(3), a3(3)
      REAL(KIND=DP),    INTENT(IN) :: b1(3), b2(3), b3(3)
      !
      CHARACTER(LEN=256) :: bravais_lattice
      !
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
      CALL iotk_write_dat( iunpun, &
                           "BRAVAIS_LATTICE", TRIM( bravais_lattice ) )
      !
      IF ( ibrav == 0 ) &
         CALL iotk_write_dat( iunpun, "CELL_SYMMETRY", symm_type )
      !
      CALL iotk_write_attr( attr, "UNIT", "Bohr", FIRST = .TRUE. )
      CALL iotk_write_dat( iunpun, "LATTICE_PARAMETER", alat, ATTR = attr )
      !
      CALL iotk_write_dat( iunpun, "CELLDM", celldm(1:6) )
      !
      CALL iotk_write_begin( iunpun, "DIRECT_LATTICE_VECTORS" )
      CALL iotk_write_dat(   iunpun, "a1", a1(:) * alat, ATTR = attr )
      CALL iotk_write_dat(   iunpun, "a2", a2(:) * alat, ATTR = attr )
      CALL iotk_write_dat(   iunpun, "a3", a3(:) * alat, ATTR = attr )
      CALL iotk_write_end(   iunpun, "DIRECT_LATTICE_VECTORS" )
      !
      CALL iotk_write_attr( attr, "UNIT", "2 pi / a", FIRST = .TRUE. )
      CALL iotk_write_begin( iunpun, "RECIPROCAL_LATTICE_VECTORS" )
      CALL iotk_write_dat(   iunpun, "b1", b1(:), ATTR = attr )
      CALL iotk_write_dat(   iunpun, "b2", b2(:), ATTR = attr )
      CALL iotk_write_dat(   iunpun, "b3", b3(:), ATTR = attr )
      CALL iotk_write_end(   iunpun, "RECIPROCAL_LATTICE_VECTORS" )
      !
      CALL iotk_write_end( iunpun, "CELL" )
      !
    END SUBROUTINE write_cell
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_ions( nsp, nat, atm, ityp, psfile, &
                           pseudo_dir, amass, tau, if_pos, dirname )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(IN) :: nsp, nat
      INTEGER,          INTENT(IN) :: ityp(:)
      CHARACTER(LEN=*), INTENT(IN) :: atm(:)
      CHARACTER(LEN=*), INTENT(IN) :: psfile(:)
      CHARACTER(LEN=*), INTENT(IN) :: pseudo_dir
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      REAL(KIND=DP),    INTENT(IN) :: amass(:)
      REAL(KIND=DP),    INTENT(IN) :: tau(:,:)
      INTEGER,          INTENT(IN) :: if_pos(:,:)
      !
      INTEGER            :: i, flen
      CHARACTER(LEN=256) :: file_pseudo
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
         CALL copy_file( TRIM( file_pseudo ), &
                         TRIM( dirname ) // "/" // TRIM( psfile( i ) ) )
         !
         CALL iotk_write_attr( attr, "UNIT", "a.m.u.", FIRST = .TRUE. )
         CALL iotk_write_dat( iunpun, TRIM( atm( i ) )//"_MASS", &
                              amass(i), ATTR = attr )
         !
         CALL iotk_link( iunpun, "PSEUDO_FOR_" // TRIM( atm( i ) ), &
                         TRIM( psfile(i) ), CREATE = .FALSE.,       &
                         BINARY = .FALSE., RAW = .TRUE. )
         !
      END DO
      !
      CALL iotk_write_attr( attr, "UNIT", "Bohr", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_ATOMIC_POSITIONS", attr )
      !
      DO i = 1, nat
         !
         CALL iotk_write_attr( attr, "SPECIES", &
                             & atm( ityp(i) ), FIRST = .TRUE. )
         CALL iotk_write_attr( attr, "INDEX",  ityp(i) )                     
         CALL iotk_write_attr( attr, "tau", tau(:,i) )
         CALL iotk_write_attr( attr, "if_pos", if_pos(:,i) )
         CALL iotk_write_empty( iunpun, &
                              & "ATOM" // TRIM( iotk_index(i) ), attr )
         !
      END DO
      !
      CALL iotk_write_end( iunpun, "IONS" )
      !
    END SUBROUTINE write_ions
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_symmetry( ibrav, symm_type, nsym, &
                               invsym, nr1, nr2, nr3, ftau, s, sname )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(IN) :: ibrav, nsym,  nr1, nr2, nr3
      CHARACTER(LEN=*), INTENT(IN) :: symm_type
      LOGICAL,          INTENT(IN) :: invsym
      INTEGER,          INTENT(IN) :: s(:,:,:), ftau(:,:)
      CHARACTER(LEN=*), INTENT(IN) :: sname(:)
      !
      INTEGER        :: i
      REAL (KIND=DP) :: tmp(3)
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
         CALL iotk_write_attr( attr, "UNIT", "Crystal", FIRST = .TRUE. )
         !
         tmp(1) = ftau(1,i) / DBLE( nr1 )
         tmp(2) = ftau(2,i) / DBLE( nr2 )
         tmp(3) = ftau(3,i) / DBLE( nr3 )
         !
         CALL iotk_write_attr( attr, "ROT",        s(:,:,i) )
         CALL iotk_write_attr( attr, "FRAC_TRANS", tmp(:) )
         CALL iotk_write_attr( attr, "NAME",       sname(i) )
         !
         CALL iotk_write_empty( iunpun, &
                                "SYMM" // TRIM( iotk_index(i) ), attr )
         !
      END DO
      !
      CALL iotk_write_end( iunpun, "SYMMETRIES" )
      !
    END SUBROUTINE write_symmetry
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_planewaves( ecutwfc, dual, npwx, gamma_only, nr1, nr2, &
                                 nr3, ngm_g, nr1s, nr2s, nr3s, ngms_g, nr1b, &
                                 nr2b, nr3b, itmp, lgvec )
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : e2
      !
      INTEGER,       INTENT(IN) :: npwx, nr1, nr2, nr3, ngm_g, &
                                   nr1s, nr2s, nr3s, ngms_g, nr1b, nr2b, nr3b
      INTEGER,       INTENT(IN) :: itmp(:,:)
      REAL(KIND=DP), INTENT(IN) :: ecutwfc, dual
      LOGICAL,       INTENT(IN) :: gamma_only, lgvec
      !
      !
      CALL iotk_write_begin( iunpun, "PLANE_WAVES" )
      !
      CALL iotk_write_attr( attr, "UNIT", "Hartree", FIRST = .TRUE. )
      !
      CALL iotk_write_dat( iunpun, "WFC_CUTOFF", ecutwfc / e2, ATTR = attr )
      !
      CALL iotk_write_dat( iunpun, "RHO_CUTOFF", &
                           ecutwfc * dual / e2, ATTR = attr )
      !
      CALL iotk_write_dat( iunpun, "MAX_NPW", npwx )
      !
      CALL iotk_write_dat( iunpun, "GAMMA_ONLY", gamma_only )
      !
      CALL iotk_write_attr( attr, "nr1", nr1, FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "nr2", nr2 )
      CALL iotk_write_attr( attr, "nr3", nr3 )
      CALL iotk_write_empty( iunpun, "FFT_GRID", attr )
      !
      CALL iotk_write_dat( iunpun, "GVECT_NUMBER", ngm_g )
      !
      CALL iotk_write_attr( attr, "nr1s", nr1s, FIRST = .TRUE. )
      CALL iotk_write_attr( attr, "nr2s", nr2s )
      CALL iotk_write_attr( attr, "nr3s", nr3s )
      CALL iotk_write_empty( iunpun, "SMOOTH_FFT_GRID", attr )
      !
      CALL iotk_write_dat( iunpun, "SMOOTH_GVECT_NUMBER", ngms_g )
      !
      IF ( lgvec ) THEN
         !
         ! ... write the G-vectors
         !
         CALL iotk_link( iunpun, "G-VECTORS_FILE", "gvectors.dat", &
                         CREATE = .TRUE., BINARY = .TRUE., RAW = .TRUE. )
         !
         CALL iotk_write_begin( iunpun, "G-VECTORS", attr = attr )
         CALL iotk_write_dat( iunpun, "g", itmp(1:3,1:ngm_g), FMT = "(3I5)" )
         CALL iotk_write_end( iunpun, "G-VECTORS" )
         !
      END IF
      !
      CALL iotk_write_attr( attr, "nr1b", nr1b , first = .TRUE. )
      CALL iotk_write_attr( attr, "nr2b", nr2b )
      CALL iotk_write_attr( attr, "nr3b", nr3b )
      CALL iotk_write_empty( iunpun, "SMALLBOX_FFT_GRID", attr )
      !
      CALL iotk_write_end( iunpun, "PLANE_WAVES" )
      !
    END SUBROUTINE write_planewaves
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_spin( lsda, noncolin, npol, lspinorb )
      !------------------------------------------------------------------------
      !
      LOGICAL, INTENT(IN) :: lsda, noncolin, lspinorb
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
      !
      CALL iotk_write_end( iunpun, "SPIN" )
      !
    END SUBROUTINE write_spin
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_xc( dft, nsp, lda_plus_u, &
                         Hubbard_lmax, Hubbard_l, Hubbard_U, Hubbard_alpha )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=*),        INTENT(IN) :: dft
      LOGICAL,                 INTENT(IN) :: lda_plus_u
      INTEGER,                 INTENT(IN) :: nsp
      INTEGER,       OPTIONAL, INTENT(IN) :: Hubbard_lmax
      INTEGER,       OPTIONAL, INTENT(IN) :: Hubbard_l(:)
      REAL(KIND=DP), OPTIONAL, INTENT(IN) :: Hubbard_U(:), Hubbard_alpha(:)
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
    END SUBROUTINE write_xc
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_occ( lgauss, ngauss, degauss, ltetra, &
                          ntetra, tfixed_occ, lsda, nelup, neldw, f_inp )
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : e2
      !
      LOGICAL,                 INTENT(IN) :: lgauss, ltetra, tfixed_occ, lsda
      INTEGER,       OPTIONAL, INTENT(IN) :: ngauss, ntetra, nelup, neldw
      REAL(KIND=DP), OPTIONAL, INTENT(IN) :: degauss, f_inp(:,:)      
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
      IF ( ltetra ) &
         CALL iotk_write_dat( iunpun, "NUMBER_OF_TETRAHEDRA", ntetra )
      !
      CALL iotk_write_dat( iunpun, "FIXED_OCCUPATIONS", tfixed_occ )
      !
      IF ( tfixed_occ ) THEN
         !
         CALL iotk_write_dat( iunpun, "INPUT_OCC_UP", f_inp(1:nelup,1) )
         !
         IF ( lsda ) &
            CALL iotk_write_dat( iunpun, "INPUT_OCC_DOWN", f_inp(1:neldw,2) )
         !
      END IF
      !
      CALL iotk_write_end( iunpun, "OCCUPATIONS" )
      !
    END SUBROUTINE write_occ
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_bz( num_k_points, xk, wk )
      !------------------------------------------------------------------------
      !
      INTEGER,       INTENT(IN) :: num_k_points
      REAL(KIND=DP), INTENT(IN) :: xk(:,:), wk(:)
      !
      INTEGER :: ik
      !
      !
      CALL iotk_write_begin( iunpun, "BRILLOUIN_ZONE" )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_K-POINTS", num_k_points )
      !
      CALL iotk_write_attr( attr, "UNIT", "2 pi / a", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_K-POINTS", attr )
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
    END SUBROUTINE write_bz
    !
END MODULE xml_io_base
