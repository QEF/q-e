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
  ! ... this module contains some common subroutines used to read and write
  ! ... in XML format the data produced by Quantum-ESPRESSO package
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
  PRIVATE
  !
  PUBLIC :: attr
  PUBLIC :: create_directory, kpoint_dir, wfc_filename, copy_file,          &
            restart_dir, check_restartfile, save_history,             &
            save_print_counter, read_print_counter, set_kpoints_vars, &
            write_cell, write_ions, write_symmetry, write_planewaves,       &
            write_spin, write_xc, write_occ, write_bz, write_phonon,        &
            write_rho_xml, write_wfc, read_wfc, read_rho_xml
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
      !
      CHARACTER(LEN=*), INTENT(IN) :: dirname
      !
      INTEGER           :: ierr
      INTEGER, EXTERNAL :: c_mkdir
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      !
      IF ( ionode ) THEN
         ierr = c_mkdir( TRIM( dirname ), LEN_TRIM( dirname ) )
      ELSE
         ierr = 0
      END IF
      !
      CALL mp_bcast ( ierr, ionode_id )
      !
      CALL errore( 'create_directory', &
                   'unable to create directory ' // TRIM( dirname ), ierr )
      !
      ! ... syncronize all jobs (not sure it is really useful)
      !
      CALL mp_barrier()
      !
      ! ... check whether the scratch directory is writable
      !
      IF ( ionode ) THEN
         OPEN( UNIT = 4, FILE = TRIM( dirname ) // '/test' // &
              TRIM( int_to_char( mpime ) ), STATUS = 'UNKNOWN', IOSTAT = ierr )
         CLOSE( UNIT = 4, STATUS = 'DELETE' )
      END IF
      !
      CALL mp_bcast ( ierr, ionode_id )
      !
      CALL errore( 'create_directory: ', &
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
    FUNCTION restart_dir( scradir, runit )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=256)           :: restart_dir
      CHARACTER(LEN=*), INTENT(IN) :: scradir
      INTEGER,          INTENT(IN) :: runit
      !
      CHARACTER(LEN=256) :: dirname
      INTEGER            :: strlen
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      !
      ! ... main restart directory
      !
      ! ... keep the line below ( this is the old style RESTARTXX ) !!!
      !
      dirname = 'RESTART' // int_to_char( runit )
      !
!      dirname = TRIM( prefix ) // '_' // TRIM( int_to_char( runit ) )// '.save'
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
      RETURN
      !
    END FUNCTION restart_dir
    !
    !------------------------------------------------------------------------
    FUNCTION check_restartfile( scradir, ndr )
      !------------------------------------------------------------------------
      !
      USE io_global, ONLY : ionode, ionode_id
      !
      IMPLICIT NONE
      !
      LOGICAL                      :: check_restartfile
      INTEGER,          INTENT(IN) :: ndr
      CHARACTER(LEN=*), INTENT(IN) :: scradir
      CHARACTER(LEN=256)           :: filename
      LOGICAL                      :: lval
      INTEGER                      :: strlen
      !
      !
      filename = restart_dir( scradir, ndr )
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( filename ) // '/' // TRIM( xmlpun )
         !
         INQUIRE( FILE = TRIM( filename ), EXIST = lval )
         !
      END IF
      !
      CALL mp_bcast( lval, ionode_id )
      !
      check_restartfile = lval
      !
      RETURN
      !
    END FUNCTION check_restartfile
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
      CHARACTER(LEN=256) :: filename
      CHARACTER(LEN=6)   :: hindex
      !
#if defined (__VERBOSE_SAVE)
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
    SUBROUTINE save_print_counter( iter, scradir, wunit )
      !------------------------------------------------------------------------
      !
      ! ... a counter indicating the last successful pritout iteration is saved
      !
      USE io_global, ONLY : ionode, ionode_id
      USE mp,        ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(IN) :: iter
      CHARACTER(LEN=*), INTENT(IN) :: scradir
      INTEGER,          INTENT(IN) :: wunit
      !
      INTEGER            :: ierr
      CHARACTER(LEN=256) :: filename, dirname
      !
      !
      dirname = restart_dir( scradir, wunit )
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
      CALL mp_bcast( ierr, ionode_id )
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
    SUBROUTINE read_print_counter( nprint_nfi, scradir, runit )
      !------------------------------------------------------------------------
      !
      ! ... the counter indicating the last successful pritout iteration 
      ! ... is read here
      !
      USE io_global, ONLY : ionode, ionode_id
      USE mp,        ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: nprint_nfi
      CHARACTER(LEN=*), INTENT(IN)  :: scradir
      INTEGER,          INTENT(IN)  :: runit
      !
      INTEGER            :: ierr
      CHARACTER(LEN=256) :: filename, dirname
      !
      !
      dirname = restart_dir( scradir, runit )
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
      CALL mp_bcast( nprint_nfi, ionode_id )
      !
      RETURN
      !
    END SUBROUTINE read_print_counter   
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_kpoints_vars( ik, nk, kunit, ngwl, igl, &
                                 npool, ikt, iks, ike, igwx, ipmask, ipsour )
      !------------------------------------------------------------------------
      !
      ! ... set working variables for k-point index (ikt) and 
      ! ... k-points number (nkt)
      !
      USE mp,         ONLY : mp_sum, mp_get, mp_max
      USE mp_global,  ONLY : mpime, nproc, me_pool, my_pool_id, &
                             nproc_pool, intra_pool_comm, root_pool
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN)  :: ik, nk, kunit
      INTEGER, INTENT(IN)  :: ngwl, igl(:)
      INTEGER, INTENT(OUT) :: npool
      INTEGER, INTENT(OUT) :: ikt, iks, ike, igwx
      INTEGER, INTENT(OUT) :: ipmask(:), ipsour
      !
      INTEGER :: ierr, i
      INTEGER :: nkl, nkr, nkbl, nkt
      !
      !
      ikt = ik
      nkt = nk
      !
      ! ... find out the number of pools
      !
      npool = nproc / nproc_pool 
      !
      ! ... find out number of k points blocks
      !
      nkbl = nkt / kunit  
      !
      ! ... k points per pool
      !
      nkl = kunit * ( nkbl / npool )
      !
      ! ... find out the reminder
      !
      nkr = ( nkt - nkl * npool ) / kunit
      !
      ! ... Assign the reminder to the first nkr pools
      !
      IF ( my_pool_id < nkr ) nkl = nkl + kunit
      !
      ! ... find out the index of the first k point in this pool
      !
      iks = nkl * my_pool_id + 1
      !
      IF ( my_pool_id >= nkr ) iks = iks + nkr * kunit
      !
      ! ... find out the index of the last k point in this pool
      !
      ike = iks + nkl - 1
      !
      ipmask = 0
      ipsour = ionode_id
      !
      ! ... find out the index of the processor which collect the data 
      ! ... in the pool of ik
      !
      IF ( npool > 1 ) THEN
         !
         IF ( ( ikt >= iks ) .AND. ( ikt <= ike ) ) THEN
            !
            IF ( me_pool == root_pool ) ipmask( mpime + 1 ) = 1
            !
         END IF
         !
         CALL mp_sum( ipmask )
         !
         DO i = 1, nproc
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
      CALL mp_max( igwx, intra_pool_comm ) 
      !
      ! ... now notify all procs if an error has been found 
      !
      CALL mp_max( ierr ) 
      !
      CALL errore( 'set_kpoint_vars ', 'wrong size ngl', ierr )
      !
      IF ( ipsour /= ionode_id ) &
         CALL mp_get( igwx, igwx, mpime, ionode_id, ipsour, 1 )
      !
      RETURN
      !
    END SUBROUTINE set_kpoints_vars
    !
    ! ... writing subroutines
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_cell( ibrav, symm_type, &
                           celldm, alat, a1, a2, a3, b1, b2, b3 )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(IN) :: ibrav
      CHARACTER(LEN=*), INTENT(IN) :: symm_type
      REAL(DP),         INTENT(IN) :: celldm(6), alat
      REAL(DP),         INTENT(IN) :: a1(3), a2(3), a3(3)
      REAL(DP),         INTENT(IN) :: b1(3), b2(3), b3(3)
      !
      CHARACTER(LEN=256) :: bravais_lattice
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
      CALL iotk_write_dat( iunpun, "CELL_SYMMETRY", symm_type )
      !
      CALL iotk_write_attr( attr, "UNITS", "Bohr", FIRST = .TRUE. )
      CALL iotk_write_dat( iunpun, "LATTICE_PARAMETER", alat, ATTR = attr )
      !
      CALL iotk_write_dat( iunpun, "CELL_DIMENSIONS", celldm(1:6) )
      !
      CALL iotk_write_begin( iunpun, "DIRECT_LATTICE_VECTORS" )
      CALL iotk_write_dat(   iunpun, "a1", a1(:) * alat, ATTR = attr )
      CALL iotk_write_dat(   iunpun, "a2", a2(:) * alat, ATTR = attr )
      CALL iotk_write_dat(   iunpun, "a3", a3(:) * alat, ATTR = attr )
      CALL iotk_write_end(   iunpun, "DIRECT_LATTICE_VECTORS" )
      !
      CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
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
    END SUBROUTINE write_ions
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_symmetry( ibrav, symm_type, nsym, &
                               invsym, nr1, nr2, nr3, ftau, s, sname, irt )
      !------------------------------------------------------------------------
      !
      INTEGER,          INTENT(IN) :: ibrav, nsym,  nr1, nr2, nr3
      CHARACTER(LEN=*), INTENT(IN) :: symm_type
      LOGICAL,          INTENT(IN) :: invsym
      INTEGER,          INTENT(IN) :: s(:,:,:), ftau(:,:)
      CHARACTER(LEN=*), INTENT(IN) :: sname(:)
      INTEGER,          INTENT(IN) :: irt(:,:)
      !
      INTEGER  :: i
      REAL(DP) :: tmp(3)
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
      INTEGER,  INTENT(IN) :: npwx, nr1, nr2, nr3, ngm_g, &
                              nr1s, nr2s, nr3s, ngms_g, nr1b, nr2b, nr3b
      INTEGER,  INTENT(IN) :: itmp(:,:)
      REAL(DP), INTENT(IN) :: ecutwfc, dual
      LOGICAL,  INTENT(IN) :: gamma_only, lgvec
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
      CALL iotk_write_dat( iunpun, "MAX_NPW", npwx )
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
         CALL iotk_link( iunpun, "G-VECTORS_FILE", &
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
      CHARACTER(LEN=*),   INTENT(IN) :: dft
      LOGICAL,            INTENT(IN) :: lda_plus_u
      INTEGER,            INTENT(IN) :: nsp
      INTEGER,  OPTIONAL, INTENT(IN) :: Hubbard_lmax
      INTEGER,  OPTIONAL, INTENT(IN) :: Hubbard_l(:)
      REAL(DP), OPTIONAL, INTENT(IN) :: Hubbard_U(:), Hubbard_alpha(:)
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
    SUBROUTINE write_occ( lgauss, ngauss, degauss, ltetra, ntetra, &
                          tetra, tfixed_occ, lsda, nelup, neldw, f_inp )
      !------------------------------------------------------------------------
      !
      USE constants, ONLY : e2
      !
      LOGICAL,            INTENT(IN) :: lgauss, ltetra, tfixed_occ, lsda
      INTEGER,  OPTIONAL, INTENT(IN) :: ngauss, ntetra, nelup, neldw
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
    SUBROUTINE write_bz( num_k_points, xk, wk, k1, k2, k3, nk1, nk2, nk3 )
      !------------------------------------------------------------------------
      !
      INTEGER,  INTENT(IN) :: num_k_points, k1, k2, k3, nk1, nk2, nk3
      REAL(DP), INTENT(IN) :: xk(:,:), wk(:)
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
    END SUBROUTINE write_bz
    !
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
      CALL iotk_write_dat( iunpun, "Q-POINT", xqq(:) )
      !
      CALL iotk_write_end( iunpun, "PHONON" )
      !
    END SUBROUTINE write_phonon
    !
    ! ... methods to write and read charge_density
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_rho_xml( rho_file_base, rho, &
                              nr1, nr2, nr3, nr1x, nr2x, ipp, npp )
      !------------------------------------------------------------------------
      !
      ! ... Writes charge density rho, one plane at a time.
      ! ... If ipp and npp are specified, planes are collected one by one from
      ! ... all processors, avoiding an overall collect of the charge density
      ! ... on a single proc.
      !
      USE io_files,  ONLY : rhounit
      USE io_global, ONLY : ionode, ionode_id
      USE mp,        ONLY : mp_sum, mp_get, mp_max
      USE mp_global, ONLY : mpime, nproc, root, me_pool, my_pool_id, &
                            nproc_pool, intra_pool_comm, root_pool, my_image_id
      !
      IMPLICIT NONE
      !
      INTEGER,           INTENT(IN) :: nr1, nr2, nr3
      INTEGER,           INTENT(IN) :: nr1x, nr2x
      CHARACTER(LEN=*),  INTENT(IN) :: rho_file_base
      REAL(DP),          INTENT(IN) :: rho(:)
      INTEGER, OPTIONAL, INTENT(IN) :: ipp(:)
      INTEGER, OPTIONAL, INTENT(IN) :: npp(:)
      !
      INTEGER               :: ierr, i, j, k, kk, ldr, ip
      CHARACTER(LEN=256)    :: rho_file
      REAL(DP), ALLOCATABLE :: rho_plane(:)
      INTEGER,  ALLOCATABLE :: kowner(:)
      !
      !
      rho_file = TRIM( rho_file_base ) // '.xml'
      !
      IF ( ionode ) &
         CALL iotk_open_write( rhounit, &
                               FILE = rho_file, BINARY = .FALSE., IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id )
      !
      CALL errore( 'write_rho_xml ', 'cannot open' // &
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
      ALLOCATE( rho_plane( nr1 * nr2 ) )
      ALLOCATE( kowner( nr3 ) )
      !
      ! ... find out the owner of each "z" plane
      !
      IF ( PRESENT( ipp ) .AND. PRESENT( npp ) ) THEN
         !
         DO ip = 1, nproc
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
      ldr = nr1x * nr2x
      !
      DO k = 1, nr3
         !
         IF( kowner(k) == mpime ) THEN
            !
            kk = k
            !
            IF ( PRESENT( ipp ) ) kk = k - ipp(mpime+1)
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
         IF ( kowner( k ) /= ionode_id ) &
            CALL mp_get( rho_plane, rho_plane, mpime, ionode_id, kowner(k), k )
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
    SUBROUTINE read_rho_xml( rho_file_base, rho, &
                             nr1, nr2, nr3, nr1x, nr2x, ipp, npp )
      !------------------------------------------------------------------------
      !
      ! ... Writes charge density rho, one plane at a time.
      ! ... If ipp and npp are specified, planes are collected one by one from
      ! ... all processors, avoiding an overall collect of the charge density
      ! ... on a single proc.
      !
      USE io_files,  ONLY : rhounit
      USE io_global, ONLY : ionode, ionode_id
      USE mp,        ONLY : mp_sum, mp_put, mp_max
      USE mp_global, ONLY : mpime, nproc, root, me_pool, my_pool_id, &
                            nproc_pool, intra_pool_comm, root_pool, my_image_id
      !
      IMPLICIT NONE
      !
      INTEGER,           INTENT(IN)  :: nr1, nr2, nr3
      INTEGER,           INTENT(IN)  :: nr1x, nr2x
      CHARACTER(LEN=*),  INTENT(IN)  :: rho_file_base
      REAL(DP),          INTENT(OUT) :: rho(:)
      INTEGER, OPTIONAL, INTENT(IN)  :: ipp(:)
      INTEGER, OPTIONAL, INTENT(IN)  :: npp(:)
      !
      INTEGER               :: ierr, i, j, k, kk, ldr, ip
      INTEGER               :: nr( 3 )
      CHARACTER(LEN=256)    :: rho_file
      REAL(DP), ALLOCATABLE :: rho_plane(:)
      INTEGER,  ALLOCATABLE :: kowner(:)
      !
      !
      rho_file = TRIM( rho_file_base ) // '.xml'
      !
      IF ( ionode ) &
         CALL iotk_open_read( rhounit, FILE = TRIM( rho_file ), &
                               BINARY = .FALSE., IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id )
      !
      CALL errore( ' read_rho_xml ', &
                   'cannot open ' // rho_file // ' file for reading', ierr )
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( rhounit, "CHARGE-DENSITY" )
         !
         CALL iotk_scan_empty( rhounit, "INFO", attr )
         !
         CALL iotk_scan_attr( attr, "nr1",   nr( 1 ) )
         CALL iotk_scan_attr( attr, "nr2",   nr( 2 ) )
         CALL iotk_scan_attr( attr, "nr3",   nr( 3 ) )
         !
      END IF
      !
      CALL mp_bcast( nr, ionode_id )
      !
      IF ( nr1 /= nr( 1 ) .OR. nr2 /= nr( 2 ) .OR. nr3 /= nr( 3 ) ) &
         CALL errore( ' read_rho_xml ', ' dimensions do not match ', 1 )
      !
      ALLOCATE( rho_plane( nr1 * nr2 ) )
      ALLOCATE( kowner( nr3 ) )
      !
      ! ... find out the owner of each "z" plane
      !
      IF ( PRESENT( ipp ) .AND. PRESENT( npp ) ) THEN
         !
         DO ip = 1, nproc
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
      ldr = nr1x * nr2x
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
         IF ( kowner( k ) /= ionode_id ) &
            CALL mp_put( rho_plane, rho_plane, mpime, ionode_id, kowner(k), k )
         !
         IF( kowner(k) == mpime ) THEN
            !
            kk = k
            !
            IF ( PRESENT( ipp ) ) kk = k - ipp(mpime+1)
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
      RETURN
      !
    END SUBROUTINE read_rho_xml
    !
    !
    ! ... methods to write and read wavefunctions
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_wfc( iuni, ik, nk, kunit, ispin, &
                          nspin, wf0, ngw, nbnd, igl, ngwl, filename, scalef )
      !------------------------------------------------------------------------
      !
      USE mp_wave,    ONLY : mergewf
      USE mp,         ONLY : mp_get
      USE mp_global,  ONLY : me_pool, nproc, nproc_pool, &
                             root_pool, intra_pool_comm, mpime
      !
      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN) :: iuni
      INTEGER,            INTENT(IN) :: ik, nk, kunit, ispin, nspin
      COMPLEX(DP),        INTENT(IN) :: wf0(:,:)
      INTEGER,            INTENT(IN) :: ngw
      INTEGER,            INTENT(IN) :: nbnd
      INTEGER,            INTENT(IN) :: ngwl
      INTEGER,            INTENT(IN) :: igl(:)
      CHARACTER(LEN=256), INTENT(IN) :: filename
      REAL(DP),           INTENT(IN) :: scalef    
        ! scale factor, usually 1.0 for pw and 1/SQRT( omega ) for CP
      !
      INTEGER                  :: i, j, ierr
      INTEGER                  :: iks, ike, ikt, igwx
      INTEGER                  :: npool, ipmask(nproc), ipsour
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
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
               CALL mp_get( wtmp, wtmp, mpime, ionode_id, ipsour, j )
            !
         ELSE
            !
            CALL mergewf( wf0(:,j), wtmp, ngwl, igl, mpime, nproc, ionode_id )
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
      RETURN
      !
    END SUBROUTINE write_wfc
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_wfc( iuni, ik, nk, kunit, ispin, &
                         nspin, wf, ngw, nbnd, igl, ngwl, filename, scalef )
      !------------------------------------------------------------------------
      !
      USE mp_wave,   ONLY : splitwf
      USE mp,        ONLY : mp_get, mp_put
      USE mp_global, ONLY : mpime, nproc, root, me_pool, my_pool_id, &
                            nproc_pool, intra_pool_comm, root_pool, my_image_id
      !
      IMPLICIT NONE
      !
      INTEGER,            INTENT(IN)    :: iuni
      COMPLEX(DP),        INTENT(OUT)   :: wf(:,:)
      INTEGER,            INTENT(IN)    :: ik, nk
      INTEGER,            INTENT(INOUT) :: kunit
      INTEGER,            INTENT(INOUT) :: ngw, nbnd, ispin, nspin
      INTEGER,            INTENT(IN)    :: ngwl
      INTEGER,            INTENT(IN)    :: igl(:)
      CHARACTER(LEN=256), INTENT(IN)    :: filename
      REAL(DP),           INTENT(OUT)   :: scalef
      !
      INTEGER                  :: i, j
      COMPLEX(DP), ALLOCATABLE :: wtmp(:)
      INTEGER                  :: ierr
      INTEGER                  :: iks, ike, nkt, ikt
      INTEGER                  :: igwx, igwx_, ik_, nk_, kunit_
      INTEGER                  :: npool, ipmask(nproc), ipdest
      !
      !
      CALL set_kpoints_vars( ik, nk, kunit, ngwl, igl, &
                             npool, ikt, iks, ike, igwx, ipmask, ipdest )
      !
      IF ( ionode ) &
         CALL iotk_open_read( iuni, FILE = TRIM( filename ), &
                              BINARY = .TRUE., IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id )
      !
      CALL errore( 'read_wfc ', &
                   'cannot open restart file for reading', ierr )
      !
      IF ( ionode ) THEN
          !
          CALL iotk_scan_empty( iuni, "INFO", attr )
          !
          CALL iotk_scan_attr( attr, "ngw",           ngw )
          CALL iotk_scan_attr( attr, "nbnd",          nbnd )
          CALL iotk_scan_attr( attr, "ik",            ik_ )
          CALL iotk_scan_attr( attr, "nk",            nk_ )
          CALL iotk_scan_attr( attr, "kunit",         kunit_ )
          CALL iotk_scan_attr( attr, "ispin",         ispin )
          CALL iotk_scan_attr( attr, "nspin",         nspin )
          CALL iotk_scan_attr( attr, "igwx",          igwx_ )
          CALL iotk_scan_attr( attr, "scale_factor",  scalef )
          !
      END IF
      !
      CALL mp_bcast( ngw,    ionode_id )
      CALL mp_bcast( nbnd,   ionode_id )
      CALL mp_bcast( ik_,    ionode_id )
      CALL mp_bcast( nk_,    ionode_id )
      CALL mp_bcast( kunit_, ionode_id )
      CALL mp_bcast( ispin,  ionode_id )
      CALL mp_bcast( nspin,  ionode_id )
      CALL mp_bcast( igwx_,  ionode_id )
      CALL mp_bcast( scalef, ionode_id )
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
               IF ( igwx > igwx_ ) wtmp((igwx_+1):igwx) = 0.D0
               !
            END IF
            !
            IF ( npool > 1 ) THEN
               
               IF ( ipdest /= ionode_id ) &
                  CALL mp_put( wtmp, wtmp, mpime, ionode_id, ipdest, j )
               !
               IF ( ( ikt >= iks ) .AND. ( ikt <= ike ) ) &
                  CALL splitwf( wf(:,j), wtmp, ngwl, igl, me_pool, &
                                nproc_pool, root_pool, intra_pool_comm )
               !
            ELSE
               !
               CALL splitwf( wf(:,j), wtmp, ngwl, igl, mpime, nproc, ionode_id )
               !
            END IF
            !
         END IF
         !
      END DO
      !
      IF ( ionode ) CALL iotk_close_read( iuni )
      !
      DEALLOCATE( wtmp )
      !
      RETURN
      !
    END SUBROUTINE read_wfc
    !
END MODULE xml_io_base
