!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
MODULE read_cards_module
  !---------------------------------------------------------------------------
  !
  ! ...  This module handles the reading of cards from standard input
  ! ...  Written by Carlo Cavazzoni and modified for "path" implementation
  ! ...  by Carlo Sbraccia
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout 
  USE constants, ONLY : angstrom_au
  USE parser,    ONLY : field_count, read_line
  USE io_global, ONLY : ionode, ionode_id
  !
  USE input_parameters
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: read_cards
  !
  ! ... end of module-scope declarations
  !
  !  ----------------------------------------------
  !
  CONTAINS
     !
     ! ... Read CARDS ....
     !
     ! ... subroutines
     !
     !----------------------------------------------------------------------
     SUBROUTINE card_default_values( prog )
       !----------------------------------------------------------------------
       !
       USE autopilot, ONLY : init_autopilot
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2) :: prog
       !
       !
       ! ... mask that control the printing of selected Kohn-Sham occupied 
       ! ... orbitals, default allocation
       !
       CALL allocate_input_iprnks( 0, nspin )
       nprnks  = 0
       !
       ! ... mask that control the printing of selected Kohn-Sham unoccupied 
       ! ... orbitals, default allocation
       !
       CALL allocate_input_iprnks_empty( 0, nspin )
       nprnks_empty  = 0
       !
       ! ... Simulation cell from standard input
       !
       trd_ht = .FALSE.
       rd_ht  = 0.0_DP
       !
       ! ... dipole
       !
       tdipole_card = .FALSE.
       !
       ! ... Constraints
       !
       nconstr_inp    = 0
       constr_tol_inp = 1.E-6_DP
       !
       ! ... ionic mass initialization
       !
       atom_mass = 0.0_DP
       !
       ! ... dimension of the real space Ewald summation
       !
       iesr_inp = 1
       !
       ! ... k-points
       !
       k_points = 'gamma'
       tk_inp   = .FALSE.
       nkstot   = 1
       nk1      = 0
       nk2      = 0
       nk3      = 0
       k1       = 0
       k2       = 0
       k3       = 0
       !
       ! ... neighbours
       !
       tneighbo   =  .FALSE.
       neighbo_radius =  0.0_DP
       !
       ! ... Turbo
       !
       tturbo_inp = .FALSE.
       nturbo_inp = 0
       !
       ! ... Grids
       !
       t2dpegrid_inp = .FALSE.
       !
       ! ... Electronic states
       !
       tf_inp = .FALSE.
       !
       ! ... Hartree planar mean
       !
       tvhmean_inp = .false.
       vhnr_inp    = 0
       vhiunit_inp = 0
       vhrmin_inp  = 0.0_DP
       vhrmax_inp  = 0.0_DP
       vhasse_inp  = 'K'
       !
       ! ... tchi
       !
       tchi2_inp = .FALSE.
       !
       ! ... ion_velocities
       !
       tavel = .FALSE.
       !
       ! ... setnfi
       !
       newnfi_card  = -1
       tnewnfi_card = .FALSE.
       !
       CALL init_autopilot()
       !
       RETURN
       !
     END SUBROUTINE card_default_values
     !
     !
     !----------------------------------------------------------------------
     SUBROUTINE read_cards( prog )
       !----------------------------------------------------------------------
       !
       USE autopilot, ONLY : card_autopilot
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2)           :: prog   ! calling program ( FP, PW, CP )
       CHARACTER(LEN=256)         :: input_line
       CHARACTER(LEN=80)          :: card
       CHARACTER(LEN=1), EXTERNAL :: capital
       LOGICAL                    :: tend
       INTEGER                    :: i
       !
       !
       CALL card_default_values( prog )
       !
 100   CALL read_line( input_line, end_of_file=tend )
       !
       IF( tend ) GO TO 120
       IF( input_line == ' ' .OR. input_line(1:1) == '#' ) GO TO 100
       !
       READ (input_line, *) card
       !
       DO i = 1, LEN_TRIM( input_line )
          input_line( i : i ) = capital( input_line( i : i ) )
       END DO
       !
       IF ( TRIM(card) == 'AUTOPILOT' ) THEN
          !
          CALL card_autopilot( input_line )
          !
       ELSE IF ( TRIM(card) == 'ATOMIC_SPECIES' ) THEN
          !
          CALL card_atomic_species( input_line, prog )
          !
       ELSE IF ( TRIM(card) == 'ATOMIC_POSITIONS' ) THEN
          !
          CALL card_atomic_positions( input_line, prog )
          !
       ELSE IF ( TRIM(card) == 'SETNFI' ) THEN
          !
          CALL card_setnfi( input_line )
          IF ( ( prog == 'PW' .OR. prog == 'CP' ) .AND. ionode ) &
             WRITE( stdout,'(A)') 'Warning: card '//trim(input_line)//' ignored'
          !
       ELSE IF ( TRIM(card) == 'CONSTRAINTS' ) THEN
          !
          CALL card_constraints( input_line )
          !
       ELSE IF ( TRIM(card) == 'COLLECTIVE_VARS' ) THEN
          !
          CALL card_collective_vars( input_line )
          !
       ELSE IF ( TRIM(card) == 'VHMEAN' ) THEN
          !
          CALL card_vhmean( input_line )
          IF ( ( prog == 'PW' .OR. prog == 'CP' ) .AND. ionode ) &
             WRITE( stdout,'(A)') 'Warning: card '//trim(input_line)//' ignored'
          !
       ELSE IF ( TRIM(card) == 'DIPOLE' ) THEN
          !
          CALL card_dipole( input_line )
          IF ( ( prog == 'PW' .OR. prog == 'CP' ) .AND. ionode ) &
             WRITE( stdout,'(A)') 'Warning: card '//trim(input_line)//' ignored'
          !
       ELSE IF ( TRIM(card) == 'ESR' ) THEN
          !
          CALL card_esr( input_line )
          IF ( ( prog == 'PW' .OR. prog == 'CP' ) .AND. ionode ) &
             WRITE( stdout,'(A)') 'Warning: card '//trim(input_line)//' ignored'
          !
       ELSE IF ( TRIM(card) == 'K_POINTS' ) THEN
          !
          IF ( prog == 'CP' ) THEN
             IF( ionode ) WRITE( stdout,'(A)') 'Warning: card '//trim(input_line)//' ignored'
          ELSE
             CALL card_kpoints( input_line )
          END IF
          !
       ELSE IF ( TRIM(card) == 'NEIGHBOURS' ) THEN
          !
          CALL card_neighbours( input_line )
          IF ( ( prog == 'PW' .OR. prog == 'CP' ) .AND. ionode ) &
             WRITE( stdout,'(A)') 'Warning: card '//trim(input_line)//' ignored'
          !
       ELSE IF ( TRIM(card) == 'OCCUPATIONS' ) THEN
          !
          CALL card_occupations( input_line )
          !
       ELSE IF ( TRIM(card) == 'CELL_PARAMETERS' ) THEN
          !
          CALL card_cell_parameters( input_line )
          !
       ELSE IF ( TRIM(card) == 'TURBO' ) THEN
          !
          CALL card_turbo( input_line )
          IF ( ( prog == 'PW' .OR. prog == 'CP' ) .AND. ionode ) &
             WRITE( stdout,'(A)') 'Warning: card '//trim(input_line)//' ignored'
          !
       ELSE IF ( TRIM(card) == 'ATOMIC_VELOCITIES' ) THEN
          !
          CALL card_ion_velocities( input_line )
          IF ( ( prog == 'PW' .OR. prog == 'CP' ) .AND. ionode ) &
             WRITE( stdout,'(A)') 'Warning: card '//trim(input_line)//' ignored'
          !
       ELSE IF ( TRIM(card) == 'KSOUT' ) THEN
          !
          CALL card_ksout( input_line )
          IF ( ( prog == 'PW' ) .AND. ionode ) &
             WRITE( stdout,'(a)') 'Warning: card '//trim(input_line)//' ignored'
          !
       ELSE IF ( TRIM(card) == 'KSOUT_EMPTY' ) THEN
          !
          CALL card_ksout_empty( input_line )
          IF ( ( prog == 'PW' ) .AND. ionode ) &
             WRITE( stdout,'(A)') 'Warning: card '//trim(input_line)//' ignored'
          !
       ELSE IF ( TRIM(card) == 'CLIMBING_IMAGES' ) THEN
          !
          CALL card_climbing_images( input_line )

       ELSE IF ( TRIM(card) == 'PLOT_WANNIER' ) THEN
          !
          CALL card_plot_wannier( input_line )

       ELSE
          !
          IF ( ionode ) &
             WRITE( stdout,'(A)') 'Warning: card '//TRIM(input_line)//' ignored'
          !
       END IF
       !
       ! ... END OF LOOP ... !
       !   
       GOTO 100
       !
120    CONTINUE
       !
       RETURN
       !
     END SUBROUTINE read_cards
     !
     !
     ! ... Description of the allowed input CARDS for FPMD code
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     !
     ! ATOMIC_SPECIES
     !
     !   set the atomic species been read and their pseudopotential file
     !
     ! Syntax:
     !
     !    ATOMIC_SPECIE
     !      label(1)    mass(1)    psfile(1)
     !       ...        ...        ...
     !      label(n)    mass(n)    psfile(n)
     !
     ! Example:
     !
     ! ATOMIC_SPECIES
     !  O 16.0 O.BLYP.UPF
     !  H 1.00 H.fpmd.UPF
     !
     ! Where:
     !
     !      label(i)  ( character(len=4) )  label of the atomic species 
     !      mass(i)   ( real )              atomic mass 
     !                                      ( in u.m.a, carbon mass is 12.0 )
     !      psfile(i) ( character(len=80) ) file name of the pseudopotential
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     SUBROUTINE card_atomic_species( input_line, prog )
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=256) :: input_line
       CHARACTER(LEN=2)   :: prog
       INTEGER            :: is, ip
       CHARACTER(LEN=4)   :: lb_pos
       CHARACTER(LEN=256) :: psfile
       LOGICAL, SAVE      :: tread = .FALSE.
       !
       !
       IF ( tread ) THEN
          CALL errore( ' card_atomic_species  ', ' two occurrence ', 2 )
       END IF
       IF ( ntyp > nsx ) THEN
         CALL errore( ' card_atomic_species ', ' nsp out of range ', ntyp )
       END IF
       !
       DO is = 1, ntyp
          !
          CALL read_line( input_line )
          READ( input_line, * ) lb_pos, atom_mass(is), psfile
          atom_pfile(is) = TRIM( psfile )
          lb_pos         = ADJUSTL( lb_pos )
          atom_label(is) = TRIM( lb_pos )
          !
!          IF ( atom_mass(is) <= 0.0_DP ) THEN
!             CALL errore( ' card_atomic_species ',' invalid  atom_mass ', is )
!          END IF
          DO ip = 1, is - 1
             IF ( atom_label(ip) == atom_label(is) ) THEN
                CALL errore( ' card_atomic_species ', &
                           & ' two occurrences of the same atomic label ', is )
             END IF
          END DO
          !
       END DO
       taspc = .TRUE.
       tread = .TRUE.
       !
       RETURN
       !
     END SUBROUTINE card_atomic_species
     !
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     !
     ! ATOMIC_POSITIONS
     !
     !   set the atomic positions in the cell
     !
     ! Syntax:
     !
     !   ATOMIC_POSITIONS (units_option)
     !     label(1) tau(1,1) tau(2,1) tau(3,1) mbl(1,1) mbl(2,1) mbl(3,1)
     !     label(2) tau(1,2) tau(2,2) tau(3,2) mbl(1,2) mbl(2,2) mbl(3,2)
     !      ...              ...               ...               ... ...
     !     label(n) tau(1,n) tau(2,n) tau(3,n) mbl(1,3) mbl(2,3) mbl(3,3)
     !
     ! Example:
     !
     ! ATOMIC_POSITIONS (bohr)
     !    O     0.0099    0.0099    0.0000  0 0 0
     !    H     1.8325   -0.2243   -0.0001  1 1 1
     !    H    -0.2243    1.8325    0.0002  1 1 1
     !
     ! Where:
     !
     !   units_option == crystal   position are given in scaled units
     !   units_option == bohr      position are given in Bohr
     !   units_option == angstrom  position are given in Angstrom
     !   units_option == alat      position are given in units of alat
     !
     !   label(k) ( character(len=4) )  atomic type
     !   tau(:,k) ( real )              coordinates  of the k-th atom
     !   mbl(:,k) ( integer )           mbl(i,k) > 0 the i-th coord. of the 
     !                                  k-th atom is allowed to be moved
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     ! ... routine modified for NEB           ( C.S. 21/10/2003 )
     ! ... routine modified for SMD           ( Y.K. 15/04/2004 )
     !
     SUBROUTINE card_atomic_positions( input_line, prog )
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=256) :: input_line
       CHARACTER(LEN=2)   :: prog
       CHARACTER(LEN=4)   :: lb_pos
       INTEGER            :: ia, k, is, nfield, idx, rep_i
       LOGICAL, EXTERNAL  :: matches
       LOGICAL            :: tend
       LOGICAL, SAVE      :: tread = .FALSE.
       !
       !
       IF ( tread ) THEN
          CALL errore( 'card_atomic_positions', 'two occurrence', 2 )
       END IF
       IF ( .NOT. taspc ) THEN
          CALL errore( 'card_atomic_positions', &
                     & 'ATOMIC_SPECIES must be present before', 2 )
       END IF
       IF ( ntyp > nsx ) THEN
          CALL errore( 'card_atomic_positions', 'nsp out of range', ntyp )
       END IF
       IF ( nat < 1 ) THEN
          CALL errore( 'card_atomic_positions', 'nat out of range', nat )
       END IF
       !
       if_pos = 1
       !
       sp_pos = 0
       rd_pos = 0.0_DP
       na_inp = 0
       !
       IF ( matches( "CRYSTAL", input_line ) ) THEN
          atomic_positions = 'crystal'
       ELSE IF ( matches( "BOHR", input_line ) ) THEN
          atomic_positions = 'bohr'
       ELSE IF ( matches( "ANGSTROM", input_line ) ) THEN
          atomic_positions = 'angstrom'
       ELSE IF ( matches( "ALAT", input_line ) ) THEN
          atomic_positions = 'alat'
       ELSE
          IF ( TRIM( ADJUSTL( input_line ) ) /= 'ATOMIC_POSITIONS' ) THEN
             CALL errore( 'read_cards ', &
                        & 'unknow unit option for ATOMIC_POSITION: '&
                        & // input_line, 1 )
          END IF
          IF ( prog == 'FP' ) atomic_positions = 'bohr'
          IF ( prog == 'CP' ) atomic_positions = 'bohr'
          IF ( prog == 'PW' ) atomic_positions = 'alat'
       END IF
       !
       IF ( full_phs_path_flag ) THEN
          !
          IF ( ALLOCATED( pos ) ) DEALLOCATE( pos )
          !
          ALLOCATE( pos( 3*nat, num_of_images ) )
          !
          pos(:,:) = 0.0_DP
          !
          IF ( calculation == 'smd' .AND. prog == 'CP' ) THEN
             !
             rep_loop : DO rep_i = 1, smd_kwnp
                !
                CALL read_line( input_line, end_of_file = tend )
                !
                IF ( tend ) &
                   CALL errore( 'read_cards', 'end of file reading ' // &
                              & 'atomic positions (smd)', rep_i )
                !
                IF ( matches( "first_image", input_line ) .OR. &
                     matches( "image",       input_line ) .OR. &
                     matches( "last_image",  input_line) ) THEN
                   !
                   CALL path_read_images( rep_i )
                   !
                ELSE
                   CALL errore( 'read_cards', 'missing or wrong image ' // &
                              & 'identifier in ATOMIC_POSITION', 1 )
                END IF
                !
             END DO rep_loop
             !
          ELSE
             !
             CALL read_line( input_line, end_of_file = tend )
             !
             IF ( tend ) &
                CALL errore( 'read_cards', &
                             'end of file reading atomic positions (path)', 1 )
             !
             IF ( matches( "first_image", input_line ) ) THEN
                !
                input_images = 1
                !
                CALL path_read_images( input_images )
                !
             ELSE
                !
                CALL errore( 'read_cards', &
                             'first_image missing in ATOMIC_POSITION', 1 )
                !
             END IF
             !
             read_conf_loop: DO
                !
                CALL read_line( input_line, end_of_file = tend )
                !
                IF ( tend ) &
                   CALL errore( 'read_cards', 'end of file reading ' // &
                              & 'atomic positions (path)', input_images + 1 )
                !
                input_images = input_images + 1
                !
                IF ( input_images > num_of_images ) &
                   CALL errore( 'read_cards', &
                              & 'too many images in ATOMIC_POSITION', 1 )
                !
                IF ( matches( "intermediate_image", input_line )  ) THEN
                   !
                   CALL path_read_images( input_images )
                   !
                ELSE
                   !
                   EXIT read_conf_loop
                   !
                END IF
                !
             END DO read_conf_loop
             !
             IF ( matches( "last_image", input_line ) ) THEN
                !
                CALL path_read_images( input_images )
                !
             ELSE
                !
                CALL errore( 'read_cards ', &
                             'last_image missing in ATOMIC_POSITION', 1 )
                !
             END IF
             !
          END IF
          !
       ELSE
          !
          DO ia = 1, nat
             !
             CALL read_line( input_line, end_of_file = tend )
             !
             IF ( tend ) &
                CALL errore( 'read_cards', &
                             'end of file reading atomic positions', ia )
             !
             CALL field_count( nfield, input_line )
             !
             IF ( sic /= 'none' .AND. nfield /= 8 ) &
                CALL errore( 'read_cards', &
                            'ATOMIC_POSITIONS with sic, 8 columns required', 1 )
             !
             IF ( nfield == 4 ) THEN
                !
                READ(input_line,*) lb_pos, ( rd_pos(k,ia), k = 1, 3 )
                !
             ELSE IF ( nfield == 7 ) THEN
                !
                READ(input_line,*) lb_pos, rd_pos(1,ia), &
                                           rd_pos(2,ia), &
                                           rd_pos(3,ia), &
                                           if_pos(1,ia), &
                                           if_pos(2,ia), &
                                           if_pos(3,ia)
                ! 
             ELSE IF ( nfield == 8 ) THEN
                !
                READ(input_line,*) lb_pos, rd_pos(1,ia), &
                                           rd_pos(2,ia), &
                                           rd_pos(3,ia), &
                                           if_pos(1,ia), &
                                           if_pos(2,ia), &
                                           if_pos(3,ia), &
                                           id_loc(ia)
                !
             ELSE
                !
                CALL errore( 'read_cards', 'wrong number of columns' // &
                           & 'in ATOMIC_POSITIONS', sp_pos(ia) )
                !
             END IF
             !
             lb_pos = ADJUSTL( lb_pos )
             !
             match_label: DO is = 1, ntyp
                !
                IF ( TRIM(lb_pos) == TRIM( atom_label(is) ) ) THEN
                   !
                   sp_pos(ia) = is
                   EXIT match_label
                   !
                END IF
                !
             END DO match_label
             !
             IF( ( sp_pos(ia) < 1 ) .OR. ( sp_pos(ia) > ntyp ) ) THEN
                !
                CALL errore( 'read_cards', &
                           & 'wrong index in ATOMIC_POSITIONS', ia )
                !
             END IF
             !
             is = sp_pos(ia)
             !
             na_inp(is) = na_inp(is) + 1
             !
          END DO
          !
       END IF
       !
       tapos = .TRUE.
       tread = .TRUE.
       !
       RETURN
       !
       CONTAINS
         !
         !-------------------------------------------------------------------
         SUBROUTINE path_read_images( image )
           !-------------------------------------------------------------------
           !
           IMPLICIT NONE
           !
           INTEGER, INTENT(IN) :: image
           !
           !
           DO ia = 1, nat
              !
              idx = 3 * ( ia - 1 )
              !
              CALL read_line( input_line, end_of_file = tend )
              !
              IF ( tend ) &
                 CALL errore( 'read_cards', &
                              'end of file reading atomic positions', ia )
              !
              CALL field_count( nfield, input_line )
              !
              IF ( nfield == 4 ) THEN
                 !
                 READ( input_line, * ) lb_pos, pos((idx+1),image), &
                                               pos((idx+2),image), &
                                               pos((idx+3),image)
                 !
              ELSE IF ( nfield == 7 ) THEN
                 !
                 IF ( image /= 1 ) THEN
                    !
                    CALL errore( 'read_cards', &
                               & 'wrong number of columns in ' // &
                               & 'ATOMIC_POSITIONS', sp_pos(ia) )
                    !
                 END IF
                 !
                 READ( input_line, * ) lb_pos, pos((idx+1),image), &
                                               pos((idx+2),image), &
                                               pos((idx+3),image), &
                                               if_pos(1,ia), &
                                               if_pos(2,ia), &
                                               if_pos(3,ia)
                 !
              ELSE
                 ! 
                 CALL errore( 'read_cards', &
                            & 'wrong number of columns in ' // &
                            & 'ATOMIC_POSITIONS', sp_pos(ia) )
                 !
              END IF
              !
              IF ( image == 1 ) THEN
                 !
                 lb_pos = ADJUSTL( lb_pos )
                 !
                 match_label_path: DO is = 1, ntyp
                    !
                    IF ( TRIM( lb_pos ) == TRIM( atom_label(is) ) ) THEN
                       !
                       sp_pos(ia) = is
                       !
                       EXIT match_label_path
                       !
                    END IF
                    !
                 END DO match_label_path
                 !
                 IF ( ( sp_pos(ia) < 1 ) .OR. ( sp_pos(ia) > ntyp ) ) THEN
                    !
                    CALL errore( 'read_cards', &
                                 'wrong index in ATOMIC_POSITIONS', ia )
                    !
                 END IF
                 !
                 is = sp_pos(ia)
                 !
                 na_inp( is ) = na_inp( is ) + 1
                 !
              END IF
              !
           END DO
           !
           RETURN
           !
         END SUBROUTINE path_read_images
         !
     END SUBROUTINE card_atomic_positions
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     !
     ! K_POINTS
     !
     !   use the specified set of k points
     !
     ! Syntax:
     !
     !   K_POINTS (mesh_option)
     !     n
     !     xk(1,1) xk(2,1) xk(3,1) wk(1)
     !     ...     ...     ...     ...
     !     xk(1,n) xk(2,n) xk(3,n) wk(n)
     !
     ! Example:
     !
     ! K_POINTS
     !   10
     !    0.1250000  0.1250000  0.1250000   1.00
     !    0.1250000  0.1250000  0.3750000   3.00
     !    0.1250000  0.1250000  0.6250000   3.00
     !    0.1250000  0.1250000  0.8750000   3.00
     !    0.1250000  0.3750000  0.3750000   3.00
     !    0.1250000  0.3750000  0.6250000   6.00
     !    0.1250000  0.3750000  0.8750000   6.00
     !    0.1250000  0.6250000  0.6250000   3.00
     !    0.3750000  0.3750000  0.3750000   1.00
     !    0.3750000  0.3750000  0.6250000   3.00
     !
     ! Where:
     !
     !   mesh_option == automatic  k points mesh is generated automatically
     !                             with Monkhorst-Pack algorithm
     !   mesh_option == crystal    k points mesh is given in stdin in scaled 
     !                             units
     !   mesh_option == tpiba      k points mesh is given in stdin in units 
     !                             of ( 2 PI / alat )
     !   mesh_option == gamma      only gamma point is used ( default in 
     !                             CPMD simulation )
     !
     !   n       ( integer )  number of k points
     !   xk(:,i) ( real )     coordinates of i-th k point
     !   wk(i)   ( real )     weights of i-th k point
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     SUBROUTINE card_kpoints( input_line )
       ! 
       IMPLICIT NONE
       !
       CHARACTER(LEN=256) :: input_line
       INTEGER            :: i
       LOGICAL, EXTERNAL  :: matches
       LOGICAL, SAVE      :: tread = .FALSE.
       LOGICAL            :: tend
       !
       !
       IF ( tread ) THEN
          CALL errore( ' card_kpoints ', ' two occurrence ', 2 )
       END IF
       !
       IF ( matches( "AUTOMATIC", input_line ) ) THEN
          !  automatic generation of k-points
          k_points = 'automatic'
       ELSE IF ( matches( "CRYSTAL", input_line ) ) THEN
          !  input k-points are in crystal (reciprocal lattice) axis
          k_points = 'crystal'
       ELSE IF ( matches( "TPIBA", input_line ) ) THEN
          !  input k-points are in 2pi/a units
          k_points = 'tpiba'
       ELSE IF ( matches( "GAMMA", input_line ) ) THEN
          !  Only Gamma (k=0) is used
          k_points = 'gamma'
       ELSE
          !  by default, input k-points are in 2pi/a units
          k_points = 'tpiba'
       END IF
       !
       IF ( k_points == 'automatic' ) THEN
          !
          ! ... automatic generation of k-points
          !
          nkstot = 0
          CALL read_line( input_line, end_of_file = tend )
          IF (tend) GO TO 10
          READ(input_line, *, END=10, ERR=10) nk1, nk2, nk3, k1, k2 ,k3
          IF ( k1 < 0 .OR. k1 > 1 .OR. &
               k2 < 0 .OR. k2 > 1 .OR. &
               k3 < 0 .OR. k3 > 1 ) CALL errore &
                  ('card_kpoints', 'invalid offsets: must be 0 or 1', 1)
          IF ( nk1 <= 0 .OR. nk2 <= 0 .OR. nk3 <= 0 ) CALL errore &
                  ('card_kpoints', 'invalid values for nk1, nk2, nk3', 1)

          !
       ELSE IF ( ( k_points == 'tpiba' ) .OR. ( k_points == 'crystal' ) ) THEN
          !
          ! ... input k-points are in 2pi/a units
          !
          CALL read_line( input_line, end_of_file = tend )
          IF (tend) GO TO 10
          READ(input_line, *, END=10, ERR=10) nkstot
          IF ( nkstot > SIZE (xk,2)  ) CALL errore &
                  ('card_kpoints', 'too many k-points',nkstot)
          !
          DO i = 1, nkstot
             CALL read_line( input_line, end_of_file = tend )
             IF (tend) GO TO 10
             READ(input_line,*, END=10, ERR=10) xk(1,i), xk(2,i), xk(3,i), wk(i)
          END DO
          !
       ELSE IF ( k_points == 'gamma' ) THEN
          !
          nkstot = 1
          xk(:,1) = 0.0_DP
          wk(1) = 1.0_DP
          !
       END IF
       !
       tread  = .TRUE.
       tk_inp = .TRUE.
       !
       RETURN
10     CALL errore ('card_kpoints', ' error or end of file while reading ' &
            & // TRIM(k_points) // ' k points', 1)
       !
     END SUBROUTINE card_kpoints
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     !
     ! SETNFI
     !
     !   Reset the step counter to the specified value
     !
     ! Syntax:
     !
     !  SETNFI
     !     nfi
     !
     ! Example:
     !
     !  SETNFI
     !     100
     !
     ! Where:
     !
     !    nfi (integer) new value for the step counter
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     SUBROUTINE card_setnfi( input_line )
       ! 
       IMPLICIT NONE
       !
       CHARACTER(LEN=256) :: input_line
       LOGICAL, SAVE      :: tread = .FALSE.
       !
       !
       IF ( tread ) THEN
          CALL errore( ' card_setnfi ', ' two occurrence ', 2 )
       END IF
       CALL read_line( input_line )
       READ(input_line,*) newnfi_card
       tnewnfi_card = .TRUE.
       tread = .TRUE.
       ! 
       RETURN
       !
     END SUBROUTINE card_setnfi
     !
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     !
     ! 2DPROCMESH
     !
     !   Distribute the Y and Z FFT dimensions across processors, 
     !   instead of Z dimension only ( default distribution )
     !
     ! Syntax:
     !
     !    2DPROCMESH
     !
     ! Where:
     !
     !    no parameters
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     !
     ! OCCUPATIONS
     !
     !   use the specified occupation numbers for electronic states.
     !   Note that you should specify 10 values per line maximum!
     !
     ! Syntax (nspin == 1):
     !
     !   OCCUPATIONS
     !      f(1)  ....   ....  f(10) 
     !      f(11) .... f(nbnd)
     !
     ! Syntax (nspin == 2):
     !
     !   OCCUPATIONS
     !      u(1)  ....   ....  u(10) 
     !      u(11) .... u(nbnd)
     !      d(1)  ....   ....  d(10)
     !      d(11) .... d(nbnd)
     !
     ! Example:
     !
     ! OCCUPATIONS
     !  2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0 2.0
     !  2.0 2.0 2.0 2.0 2.0 1.0 1.0 
     !
     ! Where:
     !
     !      f(:) (real)  these are the occupation numbers 
     !                   for LDA electronic states. 
     !
     !      u(:) (real)  these are the occupation numbers
     !                   for LSD spin == 1 electronic states 
     !      d(:) (real)  these are the occupation numbers
     !                   for LSD spin == 2 electronic states 
     !
     !      Note, maximum 10 values per line!
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     SUBROUTINE card_occupations( input_line )
       !
       IMPLICIT NONE
       ! 
       CHARACTER(LEN=256) :: input_line
       INTEGER            :: is, nx10, i, j, nspin0
       LOGICAL, SAVE      :: tread = .FALSE.
       !
       !
       IF ( tread ) THEN
          CALL errore( ' card_occupations ', ' two occurrence ', 2 )
       END IF
       nspin0=nspin
       if (nspin == 4) nspin0=1
       !
       ALLOCATE ( f_inp ( nbnd, nspin0 ) )
       DO is = 1, nspin0
          !
          nx10 = 10 * INT( nbnd / 10 )
          DO i = 1, nx10, 10
             CALL read_line( input_line )
             READ(input_line,*) ( f_inp(j,is), j = i, ( i + 9 ) )
          END DO
          IF ( MOD( nbnd, 10 ) > 0 ) THEN
             CALL read_line( input_line )
             READ(input_line,*) ( f_inp(j,is), j = ( nx10 + 1 ), nbnd)
          END IF
          !
       END DO
       !
       tf_inp = .TRUE. 
       tread = .TRUE.
       ! 
       RETURN
       !
     END SUBROUTINE card_occupations
     !
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     !
     ! VHMEAN
     !
     !   Calculation of potential average along a given axis
     !
     ! Syntax:
     !
     !   VHMEAN
     !   unit nr rmin rmax asse
     !
     ! Example:
     !
     !   ????
     !
     ! Where:
     !
     !   ????
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     SUBROUTINE card_vhmean( input_line )
       ! 
       IMPLICIT NONE
       !
       CHARACTER(LEN=256) :: input_line
       LOGICAL, SAVE :: tread = .FALSE.
       !
       !
       IF ( tread ) THEN
          CALL errore( ' card_vhmean ', ' two occurrence ', 2 )
       END IF
       !  
       tvhmean_inp = .TRUE.
       CALL read_line( input_line )
       READ(input_line,*) &
           vhiunit_inp, vhnr_inp, vhrmin_inp, vhrmax_inp, vhasse_inp
       tread = .TRUE.
       !
       RETURN
       !
     END SUBROUTINE card_vhmean
     !
     !
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     !
     ! DIPOLE
     !
     !   calculate polarizability
     !
     ! Syntax:
     !
     !   DIPOLE
     !
     ! Where:
     !
     !    no parameters
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     SUBROUTINE card_dipole( input_line )
       ! 
       IMPLICIT NONE
       !
       CHARACTER(LEN=256) :: input_line
       LOGICAL, SAVE      :: tread = .FALSE.
       !
       !
       IF ( tread ) THEN
          CALL errore( ' card_dipole ', ' two occurrence ', 2 )
       END IF
       !
       tdipole_card = .TRUE.
       tread = .TRUE.
       !
       RETURN
       !
     END SUBROUTINE card_dipole
     !
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     !
     ! IESR
     !
     !   use the specified number of neighbour cells for Ewald summations
     !
     ! Syntax:
     !
     !   ESR
     !    iesr
     !
     ! Example:
     !
     !   ESR
     !    3
     !
     ! Where:
     !
     !      iesr (integer)  determines the number of neighbour cells to be
     !                      considered:
     !                        iesr = 1 : nearest-neighbour cells (default)
     !                        iesr = 2 : next-to-nearest-neighbour cells
     !                        and so on
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     SUBROUTINE card_esr( input_line )
       ! 
       IMPLICIT NONE
       ! 
       CHARACTER(LEN=256) :: input_line
       LOGICAL, SAVE      :: tread = .FALSE.
       ! 
       IF ( tread ) THEN
          CALL errore( ' card_esr ', ' two occurrence ', 2 )
       END IF
       CALL read_line( input_line )
       READ(input_line,*) iesr_inp
       !
       tread = .TRUE.
       !
       RETURN
       !
     END SUBROUTINE card_esr
     !
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     !
     ! NEIGHBOURS
     !
     !   calculate the neighbours of (and the disance from) each atoms below 
     !   the distance specified by the parameter
     !
     ! Syntax:
     !
     !   NEIGHBOURS
     !      cut_radius
     !
     ! Example:
     !
     !   NEIGHBOURS
     !      4.0
     !
     ! Where:
     !
     !      cut_radius ( real )  radius of the region where atoms are 
     !                           considered as neighbours ( in a.u. )
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     SUBROUTINE card_neighbours( input_line )
       ! 
       IMPLICIT NONE
       ! 
       CHARACTER(LEN=256) :: input_line
       LOGICAL, SAVE      :: tread = .FALSE.
       ! 
       !
       IF ( tread ) THEN
          CALL errore( ' card_neighbours ', ' two occurrence ', 2 )
       END IF
       ! 
       CALL read_line( input_line )
       READ(input_line, *) neighbo_radius
       ! 
       tneighbo = .TRUE.
       tread = .TRUE.
       ! 
       RETURN
       !
     END SUBROUTINE card_neighbours
     !
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     !
     ! CELL_PARAMETERS
     !
     !   use the specified cell dimensions
     !
     ! Syntax:
     !
     !    CELL_PARAMETERS
     !      HT(1,1) HT(1,2) HT(1,3)
     !      HT(2,1) HT(2,2) HT(2,3)
     !      HT(3,1) HT(3,2) HT(3,3)
     !
     ! Example:
     !
     ! CELL_PARAMETERS
     !    24.50644311    0.00004215   -0.14717844
     !    -0.00211522    8.12850030    1.70624903
     !     0.16447787    0.74511792   23.07395418
     !
     ! Where:
     !
     !      HT(i,j) (real)  cell dimensions ( in a.u. ), 
     !                      note the relation with lattice vectors: 
     !                      HT(1,:) = A1, HT(2,:) = A2, HT(3,:) = A3 
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     SUBROUTINE card_cell_parameters( input_line )
       ! 
       IMPLICIT NONE
       ! 
       CHARACTER(LEN=256) :: input_line
       INTEGER            :: i, j
       LOGICAL, EXTERNAL  :: matches
       LOGICAL, SAVE      :: tread = .FALSE.
       !
       !
       IF ( tread ) THEN
          CALL errore( ' card_cell_parameters ', ' two occurrence ', 2 )
       END IF
       !
       IF ( matches( 'HEXAGONAL', input_line ) ) then
          cell_symmetry = 'hexagonal'
       ELSE
          cell_symmetry = 'cubic'
       END IF
       !
       IF ( matches( "BOHR", input_line ) ) THEN
          cell_units = 'bohr'
       ELSE IF ( matches( "ANGSTROM", input_line ) ) THEN
          cell_units = 'angstrom'
       ELSE
          cell_units = 'alat'
       END IF
       !
       DO i = 1, 3
          CALL read_line( input_line )
          READ(input_line,*) ( rd_ht( i, j ), j = 1, 3 )
       END DO
       ! 
       trd_ht = .TRUE.
       tread  = .TRUE.
       !
       RETURN
       !
     END SUBROUTINE card_cell_parameters
     !
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     !
     ! TURBO
     !
     !   allocate space to store electronic states in real space while 
     !   computing charge density, and then reuse the stored state
     !   in the calculation of forces instead of repeating the FFT
     !
     ! Syntax:
     !
     !    TURBO
     !      nturbo
     !
     ! Example:
     !
     !    TURBO
     !      64
     !
     ! Where:
     !
     !      nturbo (integer)  number of states to be stored
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     SUBROUTINE card_turbo( input_line )
       ! 
       IMPLICIT NONE
       ! 
       CHARACTER(LEN=256) :: input_line
       LOGICAL, SAVE      :: tread = .FALSE.
       ! 
       !
       IF ( tread ) THEN
          CALL errore( ' card_turbo ', ' two occurrence ', 2 )
       END IF
       !
       CALL read_line( input_line )
       READ(input_line,*) nturbo_inp
       ! 
       IF( (nturbo_inp < 0) .OR. (nturbo_inp > (nbnd/2)) ) THEN
          CALL errore( ' card_turbo ', ' NTURBO OUT OF RANGE ', nturbo_inp )
       END IF
       !
       tturbo_inp = .TRUE.
       tread = .TRUE.
       ! 
       RETURN
       !
     END SUBROUTINE
     !
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     !
     ! ATOMIC_VELOCITIES
     !
     !   read velocities (in atomic units) from standard input
     !
     ! Syntax:
     !
     !   ATOMIC_VELOCITIES
     !     label(1)  Vx(1) Vy(1) Vz(1)
     !     ....
     !     label(n)  Vx(n) Vy(n) Vz(n)
     !
     ! Example:
     !
     !   ???
     !
     ! Where:
     !
     !   label (character(len=4))       atomic label
     !   Vx(:), Vy(:) and Vz(:) (REAL)  x, y and z velocity components of
     !                                  the ions
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     SUBROUTINE card_ion_velocities( input_line )
       ! 
       IMPLICIT NONE
       ! 
       CHARACTER(LEN=256) :: input_line
       INTEGER            :: ia, k, is, nfield
       LOGICAL, SAVE      :: tread = .FALSE.
       CHARACTER(LEN=4)   :: lb_vel
       !
       !
       IF( tread ) THEN
          CALL errore( ' card_ion_velocities ', ' two occurrence ', 2 )
       END IF
       !
       IF( .NOT. taspc ) THEN
          CALL errore( ' card_ion_velocities ', &
                     & ' ATOMIC_SPECIES must be present before ', 2 )
       END IF
       !
       rd_vel = 0.0_DP
       sp_vel = 0
       !
       IF ( ion_velocities == 'from_input' ) THEN
          !
          tavel = .TRUE.
          !
          DO ia = 1, nat
             !
             CALL read_line( input_line )
             CALL field_count( nfield, input_line )
             IF ( nfield == 4 ) THEN
                READ(input_line,*) lb_vel, ( rd_vel(k,ia), k = 1, 3 )
             ELSE
                CALL errore( ' iosys ', &
                           & ' wrong entries in ION_VELOCITIES ', ia )
             END IF
             !
             match_label: DO is = 1, ntyp
                IF ( TRIM( lb_vel ) == atom_label(is) ) THEN
                   sp_vel(ia) = is
                   EXIT match_label
                END IF
             END DO match_label
             !
             IF ( sp_vel(ia) < 1 .OR. sp_vel(ia) > ntyp ) THEN
                CALL errore( ' iosys ', ' wrong LABEL in ION_VELOCITIES ', ia )
             END IF
             !
          END DO
          !
       END IF
       !
       tread = .TRUE.
       !
       RETURN
       !
     END SUBROUTINE
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     !
     ! CONSTRAINTS
     !
     !   Ionic Constraints
     !
     ! Syntax:
     !
     !    CONSTRAINTS
     !      NCONSTR CONSTR_TOL
     !      CONSTR_TYPE(.) CONSTR(1,.) CONSTR(2,.) ... { CONSTR_TARGET(.) }
     !
     ! Where:
     !
     !      NCONSTR(INTEGER)    number of constraints
     !
     !      CONSTR_TOL          tolerance for keeping the constraints 
     !                          satisfied
     !
     !      CONSTR_TYPE(.)      type of constrain: 
     !                          1: for fixed distances ( two atom indexes must
     !                             be specified )
     !                          2: for fixed planar angles ( three atom indexes
     !                             must be specified )
     !
     !      CONSTR(1,.) CONSTR(2,.) ...
     !
     !                          indices object of the constraint, as 
     !                          they appear in the 'POSITION' CARD
     !
     !      CONSTR_TARGET       target for the constrain ( in the case of 
     !                          planar angles it is the COS of the angle ).
     !                          this variable is optional.
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     SUBROUTINE card_constraints( input_line )
       ! 
       IMPLICIT NONE
       ! 
       CHARACTER(LEN=256) :: input_line
       INTEGER            :: i, nfield
       LOGICAL, SAVE      :: tread = .FALSE.
       ! 
       !
       IF ( tread ) CALL errore( 'card_constraints', 'two occurrence', 2 )
       !
       CALL read_line( input_line )
       !
       CALL field_count( nfield, input_line )
       !
       IF ( nfield == 1 ) THEN
          !
          READ( input_line, * ) nconstr_inp
          !
       ELSE IF ( nfield == 2 ) THEN
          !
          READ( input_line, * ) nconstr_inp, constr_tol_inp
          !
       ELSE
          !
          CALL errore( 'card_constraints', 'too many fields', nfield )
          !
       END IF
       !
       CALL allocate_input_constr()
       !
       DO i = 1, nconstr_inp
          !
          CALL read_line( input_line )
          !
          READ( input_line, * ) constr_type_inp(i)
          !
          CALL field_count( nfield, input_line )
          !
          IF ( nfield > nc_fields + 2 ) &
             CALL errore( 'card_constraints', &
                          'too many fields for this constraint', i )
          !
          SELECT CASE( constr_type_inp(i) )
          CASE( 'type_coord', 'atom_coord' )
             !
             IF ( nfield == 5 ) THEN
                !
                READ( input_line, * ) constr_type_inp(i), &
                                      constr_inp(1,i), &
                                      constr_inp(2,i), &
                                      constr_inp(3,i), &
                                      constr_inp(4,i)
                !
             ELSE IF ( nfield == 6 ) THEN
                !
                READ( input_line, * ) constr_type_inp(i), &
                                      constr_inp(1,i), &
                                      constr_inp(2,i), &
                                      constr_inp(3,i), &
                                      constr_inp(4,i), &
                                      constr_target(i)
                !
                constr_target_set(i) = .TRUE.
                !
             ELSE
                !
                CALL errore( 'card_constraints', 'type_coord, ' // &
                           & 'atom_coord: wrong number of fields', nfield )
                !
             END IF
             !
          CASE( 'distance' )
             !
             IF ( nfield == 3 ) THEN
                !
                READ( input_line, * ) constr_type_inp(i), &
                                      constr_inp(1,i), &
                                      constr_inp(2,i)
                !
             ELSE IF ( nfield == 4 ) THEN
                !
                READ( input_line, * ) constr_type_inp(i), &
                                      constr_inp(1,i), &
                                      constr_inp(2,i), &
                                      constr_target(i)
                !
                constr_target_set(i) = .TRUE.
                !
             ELSE
                !
                CALL errore( 'card_constraints', &
                           & 'distance: wrong number of fields', nfield )
                !
             END IF
             !
          CASE( 'planar_angle' )
             !
             IF ( nfield == 4 ) THEN
                !
                READ( input_line, * ) constr_type_inp(i), &
                                      constr_inp(1,i), &
                                      constr_inp(2,i), &
                                      constr_inp(3,i)
                !
             ELSE IF ( nfield == 5 ) THEN
                !
                READ( input_line, * ) constr_type_inp(i), &
                                      constr_inp(1,i), &
                                      constr_inp(2,i), &
                                      constr_inp(3,i), &
                                      constr_target(i)
                !
                constr_target_set(i) = .TRUE.
                !
             ELSE
                !
                CALL errore( 'card_constraints', &
                           & 'planar_angle: wrong number of fields', nfield )
                !
             END IF
             !
          CASE( 'torsional_angle' )
             !
             IF ( nfield == 5 ) THEN
                !
                READ( input_line, * ) constr_type_inp(i), &
                                      constr_inp(1,i), &
                                      constr_inp(2,i), &
                                      constr_inp(3,i), &
                                      constr_inp(4,i)
                !
             ELSE IF ( nfield == 6 ) THEN
                !
                READ( input_line, * ) constr_type_inp(i), &
                                      constr_inp(1,i), &
                                      constr_inp(2,i), &
                                      constr_inp(3,i), &
                                      constr_inp(4,i), &
                                      constr_target(i)
                !
                constr_target_set(i) = .TRUE.
                !
             ELSE
                !
                CALL errore( 'card_constraints', &
                           & 'torsional_angle: wrong number of fields', nfield )
                !
             END IF
             !
          CASE( 'bennett_proj' )
             !
             IF ( nfield == 5 ) THEN
                !
                READ( input_line, * ) constr_type_inp(i), &
                                      constr_inp(1,i), &
                                      constr_inp(2,i), &
                                      constr_inp(3,i), &
                                      constr_inp(4,i)
                !
             ELSE IF ( nfield == 6 ) THEN
                !
                READ( input_line, * ) constr_type_inp(i), &
                                      constr_inp(1,i), &
                                      constr_inp(2,i), &
                                      constr_inp(3,i), &
                                      constr_inp(4,i), &
                                      constr_target(i)
                !
                constr_target_set(i) = .TRUE.
                !
             ELSE
                !
                CALL errore( 'card_constraints', &
                           & 'bennett_proj: wrong number of fields', nfield )
                !
             END IF
             !
          CASE DEFAULT
             !
             CALL errore( 'card_constraints', 'unknown constraint ' // &
                        & 'type: ' // TRIM( constr_type_inp(i) ), 1 )
             !
          END SELECT
          !
       END DO
       !
       tread = .TRUE.
       !
       RETURN
       !
     END SUBROUTINE card_constraints
     !
     SUBROUTINE card_collective_vars( input_line )
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=256) :: input_line
       INTEGER            :: i, nfield
       LOGICAL            :: ltest
       LOGICAL, SAVE      :: tread = .FALSE.
       ! 
       !
       IF ( tread ) CALL errore( 'card_collective_vars', 'two occurrence', 2 )
       !
       CALL read_line( input_line )
       !
       CALL field_count( nfield, input_line )
       !
       IF ( nfield == 1 ) THEN
          !
          READ( input_line, * ) ncolvar_inp
          !
       ELSE IF ( nfield == 2 ) THEN
          !
          READ( input_line, * ) ncolvar_inp, colvar_tol_inp
          !
       ELSE
          !
          CALL errore( 'card_collective_vars', 'too many fields', nfield )
          !
       END IF
       !
       CALL allocate_input_colvar()
       !
       IF ( cg_phs_path_flag ) THEN
          !
          input_images = 2
          !
          IF( ALLOCATED( pos ) ) DEALLOCATE( pos )
          !
          ALLOCATE( pos( ncolvar_inp, input_images ) )
          !
          pos(:,:) = 0.0_DP
          !
       END IF
       !
       DO i = 1, ncolvar_inp
          !
          CALL read_line( input_line )
          !
          READ( input_line, * ) colvar_type_inp(i)
          !
          CALL field_count( nfield, input_line )
          !
          ltest = ( ( nfield <= nc_fields + 2 ) .OR. &
                    ( cg_phs_path_flag .AND. ( nfield <= nc_fields + 4 ) ) )
          !
          IF ( .NOT. ltest ) &
             CALL errore( 'card_collective_vars', 'too many fields for ' // &
                        & 'this constraint: ' // TRIM( constr_type_inp(i) ), i )
          !
          SELECT CASE( colvar_type_inp(i) )
          CASE( 'type_coord', 'atom_coord' )
             !
             IF ( cg_phs_path_flag ) THEN
                !
                READ( input_line, * ) colvar_type_inp(i), &
                                      colvar_inp(1,i),    &
                                      colvar_inp(2,i),    &
                                      colvar_inp(3,i),    &
                                      colvar_inp(4,i),    &
                                      pos(i,1),           &
                                      pos(i,2)
                !
             ELSE IF ( nfield == 5 ) THEN
                !
                READ( input_line, * ) colvar_type_inp(i), &
                                      colvar_inp(1,i), &
                                      colvar_inp(2,i), &
                                      colvar_inp(3,i), &
                                      colvar_inp(4,i)
                !
             ELSE
                !
                CALL errore( 'card_collective_vars', 'type_coord, ' // &
                           & 'atom_coord: wrong number of fields', nfield )
                !
             END IF
             !
          CASE( 'distance' )
             !
             IF ( cg_phs_path_flag ) THEN
                !
                READ( input_line, * ) colvar_type_inp(i), &
                                      colvar_inp(1,i),    &
                                      colvar_inp(2,i),    &
                                      pos(i,1),           &
                                      pos(i,2)
                !
             ELSE IF ( nfield == 3 ) THEN
                !
                READ( input_line, * ) colvar_type_inp(i), &
                                      colvar_inp(1,i), &
                                      colvar_inp(2,i)
                !
             ELSE
                !
                CALL errore( 'card_collective_vars', &
                           & 'distance: wrong number of fields', nfield )
                !
             END IF
             !
          CASE( 'planar_angle' )
             !
             IF ( cg_phs_path_flag ) THEN
                !
                READ( input_line, * ) colvar_type_inp(i), &
                                      colvar_inp(1,i),    &
                                      colvar_inp(2,i),    &
                                      colvar_inp(3,i),    &
                                      pos(i,1),           &
                                      pos(i,2)
                !
             ELSE IF ( nfield == 4 ) THEN
                !
                READ( input_line, * ) colvar_type_inp(i), &
                                      colvar_inp(1,i), &
                                      colvar_inp(2,i), &
                                      colvar_inp(3,i)
                !
             ELSE
                !
                CALL errore( 'card_collective_vars', &
                           & 'planar_angle: wrong number of fields', nfield )
                !
             END IF
             !
          CASE( 'torsional_angle' )
             !
             IF ( cg_phs_path_flag ) THEN
                !
                READ( input_line, * ) colvar_type_inp(i), &
                                      colvar_inp(1,i),    &
                                      colvar_inp(2,i),    &
                                      colvar_inp(3,i),    &
                                      colvar_inp(4,i),    &
                                      pos(i,1),           &
                                      pos(i,2)
                !
             ELSE IF ( nfield == 5 ) THEN
                !
                READ( input_line, * ) colvar_type_inp(i), &
                                      colvar_inp(1,i), &
                                      colvar_inp(2,i), &
                                      colvar_inp(3,i), &
                                      colvar_inp(4,i)
                !
             ELSE
                !
                CALL errore( 'card_collective_vars', &
                           & 'torsional_angle: wrong number of fields', nfield )
                !
             END IF
             !
          CASE( 'struct_fac' )
             !
             IF ( cg_phs_path_flag ) THEN
                !
                READ( input_line, * ) colvar_type_inp(i), &
                                      colvar_inp(1,i), &
                                      colvar_inp(2,i), &
                                      colvar_inp(3,i), &
                                      pos(i,1),        &
                                      pos(i,2)
                !
             ELSE IF ( nfield == 4 ) THEN
                !
                READ( input_line, * ) colvar_type_inp(i), &
                                      colvar_inp(1,i), &
                                      colvar_inp(2,i), &
                                      colvar_inp(3,i)
                !
             ELSE
                !
                CALL errore( 'card_collective_vars', &
                           & 'struct_fac: wrong number of fields', nfield )
                !
             END IF
             !
          CASE( 'sph_struct_fac' )
             !
             IF ( cg_phs_path_flag ) THEN
                !
                READ( input_line, * ) colvar_type_inp(i), &
                                      colvar_inp(1,i), &
                                      pos(i,1),        &
                                      pos(i,2)
                !
             ELSE IF ( nfield == 2 ) THEN
                !
                READ( input_line, * ) colvar_type_inp(i), &
                                      colvar_inp(1,i)
                !
             ELSE
                !
                CALL errore( 'card_collective_vars',  &
                           & 'sph_struct_fac: wrong number of fields', nfield )
                !
             END IF
             !
          CASE( 'bennett_proj' )
             !
             IF ( cg_phs_path_flag ) THEN
                !
                READ( input_line, * ) constr_type_inp(i), &
                                      constr_inp(1,i), &
                                      constr_inp(2,i), &
                                      constr_inp(3,i), &
                                      constr_inp(4,i), &
                                      pos(i,1),        &
                                      pos(i,2)
                !
             ELSE IF ( nfield == 5 ) THEN
                !
                READ( input_line, * ) constr_type_inp(i), &
                                      constr_inp(1,i), &
                                      constr_inp(2,i), &
                                      constr_inp(3,i), &
                                      constr_inp(4,i)
                !
             ELSE
                !
                CALL errore( 'card_collective_vars', &
                           & 'bennett_proj: wrong number of fields', nfield )
                !
             END IF
             !
          CASE DEFAULT
             !
             CALL errore( 'card_collective_vars', 'unknown collective ' // &
                        & 'variable: ' // TRIM( colvar_type_inp(i) ), 1 )
             !
          END SELECT
          !
       END DO
       !
       tread = .TRUE.
       !
       RETURN
       !
     END SUBROUTINE card_collective_vars
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     ! 
     ! KSOUT
     !
     !   Enable the printing of Kohn Sham states
     !
     ! Syntax ( nspin == 2 ):
     !
     !   KSOUT
     !     nu
     !     iu(1) iu(2) iu(3) .. iu(nu)
     !     nd
     !     id(1) id(2) id(3) .. id(nd)
     !
     ! Syntax ( nspin == 1 ):
     !
     !   KSOUT
     !     ns
     !     is(1) is(2) is(3) .. is(ns)
     !
     ! Example:
     !
     !   ???
     !   
     ! Where:
     !
     !   nu (integer)     number of spin=1 states to be printed 
     !   iu(:) (integer)  indexes of spin=1 states, the state iu(k) 
     !                    is saved to file KS_UP.iu(k)
     !
     !   nd (integer)     number of spin=2 states to be printed 
     !   id(:) (integer)  indexes of spin=2 states, the state id(k) 
     !                    is saved to file KS_DW.id(k)
     !
     !   ns (integer)     number of LDA states to be printed 
     !   is(:) (integer)  indexes of LDA states, the state is(k) 
     !                    is saved to file KS.is(k)
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     SUBROUTINE card_ksout( input_line )
       ! 
       IMPLICIT NONE
       ! 
       CHARACTER(LEN=256) :: input_line
       LOGICAL, SAVE      :: tread = .FALSE.
       INTEGER            :: i, s, nksx
       TYPE occupancy_type
          INTEGER, pointer :: occs(:)
       END TYPE occupancy_type
       TYPE(occupancy_type), ALLOCATABLE :: is(:)
       !
       IF ( tread ) THEN
          CALL errore( ' card_ksout ', ' two occurrence ', 2 )
       END IF
       !
       nprnks = 0 
       nksx   = 0
       !
       ALLOCATE ( is (nspin) )
       !
       DO s = 1, nspin
          !
          CALL read_line( input_line )
          READ(input_line, *) nprnks( s )
          !
          IF ( nprnks( s ) < 1 ) THEN
             CALL errore( ' card_ksout ', ' wrong number of states ', 2 )
          END IF
          !
          ALLOCATE( is(s)%occs( 1:nprnks(s) ) )
          !
          CALL read_line( input_line )
          READ(input_line, *) ( is(s)%occs(i), i = 1, nprnks( s ) )
          !
          nksx = MAX( nksx, nprnks( s ) )
          !
       END DO
       !
       CALL allocate_input_iprnks( nksx, nspin )
       !
       DO s = 1, nspin
          !
          DO i = 1, nprnks( s )
             !
             iprnks( i, s ) = is(s)%occs(i)
             !
          END DO
          !
          DEALLOCATE( is(s)%occs )
          !
       END DO
       !
       DEALLOCATE( is )
       !
       tread = .TRUE.
       !
       RETURN
       !
     END SUBROUTINE
     !
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     ! 
     ! KSOUT_EMPTY
     !
     !   Enable the printing of empty Kohn Sham states
     !
     ! Syntax ( nspin == 2 ):
     !
     !   KSOUT_EMPTY
     !     nu
     !     iu(1) iu(2) iu(3) .. iu(nu)
     !     nd
     !     id(1) id(2) id(3) .. id(nd)
     !
     ! Syntax ( nspin == 1 ):
     !
     !   KSOUT_EMPTY
     !     ns
     !     is(1) is(2) is(3) .. is(ns)
     !
     ! Example:
     !
     !   ???
     !   
     ! Where:
     !
     !   nu (integer)     number of spin=1 empty states to be printed 
     !   iu(:) (integer)  indexes of spin=1 empty states, the state iu(k) 
     !                    is saved to file KS_EMP_UP.iu(k)
     !
     !   nd (integer)     number of spin=2 empty states to be printed 
     !   id(:) (integer)  indexes of spin=2 empty states, the state id(k) 
     !                    is saved to file KS_EMP_DW.id(k)
     !
     !   ns (integer)     number of LDA empty states to be printed 
     !   is(:) (integer)  indexes of LDA empty states, the state is(k) 
     !                    is saved to file KS_EMP.is(k)
     !
     ! Note: the first empty state has index "1" !
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     SUBROUTINE card_ksout_empty( input_line )
       ! 
       IMPLICIT NONE
       ! 
       CHARACTER(LEN=256) :: input_line
       LOGICAL, SAVE      :: tread = .FALSE.
       INTEGER            :: nksx, i, s
       TYPE occupancy_type
          INTEGER, pointer :: occs(:)
       END TYPE occupancy_type
       TYPE(occupancy_type), ALLOCATABLE :: is(:)
       !
       IF ( tread ) THEN
          CALL errore( ' card_ksout_empty ', ' two occurrence ', 2 )
       END IF
       !
       ALLOCATE ( is (nspin) )
       !
       nprnks_empty = 0 
       nksx   = 0
       !
       DO s = 1, nspin
          !
          CALL read_line( input_line )
          READ(input_line,*) nprnks_empty( s )
          !
          IF ( nprnks_empty( s ) < 1 ) THEN
             CALL errore( ' card_ksout_empty ', ' wrong number of states ', 2 )
          END IF
          !
          ALLOCATE( is(s)%occs( 1:nprnks_empty( s ) ) )
          !
          CALL read_line( input_line )
          READ(input_line,*) ( is(s)%occs( i ), i = 1, nprnks_empty( s ) )
          !
          nksx = MAX( nksx, nprnks_empty( s ) )
          !
       END DO
       ! 
       CALL allocate_input_iprnks_empty( nksx, nspin )
       !
       DO s = 1, nspin
          !
          DO i = 1, nprnks_empty( s )
             !
             iprnks_empty( i, s ) = is(s)%occs( i )
             !
          END DO
          !
          DEALLOCATE( is(s)%occs )
          !
       END DO
       !
       DEALLOCATE( is )
       ! 
       tread = .TRUE.
       !
       RETURN
       !
     END SUBROUTINE
     !
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     ! 
     ! CLIMBING_IMAGES
     !
     !   Needed to explicitly specify which images have to climb 
     !
     ! Syntax:
     !
     !   CLIMBING_IMAGES
     !     index1, ..., indexN
     !
     ! Where:
     !
     !   index1, ..., indexN are indices of the images that have to climb
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     SUBROUTINE card_climbing_images( input_line )
       !
       IMPLICIT NONE
       ! 
       CHARACTER(LEN=256) :: input_line
       LOGICAL, SAVE      :: tread = .FALSE.
       LOGICAL, EXTERNAL  :: matches
       ! 
       INTEGER          :: i
       CHARACTER(LEN=5) :: i_char
       !
       CHARACTER(LEN=6), EXTERNAL :: int_to_char
       !
       !
       IF ( tread ) &
          CALL errore( ' card_climbing_images ', ' two occurrence ', 2 )
       !
       IF ( CI_scheme == 'manual' ) THEN
          !
          IF ( ALLOCATED( climbing ) ) DEALLOCATE( climbing )
          !
          ALLOCATE( climbing( num_of_images ) )   
          !
          climbing(:) = .FALSE.
          !
          CALL read_line( input_line )
          !
          DO i = 1, num_of_images
             !
             i_char = int_to_char( i ) 
             !
             IF ( matches( ' ' // TRIM( i_char ) // ',' , & 
                           ' ' // TRIM( input_line ) // ',' ) ) &
                climbing(i) = .TRUE.
             !
          END DO
          !   
       END IF
       ! 
       tread = .TRUE.
       !
       RETURN
       !
     END SUBROUTINE card_climbing_images
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     ! 
     ! PLOT WANNIER
     !
     !   Needed to specify the indices of the wannier functions that
     !   have to be plotted
     !
     ! Syntax:
     !
     !   PLOT_WANNIER
     !     index1, ..., indexN
     !
     ! Where:
     !
     !   index1, ..., indexN are indices of the wannier functions
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     SUBROUTINE card_plot_wannier( input_line )
       !
       IMPLICIT NONE
       ! 
       CHARACTER(LEN=256) :: input_line
       LOGICAL, SAVE      :: tread = .FALSE.
       LOGICAL, EXTERNAL  :: matches
       ! 
       INTEGER                    :: i, ib
       CHARACTER(LEN=5)           :: i_char
       CHARACTER(LEN=6), EXTERNAL :: int_to_char
       !
       !
       IF ( tread ) &
          CALL errore( 'card_plot_wannier', 'two occurrence', 2 )
       !
       IF ( nwf > 0 ) THEN
          !
          IF ( nwf > nwf_max ) &
             CALL errore( 'card_plot_wannier', 'too many wannier functions', 1 )
          !
          CALL read_line( input_line )
          !
          ib = 0
          !
          DO i = 1, nwf_max
             !
             i_char = int_to_char( i ) 
             !
             IF ( matches( ' ' // TRIM( i_char ) // ',', & 
                           ' ' // TRIM( input_line ) // ',' ) ) THEN
                !
                ib = ib + 1
                !
                IF ( ib > nwf ) &
                   CALL errore( 'card_plot_wannier', 'too many indices', 1 )
                !
                wannier_index(ib) = i
                !
             END IF
             !
          END DO
          !
       END IF
       ! 
       tread = .TRUE.
       !
       RETURN
       !
     END SUBROUTINE card_plot_wannier
     !
     !------------------------------------------------------------------------
     !    BEGIN manual
     !----------------------------------------------------------------------
     !
     !
     ! TEMPLATE
     !
     !      This is a template card info section
     !
     ! Syntax:
     !
     !    TEMPLATE
     !     RVALUE IVALUE
     !
     ! Example:
     !
     !    ???
     !
     ! Where:
     !
     !      RVALUE (real)     This is a real value
     !      IVALUE (integer)  This is an integer value
     !
     !----------------------------------------------------------------------
     !    END manual
     !------------------------------------------------------------------------
     !
     SUBROUTINE card_template( input_line )
       ! 
       IMPLICIT NONE
       ! 
       CHARACTER(LEN=256) :: input_line
       LOGICAL, SAVE      :: tread = .FALSE.
       ! 
       !
       IF ( tread ) THEN
          CALL errore( ' card_template ', ' two occurrence ', 2 )
       END IF
       !
       ! ....  CODE HERE
       ! 
       tread = .TRUE.
       !
       RETURN
       !
     END SUBROUTINE
     !
END MODULE read_cards_module
