!
! Copyright (C) 2002-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
MODULE read_cards_module
   !---------------------------------------------------------------------------
   !! This module handles the reading of cards from standard input.  
   !! Original version written by Carlo Cavazzoni.
   !
   USE kinds,     ONLY : DP
   USE io_global, ONLY : stdout
   USE wy_pos,    ONLY : wypos
   USE parser,    ONLY : field_count, read_line, get_field, parse_unit
   USE io_global, ONLY : ionode, ionode_id
   !
   USE input_parameters
   !
   !
   IMPLICIT NONE
   !
   SAVE
   !
   PRIVATE
   !
   PUBLIC :: read_cards, card_kpoints
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
   SUBROUTINE card_default_values( )
      !----------------------------------------------------------------------
      !
      USE autopilot, ONLY : init_autopilot
      !
      IMPLICIT NONE
      !
      !
      ! ... mask that control the printing of selected Kohn-Sham occupied
      ! ... orbitals, default allocation
      !
      CALL allocate_input_iprnks( 0, nspin )
      nprnks  = 0
      !
      ! ... Simulation cell from standard input
      !
      trd_ht = .false.
      rd_ht  = 0.0_DP
      !
      ! ... Reference Simulation cell from standard input
      !
      ref_cell = .false.
      rd_ref_ht  = 0.0_DP
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
      ! ... k-points
      !
      k_points = 'gamma'
      tk_inp   = .false.
      nkstot   = 1
      nk1      = 0
      nk2      = 0
      nk3      = 0
      k1       = 0
      k2       = 0
      k3       = 0
      !
      ! ... Electronic states
      !
      tf_inp = .false.
      !
      ! ... ion_velocities
      !
      tavel = .false.
      !
      ! ... hubbard_card
      !
      tahub = .false. 
      !
      ! ... solvent's density initialization
      !
      solv_dens1 = 0.0_DP
      solv_dens2 = 0.0_DP
      !
      CALL init_autopilot()
      !
      RETURN
      !
   END SUBROUTINE card_default_values
   !
   !
   !----------------------------------------------------------------------
   SUBROUTINE read_cards ( prog, unit )
      !----------------------------------------------------------------------
      !
      USE autopilot, ONLY : card_autopilot
      USE upf_utils, ONLY : capital
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN), optional  :: unit
      !
      CHARACTER(len=2)           :: prog   ! calling program ( PW, CP, WA )
      CHARACTER(len=256)         :: input_line
      CHARACTER(len=80)          :: card
      LOGICAL                    :: tend
      INTEGER                    :: i, ios
      !
      ! read_line reads from unit parse_unit
      !
      IF (present(unit)) THEN
         parse_unit =  unit
      ELSE
         parse_unit =  5
      END IF
      !
      CALL card_default_values( )
      !
100   CALL read_line( input_line, end_of_file=tend )
      !
      IF( tend ) GOTO 120
      IF( input_line == ' ' .OR. input_line(1:1) == '#' .OR. &
          input_line == '/' .OR. input_line(1:1) == '!' ) GOTO 100
      !
      DO i = 1, len_trim( input_line )
         input_line( i : i ) = capital( input_line( i : i ) )
      ENDDO
      !
      READ (input_line, *, iostat=ios) card
      IF(ios/=0) card=''
      !
      IF ( trim(card) == 'AUTOPILOT' ) THEN
         !
         CALL card_autopilot( input_line )
         IF ( prog == 'PW' .and. ionode ) &
            WRITE( stdout,'(A)') 'Warning: card '//trim(input_line)//' ignored'
         !
      ELSEIF ( trim(card) == 'ATOMIC_SPECIES' ) THEN
         !
         CALL card_atomic_species( input_line )
         !
      ELSEIF ( trim(card) == 'ATOMIC_POSITIONS' ) THEN
         !
         CALL card_atomic_positions( input_line, prog )
         !
      ELSEIF ( trim(card) == 'ATOMIC_FORCES' ) THEN
         !
         CALL card_atomic_forces( input_line )
         !
      ELSEIF ( trim(card) == 'CONSTRAINTS' ) THEN
         !
         CALL card_constraints( input_line )
         !
      ELSEIF ( trim(card) == 'DIPOLE' ) THEN
         !
         CALL errore('read_cards','card DIPOLE no longer existing',1)
         !
      ELSEIF ( trim(card) == 'ESR' ) THEN
         !
         CALL errore('read_cards','card ESR no longer existing',1)
         !
      ELSEIF ( trim(card) == 'K_POINTS' ) THEN
         !
         IF ( ( prog == 'CP' ) ) THEN
            IF( ionode ) &
               WRITE( stdout,'(A)') 'Warning: card '//trim(input_line)//' ignored'
         ELSE
            CALL card_kpoints( input_line )
         ENDIF
         !
      ELSEIF ( trim(card) == 'ADDITIONAL_K_POINTS' ) THEN
         !
         CALL card_add_kpoints( input_line )
 
      ELSEIF ( trim(card) == 'OCCUPATIONS' ) THEN
         !
         CALL card_occupations( input_line )
         !
      ELSEIF ( trim(card) == 'CELL_PARAMETERS' ) THEN
         !
         CALL card_cell_parameters( input_line )
         !
      ELSEIF ( trim(card) == 'REF_CELL_PARAMETERS' ) THEN
         !
         CALL card_ref_cell_parameters( input_line )
         !
      ELSEIF ( trim(card) == 'ATOMIC_VELOCITIES' ) THEN
         !
         CALL card_ion_velocities( input_line )
         !
      ELSEIF ( trim(card) == 'KSOUT' ) THEN
         !
         CALL card_ksout( input_line )
         IF ( ( prog == 'PW' ) .and. ionode ) &
            WRITE( stdout,'(a)') 'Warning: card '//trim(input_line)//' ignored'
         !
      ELSEIF ( trim(card) == 'PLOT_WANNIER' ) THEN
         !
         CALL card_plot_wannier( input_line )

      ELSEIF ( trim(card) == 'WANNIER_AC' .and. ( prog == 'WA' )) THEN
         !
         CALL card_wannier_ac( input_line )
         !
      ELSEIF ( trim(card) == 'SOLVENTS' .AND. trism ) THEN
         !
         CALL card_solvents( input_line )
         !
      ELSEIF ( trim(card) == 'TOTAL_CHARGE' ) THEN
         !
         CALL card_total_charge( input_line )
         !
      ELSEIF ( trim(card) == 'HUBBARD' ) THEN
         !
         CALL card_hubbard( input_line )
         ! 
      ELSE
         !
         IF ( ionode ) &
            WRITE( stdout,'(A)') 'Warning: card '//trim(input_line)//' ignored'
         !
      ENDIF
      !
      ! ... END OF LOOP ... !
      !
      GOTO 100
      !
120      CONTINUE
      !
      RETURN
      !
   END SUBROUTINE read_cards

   !
   ! ... Description of the allowed input CARDS
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
   SUBROUTINE card_atomic_species( input_line )
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: input_line
      INTEGER            :: is, ip, ierr
      CHARACTER(len=6)   :: lb_pos
      CHARACTER(len=256) :: psfile
      !
      !
      IF ( taspc ) THEN
         CALL errore( ' card_atomic_species  ', ' two occurrences', 2 )
      ENDIF
      IF ( ntyp > nsx ) THEN
         CALL errore( ' card_atomic_species ', ' nsp out of range ', ntyp )
      ENDIF
      !
      DO is = 1, ntyp
         !
         CALL read_line( input_line )
         READ( input_line, *, iostat=ierr ) lb_pos, atom_mass(is), psfile
            CALL errore( ' card_atomic_species ', &
                'cannot read atomic specie from: '//trim(input_line), abs(ierr))
         atom_pfile(is) = trim( psfile )
         lb_pos         = adjustl( lb_pos )
         atom_label(is) = trim( lb_pos )
         !
         DO ip = 1, is - 1
            IF ( atom_label(ip) == atom_label(is) ) THEN
               CALL errore( ' card_atomic_species ', &
                           & ' two occurrences of the same atomic label ', is )
            ENDIF
         ENDDO
         !
      ENDDO
      taspc = .true.
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
   SUBROUTINE card_atomic_positions( input_line, prog )
      !
      USE clib_wrappers, ONLY: feval_infix
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: input_line
      CHARACTER(len=2)   :: prog
      CHARACTER(len=6)   :: lb_pos
      INTEGER            :: ia, k, is, nfield, idx, rep_i
      LOGICAL, EXTERNAL  :: matches
      LOGICAL            :: tend
      REAL(DP)           :: inp(3)
      INTEGER            :: fieldused
      !
      INTEGER            :: ifield, ierr
      REAL(DP)           :: field_value
      CHARACTER(len=256) :: field_str, error_msg, wp
      !
      !
      IF ( tapos ) THEN
         CALL errore( 'card_atomic_positions', 'two occurrences', 2 )
      ENDIF
      IF ( .not. taspc ) THEN
         CALL errore( 'card_atomic_positions', &
                     & 'ATOMIC_SPECIES must be present before', 2 )
      ENDIF
      IF ( ntyp > nsx ) THEN
         CALL errore( 'card_atomic_positions', 'nsp out of range', ntyp )
      ENDIF
      IF ( nat < 1 ) THEN
         CALL errore( 'card_atomic_positions', 'nat out of range', nat )
      ENDIF
      !
      CALL allocate_input_ions(ntyp,nat)
      !
      rd_if_pos = 1
      !
      sp_pos = 0
      rd_pos = 0.0_DP
      na_inp = 0
      lsg=.FALSE.
      !
      IF ( matches( "CRYSTAL_SG", input_line ) ) THEN
         atomic_positions = 'crystal'
         lsg=.TRUE.
      ELSEIF ( matches( "CRYSTAL", input_line ) ) THEN
         atomic_positions = 'crystal'
      ELSEIF ( matches( "BOHR", input_line ) ) THEN
         atomic_positions = 'bohr'
      ELSEIF ( matches( "ANGSTROM", input_line ) ) THEN
         atomic_positions = 'angstrom'
      ELSEIF ( matches( "ALAT", input_line ) ) THEN
         atomic_positions = 'alat'
      ELSE
         IF ( trim( adjustl( input_line ) ) /= 'ATOMIC_POSITIONS' ) THEN
            CALL errore( 'read_cards ', &
                        & 'unknown option for ATOMIC_POSITION: '&
                        & // input_line, 1 )
         ENDIF
         CALL infomsg( 'read_cards ', &
            & 'DEPRECATED: no units specified in ATOMIC_POSITIONS card' )
         IF ( prog == 'CP' ) atomic_positions = 'bohr'
         IF ( prog == 'PW' ) atomic_positions = 'alat'
         CALL infomsg( 'read_cards ', &
            & 'ATOMIC_POSITIONS: units set to '//TRIM(atomic_positions) )
      ENDIF
      !
      reader_loop : DO ia = 1,nat
         !
         CALL read_line( input_line, end_of_file = tend )
         IF ( tend ) CALL errore( 'read_cards', &
                           'end of file reading atomic positions', ia )
         !
         CALL field_count( nfield, input_line )
         !
         ! read atom symbol (column 1)
         !
         CALL get_field(1, lb_pos, input_line)
         lb_pos = trim(lb_pos)
         !
         error_msg = 'Error while parsing atomic position card.'
         !
         ! read field 2 (atom X coordinate or Wyckoff position symbol)
         !
         CALL get_field(2, field_str, input_line)
         !     
         ! Check if position ia is expressed in wyckoff positions
         !
         idx = LEN_TRIM(field_str)
         IF ( lsg .AND. (idx < 4) .AND. &
              ( IACHAR(field_str(idx:idx)) > 64 .AND. &
                IACHAR(field_str(idx:idx)) < 123 ) ) THEN
            !
            ! wyckoff positions
            !
            IF ( nfield < 3 .and. nfield > 8 ) &
            CALL errore( 'read_cards', 'wrong number of columns ' // &
                           & 'in ATOMIC_POSITIONS', ia )
            wp=field_str
            inp(:)=1.d5
            !
            DO k = 3,MIN(nfield,5)
               ! read k-th field (coordinate k-2)
               CALL get_field(k, field_str, input_line)
               inp(k-2) = feval_infix(ierr, field_str )
               CALL errore('card_atomic_positions', error_msg, ierr)
            ENDDO
            !
            CALL wypos(rd_pos(1,ia),wp,inp,space_group, &
                 uniqueb,rhombohedral,origin_choice)
            !
            ! count how many fields were used to find wyckoff positions
            !
            fieldused=2
            IF ( ANY (rd_pos(1:3,ia)==inp(1)) ) fieldused=fieldused+1
            IF ( ANY (rd_pos(2:3,ia)==inp(2)) ) fieldused=fieldused+1
            IF (      rd_pos(  3,ia)==inp(3)  ) fieldused=fieldused+1
            !
         ELSE
            !
            ! no wyckoff positions 
            !
            IF ( nfield /= 4 .and. nfield /= 7 ) &
            CALL errore( 'read_cards', 'wrong number of columns ' // &
                           & 'in ATOMIC_POSITIONS', ia )
            !
            ! field just read is coordinate X
            !
            rd_pos(1,ia) = feval_infix(ierr, field_str )
            CALL errore('card_atomic_positions', error_msg, ierr)
            DO k = 3,4
               ! read fields 3 and 4 (atom Y and Z coordinate)
               CALL get_field(k, field_str, input_line)
               rd_pos(k-1,ia) = feval_infix(ierr, field_str )
               CALL errore('card_atomic_positions', error_msg, ierr)
            END DO
            !
            fieldused=4
            !
         ENDIF
         ! read constraints if present (last 3 fields)
         IF ( nfield-fieldused > 0 .AND. nfield-fieldused /= 3 ) &
            CALL errore( 'read_cards', 'unexpected number of columns ' // &
                           & 'in ATOMIC_POSITIONS', ia )
         DO k = fieldused+1, nfield
            CALL get_field(k, field_str, input_line)
            READ(field_str, *) rd_if_pos(k-fieldused,ia)
         ENDDO
         !
         match_label: DO is = 1, ntyp
            !
            IF ( trim(lb_pos) == trim( atom_label(is) ) ) THEN
               !
               sp_pos(ia) = is
               exit match_label
               !
            ENDIF
            !
         ENDDO match_label
         !
         IF( ( sp_pos(ia) < 1 ) .or. ( sp_pos(ia) > ntyp ) ) THEN
            !
            CALL errore( 'read_cards', 'species '//trim(lb_pos)// &
                           & ' in ATOMIC_POSITIONS is nonexistent', ia )
            !
         ENDIF
         !
         is = sp_pos(ia)
         !
         na_inp(is) = na_inp(is) + 1
         !
      ENDDO reader_loop
      !
      tapos = .true.
      !

      RETURN
      !
   END SUBROUTINE card_atomic_positions
   !
   !------------------------------------------------------------------------
   !    BEGIN manual
   !----------------------------------------------------------------------
   !
   ! ATOMIC_FORCES
   !
   !   read external forces (in atomic units) from standard input
   !
   ! Syntax:
   !
   !   ATOMIC_FORCES
   !     label Fx(1) Fy(1) Fz(1)
   !     .....
   !     label Fx(n) Fy(n) Fz(n)
   !
   ! Example:
   !
   !   ???
   !
   ! Where:
   !
   !   label (character(len=4))       atomic label
   !   Fx(:), Fy(:) and Fz(:) (REAL)  x, y and z component of the external force
   !                                  acting on the ions whose coordinate are given
   !                                  in the same line in card ATOMIC_POSITION
   !
   !----------------------------------------------------------------------
   !    END manual
   !------------------------------------------------------------------------
   !
   SUBROUTINE card_atomic_forces( input_line )
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: input_line
      INTEGER            :: ia, k, nfield
      CHARACTER(len=4)   :: lb
      !
      !
      IF( tforces ) THEN
         CALL errore( ' card_atomic_forces ', ' two occurrences ', 2 )
      ENDIF
      !
      IF( .not. tapos ) THEN
         CALL errore( ' card_atomic_forces ', &
                     & ' ATOMIC_SPECIES must be present before ', 2 )
      ENDIF
      !
      rd_for = 0.0_DP
      !
      DO ia = 1, nat
         !
         CALL read_line( input_line )
         CALL field_count( nfield, input_line )
         IF ( nfield == 4 ) THEN
            READ(input_line,*) lb, ( rd_for(k,ia), k = 1, 3 )
         ELSEIF( nfield == 3 ) THEN
            READ(input_line,*) ( rd_for(k,ia), k = 1, 3 )
         ELSE
            CALL errore( ' iosys ', ' wrong entries in ATOMIC_FORCES ', ia )
         ENDIF
         !
      ENDDO
      !
      tforces = .true.
      !
      RETURN
      !
   END SUBROUTINE card_atomic_forces
   !
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
   !   mesh_option == tpiba_b    as tpiba but the weights gives the
   !                             number of points between this point
   !                             and the next
   !   mesh_option == crystal_b  as crystal but the weights gives the
   !                             number of points between this point and
   !                             the next
   !   mesh_option == tpiba_c    the code expects three k points 
   !                             k_0, k_1, k_2 in tpiba units.
   !                             These points define a rectangle
   !                             in reciprocal space with vertices k_0, k_1,
   !                             k_2, k_1+k_2-k_0:  k_0 + \alpha (k_1-k_0)+
   !                             \beta (k_2-k_0) with 0<\alpha,\beta < 1. 
   !                             The code produces a uniform mesh n1 x n2 
   !                             k points in this rectangle. n1 and n2 are 
   !                             the weights of k_1 and k_2. The weight of k_0
   !                             is not used. Useful for contour plots of the 
   !                             bands.
   !   mesh_option == crystal_c  as tpiba_c but the k points are given
   !                             in crystal coordinates.
   ! 
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
      USE bz_form, ONLY : transform_label_coord
      USE cell_base, ONLY : cell_base_init, celldm_cb => celldm
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: input_line, buffer
      INTEGER            :: i, j
      INTEGER            :: nkaux, ierr
      INTEGER, ALLOCATABLE :: wkaux(:)
      REAL(DP), ALLOCATABLE :: xkaux(:,:)
      INTEGER            :: npk_label, nch
      CHARACTER(LEN=3), ALLOCATABLE :: letter(:)
      INTEGER, ALLOCATABLE :: label_list(:)
      REAL(DP) :: delta, wk0
      REAL(DP) :: dkx(3), dky(3)
      LOGICAL, EXTERNAL  :: matches
      LOGICAL            :: tend,terr
      LOGICAL            :: kband = .false.
      LOGICAL            :: kband_plane = .false.
      !
      !
      IF ( tkpoints ) THEN
         CALL errore( ' card_kpoints ', ' two occurrences', 2 )
      ENDIF
      !
      IF ( matches( "AUTOMATIC", input_line ) ) THEN
         !  automatic generation of k-points
         k_points = 'automatic'
      ELSEIF ( matches( "CRYSTAL", input_line ) ) THEN
         !  input k-points are in crystal (reciprocal lattice) axis
         k_points = 'crystal'
         IF ( matches( "_B", input_line ) ) kband=.true.
         IF ( matches( "_C", input_line ) ) kband_plane=.true.
      ELSEIF ( matches( "TPIBA", input_line ) ) THEN
         !  input k-points are in 2pi/a units
         k_points = 'tpiba'
         IF ( matches( "_B", input_line ) ) kband=.true.
         IF ( matches( "_C", input_line ) ) kband_plane=.true.
      ELSEIF ( matches( "GAMMA", input_line ) ) THEN
         !  Only Gamma (k=0) is used
         k_points = 'gamma'
      ELSE
         !  by default, input k-points are in 2pi/a units
         k_points = 'tpiba'
      ENDIF
      !
      IF ( k_points == 'automatic' ) THEN
         !
         ! ... automatic generation of k-points
         !
         nkstot = 0
         CALL read_line( input_line, end_of_file = tend, error = terr )
         IF (tend) GOTO 10
         IF (terr) GOTO 20
         READ(input_line, *, END=10, ERR=20) nk1, nk2, nk3, k1, k2 ,k3
         IF ( k1 < 0 .or. k1 > 1 .or. &
               k2 < 0 .or. k2 > 1 .or. &
               k3 < 0 .or. k3 > 1 ) CALL errore &
                  ('card_kpoints', 'invalid offsets: must be 0 or 1', 1)
         IF ( nk1 <= 0 .or. nk2 <= 0 .or. nk3 <= 0 ) CALL errore &
                  ('card_kpoints', 'invalid values for nk1, nk2, nk3', 1)
         ALLOCATE ( xk(3,1), wk(1) ) ! prevents problems with debug flags
         !                           ! when init_startk is called in iosys
      ELSEIF ( ( k_points == 'tpiba' ) .or. ( k_points == 'crystal' ) ) THEN
         !
         ! ... input k-points 
         !
         CALL read_line( input_line, end_of_file = tend, error = terr )
         IF (tend) GOTO 10
         IF (terr) GOTO 20
         READ(input_line, *, END=10, ERR=20) nkstot
         IF ( nkstot <= 0 ) GO TO 20
         !
         IF (kband) THEN
!
!        Only the initial and final k points of the lines are given
!
            nkaux=nkstot
            ALLOCATE(xkaux(3,nkstot), wkaux(nkstot))
            ALLOCATE ( letter(nkstot) )
            ALLOCATE ( label_list(nkstot) )
            npk_label=0
            DO i = 1, nkstot
               CALL read_line( input_line, end_of_file = tend, error = terr )
               IF (tend) GOTO 10
               IF (terr) GOTO 20
               DO j=1,256   ! loop over all characters of input_line
                  IF ((ICHAR(input_line(j:j)) < 58 .AND. &   ! a digit
                       ICHAR(input_line(j:j)) > 47) &
                 .OR. ICHAR(input_line(j:j)) == 43 .OR. &    ! the + sign
                      ICHAR(input_line(j:j))== 45 .OR. &     ! the - sign
                      ICHAR(input_line(j:j))== 46 ) THEN     ! a dot .
!
!   This is a digit, therefore this line contains the coordinates of the
!   k point. We read it and exit from the loop on the characters
!
                     READ(input_line,*, END=10, ERR=20) xkaux(1,i), &
                                           xkaux(2,i), xkaux(3,i), wk0
                     wkaux(i) = NINT ( wk0 ) ! beware: wkaux is integer
                     EXIT
                  ELSEIF ((ICHAR(input_line(j:j)) < 123 .AND. &
                           ICHAR(input_line(j:j)) > 64))  THEN
!
!   This is a letter, not a space character. We read the next three 
!   characters and save them in the letter array, save also which k point
!   it is
!
                     npk_label=npk_label+1
                     READ(input_line(j:),'(a3)') letter(npk_label)
                     label_list(npk_label)=i
!
!  now we remove the letters from input_line and read the number of points
!  of the line. The next two line should account for the case in which
!  there is only one space between the letter and the number of points.
!
                     nch=3
                     IF ( ICHAR(input_line(j+1:j+1))==32 .OR. &
                          ICHAR(input_line(j+2:j+2))==32 ) nch=2
                     buffer=input_line(j+nch:)
                     READ(buffer,*,err=20) wkaux(i)
                     EXIT
                  ENDIF
               ENDDO
            ENDDO
            IF ( npk_label > 0 ) THEN
               CALL cell_base_init ( ibrav, celldm, a, b, c, cosab, &
                              cosac, cosbc, trd_ht, rd_ht, cell_units )
               CALL transform_label_coord(ibrav, celldm_cb, xkaux, letter, &
                    label_list, npk_label, nkstot, k_points, point_label_type )
            END IF

            DEALLOCATE(letter)
            DEALLOCATE(label_list)
            ! Count k-points first
            nkstot=SUM(wkaux(1:nkaux-1))+1
            DO i=1,nkaux-1
              IF (wkaux(i)==0) nkstot=nkstot+1
            ENDDO
            ALLOCATE ( xk(3,nkstot), wk(nkstot) )
            !
            !  generate the points along the lines
            !
            CALL generate_k_along_lines(nkaux, xkaux, wkaux, xk, wk, nkstot)
            !
            !  workaround: discard current wk (contains the length of k-path, 
            !  never used), replace with wk=1 so that band occupations (wg)
            !  are correctly written to file - needed by BerkeleyGW interface
            !
            wk(:) = 1.0_dp
            DEALLOCATE(xkaux)
            DEALLOCATE(wkaux)
            !
         ELSEIF (kband_plane) THEN
!
!        Generate a uniform mesh of k points on the plane defined by
!        the origin k_0, and two vectors applied in k_0, k_1 and k_2.
!
            IF (nkstot /= 3) CALL errore ('card_kpoints', &
                                'option _c requires 3 k points',i)
            nkaux=nkstot
            ALLOCATE(xkaux(3,nkstot), wkaux(nkstot))
            DO i = 1, nkstot
               CALL read_line( input_line, end_of_file = tend, error = terr )
               IF (tend) GOTO 10
               IF (terr) GOTO 20
               READ(input_line,*, END=10, ERR=20) xkaux(1,i), xkaux(2,i), &
                                                  xkaux(3,i), wk0
               wkaux(i) = NINT ( wk0 ) ! beware: wkaux is integer
            ENDDO
            ! Count k-points first
            nkstot = wkaux(2) * wkaux(3)
            ALLOCATE ( xk(3,nkstot), wk(nkstot) )
            CALL generate_k_in_plane(nkaux, xkaux, wkaux, xk, wk, nkstot)
            DEALLOCATE(xkaux)
            DEALLOCATE(wkaux)
         ELSE
!
!    Reads on input the k points
!
            ALLOCATE ( xk(3, nkstot), wk(nkstot), labelk(nkstot) )
            DO i = 1, nkstot
               labelk(i) = ''
               !
               CALL read_line( input_line, end_of_file = tend, error = terr )
               IF (tend) GOTO 10
               IF (terr) GOTO 20
               !
               ! Try to read with optional label
               READ(input_line,*, IOSTAT=ierr) xk(1,i),xk(2,i),xk(3,i),wk(i),labelk(i)
               !
               IF (ierr /= 0 ) &
                  READ(input_line,*, END=10, ERR=20) xk(1,i),xk(2,i),xk(3,i),wk(i)
               !
            ENDDO
         ENDIF
         !
      ELSEIF ( k_points == 'gamma' ) THEN
         !
         nkstot = 1
         ALLOCATE ( xk(3,1), wk(1) )
         xk(:,1) = 0.0_DP
         wk(1) = 1.0_DP
         !
      ENDIF
      !
      tkpoints  = .true.
      tk_inp = .true.
      !
      RETURN
10     CALL errore ('card_kpoints', ' end of file while reading ' &
            & // trim(k_points) // ' k points', 1)
20     CALL errore ('card_kpoints', ' error while reading ' &
            & // trim(k_points) // ' k points', 1)
      !
   END SUBROUTINE card_kpoints

   SUBROUTINE card_add_kpoints( input_line )
     USE additional_kpoints, ONLY : nkstot_add, xk_add, k_points_add
     IMPLICIT NONE
     CHARACTER(len=*),INTENT(in) :: input_line
     CHARACTER(len=256) :: input_line_aux
     REAL(DP),ALLOCATABLE :: xk_old(:,:), wk_old(:)
     INTEGER :: nk1_old, nk2_old, nk3_old, nkstot_old
     INTEGER :: k1_old,  k2_old,  k3_old
     LOGICAL, EXTERNAL  :: matches
     CHARACTER(len=80) :: k_points_old
     !
     IF(.not.allocated(xk) .or. .not.allocated(wk))&
       CALL errore("add_kpoints", "ADDITIONAL_K_POINTS must appear after K_POINTS",1)
     IF(.not.tkpoints) &
       CALL errore("add_kpoints", "ADDITIONAL_K_POINTS must appear after K_POINTS",2)
     IF(matches( "AUTOMATIC", input_line )) &
       CALL errore("add_kpoints", "ADDITIONAL_K_POINTS cannot be 'automatic'", 3)

     ! Back-up existing points
     nkstot_old = nkstot
     ALLOCATE(xk_old(3,nkstot_old))
     ALLOCATE(wk_old(nkstot_old))
     k_points_old = k_points
     xk_old  = xk
     wk_old  = wk
     nk1_old = nk1
     nk2_old = nk2
     nk3_old = nk3
     k1_old  = k1
     k2_old  = k2
     k3_old  = k3
     DEALLOCATE(xk,wk)

     ! Prepare to read k-points again
     nkstot = 0
     input_line_aux = TRIM(ADJUSTL(input_line))
     input_line_aux = input_line_aux(12:)
     tkpoints = .false.
     CALL card_kpoints(input_line_aux)
     !
     ! Backup new points to module
     nkstot_add = nkstot
     IF(nkstot_add==0) CALL errore("add_kpoints", "No new k_points?",1)
     ALLOCATE(xk_add(3,nkstot_add))
     xk_add = xk
     k_points_add = k_points

     ! Put back previous stuff
     DEALLOCATE(xk, wk)
     nkstot = nkstot_old
     ALLOCATE(xk(3,nkstot))
     ALLOCATE(wk(nkstot))
     k_points = k_points_old
     xk  = xk_old
     wk  = wk_old
     nk1 = nk1_old
     nk2 = nk2_old
     nk3 = nk3_old
     k1  = k1_old
     k2  = k2_old
     k3  = k3_old
     DEALLOCATE(xk_old,wk_old)

     RETURN 
   END SUBROUTINE card_add_kpoints
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
      USE clib_wrappers, ONLY: feval_infix
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: input_line, field_str
      INTEGER            :: is, nx10, i, j, nspin0
      INTEGER            :: nfield, nbnd_read, nf, ierr
      LOGICAL :: tef
      !
      !
      IF ( tocc ) THEN
         CALL errore( ' card_occupations ', ' two occurrences', 2 )
      ENDIF
      nspin0=nspin
      IF (nspin == 4) nspin0=1
      !
      ALLOCATE ( f_inp ( nbnd, nspin0 ) )
      DO is = 1, nspin0
         !
         nbnd_read = 0
         DO WHILE ( nbnd_read < nbnd)
            CALL read_line( input_line, end_of_file=tef )
            IF (tef) CALL errore('card_occupations',&
                        'Missing occupations, end of file reached',1)
            CALL field_count( nfield, input_line )
            !
            DO nf = 1,nfield
               nbnd_read = nbnd_read+1
               IF (nbnd_read > nbnd ) EXIT
               CALL get_field(nf, field_str, input_line)
               !
               f_inp(nbnd_read,is) = feval_infix(ierr, field_str )
               CALL errore('card_occupations',&
                  'Error parsing occupation: '//trim(field_str), nbnd_read*ierr)
            ENDDO
         ENDDO
         !
      ENDDO
      !
      tf_inp = .true.
      tocc = .true.
      !
      RETURN
      !
   END SUBROUTINE card_occupations
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
   !    CELL_PARAMETERS (cell_option)
   !      HT(1,1) HT(1,2) HT(1,3)
   !      HT(2,1) HT(2,2) HT(2,3)
   !      HT(3,1) HT(3,2) HT(3,3)
   !
   !   cell_option == alat      lattice vectors in units of alat
   !   cell_option == bohr      lattice vectors in Bohr
   !   cell_option == angstrom  lattice vectors in Angstrom
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
      CHARACTER(len=256) :: input_line
      INTEGER            :: i, j
      LOGICAL, EXTERNAL  :: matches
      !
      !
      IF ( tcell ) THEN
         CALL errore( ' card_cell_parameters ', ' two occurrences', 2 )
      ENDIF
      !
      IF ( matches( "BOHR", input_line ) ) THEN
         cell_units = 'bohr'
      ELSEIF ( matches( "ANGSTROM", input_line ) ) THEN
         cell_units = 'angstrom'
      ELSEIF ( matches( "ALAT", input_line ) ) THEN
         cell_units = 'alat'
      ELSE
         cell_units = 'none'
         CALL infomsg( 'read_cards ', &
            & 'DEPRECATED: no units specified in CELL_PARAMETERS card' )
         ! Cell parameters are set in cell_base_init
      ENDIF
      !
      DO i = 1, 3
         CALL read_line( input_line )
         READ(input_line,*) ( rd_ht( i, j ), j = 1, 3 )
      ENDDO
      !
      trd_ht = .true.
      tcell  = .true.
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
   ! REF_CELL_PARAMETERS
   !
   !   use the specified cell dimensions
   !
   ! Syntax:
   !
   !    REF_CELL_PARAMETERS (cell_option)
   !      rd_ref_HT(1,1) rd_ref_HT(1,2) rd_ref_HT(1,3)
   !      rd_ref_HT(2,1) rd_ref_HT(2,2) rd_ref_HT(2,3)
   !      rd_ref_HT(3,1) rd_ref_HT(3,2) rd_ref_HT(3,3)
   !
   !   cell_option == alat      lattice vectors in units of alat set by ref_alat keyword (default)
   !   cell_option == bohr      lattice vectors in Bohr
   !   cell_option == angstrom  lattice vectors in Angstrom
   !
   ! Example:
   !
   ! REF_CELL_PARAMETERS
   !    24.50644311    0.00004215   -0.14717844
   !    -0.00211522    8.12850030    1.70624903
   !     0.16447787    0.74511792   23.07395418
   !
   ! Where:
   !
   !      rd_ref_HT(i,j) (real)  cell dimensions ( in a.u. ),
   !        note the relation with reference lattice vectors:
   !        rd_ref_HT(1,:) = ref_A1, rd_ref_HT(2,:) = ref_A2, rd_ref_HT(3,:) = re_A3
   !
   !----------------------------------------------------------------------
   !    END manual
   !------------------------------------------------------------------------
   !
   SUBROUTINE card_ref_cell_parameters( input_line )
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: input_line
      INTEGER            :: i, j
      LOGICAL, EXTERNAL  :: matches
      !
      !
      IF ( ref_cell ) THEN
         CALL errore( ' card_reference_cell_parameters ', ' two occurrences', 2 )
      ENDIF
      !
      IF ( matches( "BOHR", input_line ) ) THEN
         ref_cell_units = 'bohr'
      ELSEIF ( matches( "ANGSTROM", input_line ) ) THEN
         ref_cell_units = 'angstrom'
      ELSE
         ref_cell_units = 'alat'
      ENDIF
      !
      DO i = 1, 3
         CALL read_line( input_line )
         READ(input_line,*) ( rd_ref_ht( i, j ), j = 1, 3 )
      ENDDO
      !
      ref_cell = .true.
      !
      RETURN
      !
   END SUBROUTINE card_ref_cell_parameters
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
      CHARACTER(len=256) :: input_line
      INTEGER            :: ia, k, is, nfield
      CHARACTER(len=4)   :: lb_vel
      !
      !
      IF( tionvel ) THEN
         CALL errore( ' card_ion_velocities ', ' two occurrences', 2 )
      ENDIF
      !
      IF( .not. tapos ) THEN
         CALL errore( ' card_ion_velocities ', &
                     & ' ATOMIC_SPECIES must be present before ', 2 )
      ENDIF
      !
      rd_vel = 0.0_DP
      sp_vel = 0
      !
      IF ( ion_velocities == 'from_input' ) THEN
         !
         tavel = .true.
         !
         DO ia = 1, nat
            !
            CALL read_line( input_line )
            CALL field_count( nfield, input_line )
            IF ( nfield == 4 ) THEN
               READ(input_line,*) lb_vel, ( rd_vel(k,ia), k = 1, 3 )
            ELSE
               CALL errore( ' iosys ', &
                           & ' wrong entries in ATOMIC_VELOCITIES ', ia )
            ENDIF
            !
            match_label: DO is = 1, ntyp
               IF ( trim( lb_vel ) == atom_label(is) ) THEN
                  sp_vel(ia) = is
                  exit match_label
               ENDIF
            ENDDO match_label
            !
            IF ( sp_vel(ia) < 1 .or. sp_vel(ia) > ntyp ) THEN
               CALL errore( ' iosys ', ' wrong LABEL in ION_VELOCITIES ', ia )
            ENDIF
            !
         ENDDO
         !
      ENDIF
      !
      tionvel = .true.
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
      CHARACTER(len=256) :: input_line
      INTEGER            :: i, nfield
      !
      !
      IF ( tconstr ) CALL errore( 'card_constraints', 'two occurrences', 2 )
      !
      CALL read_line( input_line )
      !
      CALL field_count( nfield, input_line )
      !
      IF ( nfield == 1 ) THEN
         !
         READ( input_line, * ) nconstr_inp
         !
      ELSEIF ( nfield == 2 ) THEN
         !
         READ( input_line, * ) nconstr_inp, constr_tol_inp
         !
      ELSE
         !
         CALL errore( 'card_constraints', 'too many fields', nfield )
         !
      ENDIF
      WRITE(stdout,'(5x,a,i4,a,f12.6)') &
         'Reading',nconstr_inp,' constraints; tolerance:', constr_tol_inp
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
               WRITE(stdout,'(7x,i3,a,i3,a,i2,a,2f12.6)') i, &
                  ') '//constr_type_inp(i)(1:4),int(constr_inp(1,i)) ,&
                  ' coordination wrt type:', int(constr_inp(2,i)), &
                  ' cutoff distance and smoothing:',  constr_inp(3:4,i)
            ELSEIF ( nfield == 6 ) THEN
               !
               READ( input_line, * ) constr_type_inp(i), &
                                    constr_inp(1,i), &
                                    constr_inp(2,i), &
                                    constr_inp(3,i), &
                                    constr_inp(4,i), &
                                    constr_target_inp(i)
               !
               constr_target_set(i) = .true.
               !
               WRITE(stdout,'(7x,i3,a,i3,a,i2,a,2f12.6,a,f12.6)') i, &
                  ') '//constr_type_inp(i)(1:4),int(constr_inp(1,i)) , &
                  ' coordination wrt type:', int(constr_inp(2,i)), &
                  ' cutoff distance and smoothing:',  constr_inp(3:4,i), &
                  '; target:', constr_target_inp(i)
            ELSE
               !
               CALL errore( 'card_constraints', 'type_coord, ' // &
                           & 'atom_coord: wrong number of fields', nfield )
               !
            ENDIF
            !
         CASE( 'distance' )
            !
            IF ( nfield == 3 ) THEN
               !
               READ( input_line, * ) constr_type_inp(i), &
                                    constr_inp(1,i), &
                                    constr_inp(2,i)
               !
               WRITE(stdout,'(7x,i3,a,2i3)') &
                  i,') distance between atoms: ', int(constr_inp(1:2,i))
            ELSEIF ( nfield == 4 ) THEN
               !
               READ( input_line, * ) constr_type_inp(i), &
                                    constr_inp(1,i), &
                                    constr_inp(2,i), &
                                    constr_target_inp(i)
               !
               constr_target_set(i) = .true.
               !
               WRITE(stdout,'(7x,i3,a,2i3,a,f12.6)') i, &
                  ') distance between atoms: ', int(constr_inp(1:2,i)), &
                  '; target:',  constr_target_inp(i)
            ELSE
               !
               CALL errore( 'card_constraints', &
                           & 'distance: wrong number of fields', nfield )
               !
            ENDIF
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
               WRITE(stdout, '(7x,i3,a,3i3)') &
                  i,') planar angle between atoms: ', int(constr_inp(1:3,i))
            ELSEIF ( nfield == 5 ) THEN
               !
               READ( input_line, * ) constr_type_inp(i), &
                                    constr_inp(1,i), &
                                    constr_inp(2,i), &
                                    constr_inp(3,i), &
                                    constr_target_inp(i)
               !
               constr_target_set(i) = .true.
               !
               WRITE(stdout, '(7x,i3,a,3i3,a,f12.6)') i, &
                  ') planar angle between atoms: ', int(constr_inp(1:3,i)), &
                  '; target:', constr_target_inp(i)
            ELSE
               !
               CALL errore( 'card_constraints', &
                           & 'planar_angle: wrong number of fields', nfield )
               !
            ENDIF
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
               WRITE(stdout, '(7x,i3,a,4i3)') &
                  i,') torsional angle between atoms: ', int(constr_inp(1:4,i))
            ELSEIF ( nfield == 6 ) THEN
               !
               READ( input_line, * ) constr_type_inp(i), &
                                    constr_inp(1,i), &
                                    constr_inp(2,i), &
                                    constr_inp(3,i), &
                                    constr_inp(4,i), &
                                    constr_target_inp(i)
               !
               constr_target_set(i) = .true.
               !
               WRITE(stdout, '(7x,i3,a,4i3,a,f12.6)') i, &
                  ') torsional angle between atoms: ', int(constr_inp(1:4,i)),&
                  '; target:', constr_target_inp(i)
            ELSE
               !
               CALL errore( 'card_constraints', &
                           & 'torsional_angle: wrong number of fields', nfield )
               !
            ENDIF
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
               WRITE(stdout, '(7x,i3,a,i3,a,3f12.6)') i, &
                  ') bennet projection of atom ', int(constr_inp(1,i)), &
                  ' along vector:', constr_inp(2:4,i)
            ELSEIF ( nfield == 6 ) THEN
               !
               READ( input_line, * ) constr_type_inp(i), &
                                    constr_inp(1,i), &
                                    constr_inp(2,i), &
                                    constr_inp(3,i), &
                                    constr_inp(4,i), &
                                    constr_target_inp(i)
               !
               constr_target_set(i) = .true.
               !
               WRITE(stdout, '(7x,i3,a,i3,a,3f12.6,a,f12.6)') i, &
                  ') bennet projection of atom ', int(constr_inp(1,i)), &
                  ' along vector:', constr_inp(2:4,i), &
                  '; target:', constr_target_inp(i)
            ELSE
               !
               CALL errore( 'card_constraints', &
                           & 'bennett_proj: wrong number of fields', nfield )
               !
            ENDIF
            !
         CASE( 'potential_wall' )
            !
            IF ( nfield == 4 ) THEN
               !
               READ( input_line, * ) constr_type_inp(i), &
                                    constr_inp(1,i), &
                                    constr_inp(2,i), &
                                    constr_inp(3,i)
               !
               WRITE(stdout, '(7x,i3,a)') &
                  i,') potential wall at origin normal to z-axis is requested'
               WRITE(stdout, '(9x,a)') 'External force is proportional to:'
               WRITE(stdout, '(11x,f12.6,a,f12.6,a,f12.6,a)') constr_inp(1,i), &
                  ' (a.u.) * ', constr_inp(2,i), ' (a.u.) * exp(', &
                  (-1._dp) * constr_inp(2,i), ').'
               WRITE(stdout, '(9x,a,f12.6,a)') 'Force is applied when atom is within ',&
                  constr_inp(3,i), ' (a.u.) from the wall.'
            !
            ELSE
               !
               CALL errore( 'card_constraints', &
                           & 'potential_wall: wrong number of fields', nfield )
               !
            ENDIF
            !
         CASE DEFAULT
            !
            CALL errore( 'card_constraints', 'unknown constraint ' // &
                        & 'type: ' // trim( constr_type_inp(i) ), 1 )
            !
         END SELECT
         !
      ENDDO
      !
      tconstr = .true.
      !
      RETURN
      !
   END SUBROUTINE card_constraints
   !
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
      CHARACTER(len=256) :: input_line
      INTEGER            :: i, s, nksx
      TYPE occupancy_type
         INTEGER, POINTER :: occs(:)
      END TYPE occupancy_type
      TYPE(occupancy_type), ALLOCATABLE :: is(:)
      !
      IF ( tksout ) THEN
         CALL errore( ' card_ksout ', ' two occurrences', 2 )
      ENDIF
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
         ENDIF
         !
         ALLOCATE( is(s)%occs( 1:nprnks(s) ) )
         !
         CALL read_line( input_line )
         READ(input_line, *) ( is(s)%occs(i), i = 1, nprnks( s ) )
         !
         nksx = max( nksx, nprnks( s ) )
         !
      ENDDO
      !
      CALL allocate_input_iprnks( nksx, nspin )
      !
      DO s = 1, nspin
         !
         DO i = 1, nprnks( s )
            !
            iprnks( i, s ) = is(s)%occs(i)
            !
         ENDDO
         !
         DEALLOCATE( is(s)%occs )
         !
      ENDDO
      !
      DEALLOCATE( is )
      !
      tksout = .true.
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
      CHARACTER(len=256) :: input_line
      LOGICAL, EXTERNAL  :: matches
      !
      INTEGER                    :: i, ib
      CHARACTER(len=6)           :: i_char
      CHARACTER(len=6), EXTERNAL :: int_to_char
      !
      !
      IF ( twannier ) &
         CALL errore( 'card_plot_wannier', 'two occurrences', 2 )
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
            IF ( matches( ' ' // trim( i_char ) // ',', &
                           ' ' // trim( input_line ) // ',' ) ) THEN
               !
               ib = ib + 1
               !
               IF ( ib > nwf ) &
                  CALL errore( 'card_plot_wannier', 'too many indices', 1 )
               !
               wannier_index(ib) = i
               !
            ENDIF
            !
         ENDDO
         !
      ENDIF
      !
      twannier = .true.
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
      CHARACTER(len=256) :: input_line
      !
      !
      IF ( ttemplate ) THEN
         CALL errore( ' card_template ', ' two occurrences', 2 )
      ENDIF
      !
      ! ....  CODE HERE
      !
      ttemplate = .true.
      !
      RETURN
      !
   END SUBROUTINE
   !
   !
   !------------------------------------------------------------------------
   !    BEGIN manual
   !----------------------------------------------------------------------
   !WANNIER_AC
   !Wannier# 1 10.5 15.7 2
   !atom 1
   !d 1 0.45
   !p 3 0.55
   !Wannier# 2 10.5 15.7 1
   !atom 3
   !p 1 0.8
   !Spin#2:
   !Wannier# 1 10.5 15.7 2
   !atom 1
   !d 1 0.45
   !p 3 0.55
   !Wannier# 2 10.5 15.7 1
   !atom 3
   !p 1 0.8
   !----------------------------------------------------------------------
   !    END manual
   !------------------------------------------------------------------------
   !
   SUBROUTINE card_wannier_ac( input_line )
      !
      USE wannier_new, ONLY: nwan

      IMPLICIT NONE
      !
      CHARACTER(len=256) :: input_line
      INTEGER :: i,j,k, nfield, iwan, ning, iatom,il,im,ispin
      LOGICAL :: tend
      REAL :: c, b_from, b_to
      CHARACTER(len=10) :: text, lo

      ispin = 1
      !
      DO i = 1, nwan
         !
         CALL read_line( input_line, end_of_file = tend )
         !
         IF ( tend ) &
            CALL errore( 'read_cards', &
                        'end of file reading trial wfc composition', i )
         !
         CALL field_count( nfield, input_line )
         !
         IF ( nfield == 4 ) THEN
            READ(input_line,*) text, iwan, b_from, b_to
            ning = 1
         ELSEIF ( nfield == 5 ) THEN
            READ(input_line,*) text, iwan, b_from, b_to, ning
         ELSE
            CALL errore( 'read_cards', &
                        'wrong format', nfield )
         ENDIF
         IF(iwan/=i) CALL errore( 'read_cards', 'wrong wannier order', iwan)

         ! Read atom number
         CALL read_line( input_line, end_of_file = tend )
         READ(input_line,*) text, iatom
         !
         wan_data(iwan,ispin)%iatom = iatom
         wan_data(iwan,ispin)%ning = ning
         wan_data(iwan,ispin)%bands_from = b_from
         wan_data(iwan,ispin)%bands_to = b_to
         !
         DO j=1, ning
            CALL read_line( input_line, end_of_file = tend )
            !
            IF ( tend ) &
               CALL errore( 'read_cards', &
                           'not enough wavefunctions', j )
            IF (ning==1) THEN
               READ(input_line,*) lo,im
               c = 1.d0
            ELSE
               READ(input_line,*) lo,im,c
            ENDIF

            SELECT CASE(trim(lo))
            CASE('s')
               il = 0
            CASE('p')
               il = 1
            CASE('d')
               il = 2
            CASE('f')
               il = 3
            CASE DEFAULT
               CALL errore( 'read_cards', &
                           'wrong l-label', 1 )
            END SELECT

            wan_data(iwan,ispin)%ing(j)%l = il
            wan_data(iwan,ispin)%ing(j)%m = im
            wan_data(iwan,ispin)%ing(j)%c = c
         ENDDO
      ENDDO

      !Is there spin 2 information?
      CALL read_line( input_line, end_of_file = tend )
      !
      IF ( .not. tend ) THEN
         READ(input_line,*) text
         IF ( trim(text) == 'Spin#2:') THEN ! ok, there is spin 2 data
            ispin = 2
            !
            DO i = 1, nwan
               !
               CALL read_line( input_line, end_of_file = tend )
               !
               IF ( tend ) &
                  CALL errore( 'read_cards', &
                              'end of file reading trial wfc composition', i )
               !
               CALL field_count( nfield, input_line )
               !
               IF ( nfield == 4 ) THEN
                  READ(input_line,*) text, iwan, b_from, b_to
                  ning = 1
               ELSEIF ( nfield == 5 ) THEN
                  READ(input_line,*) text, iwan, b_from, b_to, ning
               ELSE
                  CALL errore( 'read_cards', &
                              'wrong format', nfield )
               ENDIF
               IF(iwan/=i) CALL errore( 'read_cards', 'wrong wannier order', iwan)

               ! Read atom number
               CALL read_line( input_line, end_of_file = tend )
               READ(input_line,*) text, iatom
               !
               wan_data(iwan,ispin)%iatom = iatom
               wan_data(iwan,ispin)%ning = ning
               wan_data(iwan,ispin)%bands_from = b_from
               wan_data(iwan,ispin)%bands_to = b_to
               !
               DO j=1, ning
                  CALL read_line( input_line, end_of_file = tend )
                  !
                  IF ( tend ) &
                     CALL errore( 'read_cards', &
                                 'not enough wavefunctions', j )
                  IF (ning==1) THEN
                     READ(input_line,*) lo,im
                     c = 1.d0
                  ELSE
                     READ(input_line,*) lo,im,c
                  ENDIF

                  SELECT CASE(trim(lo))
                  CASE('s')
                     il = 0
                  CASE('p')
                     il = 1
                  CASE('d')
                     il = 2
                  CASE('f')
                     il = 3
                  CASE DEFAULT
                     CALL errore( 'read_cards', &
                                 'wrong l-label', 1 )
                  END SELECT

                  wan_data(iwan,ispin)%ing(j)%l = il
                  wan_data(iwan,ispin)%ing(j)%m = im
                  wan_data(iwan,ispin)%ing(j)%c = c
               ENDDO
            ENDDO
         ELSE
         ! oups - that is not our data - let's move one line up in input file
         ! not sure that a direct access to the parce_unit is safe enougth
            IF (ionode) BACKSPACE(parse_unit)
         ENDIF
      ELSE
         ! ok, that's the end of file. But I will move one line up
         ! for a correct handling of EOF in the parent read_cards subroutine
         ! otherwise (at least with gfortran on Mac) there will be the read error
         IF (ionode) BACKSPACE(parse_unit)
      ENDIF
      !
      RETURN
      !
   END SUBROUTINE card_wannier_ac
   !
   !
   !------------------------------------------------------------------------
   !    BEGIN manual
   !----------------------------------------------------------------------
   !
   ! SOLVENTS
   !
   !   set the solvents been read and their molecular files
   !
   ! Syntax:
   !
   !   SOLVENTS (units_option)
   !      label(1)    density(1)    molfile(1)
   !       ...        ...           ...
   !      label(n)    density(n)    molfile(n)
   !
   ! Example:
   !
   ! SOLVENTS (mol/L)
   !   H2O  55.3   H2O.spc.MOL
   !   Na+   0.1   Na.aq.MOL
   !   Cl-   0.1   Cl.aq.MOL
   !
   ! Where:
   !
   !   units_option == 1/cell  densities are given in particle's numbers in a cell
   !   units_option == mol/L   densities are given in mol/L
   !   units_option == g/cm^3  densities are given in g/cm^3
   !
   !   label(i)   ( character(len=10) )  label of the solvent
   !   density(i) ( real )               solvent's density
   !   molfile(i) ( character(len=80) )  file name of the pseudopotential
   !
   !----------------------------------------------------------------------
   !    END manual
   !------------------------------------------------------------------------
   !
   SUBROUTINE card_solvents( input_line )
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: input_line
      LOGICAL, EXTERNAL  :: matches
      INTEGER            :: iv, ip, ierr
      CHARACTER(len=10)  :: lb_mol
      CHARACTER(len=256) :: molfile
      !
      !
      IF ( tsolvents ) THEN
         CALL errore( ' card_solvents ', 'two occurrences', 2 )
      ENDIF
      IF ( nsolv > nsolx ) THEN
         CALL errore( ' card_solvents ', 'nsolv out of range', nsolv )
      ENDIF
      !
      IF ( matches( "1/CELL", input_line ) ) THEN
         solvents_unit = '1/cell'
      ELSEIF ( matches( "MOL/L", input_line ) ) THEN
         solvents_unit = 'mol/L'
      ELSEIF ( matches( "G/CM^3", input_line ) ) THEN
         solvents_unit = 'g/cm^3'
      ELSE
         IF ( trim( adjustl( input_line ) ) /= 'SOLVENTS' ) THEN
            CALL errore( 'read_cards ', &
                        & 'unknown option for SOLVENTS: '&
                        & // input_line, 1 )
         ENDIF
         CALL infomsg( 'read_cards ', &
            & 'DEPRECATED: no units specified in SOLVENTS card' )
         solvents_unit = '1/cell'
         CALL infomsg( 'read_cards ', &
            & 'SOLVENTS: units set to '//trim(solvents_unit) )
      ENDIF
      !
      DO iv = 1, nsolv
         !
         CALL read_line( input_line )
         IF (.NOT. laue_both_hands) THEN
           READ( input_line, *, iostat=ierr ) lb_mol, solv_dens1(iv), molfile
         ELSE
           READ( input_line, *, iostat=ierr ) lb_mol, solv_dens2(iv), solv_dens1(iv), molfile
         END IF
         CALL errore( ' card_solvents ', &
            & 'cannot read solvents from: '//trim(input_line), abs(ierr))
         solv_mfile(iv) = trim( molfile )
         lb_mol         = adjustl( lb_mol )
         solv_label(iv) = trim( lb_mol )
         !
         DO ip = 1, iv - 1
            IF ( solv_label(ip) == solv_label(iv) ) THEN
               CALL errore( ' card_solvents ', &
                           & " two occurrences of the same solvent's label ", iv )
            ENDIF
         ENDDO
         !
      ENDDO
      tsolvents = .true.
      !
      RETURN
      !
   END SUBROUTINE card_solvents
   !
   !
   !------------------------------------------------------------------------
   !    BEGIN manual
   !----------------------------------------------------------------------
   !
   ! TOTAL_CHARGE
   !
   !   set the total charge
   !
   ! Syntax:
   !
   !   TOTAL_CHARGE
   !      tot_charge
   !
   ! Example:
   !
   ! TOTAL_CHARGE
   !   0.1
   !
   !----------------------------------------------------------------------
   !    END manual
   !------------------------------------------------------------------------
   !
   SUBROUTINE card_total_charge( input_line )
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: input_line
      LOGICAL, EXTERNAL  :: matches
      INTEGER            :: iv, ip, ierr
      CHARACTER(len=10)  :: lb_mol
      CHARACTER(len=256) :: molfile
      !
      !
      IF ( ttotcharge ) THEN
         CALL errore( ' card_total_charge ', 'two occurrences', 2 )
      ENDIF
      !
      CALL read_line( input_line )
      READ( input_line, *, iostat=ierr ) tot_charge
      !
      CALL errore( ' card_total_charge ', &
         & 'cannot read total_charge from: '//trim(input_line), abs(ierr))
      !
      ttotcharge = .true.
      !
      RETURN
      !
   END SUBROUTINE card_total_charge
   !
   !------------------------------------------------------------------------
   !    BEGIN manual
   !----------------------------------------------------------------------
   !
   !    Needed to read Hubbard parameters
   !
   !  Syntax:
   !
   !    HUBBARD (projectors_option)
   !    hu_param hu_type-hu_manifold {case-specific-variables} hu_value 
   !
   !  Example:
   !
   !    HUBBARD (ortho-atomic) 
   !      1
   !    U  Fe-3d 5.0
   !    J0 Fe-3d 1.0
   !    J  Fe-3d 1 1.0
   !    V  Fe-3d O-2p 1 2 0.8
   !
   !  Where:
   !
   !    hu_param                          Hubbard parameter (U, J0, J, V, ...)
   !    hu_type                           Hubbard atom type 
   !    hu_manifold                       Hubbard manifold (2p, 3d, 4f, ...)
   !    hu_value                          value of the Hubbard parameter  
   !    projectors_option = ortho-atomic  ortho-atomic Hubbard projectors
   !    projectors_option = norm-atomic   norm-atomic Hubbard projectors
   !    projectors_option = atomic        atomic Hubbard projectors
   !    projectors_option = wf            Hubbard projectors read from file produced by pmw.x
   !    projectors_option = pseudo        Hubbard projectors are beta projectors from PP
   ! 
   !    'Fe' is the atomic type, '3d' is the Hubbard manifold, and
   !    '5.0' is the value of the Hubbard U parameter
   !
   !----------------------------------------------------------------------
   !    END manual
   !------------------------------------------------------------------------
   !
   SUBROUTINE card_hubbard ( input_line )
      !
      USE parameters,  ONLY : natx, sc_size
      USE constants,   ONLY : eps16
      USE upf_utils,   ONLY : spdf_to_l
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=256), INTENT(INOUT) :: input_line
      !
      CHARACTER(len=256) :: aux
      LOGICAL, EXTERNAL  :: imatches
      LOGICAL            :: tend, terr
      !
      ! Internal variables
      INTEGER :: i, nt, hu_nt, hu_nt2, nfield, na, nb, nc,&
                 nx, ny, nz, ldim, neigvals, eigval_index,&
                 io_stat, field, j, is, m
      REAL(DP):: hu_um_temp, hu_alpha_m_temp
      LOGICAL :: is_u, is_um, is_j0, is_v, is_j, is_b, is_e2, is_e3, is_alpha, is_alpha_m
      CHARACTER(LEN=20)  :: hu_param, field_str, hu_val, &
                            hu_at, hu_wfc, hu_at2, hu_wfc2, string, str, &
                            temp, hu_wfc_, hu_wfc2_, eigval_str
      INTEGER, ALLOCATABLE :: counter_u(:), counter_j0(:), counter_j(:), counter_b(:), counter_alpha(:), &
                              counter_e2(:), counter_e3(:), counter_v(:,:), ityp(:), target_indices(:)
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      !
      ! Output variables
      REAL(DP) :: hu_u,  &   ! Hubbard U (on-site)  
                  hu_j0, &   ! Hund's J0 (on-site)
                  hu_j,  &   ! Hund's J  (on-site)
                  hu_b,  &   ! Hund's B  (on-site) - only for d shell
                  hu_e2, &   ! Hund's E2 (on-site) - only for f shell
                  hu_e3, &   ! Hund's E3 (on-site) - only for f shell
                  hu_v,  &   ! Hubbard V (inter-site)
                  hu_alpha   ! Hubbard alpha (on-site)
      REAL(DP), ALLOCATABLE :: hu_um(:,:), &    ! Hubbard U (orbital-resolved, on-site)
                               hu_alpha_m(:,:)  ! Hubbard alpha (orbital-resolved, on-site)
      !
      INTEGER  :: hu_l,   &  ! orbital quantum number
                  hu_n,   &  ! principal quantum number
                  hu_l2,  &  ! orbital quantum number
                  hu_n2,  &  ! principal quantum number
                  hu_l_,  &  ! orbital quantum number
                  hu_n_,  &  ! principal quantum number
                  hu_l2_, &  ! orbital quantum number
                  hu_n2_     ! principal quantum number
      !
      IF ( tahub ) THEN
         CALL errore( 'card_hubbard', 'two occurrences', 2 )
      ENDIF
      IF ( .not. taspc ) THEN
         CALL errore( 'card_hubbard', 'ATOMIC_SPECIES must be present before HUBBARD', 2 )
      ENDIF
      !
      ! Read the type of Hubbard projectors
      IF ( imatches( "ORTHO-ATOMIC", input_line ) ) THEN
         Hubbard_projectors = 'ortho-atomic'
      ELSEIF ( imatches( "NORM-ATOMIC", input_line ) ) THEN
         Hubbard_projectors = 'norm-atomic'
      ELSEIF ( imatches( "-ATOMIC", input_line ) ) THEN
         ! Sanity check
         ! This is the case when the first part of the name was misspelled 
         CALL errore( 'card_hubbard', 'Wrong name of the Hubbard projectors',1)
      ELSEIF ( imatches( "ORTHOATOMIC", input_line ) .OR. &
               imatches( "NORMATOMIC", input_line ) ) THEN
         ! Sanity check
         ! This is the case when the dash was forgotten in the name
         CALL errore( 'card_hubbard', 'Wrong name of the Hubbard projectors',1)
      ELSEIF ( imatches( "ATOMIC", input_line ) ) THEN
         Hubbard_projectors = 'atomic'
      ELSEIF ( imatches( "WF", input_line ) ) THEN 
         Hubbard_projectors = 'wf'
      ELSEIF ( imatches( "PSEUDO", input_line ) ) THEN
         Hubbard_projectors = 'pseudo'
      ELSE
         IF ( trim( adjustl( input_line ) ) /= 'HUBBARD' ) THEN
            CALL errore( 'card_hubbard', &
                        & 'unknown option for HUBBARD: '&
                        & // input_line, 1 )
         ELSE
           CALL errore( 'card_hubbard', &
                        & 'None or wrong Hubbard projectors specified in the HUBBARD card: ',1)
         ENDIF
      ENDIF
      !
      Hubbard_projectors = TRIM(ADJUSTL(Hubbard_projectors))
      !
      ALLOCATE(counter_u(ntyp))
      counter_u(:) = 0
      ALLOCATE(counter_j0(ntyp))
      counter_j0(:) = 0
      ALLOCATE(counter_j(ntyp))
      counter_j(:) = 0
      ALLOCATE(counter_b(ntyp))
      counter_b(:) = 0
      ALLOCATE(counter_e2(ntyp))
      counter_e2(:) = 0
      ALLOCATE(counter_e3(ntyp))
      counter_e3(:) = 0
      ALLOCATE(counter_alpha(ntyp))
      counter_alpha(:) = 0
      !
      ! Read Hubbard parameters, principal and orbital quantum numbers
      !
      i = 0
      DO WHILE (.TRUE.)
         !
         ! We can exit this loop in two cases:
         ! 1) HUBBARD card is the last in the input and we reached the end of the input file, or
         ! 2) HUBBARD card is NOT the last in the input and we reached the next card. 
         !
         i = i + 1
         !
         ! Initialize different parameters
         hu_l=-1;  hu_n=-1;  hu_l2=-1;  hu_n2=-1 
         hu_l_=-1; hu_n_=-1; hu_l2_=-1; hu_n2_=-1
         hu_nt=-1; hu_nt2=-1
         hu_u=0.0; hu_j0=0.0; hu_j=0.0; hu_b=0.0
         hu_e2=0.0; hu_e3=0.0; hu_v=0.0
         hu_alpha_m_temp=0.0; hu_um_temp=0.0
         hu_wfc=''; hu_wfc_=''; hu_wfc2=''; hu_wfc2_=''
         !
         ! Read the i-th input line
         CALL read_line( input_line, end_of_file = tend, error = terr )
         IF ( tend ) THEN
            IF (i==1) THEN
               ! In this case nothing was read so far, so stopping
               CALL errore ('card_hubbard', ' Nothing was read from the HUBBARD card. Stopping... ' &
                          & // TRIM(int_to_char(i)), i)
            ELSE
               ! All lines were read successfully and we reached the end of the HUBBARD card. Exit smoothly.
               ! However, before exiting, we move one line up in the input file for a correct handling of 
               ! EOF in the parent read_cards subroutine, because otherwise (with gfortran) there will be 
               ! the read error
               IF (ionode) BACKSPACE (parse_unit)
               GO TO  11
            ENDIF
         ENDIF
         IF ( terr ) CALL errore ('card_hubbard', ' Error while reading the HUBBARD card on line ' &
                          & // TRIM(int_to_char(i)), i)
         !
         ! Determine how many columns in the i-th row
         CALL field_count( nfield, input_line )
         !
         ! Column 1: Read the Hubbard parameter name (e.g. U, J0, J, V)
         CALL get_field(1, hu_param, input_line)
         hu_param = TRIM(hu_param)
         IF ( LEN_TRIM(hu_param) > 5 ) THEN
            ! This is the case when most likely we reached the end of the HUBBARD card 
            ! and started reading the next card in the input. So we need to exit smoothly.
            ! This case will not happen if the HUBBARD card is the last in the input file.
            ! Let's move one line up in the input file
            IF(ionode) BACKSPACE (parse_unit)
            GO TO 11
         ENDIF
         ! Check whether the length of the Hubbard parameter is within the allowed ranges
         IF ( LEN_TRIM(hu_param) < 1 .or. LEN_TRIM(hu_param) > 5) &
           CALL errore( 'card_hubbard', &
                      'Hubbard parameter name missing or too long', i )
         !
         is_u  =     ( hu_param == 'u' .OR. hu_param == 'U' )
         is_um =     ( is_u            .AND. nfield   >  3  ) ! error handling below
         is_j0 =     ( hu_param == 'j0'.OR. hu_param == 'J0')
         is_j  =     ( hu_param == 'j' .OR. hu_param == 'J' )
         is_b  =     ( hu_param == 'b' .OR. hu_param == 'B' ) ! for d shell
         is_e2 =     ( hu_param == 'e2'.OR. hu_param == 'E2') ! for f shell
         is_e3 =     ( hu_param == 'e3'.OR. hu_param == 'E3') ! for f shell
         is_v  =     ( hu_param == 'v' .OR. hu_param == 'V' )
         is_alpha =  ( hu_param == 'alpha'.OR. &
                hu_param == 'ALPHA' .OR. hu_param == 'Alpha')
         is_alpha_m= ( is_alpha        .AND. nfield   >  3  )
         !
         IF (is_um .OR. is_alpha_m) THEN
            is_u  = .FALSE.    ! orbital-resolved U switches off standard U (more checks below)
            is_alpha = .FALSE. ! same for orbital-resolved Hubbard_alpha vs. standard alpha
         ENDIF
         !
         IF (.NOT.is_u .AND. .NOT.is_um .AND. .NOT.is_j0 .AND. .NOT.is_j .AND. &
             .NOT.is_b .AND. .NOT.is_e2 .AND. .NOT.is_e3 .AND. .NOT.is_v .AND. &
             .NOT.is_alpha .AND. .NOT.is_alpha_m) THEN
            WRITE(stdout,'(/5x,"Problem in the HUBBARD card on line ",i5)') i
            CALL errore('card_hubbard', 'Unknown label of the Hubbard parameter', i)
         ENDIF
         !
         IF ( is_u .OR. is_um .OR. is_j0 .OR. is_j .OR. is_b .OR. is_e2 .OR. &
              is_e3 .OR. is_alpha .OR. is_alpha_m) THEN
            !     
            ! Column 2: Read the atomic type name and the Hubbard manifold (e.g. Fe-3d)
            CALL get_field(2, field_str, input_line)
            field_str = TRIM(field_str)
            !
            ! Read the Hubbard atom name (e.g. Fe)
            hu_at = between( field_str, '', '-' )
            IF ( LEN_TRIM(ADJUSTL(hu_at)) < 1 .OR. LEN_TRIM(ADJUSTL(hu_at)) > 4 ) &
               CALL errore( 'card_hubbard', &
                      'Hubbard atom name missing or wrong or too long', i )
            !
            ! Determine the index of the atomic type of the current Hubbard atom
            DO nt = 1, ntyp
               IF (TRIM(ADJUSTL(hu_at)) == atom_label(nt)) THEN
                  hu_nt = nt
                  GO TO 16
               ENDIF
            ENDDO
            IF (hu_nt == -1) CALL errore( 'card_hubbard', &
                   'Name of the Hubbard atom does not match with any type in ATOMIC_SPECIES', i )
16          CONTINUE
            !
            ! Setup the counter to monitor how many Hubbard manifolds per atomic type do we have
            IF (is_u) THEN
               counter_u(hu_nt) = counter_u(hu_nt) + 1
               IF (counter_u(hu_nt)>3) CALL errore( 'card_hubbard', &
                  'Too many entries for U for the same atomic type', i )
            ELSEIF (is_j0) THEN
               counter_j0(hu_nt) = counter_j0(hu_nt) + 1
               IF (counter_j0(hu_nt)>3) CALL errore( 'card_hubbard', &
                  'Too many entries for J0 for the same atomic type', i )
            ELSEIF (is_j) THEN
               counter_j(hu_nt) = counter_j(hu_nt) + 1
               IF (counter_j(hu_nt) > 1) CALL errore( 'card_hubbard', &
                       'More than 1 entry for J of the same atomic type is not allowed', i )
            ELSEIF (is_b) THEN
               counter_b(hu_nt) = counter_b(hu_nt) + 1
               IF (counter_b(hu_nt) > 1) CALL errore( 'card_hubbard', &
                       'More than 1 entry for B of the same atomic type is not allowed', i )
            ELSEIF (is_e2) THEN
               counter_e2(hu_nt) = counter_e2(hu_nt) + 1
               IF (counter_e2(hu_nt) > 1) CALL errore( 'card_hubbard', &
                       'More than 1 entry for E2 of the same atomic type is not allowed', i )
            ELSEIF (is_e3) THEN
               counter_e3(hu_nt) = counter_e3(hu_nt) + 1
               IF (counter_e3(hu_nt) > 1) CALL errore( 'card_hubbard', &
                       'More than 1 entry for E3 of the same atomic type is not allowed', i )
            ELSEIF (is_alpha) THEN
               counter_alpha(hu_nt) = counter_alpha(hu_nt) + 1
               IF (counter_alpha(hu_nt) > 1) CALL errore( 'card_hubbard', &
                       'More than 1 entry for ALPHA of the same atomic type is not allowed', i )
            ELSEIF (is_alpha_m) THEN
               counter_alpha(hu_nt) = counter_alpha(hu_nt) + 1
               IF (counter_alpha(hu_nt) > 1) CALL errore( 'card_hubbard', &
                       'More than 1 entry for orbital-resolved ALPHA of the same atomic type is not allowed', i )
            ENDIF
            !
            ! Read the Hubbard manifold(s)
            ! Note: There may be two manifolds at the same time, though this is not 
            ! allowed for the first (main) Hubbard channel.
            IF ( (is_u  .AND. counter_u(hu_nt)==1)  .OR. &
                 (is_j0 .AND. counter_j0(hu_nt)==1) .OR. &
                  is_j .OR. is_um .OR. is_b .OR. is_e2 .OR. &
                  is_e3 .OR. is_alpha .OR. is_alpha_m ) THEN
               ! e.g. Fe-3d
               hu_wfc = between( field_str, '-', '' )
            ELSEIF ((is_u  .AND. counter_u(hu_nt)==2) .OR. &
                    (is_j0 .AND. counter_j0(hu_nt)==2) ) THEN
               ! e.g. Fe-3p or Fe-3p-3s
               temp = between( field_str, '-', '' ) 
               hu_wfc = between( temp, '', '-' )
               IF (hu_wfc/='') THEN     
                  ! two manifolds found
                  hu_wfc_ = between( temp, '-', '' )
               ELSE
                  ! one manifold found
                  hu_wfc = temp
               ENDIF
            ENDIF
            IF ( LEN_TRIM(hu_wfc) /= 2 ) THEN
               WRITE(stdout,'(/5x,"Hubbard manifold ",a6," in the HUBBARD card on line ",i5)') TRIM(hu_wfc), i
               CALL errore( 'card_hubbard', &
                     'Hubbard manifold name missing or wrong or too long or in wrong order', i )
            ENDIF 
            IF ( hu_wfc_/='' .AND. LEN_TRIM(hu_wfc_) /= 2 ) THEN
               WRITE(stdout,'(/5x,"Hubbard manifold ",a6," in the HUBBARD card on line ",i5)') TRIM(hu_wfc_), i
               CALL errore( 'card_hubbard', &
                     'Hubbard manifold name missing or wrong or too long', i )
            ENDIF
            !
            ! Determine the principal and orbital quantum numbers of the Hubbard manifold
            READ (hu_wfc(1:1),'(i1)', END=14, ERR=15) hu_n
            hu_l = spdf_to_l( hu_wfc(2:2) )
            IF ( hu_n == -1 ) CALL errore( 'card_hubbard', 'Hubbard n is wrong', i)
            IF ( hu_l == -1 ) CALL errore( 'card_hubbard', 'Hubbard l is wrong', i)
            !
            ! Determine the principal and orbital quantum numbers of the Hubbard manifold
            IF (hu_wfc_/='') THEN
               READ (hu_wfc_(1:1),'(i1)', END=14, ERR=15) hu_n_
               hu_l_ = spdf_to_l( hu_wfc_(2:2) )
               IF ( hu_n_ == -1 ) CALL errore( 'card_hubbard', 'Hubbard n is wrong', i)
               IF ( hu_l_ == -1 ) CALL errore( 'card_hubbard', 'Hubbard l is wrong', i)
            ENDIF
            !
            ldim =  2*hu_l+1 ! Compute number of magnetic quantum orbitals
            !
            ! Set up array to store orbital-resolved DFT+U parameters
            IF (is_um .OR. is_alpha_m) THEN
               IF (noncolin) THEN
                  ! only need one dimension for noncollinear case
                  ! retain the extra dimension of spin for compatibility
                  ALLOCATE(hu_um(2*ldim,1))
                  ALLOCATE(hu_alpha_m(2*ldim,1)) 
               ELSE
                  ALLOCATE(hu_um(ldim,nspin))
                  ALLOCATE(hu_alpha_m(ldim,nspin))
               ENDIF
               hu_alpha_m(:,:) = 0.0
               hu_um(:,:) = 0.0
            ENDIF
            !
            ! Sanity check
            IF (is_b .AND. hu_l/=2) THEN
               ! allowed only for d electrons
               WRITE(stdout,'(/5x,"Problem in the HUBBARD card on line ",i5)') i
               CALL errore( 'card_hubbard', &
                       & 'B parameter can be specified only for d electrons (see documentation)', i )
            ELSEIF (is_e2 .AND. hu_l/=3) THEN
               ! allowed only for f electrons
               WRITE(stdout,'(/5x,"Problem in the HUBBARD card on line ",i5)') i
               CALL errore( 'card_hubbard', &
                       & 'E2 parameter can be specified only for f electrons (see documentation)', i )
            ELSEIF (is_e3 .AND. hu_l/=3) THEN
               ! allowed only for f electrons
               WRITE(stdout,'(/5x,"Problem in the HUBBARD card on line ",i5)') i
               CALL errore( 'card_hubbard', &
                       & 'E3 parameter can be specified only for f electrons (see documentation)', i )
            ENDIF
            !
            ! Assign the principal and orbital quantum numbers
            IF ( (is_u  .AND. counter_u(hu_nt)==1)  .OR. &
                 (is_j0 .AND. counter_j0(hu_nt)==1) .OR. &
                  is_j .OR. is_um .OR. is_b .OR. is_e2 .OR. &
                  is_e3 .OR. is_alpha .OR. is_alpha_m ) THEN
               ! First Hubbard manifold
               IF (Hubbard_n(hu_nt)<0 .AND. Hubbard_l(hu_nt)<0) THEN
                  ! initialization
                  Hubbard_n(hu_nt) = hu_n
                  Hubbard_l(hu_nt) = hu_l
               ELSE
                  ! sanity checks
                  IF (hu_n/=Hubbard_n(hu_nt) .OR. hu_l/=Hubbard_l(hu_nt)) THEN
                     WRITE(stdout,'(/5x,"Problem in the HUBBARD card on line ",i5)') i
                     IF (is_j0) CALL errore( 'card_hubbard', &
                          & 'Mismatch in the quantum numbers for U and J0 for the same atomic type', i )
                     IF (is_j)  WRITE(stdout,'(/5x,"Only one manifold is allowed for Hund J")')
                     IF (is_um) CALL errore( 'card_hubbard', 'Currently, orbital-resolved Hubbard U params &
                           &can only be assigned to the orbitals of one shell per atom', i)
                     IF (is_b)  WRITE(stdout,'(/5x,"Only one manifold is allowed for Hund B")')
                     IF (is_e2) WRITE(stdout,'(/5x,"Only one manifold is allowed for Hund E2")')
                     IF (is_e3) WRITE(stdout,'(/5x,"Only one manifold is allowed for Hund E3")')
                     IF (is_alpha) WRITE(stdout,'(/5x,"Only one manifold is allowed for Hubbard ALPHA")')
                     IF (is_alpha_m) CALL errore( 'card_hubbard', 'Orbital-resolved Hubbard ALPHA params &
                           &can only be assigned to the orbitals of one shell per atom', i)
                     CALL errore( 'card_hubbard', &
                          & 'Mismatch in the quantum numbers for the same atomic type', i )
                  ENDIF
               ENDIF
            ELSEIF ((is_u  .AND. counter_u(hu_nt)==2) .OR. &
                    (is_j0 .AND. counter_j0(hu_nt)==2)) THEN
               ! Second Hubbard manifold
               ! Check whether we have different Hubbard manifolds for the same atomic type
               IF ( hu_n==Hubbard_n(hu_nt) .AND. hu_l==Hubbard_l(hu_nt) ) THEN
                  WRITE(stdout,'(/5x,"Problem in the HUBBARD card on line ",i5)') i
                  WRITE(stdout,'(5x,"Two Hubbard channels are the same for ",a)') atom_label(hu_nt)
                  CALL errore( 'card_hubbard', &
                      'Not allowed to specify two Hubbard channels that are the same for the same atom', i )
               ENDIF
               IF (Hubbard_n2(hu_nt)<0 .AND. Hubbard_l2(hu_nt)<0) THEN
                  ! initialization
                  Hubbard_n2(hu_nt) = hu_n
                  Hubbard_l2(hu_nt) = hu_l
               ELSE
                  ! sanity check (needed for DFT+U+V)
                  IF (hu_n/=Hubbard_n2(hu_nt) .OR. hu_l/=Hubbard_l2(hu_nt)) THEN
                     WRITE(stdout,'(/5x,"Problem in the HUBBARD card for U on line ",i5)') i
                     IF (is_j0) CALL errore( 'card_hubbard', &
                          & 'Mismatch in the quantum numbers for U and J0 for the same atomic type (2nd channel)', i )
                     CALL errore( 'card_hubbard', &
                          & 'Mismatch in the quantum numbers for the same atomic type', i )
                  ENDIF
               ENDIF
               IF (hu_n_>-1 .AND. hu_l_>-1) THEN
                  ! Third Hubbard manifold
                  ! Check whether we have different Hubbard manifolds for the same atomic type
                  IF ( hu_n_==Hubbard_n(hu_nt)  .AND. hu_l_==Hubbard_l(hu_nt) .OR. &
                       hu_n_==Hubbard_n2(hu_nt) .AND. hu_l_==Hubbard_l2(hu_nt)) THEN
                     WRITE(stdout,'(/5x,"Problem in the HUBBARD card on line ",i5)') i
                     WRITE(stdout,'(5x,"Two Hubbard channels are the same for ",a)') atom_label(hu_nt)
                     CALL errore( 'card_hubbard', &
                         'Not allowed to specify two Hubbard channels that are the same for the same atom', i )
                  ENDIF
                  backall(hu_nt) = .TRUE.
                  IF (Hubbard_n3(hu_nt)<0 .AND. Hubbard_l3(hu_nt)<0) THEN
                     ! initialization
                     Hubbard_n3(hu_nt) = hu_n_
                     Hubbard_l3(hu_nt) = hu_l_
                  ELSE
                     ! sanity check (needed for DFT+U+V)
                     IF (hu_n_/=Hubbard_n3(hu_nt) .OR. hu_l_/=Hubbard_l3(hu_nt)) THEN
                        WRITE(stdout,'(/5x,"Problem in the HUBBARD card for U on line ",i5)') i
                        IF (is_j0) CALL errore( 'card_hubbard', &
                          & 'Mismatch in the quantum numbers for U and J0 for the same atomic type (3rd channel)', i )
                        CALL errore( 'card_hubbard', &
                             & 'Mismatch in the quantum numbers for the same atomic type', i )
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            !
         ELSEIF ( is_v ) THEN
            !
            ! Here is the case of V
            !
            ! Sanity check
            IF (nat>natx) CALL errore('card_hubbard', 'Too many atoms. &
                Increase the value of natx in Modules/parameters.f90 and recompile the code.',1)
            !
            ! Initialize the atomic types for 
            ! the virtual atoms in the same way as it is done in
            ! PW/src/intersite_V.f90
            ! sp_pos(na) is the atomic type of the atom na
            IF (.NOT.ALLOCATED(ityp)) THEN
               IF (.NOT.ALLOCATED(sp_pos)) CALL errore ('card_hubbard', &
                       'card HUBBARD must follow card ATOMIC_SPECIES',1)
               ALLOCATE(ityp(natx*(2*sc_size+1)**3))
               ityp(1:nat) = sp_pos(1:nat)
               i = nat
               DO nx = -sc_size, sc_size
                  DO ny = -sc_size, sc_size
                     DO nz = -sc_size, sc_size
                        IF ( nx.NE.0 .OR. ny.NE.0 .OR. nz.NE.0 ) THEN
                           DO na = 1, nat
                              i = i + 1
                              ityp(i) = sp_pos(na)
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF
            !
            IF (.NOT.ALLOCATED(counter_v)) THEN
               ALLOCATE(counter_v(natx,natx*(2*sc_size+1)**3))
               counter_v(:,:) = 0
            ENDIF
            !
            ! First of all, we read the indices na and nb that correspond to the location 
            ! of atoms Fe and O in the ATOMIC_POSITIONS card (columns 4 and 5)
            CALL get_field(4, field_str, input_line)
            READ(field_str,'(i8)', END=14, ERR=15) na
            IF ( na < 0 .or. na > nat ) &
               CALL errore( 'card_hubbard', 'Not allowed value of the atomic index na', i)
            CALL get_field(5, field_str, input_line)
            READ(field_str,'(i8)', END=14, ERR=15) nb
            IF ( nb < 0 .or. nb > nat*natx ) &
               CALL errore( 'card_hubbard', 'Not allowed value of the atomic index nb', i)
            !
            ! In the DFT+U+V case there are maximum 4 Hubbard_V parameters per couple (na,nb)
            ! that cover interactions: standard-standard (nc=1), standard-background (nc=2),
            ! background-background (nc=3), and background-standard (nc=4).
            !
            IF (na==nb) THEN
               nt = ityp(na)
               IF (counter_u(nt)>0) THEN
                  IF (counter_u(nt)==1 .AND. (counter_v(na,na)==0)) THEN
                     ! In this case, the on-site term was already read using the
                     ! syntax (e.g. U Fe-3d instead of V Fe-3d Fe-3d). Therefore, 
                     ! now we are reading the next occurrence (cross terms between
                     ! the first and the second Hubbard manifold of Fe) using 
                     ! the V syntax.
                     counter_v(na,na) = counter_u(nt) + 1
                  ELSEIF (counter_u(nt)>1) THEN
                     CALL errore( 'card_hubbard', 'Use V instead of U to specify Hubbard &
                             manifolds other than the first in the DFT+U+V scheme', i)
                  ELSE
                     counter_v(na,nb) = counter_v(na,nb) + 1
                  ENDIF  
               ELSE
                  counter_v(na,nb) = counter_v(na,nb) + 1
               ENDIF
            ELSE
               counter_v(na,nb) = counter_v(na,nb) + 1
            ENDIF
            !
            ! Setup the value of "nc"
            IF (counter_v(na,nb)==1) THEN
               nc = 1
            ELSEIF (counter_v(na,nb)==2) THEN
               nc = 2
            ELSEIF (counter_v(na,nb)==3) THEN
               nc = 3
            ELSEIF (counter_v(na,nb)==4) THEN
               nc = 4
            ELSE
               WRITE(stdout,'(/5x,"Problem in the HUBBARD card for V on line ",i5)') i
               CALL errore( 'card_hubbard', 'Too many occurrences of V for the same couple of atoms', i)
            ENDIF
            !
            !**********************************************************************************!
            !*                   Read the data for the first atom                             *!
            !**********************************************************************************!
            !
            ! Column 3: Read the atomic type name and the Hubbard manifold (e.g. Fe-3d)
            CALL get_field(2, field_str, input_line)
            field_str = TRIM(field_str)
            !
            ! Read the Hubbard atom name (e.g. Fe)
            hu_at = between( field_str, '', '-' )
            IF ( LEN_TRIM(ADJUSTL(hu_at)) < 1 .OR. LEN_TRIM(ADJUSTL(hu_at)) > 4 ) &
               CALL errore( 'card_hubbard', &
                      'Hubbard V: 1st atom name missing or wrong or too long', i )
            !
            ! Determine the index of the atomic type of the first Hubbard atom
            DO nt = 1, ntyp
               IF (TRIM(ADJUSTL(hu_at)) == atom_label(nt)) THEN
                  hu_nt = nt
                  GO TO 12
               ENDIF
            ENDDO
            IF (hu_nt == -1) CALL errore( 'card_hubbard', &
                 'Name of the Hubbard atom does not match with any type in ATOMIC_SPECIES', i )
12          CONTINUE
            !
            ! Read the Hubbard manifold(s)
            ! Note: There may be two manifolds at the same time, though this is
            ! allowed only for counter_v(na,nb)=3 and counter_v(na,nb)=4 for the first atom.
            IF (counter_v(na,nb)==1 .OR. counter_v(na,nb)==2) THEN
               ! e.g. Fe-3d
               hu_wfc = between( field_str, '-', '' )
               IF ( LEN_TRIM(hu_wfc) /= 2 ) &
                  CALL errore( 'card_hubbard', &
                        'Hubbard V: manifold of the 1st atom missing or wrong or too long', i )
            ELSEIF (counter_v(na,nb)==3 .OR. counter_v(na,nb)==4) THEN
               ! e.g. Fe-3p or Fe-3p-3s
               temp = between( field_str, '-', '' )
               hu_wfc = between( temp, '', '-' )
               IF (hu_wfc/='') THEN
                  ! two manifolds found
                  hu_wfc_ = between( temp, '-', '' )
               ELSE
                  ! one manifold found
                  hu_wfc = temp
                  IF (Hubbard_n3(hu_nt)>-1 .OR. Hubbard_l3(hu_nt)>-1) &
                     CALL errore( 'card_hubbard', &
                     'Three Hubbard manifolds were previously found for this atomic type', i )
               ENDIF
            ENDIF
            IF ( LEN_TRIM(hu_wfc) /= 2 ) THEN
               WRITE(stdout,'(/5x,"Hubbard manifold ",a6," in the HUBBARD card on line ",i5)') TRIM(hu_wfc), i
               CALL errore( 'card_hubbard', &
                     'Hubbard manifold name missing or wrong or too long or in wrong order', i )
            ENDIF
            IF ( hu_wfc_/='' .AND. LEN_TRIM(hu_wfc_) /= 2 ) THEN
               WRITE(stdout,'(/5x,"Hubbard manifold ",a6," in the HUBBARD card on line ",i5)') TRIM(hu_wfc_), i
               CALL errore( 'card_hubbard', &
                     'Hubbard manifold name missing or wrong or too long', i )
            ENDIF 
            !
            ! Determine the principal and orbital quantum numbers of the Hubbard manifold
            READ(hu_wfc(1:1),'(i1)', END=14, ERR=15) hu_n
            hu_l = spdf_to_l( hu_wfc(2:2) )
            IF ( hu_n == -1 ) CALL errore( 'card_hubbard', 'Hubbard n is wrong (1st atom)', i)
            IF ( hu_l == -1 ) CALL errore( 'card_hubbard', 'Hubbard l is wrong (1st atom)', i)
            !
            ! Determine the principal and orbital quantum numbers of the Hubbard manifold
            IF (hu_wfc_/='') THEN
               READ (hu_wfc_(1:1),'(i1)', END=14, ERR=15) hu_n_
               hu_l_ = spdf_to_l( hu_wfc_(2:2) )
               IF ( hu_n_ == -1 ) CALL errore( 'card_hubbard', 'Hubbard n is wrong', i)
               IF ( hu_l_ == -1 ) CALL errore( 'card_hubbard', 'Hubbard l is wrong', i)
            ENDIF
            !
            ! Assign the principal and orbital quantum numbers
            IF (counter_v(na,nb)==1 .OR. counter_v(na,nb)==2) THEN
               ! First Hubbard manifold of the first atom
               IF (Hubbard_n(hu_nt)<0 .AND. Hubbard_l(hu_nt)<0) THEN
                  ! initialization
                  Hubbard_n(hu_nt) = hu_n
                  Hubbard_l(hu_nt) = hu_l
               ELSE
                  ! sanity check
                  IF (hu_n/=Hubbard_n(hu_nt) .OR. hu_l/=Hubbard_l(hu_nt)) THEN
                     WRITE(stdout,'(/5x,"Problem in the HUBBARD card on line ",i5)') i
                     CALL errore( 'card_hubbard', &
                          & 'Mismatch in the quantum numbers for the same atomic type', i )
                  ENDIF
               ENDIF
            ELSEIF (counter_v(na,nb)==3 .OR. counter_v(na,nb)==4) THEN
               ! Second Hubbard manifold of the first atom
               ! Check whether we have different Hubbard manifolds for the same atomic type
               IF ( hu_n==Hubbard_n(hu_nt) .AND. hu_l==Hubbard_l(hu_nt) ) THEN
                  WRITE(stdout,'(/5x,"Problem in the HUBBARD card on line ",i5)') i
                  WRITE(stdout,'(5x,"Two Hubbard channels are the same for ",a)') atom_label(hu_nt)
                  CALL errore( 'card_hubbard', &
                      'Not allowed to specify two Hubbard channels that are the same for the same atom', i )
               ENDIF
               IF (Hubbard_n2(hu_nt)<0 .AND. Hubbard_l2(hu_nt)<0) THEN
                  ! initialization
                  Hubbard_n2(hu_nt) = hu_n
                  Hubbard_l2(hu_nt) = hu_l
               ELSE
                  ! sanity check (needed for DFT+U+V)
                  IF (hu_n/=Hubbard_n2(hu_nt) .OR. hu_l/=Hubbard_l2(hu_nt)) THEN
                     WRITE(stdout,'(/5x,"Problem in the HUBBARD card for U on line ",i5)') i
                     CALL errore( 'card_hubbard', &
                          & 'Mismatch in the quantum numbers for the same atomic type', i )
                  ENDIF
               ENDIF
               IF (hu_n_>-1 .AND. hu_l_>-1) THEN
                  ! Third Hubbard manifold of the first atom
                  ! Check whether we have different Hubbard manifolds for the same atomic type
                  IF ( hu_n_==Hubbard_n(hu_nt)  .AND. hu_l_==Hubbard_l(hu_nt) .OR. &
                       hu_n_==Hubbard_n2(hu_nt) .AND. hu_l_==Hubbard_l2(hu_nt)) THEN
                     WRITE(stdout,'(/5x,"Problem in the HUBBARD card on line ",i5)') i
                     WRITE(stdout,'(5x,"Two Hubbard channels are the same for ",a)') atom_label(hu_nt)
                     CALL errore( 'card_hubbard', &
                         'Not allowed to specify two Hubbard channels that are the same for the same atom', i )
                  ENDIF
                  backall(hu_nt) = .TRUE.
                  IF (Hubbard_n3(hu_nt)<0 .AND. Hubbard_l3(hu_nt)<0) THEN
                     ! initialization
                     Hubbard_n3(hu_nt) = hu_n_
                     Hubbard_l3(hu_nt) = hu_l_
                  ELSE
                     ! sanity check (needed for DFT+U+V)
                     IF (hu_n_/=Hubbard_n3(hu_nt) .OR. hu_l_/=Hubbard_l3(hu_nt)) THEN
                        WRITE(stdout,'(/5x,"Problem in the HUBBARD card for U on line ",i5)') i
                        CALL errore( 'card_hubbard', &
                             & 'Mismatch in the quantum numbers for the same atomic type', i )
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            !
            !**********************************************************************************!
            !*                   Read the data for the second atom                            *!
            !**********************************************************************************!
            !
            ! Column 3: Read the atomic type name and the Hubbard manifold (e.g. O-2p)
            CALL get_field(3, field_str, input_line)
            field_str = TRIM(field_str)
            !
            ! Read the Hubbard atom name (e.g. O)
            hu_at2 = between( field_str, '', '-' )
            IF ( LEN_TRIM(ADJUSTL(hu_at2)) < 1 .OR. LEN_TRIM(ADJUSTL(hu_at2)) > 4 ) &
               CALL errore( 'card_hubbard', &
                      'Hubbard V: 2nd atom name missing or wrong or too long', i )
            !
            ! Determine the index of the atomic type of the second Hubbard atom
            ! (i.e. of the neighbor atom)
            DO nt = 1, ntyp
               IF (TRIM(ADJUSTL(hu_at2)) == atom_label(nt)) THEN
                  hu_nt2 = nt
                  GO TO 13
               ENDIF
            ENDDO
            IF (hu_nt2 == -1) CALL errore( 'card_hubbard', &
                 'Name of the Hubbard atom does not match with any type in ATOMIC_SPECIES', i )
13          CONTINUE
            !
            ! Read the Hubbard manifold(s)
            ! Note: There may be two manifolds at the same time, though this is
            ! allowed only for counter_v(na,nb)=2 and counter_v(na,nb)=3 for the second atom.
            IF (counter_v(na,nb)==1 .OR. counter_v(na,nb)==4) THEN
               ! e.g. O-2p
               hu_wfc2 = between( field_str, '-', '' )
               IF ( len_trim(hu_wfc2) /= 2 ) &
                  CALL errore( 'card_hubbard', &
                        'Hubbard V: manifold of the 2nd atom missing or wrong or too long', i )
            ELSEIF (counter_v(na,nb)==2 .OR. counter_v(na,nb)==3) THEN
               ! e.g. O-2s or O-2s-1s
               temp = between( field_str, '-', '' )
               hu_wfc2 = between( temp, '', '-' )
               IF (hu_wfc2/='') THEN
                  ! two manifolds found
                  hu_wfc2_ = between( temp, '-', '' )
               ELSE
                  ! one manifold found
                  hu_wfc2 = temp
                  IF (Hubbard_n3(hu_nt2)>-1 .OR. Hubbard_l3(hu_nt2)>-1) &
                     CALL errore( 'card_hubbard', &
                     'Three Hubbard manifolds were previously found for this atomic type', i )
               ENDIF
            ENDIF
            IF ( LEN_TRIM(hu_wfc2) /= 2 ) THEN
               WRITE(stdout,'(/5x,"Hubbard manifold ",a6," in the HUBBARD card on line ",i5)') TRIM(hu_wfc2), i
               CALL errore( 'card_hubbard', &
                     'Hubbard manifold name missing or wrong or too long or in wrong order', i )
            ENDIF
            IF ( hu_wfc2_/='' .AND. LEN_TRIM(hu_wfc2_) /= 2 ) THEN
               WRITE(stdout,'(/5x,"Hubbard manifold ",a6," in the HUBBARD card on line ",i5)') TRIM(hu_wfc2_), i
               CALL errore( 'card_hubbard', &
                     'Hubbard manifold name missing or wrong or too long', i )
            ENDIF
            !
            ! Determine the principal and orbital quantum numbers
            ! of the 2nd Hubbard manifold (i.e. of the neighbor atom)
            READ(hu_wfc2(1:1),'(i1)', END=14, ERR=15) hu_n2
            hu_l2 = spdf_to_l( hu_wfc2(2:2) )
            IF ( hu_n2 == -1 ) CALL errore( 'card_hubbard', 'Hubbard n is wrong (2nd atom)', i)
            IF ( hu_l2 == -1 ) CALL errore( 'card_hubbard', 'Hubbard l is wrong (2nd atom)', i)
            !
            ! Determine the principal and orbital quantum numbers of the Hubbard manifold
            IF (hu_wfc2_/='') THEN
               READ (hu_wfc2_(1:1),'(i1)', END=14, ERR=15) hu_n2_
               hu_l2_ = spdf_to_l( hu_wfc2_(2:2) )
               IF ( hu_n2_ == -1 ) CALL errore( 'card_hubbard', 'Hubbard n is wrong', i)
               IF ( hu_l2_ == -1 ) CALL errore( 'card_hubbard', 'Hubbard l is wrong', i)
            ENDIF
            !
            ! Assign the principal and orbital quantum numbers
            IF (counter_v(na,nb)==1 .OR. counter_v(na,nb)==4) THEN
               ! First Hubbard manifold of the second atom
               IF (Hubbard_n(hu_nt2)<0 .AND. Hubbard_l(hu_nt2)<0) THEN
                  ! initialization
                  Hubbard_n(hu_nt2) = hu_n2
                  Hubbard_l(hu_nt2) = hu_l2
               ELSE
                  ! sanity check
                  IF (hu_n2/=Hubbard_n(hu_nt2) .OR. hu_l2/=Hubbard_l(hu_nt2)) THEN
                     WRITE(stdout,'(/5x,"Problem in the HUBBARD card on line ",i5)') i
                     CALL errore( 'card_hubbard', &
                          & 'Mismatch in the quantum numbers for the same atomic type', i )
                  ENDIF
               ENDIF
            ELSEIF (counter_v(na,nb)==2 .OR. counter_v(na,nb)==3) THEN
               ! Second Hubbard manifold of the second atom
               ! Check whether we have different Hubbard manifolds for the same atomic type
               IF ( hu_n2==Hubbard_n(hu_nt2) .AND. hu_l2==Hubbard_l(hu_nt2) ) THEN
                  WRITE(stdout,'(/5x,"Problem in the HUBBARD card on line ",i5)') i
                  WRITE(stdout,'(5x,"Two Hubbard channels are the same for ",a)') atom_label(hu_nt2)
                  CALL errore( 'card_hubbard', &
                      'Not allowed to specify two Hubbard channels that are the same for the same atom', i )
               ENDIF
               IF (Hubbard_n2(hu_nt2)<0 .AND. Hubbard_l2(hu_nt2)<0) THEN
                  ! initialization
                  Hubbard_n2(hu_nt2) = hu_n2
                  Hubbard_l2(hu_nt2) = hu_l2
               ELSE
                  ! sanity check
                  IF (hu_n2/=Hubbard_n2(hu_nt2) .OR. hu_l2/=Hubbard_l2(hu_nt2)) THEN
                     WRITE(stdout,'(/5x,"Problem in the HUBBARD card for U on line ",i5)') i
                     CALL errore( 'card_hubbard', &
                          & 'Mismatch in the quantum numbers for the same atomic type', i )
                  ENDIF
               ENDIF
               IF (hu_n2_>-1 .AND. hu_l2_>-1) THEN
                  ! Third Hubbard manifold of the second atom
                  ! Check whether we have different Hubbard manifolds for the same atomic type
                  IF ( hu_n2_==Hubbard_n(hu_nt2)  .AND. hu_l2_==Hubbard_l(hu_nt2) .OR. &
                       hu_n2_==Hubbard_n2(hu_nt2) .AND. hu_l2_==Hubbard_l2(hu_nt2)) THEN
                     WRITE(stdout,'(/5x,"Problem in the HUBBARD card on line ",i5)') i
                     WRITE(stdout,'(5x,"Two Hubbard channels are the same for ",a)') atom_label(hu_nt)
                     CALL errore( 'card_hubbard', &
                         'Not allowed to specify two Hubbard channels that are the same for the same atom', i )
                  ENDIF
                  backall(hu_nt2) = .TRUE.
                  IF (Hubbard_n3(hu_nt2)<0 .AND. Hubbard_l3(hu_nt2)<0) THEN
                     ! initialization
                     Hubbard_n3(hu_nt2) = hu_n2_
                     Hubbard_l3(hu_nt2) = hu_l2_
                  ELSE
                     ! sanity check
                     IF (hu_n2_/=Hubbard_n3(hu_nt2) .OR. hu_l2_/=Hubbard_l3(hu_nt2)) THEN
                        WRITE(stdout,'(/5x,"Problem in the HUBBARD card for U on line ",i5)') i
                        CALL errore( 'card_hubbard', &
                             & 'Mismatch in the quantum numbers for the same atomic type', i )
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            !
            ! Sanity check 1
            IF (ityp(na)/=hu_nt .OR. ityp(nb)/=hu_nt2) THEN
               WRITE(stdout,'(/5x,"Problem in the HUBBARD card for V on line ",i5)') i
               CALL errore( 'card_hubbard', 'Mismatch between the atomic types and atomic indices', i)
            ENDIF
            ! Sanity check 2
            IF ((na==nb) .AND. (hu_nt/=hu_nt2)) THEN
               WRITE(stdout,'(/5x,"Problem in the HUBBARD card for V on line ",i5)') i
               CALL errore( 'card_hubbard', 'Atomic indices are equal but the atomic types are different', i)
            ENDIF
            !
         ENDIF
         !
         ! Read the value of the Hubbard parameter
         IF ( is_v ) THEN
            ! Column 6
            CALL get_field(6, hu_val, input_line)
         ELSE
            ! Column 3
            CALL get_field(3, hu_val, input_line)
         ENDIF
         hu_val = TRIM(hu_val)
         !
         IF ( LEN_TRIM(hu_val) < 1 ) &
            CALL errore( 'card_hubbard', &
                      'Value for the Hubbard parameter is missing', i )
         IF ( is_u ) THEN
            READ(hu_val,*, END=14, ERR=15) hu_u
         ELSEIF ( is_um ) THEN
            READ(hu_val,*, END=14, ERR=15) hu_um_temp
            ! cannot store hu_val directly in hu_um because we 
            ! need to read the indices of the eigenstates 
            ! to which hu_val will be applied
         ELSEIF ( is_j0 ) THEN
            READ(hu_val,*, END=14, ERR=15) hu_j0
         ELSEIF ( is_j ) THEN
            READ(hu_val,*, END=14, ERR=15) hu_j
         ELSEIF ( is_b ) THEN
            READ(hu_val,*, END=14, ERR=15) hu_b
         ELSEIF ( is_e2 ) THEN
            READ(hu_val,*, END=14, ERR=15) hu_e2
         ELSEIF ( is_e3 ) THEN
            READ(hu_val,*, END=14, ERR=15) hu_e3
         ELSEIF ( is_v ) THEN
            READ(hu_val,*, END=14, ERR=15) hu_v
         ELSEIF ( is_alpha ) THEN
            READ(hu_val,*, END=14, ERR=15) hu_alpha
         ELSEIF ( is_alpha_m ) THEN
            READ(hu_val,*, END=14, ERR=15) hu_alpha_m_temp 
         ELSE
            CALL errore( 'card_hubbard', &
                      'Incorrect name of the Hubbard parameter', i )
         ENDIF
         !
         IF ( is_um .OR. is_alpha_m ) THEN
            !
            ! ... orbital-resolved DFT+U: read and validate eigenvalue numbers:
            !     read the fields that follow the Hubbard U value (nfield>3)
            !     and check if these are integers between 1 and 2*l+1 (if nspin=1)
            !     or 2*(2*l+1) if nspin==2 or nspin==4 (noncolinear).
            !     If so, store all those indices in a vector, which we later
            !     use to transfer the hu_um_temp into the actual Hubbard_Um array.
            neigvals = nfield - 3
            ! sanity check: there can't be more eigenvalues than 2*(2l+1)
            IF ( neigvals > (2*ldim) ) THEN
            ! ... or more than (2l+1) in the closed-shell case
               IF ( nspin == 1 .AND. (.NOT. noncolin)) &
               CALL errore( 'card_hubbard', &
               'Too many target orbitals selected for orbital-resolved DFT+U', i )
            ENDIF
            !
            ALLOCATE(target_indices(neigvals)) ! allocate a vector to store the eigenvalue indices
            !
            DO field = 4, nfield
               ! start from 4 because the first three fields don't contain indices
               ! read eigenvalue indices and check if they are (allowed) integers
               CALL get_field(field, eigval_str, input_line)
               eigval_str = TRIM(eigval_str)
               READ(eigval_str, *, iostat=io_stat) eigval_index
               !
               IF (io_stat /= 0 )  CALL errore( 'card_hubbard', &
                     'Must use integer values to specify eigenvalue indices &
                         in orbital-resolved DFT+U.', i )
               IF ( eigval_index < 1 .OR. eigval_index > (2*ldim) )  CALL errore( 'card_hubbard', &
                     'The eigenvalues targeted by orbital-resolved U or ALPHA must range &
                      &from 1 to 2*(2*l+1)', i )
               IF ( (.NOT. noncolin) .AND. eigval_index > (nspin*ldim) )  CALL errore( 'card_hubbard', &
                     'The eigenvalues targeted by orbital-resolved U or ALPHA cannot be larger &
                      &than nspin*(2*l+1) in the colinear case.', i )
               !
               ! everything is alright: store the eigenvalue indices in the vector
               target_indices(field-3) = eigval_index
            ENDDO
            !
            ! Assign the orbital-resolved Hubbard parameter to the
            ! temporary array using the vector of target eigvals
            DO j = 1, SIZE(target_indices)
               IF ( (target_indices(j) <= ldim) .OR. noncolin ) THEN 
                  ! We are EITHER 1) in the spin-up channel with nspin==2,
                  ! 2) nspin==1, or 3) nspin==4 (noncolinear magnetism).
                  ! In this case, we also want to store indices 
                  ! >ldim in the "first" spin channel
                  IF ( is_um ) THEN
                     hu_um(target_indices(j),1) = hu_um_temp
                  ELSE
                     hu_alpha_m(target_indices(j),1) = hu_alpha_m_temp
                  ENDIF
               ELSE
                  ! If nspin == 2 and we're targeting a spin-down eigenvalue.
                  ! Example: if U=4.0 is assigned to eigval 6 in a 3d manifold,
                  ! store hu_um((6-5),2)=4.0
                  IF ( is_um ) THEN
                     hu_um(target_indices(j)-ldim, 2) = hu_um_temp 
                  ELSE
                     hu_alpha_m(target_indices(j)-ldim, 2) = hu_alpha_m_temp
                  ENDIF
               ENDIF 
            ENDDO
            DEALLOCATE(target_indices)
         ENDIF
         !
         ! Assign Hubbard parameters to the corresponding Hubbard arrays
         ! Allow positive and negative values of Hubbard parameters
         ! (in case users want to experiment with negative values)
         IF (is_u) THEN
            IF (counter_u(hu_nt)==1) THEN
               ! Hubbard parameter for the first (main) channel of the atomic type hu_nt
               IF (ABS(Hubbard_U(hu_nt))<eps16) THEN
                   Hubbard_U(hu_nt) = hu_u
               ELSE
                   WRITE(stdout,'(/5x,"Problem in the HUBBARD card for U on line ",i5)') i
                   CALL errore( 'card_hubbard', &
                           & 'U for this atomic type was already set', i )
                ENDIF
            ELSEIF (counter_u(hu_nt)==2) THEN
               ! Hubbard parameter for the second (and third) channel(s) of the atomic type hu_nt
               IF (ABS(Hubbard_U2(hu_nt))<eps16) THEN
                   Hubbard_U2(hu_nt) = hu_u
               ELSE
                   WRITE(stdout,'(/5x,"Problem in the HUBBARD card for U on line ",i5)') i
                   CALL errore( 'card_hubbard', &
                           & 'U for this atomic type was already set', i )
               ENDIF
            ENDIF
         ELSEIF (is_um .OR. is_alpha_m) THEN
            IF (noncolin) THEN
               DO m = 1, 2*ldim
                  IF ( ABS(hu_um(m,1)) > eps16 ) THEN 
                     !
                     ! make sure we're not overwriting values that were already set
                     IF ( ABS(Hubbard_Um_nc(m,hu_nt)) < eps16 ) THEN
                        !
                        Hubbard_Um_nc(m,hu_nt) = hu_um(m,1)
                     ELSE
                        WRITE(stdout,'(/5x,"Problem in the HUBBARD card for U on line ",i5)') i
                        CALL errore( 'card_hubbard', &
                                 & 'A Hubbard U parameter for this eigenvalue index was already set.', i )
                     ENDIF
                  !
                  ELSEIF ( ABS(hu_alpha_m(m,1)) > eps16 ) THEN
                     !
                     IF ( ABS(Hubbard_alpha_m_nc(m,hu_nt)) < eps16 ) THEN
                        Hubbard_alpha_m_nc(m,hu_nt) = hu_alpha_m(m,1)
                     ELSE
                        WRITE(stdout,'(/5x,"Problem in the HUBBARD card for ALPHA on line ",i5)') i
                        CALL errore( 'card_hubbard', &
                                 & 'A Hubbard ALPHA parameter for this eigenvalue index was already set.', i )
                     ENDIF
                  ENDIF
               ENDDO
            ELSE !colinear or closed-shell case
               DO is = 1, nspin
                  DO m = 1, ldim
                     IF ( ABS(hu_um(m,is)) > eps16 ) THEN 
                        !
                        IF ( ABS(Hubbard_Um(m,is,hu_nt)) < eps16 ) THEN
                           ! make sure we're not overwriting values that were already set
                           Hubbard_Um(m,is,hu_nt) = hu_um(m,is)
                        ELSE
                           WRITE(stdout,'(/5x,"Problem in the HUBBARD card for U on line ",i5)') i
                           CALL errore( 'card_hubbard', &
                                    & 'A Hubbard U parameter for this eigenvalue index was already set.', i )
                        ENDIF
                     !
                     ELSEIF ( ABS(hu_alpha_m(m,is)) > eps16 ) THEN
                        IF ( ABS(Hubbard_alpha_m(m,is,hu_nt)) < eps16 ) THEN
                           Hubbard_alpha_m(m,is,hu_nt) = hu_alpha_m(m,is)
                        ELSE
                           WRITE(stdout,'(/5x,"Problem in the HUBBARD card for ALPHA on line ",i5)') i
                           CALL errore( 'card_hubbard', &
                                    & 'A Hubbard ALPHA parameter for this eigenvalue index was already set.', i )
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
            DEALLOCATE(hu_um)
            DEALLOCATE(hu_alpha_m)
         !
         ELSEIF (is_j0) THEN
            IF (counter_j0(hu_nt)==1) THEN
               ! Hubbard parameter for the first (main) channel of the atomic type hu_nt
               IF (ABS(Hubbard_J0(hu_nt))<eps16) THEN
                   Hubbard_J0(hu_nt) = hu_j0
               ELSE
                   WRITE(stdout,'(/5x,"Problem in the HUBBARD card for J0 on line ",i5)') i
                   CALL errore( 'card_hubbard', &
                           & 'J0 for this atomic type was already set', i )
               ENDIF
            ELSEIF (counter_j0(hu_nt)==2) THEN
               CALL errore( 'card_hubbard', &
                       & 'Two channels for J0 for the same atomic type is not implemented', i )
            ENDIF
         ELSEIF (is_alpha) THEN
            IF (counter_alpha(hu_nt)==1) THEN
               IF (ABS(Hubbard_alpha(hu_nt)) < eps16) THEN
                   Hubbard_alpha(hu_nt) = hu_alpha
               ELSE
                   WRITE(stdout,'(/5x,"Problem in the HUBBARD card for ALPHA on line ",i5)') i
                   CALL errore( 'card_hubbard', &
                           & 'ALPHA for this atomic type was already set', i )
               ENDIF
            ELSEIF (counter_alpha(hu_nt)==2) THEN
               CALL errore( 'card_hubbard', &
                       & 'Two channels for ALPHA for the same atomic type is not implemented', i )
            ENDIF
         ELSEIF (is_j) THEN
            IF (ABS(Hubbard_J(1,hu_nt))<eps16) THEN
                Hubbard_J(1,hu_nt) = hu_j
            ELSE
                WRITE(stdout,'(/5x,"Problem in the HUBBARD card for J on line ",i5)') i
                CALL errore( 'card_hubbard', &
                        & 'J for this atomic type was already set', i )
            ENDIF
         ELSEIF (is_b) THEN
            IF (ABS(Hubbard_J(2,hu_nt))<eps16) THEN
                Hubbard_J(2,hu_nt) = hu_b
            ELSE
                WRITE(stdout,'(/5x,"Problem in the HUBBARD card for B on line ",i5)') i
                CALL errore( 'card_hubbard', &
                        & 'B for this atomic type was already set', i )
            ENDIF
         ELSEIF (is_e2) THEN
            IF (ABS(Hubbard_J(2,hu_nt))<eps16) THEN
                Hubbard_J(2,hu_nt) = hu_e2
            ELSE
                WRITE(stdout,'(/5x,"Problem in the HUBBARD card for E2 on line ",i5)') i
                CALL errore( 'card_hubbard', &
                        & 'E2 for this atomic type was already set', i )
            ENDIF
         ELSEIF (is_e3) THEN
            IF (ABS(Hubbard_J(3,hu_nt))<eps16) THEN
                Hubbard_J(3,hu_nt) = hu_e3
            ELSE
                WRITE(stdout,'(/5x,"Problem in the HUBBARD card for E3 on line ",i5)') i
                CALL errore( 'card_hubbard', &
                        & 'E3 for this atomic type was already set', i )
            ENDIF
         ELSEIF (is_v) THEN
            IF (ABS(Hubbard_V(na,nb,nc))<eps16) THEN
                Hubbard_V(na,nb,nc) = hu_v
            ELSE
                WRITE(stdout,'(/5x,"Problem in the HUBBARD card for V on line ",i5)') i
                CALL errore( 'card_hubbard', &
                     & 'Hubbard V for this couple was already set', i )
            ENDIF
         ENDIF
         !
      ENDDO
      !
11    CONTINUE
      !
      IF ( i > 0 ) THEN
         !     
         lda_plus_u = .TRUE.
         !
         ! We need to determine automatically which case we are dealing with,
         ! based on what Hubbard parameters we found in the HUBBARD card.
         ! Allow positive and negative values of Hubbard parameters
         ! (just in case if users want to experiment with negative values)
         !
         IF (ANY(ABS(Hubbard_J(:,:))>eps16)) THEN
            ! DFT+U+J
            lda_plus_u_kind = 1
            IF (ANY(ABS(Hubbard_J0(:))>eps16)) CALL errore('card_hubbard', &
                    & 'Hund J is not compatible with Hund J0', i)
            IF (ANY(ABS(Hubbard_V(:,:,:))>eps16)) CALL errore('card_hubbard', &
                    & 'Currently Hund J is not compatible with Hubbard V', i)
            IF (ANY(ABS(Hubbard_Um(:,:,:))>eps16)) CALL errore('card_hubbard', &
                    & 'Currently Hund J is not compatible with orbital-resolved Hubbard U', i)                    
            IF (ANY(ABS(Hubbard_alpha_m(:,:,:))>eps16)) CALL errore('card_hubbard', &
                    & 'Currently Hund J is not compatible with orbital-resolved Hubbard ALPHA', i)   
         ELSEIF (ANY(ABS(Hubbard_V(:,:,:))>eps16)) THEN
            ! DFT+U+V(+J0)
            lda_plus_u_kind = 2
            ! 
            IF (noncolin .and. ANY(Hubbard_J0(:)>eps16)) CALL errore('card_hubbard', &
                    & 'Currently Hund J0 is not compatible with noncolin=.true.', i)
            IF (ANY(ABS(Hubbard_Um(:,:,:))>eps16)) CALL errore('card_hubbard', &
                    & 'Currently DFT+U+V does not support orbital-resolved Hubbard U parameters', i)
            IF (ANY(ABS(Hubbard_alpha_m(:,:,:))>eps16)) CALL errore('card_hubbard', &
                    & 'Currently DFT+U+V does not support orbital-resolved Hubbard ALPHA parameters', i)
            !
         ELSEIF (ANY(ABS(Hubbard_Um(:,:,:))>eps16) .OR. ANY(ABS(Hubbard_alpha_m(:,:,:))>eps16) .OR. &
            & ANY(ABS(Hubbard_Um_nc(:,:))>eps16) .OR. ANY(ABS(Hubbard_alpha_m_nc(:,:))>eps16)) THEN
            ! Orbital-resolved DFT+U
            lda_plus_u_kind = 0
            orbital_resolved = .true.
            !
            IF (ANY(Hubbard_U(:)>eps16)) CALL errore('card_hubbard', &
                    & 'Cannot use shell-averaged Hubbard U parameters when using orbital-resolved DFT+U', i)
            IF (ANY(Hubbard_alpha(:)>eps16)) CALL errore('card_hubbard', &
                    & 'Cannot use shell-averaged Hubbard ALPHA parameters when using orbital-resolved DFT+U', i)
            IF (ANY(Hubbard_J0(:)>eps16)) CALL errore('card_hubbard', &
                    & 'Orbital-resolved DFT+U does not (yet) support Hund J0 parameters', i)
         ELSEIF (ANY(Hubbard_U(:)>eps16) .OR. ANY(Hubbard_J0(:)>eps16) .OR. &
                  ANY(ABS(Hubbard_alpha(:)) >eps16)) THEN
            ! DFT+U(+J0)
            lda_plus_u_kind = 0
            IF (noncolin .and. ANY(Hubbard_J0(:)>eps16)) CALL errore('card_hubbard', &
                    & 'Currently Hund J0 is not compatible with noncolin=.true.', i)
         ELSE
            CALL errore('card_hubbard', 'Unknown case for lda_plus_u_kind...', i)
         ENDIF
         !
         ! DFT+U+V: copy Hubbard_U to Hubbard_V
         IF (lda_plus_u_kind==2) THEN
            DO na = 1, nat
               nt = ityp(na)
               IF (counter_u(nt)==1) THEN
                  IF (Hubbard_V(na,na,1)<eps16 .AND. Hubbard_U(nt)>eps16) THEN
                      Hubbard_V(na,na,1) = Hubbard_U(nt)
                  ELSEIF (Hubbard_V(na,na,1)>eps16 .AND. Hubbard_U(nt)>eps16) THEN
                      CALL errore('card_hubbard', 'On-site Hubbard parameter for ' // &
                              & atom_label(nt) // ' was specified more than once', na)
                  ENDIF
               ELSEIF (counter_u(nt)==2) THEN
                  IF (Hubbard_V(na,na,3)<eps16 .AND. Hubbard_U(nt)>eps16) THEN
                      Hubbard_V(na,na,3) = Hubbard_U2(nt)
                  ELSEIF (Hubbard_V(na,na,3)>eps16 .AND. Hubbard_U(nt)>eps16) THEN
                      CALL errore('card_hubbard', 'On-site Hubbard parameter for ' // &
                              & atom_label(nt) // ' was specified more than once', na)
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
         !
      ENDIF
      !
      DEALLOCATE(counter_u)
      DEALLOCATE(counter_j0)
      DEALLOCATE(counter_j)
      DEALLOCATE(counter_b)
      DEALLOCATE(counter_e2)
      DEALLOCATE(counter_e3)
      DEALLOCATE(counter_alpha)
      IF (ALLOCATED(counter_v)) DEALLOCATE(counter_v)
      IF (ALLOCATED(ityp)) DEALLOCATE(ityp)
      tahub = .true.
      !
      RETURN
      !
14    CALL errore ('card_hubbard', ' End of file while parsing Hubbard parameters', 1)
15    CALL errore ('card_hubbard', ' Error while parsing Hubbard parameters', 1)
      !
   END SUBROUTINE card_hubbard
   !
   FUNCTION between( string, delimiter1, delimiter2, n )
      !
      ! Return what is found between characters delimiter1 and delimiter2
      ! if delimiter1 = '' , use beginning of string
      ! if delimiter2 = '' , use end of string
      ! if n >= 1 is present, use the n-th occurrence of delimiter1 and
      ! the first occurence of delimiter2 at the right of it
      !
      IMPLICIT NONE
      CHARACTER(len=:), ALLOCATABLE :: between
      CHARACTER(len=*), INTENT(IN)  :: string
      CHARACTER(len=*), INTENT(IN)  :: delimiter1
      CHARACTER(len=*), INTENT(IN)  :: delimiter2
      INTEGER, INTENT(IN), OPTIONAL :: n
      INTEGER :: n_, i_, i1,i2
      !
      between = ''
      n_ = 1
      IF ( PRESENT(n) ) n_ = n
      IF ( n_ < 1 ) RETURN
      !
      IF ( LEN(delimiter1) == 0 ) THEN
         IF ( n_ > 1 ) RETURN
         i1 = 1
      ELSE
         i1 = 1
         DO i_= 1, n_
            i1 = INDEX(string(i1:),delimiter1(1:1)) + i1
         ENDDO
         IF ( i1 < 2 ) RETURN
      ENDIF
      !
      IF ( LEN(delimiter2) == 0 ) THEN
         i2 = LEN_TRIM(string(i1:))
      ELSE
         i2 = INDEX(string(i1:),delimiter2(1:1)) - 1
         IF ( i2 < 1 ) RETURN
      ENDIF
      !
      between = ADJUSTL(TRIM(string(i1:i1+i2-1)))
      !
  END FUNCTION between
  !
END MODULE read_cards_module
