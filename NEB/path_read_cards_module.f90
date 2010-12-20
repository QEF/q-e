!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
MODULE path_read_cards_module
   !---------------------------------------------------------------------------
   !
   ! ...  This module handles the reading of cards from standard input
   ! ...  Written by Carlo Cavazzoni and modified for "path" implementation
   ! ...  by Carlo Sbraccia
   !
   USE kinds,     ONLY : DP
   USE io_global, ONLY : stdout
   USE constants, ONLY : angstrom_au
   USE parser,    ONLY : parse_unit,field_count, read_line, get_field
   USE io_global, ONLY : ionode, ionode_id
   !
   USE path_input_parameters_module
   USE input_parameters, ONLY : atom_mass, atom_label, ntyp, nat, taspc, &
                                if_pos, sp_pos, na_inp, allocate_input_ions, &
                                atomic_positions, tapos
   !
   IMPLICIT NONE
   !
   SAVE
   !
   PRIVATE
   !
   PUBLIC :: path_read_cards
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
   !
   !----------------------------------------------------------------------
   SUBROUTINE path_read_cards(unit)
      !----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: unit
      !
      CHARACTER(len=256)         :: input_line
      CHARACTER(len=80)          :: card
      CHARACTER(len=1), EXTERNAL :: capital
      LOGICAL                    :: tend
      INTEGER                    :: i
      !
      !
      parse_unit = unit
      !
100   CALL read_line( input_line, end_of_file=tend )
      !
      IF( tend ) GOTO 120
      IF( input_line == ' ' .or. input_line(1:1) == '#' ) GOTO 100
      !
      READ (input_line, *) card
      !
      DO i = 1, len_trim( input_line )
         input_line( i : i ) = capital( input_line( i : i ) )
      ENDDO
      !
      IF ( trim(card) == 'PATH_ATOMIC_POSITIONS' ) THEN
         !
         CALL card_path_atomic_positions( input_line )
         !
      ELSEIF( trim(card) =='CLIMBING_IMAGES' ) THEN
         !
         CALL card_climbing_images( input_line )
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
   END SUBROUTINE path_read_cards

   !
   ! ... Description of the allowed input CARDS
   !
   !------------------------------------------------------------------------
   !    BEGIN manual
   !----------------------------------------------------------------------
   !
   !
   ! ATOMIC_POSITIONS
   !
   !
   !----------------------------------------------------------------------
   !    END manual
   !------------------------------------------------------------------------
   !
   ! ... routine modified for NEB           ( C.S. 21/10/2003 )
   ! ... routine modified for SMD           ( Y.K. 15/04/2004 )
   !
   SUBROUTINE card_path_atomic_positions( input_line )
      !
      !
      IMPLICIT NONE
      !
      CHARACTER(len=256) :: input_line
      CHARACTER(len=4)   :: lb_pos
      INTEGER            :: ia, k, is, nfield, idx, rep_i
      LOGICAL, EXTERNAL  :: matches
      LOGICAL            :: tend
      LOGICAL, SAVE      :: tread = .false.
      !
      INTEGER            :: ifield, ierr
      REAL(DP)           :: field_value
      CHARACTER(len=256) :: field_str, error_msg
      !
      !
      IF ( tread ) THEN
         CALL errore( 'card_path_atomic_positions', 'two occurrences', 2 )
      ENDIF
      IF ( .not. taspc ) THEN
         CALL errore( 'card_atomic_positions', &
                     & 'ATOMIC_SPECIES must be present before', 2 )
      ENDIF
      IF ( nat < 1 ) THEN
         CALL errore( 'card_atomic_positions', 'nat out of range', nat )
      ENDIF
      !
      !
!new
      CALL allocate_input_ions(ntyp,nat)
      !
      !
      if_pos = 1
      !
      sp_pos = 0
      na_inp = 0
      !
      IF ( matches( "CRYSTAL", input_line ) ) THEN
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
         atomic_positions = 'alat'
      ENDIF
      !

         !
         IF ( allocated( pos ) ) DEALLOCATE( pos )
         ALLOCATE( pos( 3*nat, num_of_images ) )
         pos(:,:) = 0.0_DP
         !
            !
            CALL read_line( input_line, end_of_file = tend )
            IF ( tend ) &
               CALL errore( 'read_cards', &
                           'end of file reading atomic positions (path)', 1 )
            !
            IF ( matches( "first_image", input_line ) ) THEN
               !
               input_images = 1
               CALL path_read_images( input_images )
               !
            ELSE
               !
               CALL errore( 'read_cards', &
                           'first_image missing in ATOMIC_POSITION', 1 )
               !
            ENDIF
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
                  exit read_conf_loop
                  !
               ENDIF
               !
            ENDDO read_conf_loop
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
            ENDIF
            !
         !
      tread = .true.
      tapos = .true.
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
         INTEGER, INTENT(in) :: image
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
            ELSEIF ( nfield == 7 ) THEN
               !
               IF ( image /= 1 ) THEN
                  !
                  CALL errore( 'read_cards', &
                              & 'wrong number of columns in ' // &
                              & 'ATOMIC_POSITIONS', sp_pos(ia) )
                  !
               ENDIF
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
            ENDIF
            !
            IF ( image == 1 ) THEN
               !
               lb_pos = adjustl( lb_pos )
               !
               match_label_path: DO is = 1, ntyp
                  !
                  IF ( trim( lb_pos ) == trim( atom_label(is) ) ) THEN
                     !
                     sp_pos(ia) = is
                     !
                     exit match_label_path
                     !
                  ENDIF
                  !
               ENDDO match_label_path
               !
               IF ( ( sp_pos(ia) < 1 ) .or. ( sp_pos(ia) > ntyp ) ) THEN
                  !
                  CALL errore( 'read_cards', &
                                 'wrong index in ATOMIC_POSITIONS', ia )
                  !
               ENDIF
               !
               is = sp_pos(ia)
               !
               na_inp( is ) = na_inp( is ) + 1
               !
            ENDIF
            !
         ENDDO
         !
         RETURN
         !
         END SUBROUTINE path_read_images
         !
   END SUBROUTINE card_path_atomic_positions
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
      CHARACTER(len=256) :: input_line
      LOGICAL, SAVE      :: tread = .false.
      LOGICAL, EXTERNAL  :: matches
      !
      INTEGER          :: i
      CHARACTER(len=5) :: i_char
      !
      CHARACTER(len=6), EXTERNAL :: int_to_char
      !
      !
      IF ( tread ) &
         CALL errore( ' card_climbing_images ', ' two occurrences', 2 )
      !
      IF ( CI_scheme == 'manual' ) THEN
         !
         IF ( allocated( climbing ) ) DEALLOCATE( climbing )
         !
         ALLOCATE( climbing( num_of_images ) )
         !
         climbing(:) = .false.
         !
         CALL read_line( input_line )
         !
         DO i = 1, num_of_images
            !
            i_char = int_to_char( i )
            !
            IF ( matches( ' ' // trim( i_char ) // ',' , &
                           ' ' // trim( input_line ) // ',' ) ) &
               climbing(i) = .true.
            !
         ENDDO
         !
      ENDIF
      !
      tread = .true.
      !
      RETURN
      !
   END SUBROUTINE card_climbing_images
   !
END MODULE path_read_cards_module
