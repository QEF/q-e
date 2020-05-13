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
   USE constants, ONLY : angstrom_au
   USE parser,    ONLY : parse_unit,field_count, read_line, get_field
   USE io_global, ONLY : meta_ionode
   !
   USE path_input_parameters_module
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
      IF( trim(card) =='CLIMBING_IMAGES' ) THEN
         !
         CALL card_climbing_images( input_line )
         !
      ELSE
         !
         IF ( meta_ionode ) CALL infomsg ('read_cards_module',&
            'card '//trim(input_line)//' ignored' )
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
