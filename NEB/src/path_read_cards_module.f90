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
   USE io_global, ONLY : meta_ionode, meta_ionode_id
   USE mp,        ONLY : mp_bcast
   USE mp_world,  ONLY : world_comm
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
      CHARACTER(len=80)          :: climbing_images
      CHARACTER(len=1), EXTERNAL :: capital
      INTEGER                    :: i
      !
      climbing_images = ' '
      !
      IF ( meta_ionode) THEN
         !
         read_loop:  DO
            !
            READ (unit, fmt='(A256)', END=100, ERR=101) input_line
            !
            IF( input_line == ' ' .or. input_line(1:1) == '#' ) CYCLE
            !
            DO i = 1, len_trim( input_line )
               input_line( i : i ) = capital( input_line( i : i ) )
            ENDDO
            !
            IF( trim(adjustl(input_line)) == 'CLIMBING_IMAGES' ) THEN
               !
               READ (unit, fmt='(A80)', END=101, ERR=101) climbing_images
               !
            ELSE
               !
               CALL infomsg ('read_cards_module',&
                    'card '//trim(input_line)//' ignored' )
               !
            ENDIF
            !
         END DO read_loop
         !
      END IF
      !
100   CONTINUE
      !
      CALL mp_bcast(climbing_images, meta_ionode_id, world_comm)
      CALL card_climbing_images( climbing_images )
      !
      RETURN
      !
101   CALL errore ('read_cards_module','error reading neb.dat file', 1 )
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
   SUBROUTINE card_climbing_images( climbing_images )
      !
      IMPLICIT NONE
      !
      CHARACTER(len=80)  :: climbing_images
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
         DO i = 1, num_of_images
            !
            i_char = int_to_char( i )
            !
            IF ( matches( ' ' // trim( i_char ) // ',' , &
                          ' ' // trim( climbing_images ) // ',' ) ) &
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
