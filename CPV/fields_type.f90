!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
      MODULE fields_type

        USE kinds, ONLY: DP
        USE parallel_types, ONLY: descriptor, processors_grid
        IMPLICIT NONE
        PRIVATE
        SAVE

        TYPE scalar_field
          TYPE (descriptor) :: desc
          REAL(DP), POINTER :: r(:,:,:)
        END TYPE

        TYPE chden
          TYPE (descriptor) :: desc
          REAL(DP), POINTER :: r(:,:,:)
        END TYPE

        INTERFACE allocate_chden
          MODULE PROCEDURE allocate_chden_s, allocate_chden_v
        END INTERFACE
        INTERFACE deallocate_chden
          MODULE PROCEDURE deallocate_chden_s, deallocate_chden_v
        END INTERFACE
        INTERFACE ASSIGNMENT(=)
          MODULE PROCEDURE copy_chden_s, copy_chden_v, copy_chden_DP_s
        END INTERFACE

        !PUBLIC :: scalar_field, chden, allocate_chden, deallocate_chden, &
        !  ASSIGNMENT(=)

        CONTAINS

        SUBROUTINE ALLOCATE_CHDEN_S(f,desc) 
          TYPE (chden) :: f
          TYPE (descriptor) :: desc
          f%desc = desc
          ALLOCATE(f%r(desc%ldx,desc%ldy,desc%nzl))
          RETURN
        END SUBROUTINE ALLOCATE_CHDEN_S

        SUBROUTINE DEALLOCATE_CHDEN_S(f)
          TYPE (chden) :: f
          IF(ASSOCIATED(f%r)) DEALLOCATE(f%r)
          RETURN
        END SUBROUTINE DEALLOCATE_CHDEN_S

        SUBROUTINE ALLOCATE_CHDEN_V(f,desc) 
          TYPE (chden) :: f(:)
          TYPE (descriptor) :: desc
          INTEGER :: i
          DO i = 1, SIZE(f)
            CALL ALLOCATE_CHDEN_S(f(i),desc)
          END DO
          RETURN
        END SUBROUTINE ALLOCATE_CHDEN_V

        SUBROUTINE DEALLOCATE_CHDEN_V(f)
          TYPE (chden) :: f(:)
          INTEGER :: i
          DO i = 1, SIZE(f)
            CALL DEALLOCATE_CHDEN_S(f(i))
          END DO
          RETURN
        END SUBROUTINE DEALLOCATE_CHDEN_V

        SUBROUTINE ALLOCATE_SCALAR_FIELD(f,desc) 
          TYPE (scalar_field) :: f
          TYPE (descriptor) :: desc
          f%desc = desc
          ALLOCATE(f%r(desc%ldx,desc%ldy,desc%nzl))
          RETURN
        END SUBROUTINE ALLOCATE_SCALAR_FIELD

        SUBROUTINE DEALLOCATE_SCALAR_FIELD(f)
          TYPE (scalar_field) :: f
          IF(ASSOCIATED(f%r)) DEALLOCATE(f%r)
          RETURN
        END SUBROUTINE DEALLOCATE_SCALAR_FIELD

        SUBROUTINE copy_chden_s( rhoa, rhob )
          TYPE (chden), INTENT(INOUT) :: rhoa
          TYPE (chden), INTENT(IN)    :: rhob
          IF( SIZE( rhoa%r, 1) /= SIZE( rhob%r, 1) .OR. &
              SIZE( rhoa%r, 2) /= SIZE( rhob%r, 2) .OR. &
              SIZE( rhoa%r, 3) /= SIZE( rhob%r, 3) ) THEN
            CALL DEALLOCATE_CHDEN_S( rhoa )
            CALL ALLOCATE_CHDEN_S( rhoa, rhob%desc )
          END IF
          rhoa%r = rhob%r
          RETURN
        END SUBROUTINE copy_chden_s

        SUBROUTINE copy_chden_v( rhoa, rhob )
          TYPE (chden), INTENT(INOUT) :: rhoa(:)
          TYPE (chden), INTENT(IN)    :: rhob(:)
          INTEGER :: i
          DO i = 1, MIN( SIZE(rhoa), SIZE(rhob) )
            CALL copy_chden_s( rhoa(i), rhob(i) )
          END DO
          RETURN
        END SUBROUTINE copy_chden_v

        SUBROUTINE copy_chden_DP_s( rhoa, rb )
          TYPE (chden), INTENT(INOUT) :: rhoa
          REAL (DP), INTENT(IN)    :: rb(:,:,:)
          IF( SIZE( rhoa%r, 1) /= SIZE( rb, 1) .OR. &
              SIZE( rhoa%r, 2) /= SIZE( rb, 2) .OR. &
              SIZE( rhoa%r, 3) /= SIZE( rb, 3) ) THEN
            CALL errore(' charge_density: copy_chden_DP_s ', ' sizes do not match ', 0)
          END IF
          rhoa%r = rb
          RETURN
        END SUBROUTINE copy_chden_DP_s



      END MODULE fields_type
