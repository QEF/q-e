!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!



      MODULE pseudo_projector

        USE kinds

        IMPLICIT NONE
        SAVE

        PRIVATE

        TYPE projector
          LOGICAL :: gamma_symmetry
          REAL(dbl), POINTER    :: r(:,:,:)
          COMPLEX(dbl), POINTER :: c(:,:,:)
        END TYPE

        INTERFACE allocate_projector
          MODULE PROCEDURE allocate_projector_s, allocate_projector_v, allocate_projector_m
        END INTERFACE
        INTERFACE deallocate_projector 
          MODULE PROCEDURE deallocate_projector_s, deallocate_projector_v, &
            deallocate_projector_m
        END INTERFACE

        PUBLIC :: projector, allocate_projector, deallocate_projector
          
      CONTAINS
!
!----------------------------------------------------------------------
!

        SUBROUTINE allocate_projector_s(fnl, nsanl, nx, ngh, gamma)
          IMPLICIT NONE
          TYPE (projector) :: fnl
          INTEGER, INTENT(IN) :: nsanl, nx, ngh
          LOGICAL, INTENT(IN) :: gamma
          INTEGER :: ierr
          fnl%gamma_symmetry = gamma
          IF(gamma) THEN
            ALLOCATE(fnl%r(MAX(nsanl,1), MAX(ngh,1), nx), STAT=ierr )
            IF( ierr /= 0 ) CALL errore(' allocate_projector ',' allocating fnl%r ',ierr)
            NULLIFY(fnl%c)
            fnl%r = 0.0d0
          ELSE
            ALLOCATE(fnl%c(MAX(nsanl,1), MAX(ngh,1), nx), STAT=ierr )
            IF( ierr /= 0 ) CALL errore(' allocate_projector ',' allocating fnl%c ',ierr)
            NULLIFY(fnl%r)
            fnl%c = 0.0d0
          END IF
          RETURN
        END SUBROUTINE

        SUBROUTINE allocate_projector_v(fnl, nsanl, nx, ngh, gamma)
          IMPLICIT NONE
          TYPE (projector) :: fnl(:)
          INTEGER, INTENT(IN) :: nsanl, nx, ngh
          LOGICAL, INTENT(IN) :: gamma
          INTEGER :: i
          DO i = 1, SIZE(fnl)
            CALL allocate_projector_s(fnl(i), nsanl, nx, ngh, gamma)
          END DO
          RETURN
        END SUBROUTINE

        SUBROUTINE allocate_projector_m(fnl, nsanl, nx, ngh, gamma)
          IMPLICIT NONE
          TYPE (projector) :: fnl(:,:)
          INTEGER, INTENT(IN) :: nsanl, nx, ngh
          LOGICAL, INTENT(IN) :: gamma
          INTEGER :: i, j
          DO j = 1, SIZE(fnl, 2)
            DO i = 1, SIZE(fnl, 1)
              CALL allocate_projector_s(fnl(i,j), nsanl, nx, ngh, gamma)
            END DO
          END DO
          RETURN
        END SUBROUTINE


        SUBROUTINE deallocate_projector_s(fnl)
          IMPLICIT NONE
          TYPE (projector) :: fnl
          INTEGER :: ierr
          IF( ASSOCIATED(fnl%r) ) THEN
            DEALLOCATE(fnl%r, STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_projector ',' deallocating fnl%r ',ierr)
          END IF
          IF( ASSOCIATED(fnl%c) ) THEN
            DEALLOCATE(fnl%c, STAT=ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_projector ',' deallocating fnl%c ',ierr)
          END IF
          RETURN
        END SUBROUTINE

        SUBROUTINE deallocate_projector_v(fnl)
          IMPLICIT NONE
          TYPE (projector) :: fnl(:)
          INTEGER :: i
          DO i = 1, SIZE(fnl)
            CALL deallocate_projector_s(fnl(i))
          END DO
          RETURN
        END SUBROUTINE

        SUBROUTINE deallocate_projector_m(fnl)
          IMPLICIT NONE
          TYPE (projector) :: fnl(:,:)
          INTEGER :: i, j
          DO j = 1, SIZE(fnl, 2)
            DO i = 1, SIZE(fnl, 1)
              CALL deallocate_projector_s(fnl(i,j))
            END DO
          END DO
          RETURN
        END SUBROUTINE


      END MODULE pseudo_projector
