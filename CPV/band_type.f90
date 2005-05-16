!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!  ----------------------------------------------
!  BEGIN manual

!=----------------------------------------------------------------------------=!
   MODULE band_type
!=----------------------------------------------------------------------------=!

!  this module contains the definitions of several TYPE structures,
!  together with their allocation/deallocation routines
!  ----------------------------------------------
!  routines in this module:
!  SUBROUTINE allocate_band(f,nel,nb_l,nb_g)
!  SUBROUTINE deallocate_band(f)
!  ----------------------------------------------
!  END manual

        USE kinds

        IMPLICIT NONE
        SAVE

        PRIVATE

!  ----------------------------------------------
!  BEGIN manual
!  TYPE DEFINITIONS

!! ...  occupation numbers
        TYPE band
          INTEGER nel   ! number of electrons per unit cell
          INTEGER nb_l  ! local number of bands
          INTEGER nb_g  ! global number of bands
          REAL(dbl), POINTER :: s(:) ! occupation numbers
        END TYPE band

!  END manual
!  ----------------------------------------------

        INTERFACE allocate_band
          MODULE PROCEDURE allocate_band_s, allocate_band_v
        END INTERFACE

        INTERFACE deallocate_band
          MODULE PROCEDURE deallocate_band_s, deallocate_band_v, deallocate_band_m
        END INTERFACE

!  end of module-scope declarations
!  ----------------------------------------------

!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!

!  ... subroutines

        SUBROUTINE allocate_band_s(f, nel, nb_l, nb_g)

!  (describe briefly what this routine does...)
!  ----------------------------------------------
          TYPE (band) :: f
          INTEGER, INTENT(IN) :: nel,nb_l,nb_g
          INTEGER :: ierr
          f%nel  = nel
          f%nb_l = nb_l
          f%nb_g = nb_g
          ALLOCATE( f%s( MAX(nb_l, 1) ), STAT = ierr )
          IF( ierr /= 0 ) CALL errore(' allocate_band ', ' cannot allocate ', ABS(ierr) )
          RETURN
        END SUBROUTINE allocate_band_s

!=----------------------------------------------------------------------------=!

        SUBROUTINE allocate_band_v(f, nel, nb_l, nb_g)
          TYPE (band) :: f(:)
          INTEGER, INTENT(IN) :: nel,nb_l,nb_g
          INTEGER :: i
          DO i = 1, SIZE( f )
            CALL allocate_band_s(f(i), nel, nb_l, nb_g)
          END DO
          RETURN
        END SUBROUTINE allocate_band_v

!=----------------------------------------------------------------------------=!

        SUBROUTINE deallocate_band_s(f)
          TYPE (band) :: f
          INTEGER :: ierr
          f%nel  = 0
          f%nb_l = 0
          f%nb_g = 0
          IF(ASSOCIATED(f%s)) THEN
            DEALLOCATE(f%s, STAT = ierr)
            IF( ierr /= 0 ) CALL errore(' deallocate_band ', ' deallocating ', ABS(ierr) )
          END IF
          RETURN
        END SUBROUTINE deallocate_band_s

!=----------------------------------------------------------------------------=!

        SUBROUTINE deallocate_band_v(f)
          TYPE (band) :: f(:)
          INTEGER :: i
          DO i = 1, SIZE( f )
            CALL deallocate_band_s(f(i))
          END DO
          RETURN
        END SUBROUTINE deallocate_band_v

!=----------------------------------------------------------------------------=!

        SUBROUTINE deallocate_band_m(f)
          TYPE (band) :: f(:,:)
          INTEGER :: i, j
          DO j = 1, SIZE( f, 2 )
            DO i = 1, SIZE( f, 1 )
              CALL deallocate_band_s(f(i,j))
            END DO
          END DO
          RETURN
        END SUBROUTINE deallocate_band_m   

!=----------------------------------------------------------------------------=!
   END MODULE band_type
!=----------------------------------------------------------------------------=!

