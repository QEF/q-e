!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

    MODULE util

      USE kinds, ONLY : dbl
      IMPLICIT NONE
      SAVE

      PRIVATE

      INTEGER, PARAMETER :: i4b = kind(0)
      INTEGER, PARAMETER :: sp = kind(1._dbl)
      INTEGER, PARAMETER :: spc = kind((1._dbl,1._dbl))
      INTEGER, PARAMETER :: lgt = kind(.TRUE.)

      INTERFACE swap
        MODULE PROCEDURE swap_i, swap_r, swap_rv, swap_c, swap_cv, swap_cm, &
          masked_swap_rs, masked_swap_rv, masked_swap_rm
      END INTERFACE

      PUBLIC :: swap

!!*****
    CONTAINS
!BL
      SUBROUTINE swap_i(a,b)
        INTEGER (i4b), INTENT (INOUT) :: a, b
        INTEGER (i4b) :: dum

        dum = a
        a = b
        b = dum
      END SUBROUTINE swap_i
!BL
      SUBROUTINE swap_r(a,b)
        REAL (sp), INTENT (INOUT) :: a, b
        REAL (sp) :: dum

        dum = a
        a = b
        b = dum
      END SUBROUTINE swap_r
!BL
      SUBROUTINE swap_rv(a,b)
        REAL (sp), DIMENSION (:), INTENT (INOUT) :: a, b
        REAL (sp), DIMENSION (size(a)) :: dum

        dum = a
        a = b
        b = dum
      END SUBROUTINE swap_rv
!BL
      SUBROUTINE swap_c(a,b)
        COMPLEX (spc), INTENT (INOUT) :: a, b
        COMPLEX (spc) :: dum

        dum = a
        a = b
        b = dum
      END SUBROUTINE swap_c
!BL
      SUBROUTINE swap_cv(a,b)
        COMPLEX (spc), DIMENSION (:), INTENT (INOUT) :: a, b
        COMPLEX (spc), DIMENSION (size(a)) :: dum

        dum = a
        a = b
        b = dum
      END SUBROUTINE swap_cv
!BL
      SUBROUTINE swap_cm(a,b)
        COMPLEX (spc), DIMENSION (:,:), INTENT (INOUT) :: a, b
        COMPLEX (spc), DIMENSION (size(a,1),size(a,2)) :: dum

        dum = a
        a = b
        b = dum
      END SUBROUTINE swap_cm
!BL
      SUBROUTINE masked_swap_rs(a,b,mask)
        REAL (sp), INTENT (INOUT) :: a, b
        LOGICAL (lgt), INTENT (IN) :: mask
        REAL (sp) :: swp

        IF (mask) THEN
          swp = a
          a = b
          b = swp
        END IF
      END SUBROUTINE masked_swap_rs
!BL
      SUBROUTINE masked_swap_rv(a,b,mask)
        REAL (sp), DIMENSION (:), INTENT (INOUT) :: a, b
        LOGICAL (lgt), DIMENSION (:), INTENT (IN) :: mask
        REAL (sp), DIMENSION (size(a)) :: swp

        WHERE (mask)
          swp = a
          a = b
          b = swp
        END WHERE
      END SUBROUTINE masked_swap_rv
!BL
      SUBROUTINE masked_swap_rm(a,b,mask)
        REAL (sp), DIMENSION (:,:), INTENT (INOUT) :: a, b
        LOGICAL (lgt), DIMENSION (:,:), INTENT (IN) :: mask
        REAL (sp), DIMENSION (size(a,1),size(a,2)) :: swp

        WHERE (mask)
          swp = a
          a = b
          b = swp
        END WHERE
      END SUBROUTINE masked_swap_rm
!BL
!BL
    END MODULE util
