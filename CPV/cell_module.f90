!------------------------------------------------------------------------------!
    MODULE cell_module
!------------------------------------------------------------------------------!

      USE kinds, ONLY : dbl
      USE cell_base
!
      IMPLICIT NONE
      SAVE
!
      PRIVATE
      PUBLIC :: boxdimensions, r_to_s, s_to_r, cell_init

!
!     Cell variables read from stdin
!     They should never be changed.

      REAL(KIND=8), PUBLIC :: h(3,3), deth, hold(3,3), wmass
!
!
!------------------------------------------------------------------------------!
    END MODULE cell_module
!------------------------------------------------------------------------------!
