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
        REAL(dbl), PRIVATE  :: A1(3)
        REAL(dbl), PRIVATE  :: A2(3)
        REAL(dbl), PRIVATE  :: A3(3)

        REAL(dbl), PRIVATE :: celldm(6)
        INTEGER, PRIVATE :: IBRAV

        REAL(dbl), PRIVATE  :: PRESS
        REAL(dbl), PRIVATE  :: WC

        LOGICAL, PRIVATE :: taxis(3) = (/ .FALSE., .FALSE., .FALSE. /)

        DATA IBRAV/-1/
        DATA WC/0.0d0/


      contains

!------------------------------------------------------------------------------!


        SUBROUTINE get_lattice_vectors(a1_out,a2_out,a3_out)
          REAL(dbl), intent(out) :: a1_out(3), a2_out(3), a3_out(3)
            a1_out   = a1
            a2_out   = a2
            a3_out   = a3
          RETURN
        END SUBROUTINE get_lattice_vectors

!------------------------------------------------------------------------------!

        SUBROUTINE get_celldm( ibrav_out, celldm_out)
          REAL(dbl), intent(out) :: celldm_out(6)
          INTEGER, intent(out) :: ibrav_out
            ibrav_out = ibrav
            celldm_out = celldm
          RETURN
        END SUBROUTINE get_celldm

!------------------------------------------------------------------------------!

!
!------------------------------------------------------------------------------!
    END MODULE cell_module
!------------------------------------------------------------------------------!
