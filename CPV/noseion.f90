!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE nose_ions
!=----------------------------------------------------------------------------=!

        USE kinds
        USE control_flags, ONLY: tnosep
        USE ions_nose, ONLY: &
          xnos0 => xnhp0, &
          xnosm => xnhpm, &
          xnos2m => xnhpm2, &
          xnosp => xnhpp, &
          vnhp => vnhp, &
          qnp, &
          tempw, & 
          fnosep, &
          ndega, &
          gkbt

        IMPLICIT NONE
        SAVE

        PRIVATE

        REAL(dbl) :: ekincm = 0.0d0
        REAL(dbl) :: enosep = 0.0d0

        PUBLIC :: enosep
        PUBLIC :: nosep_velocity
        PUBLIC :: movenosep

!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!


      REAL(dbl) FUNCTION nosep_velocity()
        USE time_step, ONLY: delt
        nosep_velocity = ( 3.D0 * xnos0(1) - 4.D0 * xnosm(1) + xnos2m(1) ) / ( 2.0d0 * delt )
        RETURN
      END FUNCTION nosep_velocity

      SUBROUTINE movenosep( ekinp )
        USE time_step, only: delt
        IMPLICIT NONE
        REAL(dbl), INTENT(IN) :: ekinp
        xnosp(1)  = 2.D0 * xnos0(1) - xnosm(1) + 2.D0 * delt * delt / qnp(1) * ( ekinp -  0.5D0 * gkbt )
        vnhp(1) = ( xnosp(1) - xnosm(1) ) / ( 2.0d0 * delt ) 
        RETURN
      END SUBROUTINE movenosep

!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
   END MODULE nose_ions
!=----------------------------------------------------------------------------=!
