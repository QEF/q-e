!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Sat Feb 12 11:43:48 MET 2000
!  ----------------------------------------------
!  BEGIN manual

      MODULE time_step

!  (describe briefly what this module does...)
!  ----------------------------------------------
!  routines in this module:
!  SUBROUTINE set_time_step(dt)
!  ----------------------------------------------
!  END manual
!  ----------------------------------------------

        USE kinds
        IMPLICIT NONE
        SAVE

        PRIVATE

!  ...  declare module-scope variables
        REAL(DP)  :: delthal, twodelt, fordt2, dt2, dt2by2, delt
        REAL(DP)  :: tps ! elapsed simulated time in picoseconds

        PUBLIC :: set_time_step, tps, delt, twodelt, dt2, dt2by2

!  end of module-scope declarations
!  ----------------------------------------------

      CONTAINS

!  subroutines
!  ----------------------------------------------
!  ----------------------------------------------
        SUBROUTINE set_time_step(dt)

           REAL(DP), INTENT(IN) :: dt

           delt    = dt
           dt2     = dt ** 2
           fordt2  = 4.0_DP * dt2
           delthal = 0.5_DP * delt
           twodelt = 2.0_DP * delt
           dt2by2  = 0.5_DP * dt2
           tps     = 0.0_DP 

           RETURN
        END SUBROUTINE set_time_step

!  ----------------------------------------------
!  ----------------------------------------------

      END MODULE time_step

