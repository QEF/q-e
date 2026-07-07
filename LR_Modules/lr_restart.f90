!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
MODULE lr_restart
   !! Define abstract interface for write_rec subroutine that can be passed to dfpt_kernels.
   !
   ABSTRACT INTERFACE
      !--------------------------------------------------------------------------
      SUBROUTINE write_rec_interface( where, irr, dr2, iter, convt, dfpt_data )
         !! This routine saves the information needed to recover the DFPT calculation.
         USE kinds,     ONLY : DP
         USE dfpt_type, ONLY : dfpt_data_type
         IMPLICIT NONE
         CHARACTER(LEN=10), INTENT(IN) :: where
         INTEGER, INTENT(IN) :: irr, iter
         LOGICAL, INTENT(IN) :: convt
         REAL(DP), INTENT(IN) :: dr2
         TYPE(dfpt_data_type), INTENT(INOUT) :: dfpt_data
      END SUBROUTINE write_rec_interface
      !--------------------------------------------------------------------------
      SUBROUTINE stop_callback_interface( flag )
         !! This routine is called to stop the calculation gracefully.
         IMPLICIT NONE
         LOGICAL, INTENT(IN) :: flag
      END SUBROUTINE stop_callback_interface
      !--------------------------------------------------------------------------
   END INTERFACE
END MODULE lr_restart
