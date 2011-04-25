!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
MODULE step_penalty
!-------------------------------------------------------------------------
  !
  ! LDA+U with occupation constraint 
  !
  USE kinds
  implicit none
  integer :: natx
  real(DP) :: E_pen = 0.d0
  real(DP), allocatable :: A_pen(:,:), sigma_pen(:), alpha_pen(:)
  logical :: step_pen

CONTAINS
  !
  subroutine ldaUpen_init ( natx_, step_pen_, sigma_pen_, alpha_pen_, A_pen_ )
  !-----------------------------------------------------------------------
  !
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: natx_
      LOGICAL, INTENT(IN) :: step_pen_
      REAL(DP),INTENT(IN) :: sigma_pen_(natx_), alpha_pen_(natx_), A_pen_(natx_,2)
      
      step_pen=step_pen_
      IF ( step_pen ) THEN
         allocate (A_pen(natx,2), sigma_pen(natx), alpha_pen(natx) )
         sigma_pen=sigma_pen_
         alpha_pen=alpha_pen_
         A_pen=A_pen_
      END IF
  END SUBROUTINE ldaUpen_init
  !
  subroutine deallocate_step_pen()
  !-----------------------------------------------------------------------
  !
     IF( ALLOCATED( alpha_pen ) ) DEALLOCATE( alpha_pen )
     IF( ALLOCATED( sigma_pen ) ) DEALLOCATE( sigma_pen )
     IF( ALLOCATED( A_pen ) ) DEALLOCATE( A_pen )
     !
  end subroutine
  !
end module step_penalty
