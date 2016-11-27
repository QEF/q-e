!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=---------------------------------------------------------------------------=!
MODULE fft_interfaces

  IMPLICIT NONE
  PRIVATE


  PUBLIC :: fwfft, invfft

  
  INTERFACE invfft
     !! invfft is the interface to both the standard fft **invfft_x**,
     !! and to the "box-grid" version **invfft_b**, used only in CP 
     !! (the latter has an additional argument)
     
     SUBROUTINE invfft_x( grid_type, f, dfft, dtgs, howmany )
       USE fft_types,  ONLY: fft_type_descriptor
       USE task_groups,   ONLY: task_groups_descriptor
       USE fft_param,  ONLY :DP
       IMPLICIT NONE
       CHARACTER(LEN=*),  INTENT(IN) :: grid_type
       TYPE(fft_type_descriptor), INTENT(IN) :: dfft
       TYPE(task_groups_descriptor), OPTIONAL, INTENT(IN) :: dtgs
       INTEGER, OPTIONAL, INTENT(IN) :: howmany
       COMPLEX(DP) :: f(:)
     END SUBROUTINE invfft_x
     !
     SUBROUTINE invfft_b( f, dfft, ia )
       USE fft_smallbox_type,  ONLY: fft_box_descriptor
       USE fft_param,  ONLY :DP
       IMPLICIT NONE
       INTEGER, INTENT(IN) :: ia
       TYPE(fft_box_descriptor), INTENT(IN) :: dfft
       COMPLEX(DP) :: f(:)
     END SUBROUTINE invfft_b
  END INTERFACE

  INTERFACE fwfft
     SUBROUTINE fwfft_x( grid_type, f, dfft, dtgs, howmany )
       USE fft_types,  ONLY: fft_type_descriptor
       USE task_groups,   ONLY: task_groups_descriptor
       USE fft_param,  ONLY :DP
       IMPLICIT NONE
       CHARACTER(LEN=*), INTENT(IN) :: grid_type
       TYPE(fft_type_descriptor), INTENT(IN) :: dfft
       TYPE(task_groups_descriptor), OPTIONAL, INTENT(IN) :: dtgs
       INTEGER, OPTIONAL, INTENT(IN) :: howmany
       COMPLEX(DP) :: f(:)
     END SUBROUTINE fwfft_x
  END INTERFACE

END MODULE fft_interfaces
!=---------------------------------------------------------------------------=!
