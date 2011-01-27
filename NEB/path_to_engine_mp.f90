!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE engine_mp_start()
  !-----------------------------------------------------------------------------
  !
  ! fill engine modules for parallel calculations
  !
  USE kinds,         ONLY : DP
  !
  USE mp_image_global_module, ONLY : intra_image_comm
  USE mp_image_global_module, ONLY : root_image, me_image, nproc_image_file, &
                                    nproc_image, my_image_id, nimage
  !
  USE mp_global, ONLY : mp_startup_new, &
                        intra_image_comm_ => intra_image_comm, &
                        me_image_ => me_image, &
                        root_image_ => root_image, &
                        nproc_image_file_ => nproc_image_file, &
                        nproc_image_ => nproc_image, &
                        my_image_id_ => my_image_id, &
                        nimage_ => nimage  
  !
  ! ... "path" specific
  !
  IMPLICIT NONE
  !
  intra_image_comm_ = intra_image_comm
  me_image_ = me_image
  root_image_ = root_image
  nproc_image_ = nproc_image
  nproc_image_file = nproc_image_file
  my_image_id_ = my_image_id
  nimage_ = nimage
  !
  CALL mp_startup_new(root_image, intra_image_comm) 
  !
  !
  RETURN
  !
END SUBROUTINE engine_mp_start
!
