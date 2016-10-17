  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  !
  ! Copyright (C) 2001 PWSCF group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  !
  ! Modified from stop_ph
  !
  !-----------------------------------------------------------------------
  subroutine stop_epw
  !-----------------------------------------------------------------------
  !!
  !! Close all files and synchronize processes before stopping.
  !! Called at the end of the run (removes 'recover')
  !!
  use mp,            ONLY : mp_end, mp_barrier
  USE mp_global,     ONLY : inter_pool_comm, mp_global_end
  USE io_global,     ONLY : stdout
  ! 
  implicit none
  !
  CALL print_clock_epw
write(stdout,'(a)') "                                                                                          "
write(stdout,'(a)') " Please consider citing:                                                                  "
write(stdout,'(a)') " S. Ponce, E. R. Margine, C. Verdi and F. Giustino, Comput. Phys. Commun. 209, 116 (2016) " 
write(stdout,'(a)') "                                                                                          "
  CALL mp_barrier(inter_pool_comm)
  CALL mp_end(inter_pool_comm)
  !
  CALL mp_global_end()
!#ifdef __T3E
!  !
!  ! set streambuffers off
!  !
!  CALL set_d_stream (0)
!#endif

!  CALL deallocate_part
  ! 
  STOP
  ! 
  RETURN
end subroutine stop_epw
