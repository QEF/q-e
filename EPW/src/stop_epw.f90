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
  subroutine stop_epw (flag)
  !-----------------------------------------------------------------------
  !
  ! Close all files and synchronize processes before stopping.
  ! Called at the end of the run with flag=.true. (removes 'recover')
  ! or during execution with flag=.false. (does not remove 'recover')
  !
#ifdef __PARA
  use mp,            ONLY: mp_end, mp_barrier
  USE mp_global,     ONLY :  inter_pool_comm, mp_global_end
#endif
  implicit none
  logical :: flag
  INTEGER :: ierr
  !
  CALL print_clock_epw
#ifdef __PARA
  CALL mp_barrier(inter_pool_comm)
  CALL mp_end(inter_pool_comm)
  !
  CALL mp_global_end()
#endif
#ifdef __T3E
  !
  ! set streambuffers off
  !
  CALL set_d_stream (0)
#endif

  CALL deallocate_part
  IF (flag) then
     stop
  ELSE
     stop 1
  ENDIF

  RETURN
end subroutine stop_epw
