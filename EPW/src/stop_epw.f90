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
  SUBROUTINE stop_epw
  !-----------------------------------------------------------------------
  !!
  !! Close all files and synchronize processes before stopping.
  !! Called at the end of the run
  !!
  use mp,            ONLY : mp_end, mp_barrier
  USE mp_global,     ONLY : inter_pool_comm, mp_global_end
  USE io_global,     ONLY : stdout
  USE printing,      ONLY : print_clock_epw
  ! 
  implicit none
  !
  CALL print_clock_epw
  write(stdout,'(a)') "                                                                                          "
  write(stdout,'(a)') " Please consider citing:                                                                  "
  write(stdout,'(a)') " S. Ponce, E. R. Margine, C. Verdi and F. Giustino, Comput. Phys. Commun. 209, 116 (2016) " 
  write(stdout,'(a)') "                                                                                          "
  write(stdout,'(a)') " In addition, if you used anisotropic Eliashberg superconductivity please cite:           "
  write(stdout,'(a)') "              E. R. Margine and F. Giustino, Phys. Rev. B 87, 024505 (2013)               "
  write(stdout,'(a)') "              if you used transport properties (scattering rates, mobility) please cite:  "
  write(stdout,'(a)') "              S. Ponce, E. R. Margine and F. Giustino, Phys. Rev. B 97, 121201 (2018)     "
  CALL mp_barrier(inter_pool_comm)
  CALL mp_end(inter_pool_comm)
  !
  CALL mp_global_end()
  ! 
  STOP
  ! 
  RETURN
END SUBROUTINE stop_epw
