!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!
subroutine start_pw4gww 
  !
  !  Usage: [mpirun, mpprun, whatever] postproc [-npool N]
  !
  !  Wrapper routine for postprocessing initialization
  !
  USE mp_global,     ONLY: mp_startup
  USE environment,   ONLY: environment_start
  implicit none
  character(len=9) :: code = 'PW4GWW'
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( code )
  ! 
  return
end subroutine start_pw4gww
