!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Author: L. Martin-Samos 
!
subroutine start_pw4gww (nodenumber)
  !
  !  Usage: [mpirun, mpprun, whatever] postproc [-npool N]
  !
  !  Wrapper routine for postprocessing initialization
  !
  USE global_version, ONLY: version_number
  implicit none
  character(len=3) :: nodenumber
  character(len=9) :: code = 'PW4GWW'
  !
#if defined __INTEL
  ! ... Intel compilers v .ge.8 allocate a lot of stack space
  ! ... Stack limit is often small, thus causing SIGSEGV and crash
  CALL remove_stack_limit ( )
#endif
  call startup (nodenumber, code, version_number)
  return
end subroutine start_pw4gww
