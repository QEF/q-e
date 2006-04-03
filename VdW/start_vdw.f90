!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine start_vdw (nodenumber)
  !
  !  Usage: [mpirun, mpprun, whatever] postproc [-npool N]
  !
  !  Wrapper routine for postprocessing initialization
  !
  USE global_version, ONLY: version_number
  USE wvfct,  ONLY: gamma_only
  implicit none
  character(len=3) :: nodenumber
  character(len=9) :: code = 'VdW'
  !
  call startup (nodenumber, code, version_number)
  return
end subroutine start_vdw
