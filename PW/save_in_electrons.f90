!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine save_in_electrons (iter, dr2)
  !-----------------------------------------------------------------------
  use pwcom
  use io_files, only: iunres
  implicit none
  character :: where * 20
  ! are we in the right place?
  integer :: ik, ibnd, ik_, iter
  ! counters
  ! last completed kpoint
  ! last completed iteration
  logical :: exst


  real(kind=DP) :: dr2
  if (reduce_io) return
  !
  ! open recover file
  !
  call seqopn (iunres, 'restart', 'unformatted', exst)
  !
  ! save restart information
  !
  if (conv_elec) then
     !
     ! step on electrons has been completed. restart from ions
     !
     where = 'IONS'
     write (iunres) where
     write (iunres) ( (et (ibnd, ik), ibnd = 1, nbnd), ik = 1, nks)
     write (iunres) etot, tr2
     ! vnew = V(in)-V(out) is needed in the scf correction term to forces
     write (iunres) vnew
  else
     !
     ! save iteration number
     !
     ! iteration iter has been completed
     ik_ = 0
     where = 'ELECTRONS'
     write (iunres) where
     write (iunres) ( (et (ibnd, ik), ibnd = 1, nbnd), ik = 1, nks)
     write (iunres) iter, ik_, dr2, tr2, ethr
  endif
  !

  close (unit = iunres, status = 'keep')
  return

end subroutine save_in_electrons
