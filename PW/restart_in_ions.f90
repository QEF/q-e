!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine restart_in_ions (iter, ik_, dr2)
  !-----------------------------------------------------------------------
  use pwcom
  USE wavefunctions,    ONLY : evc, psic
  implicit none
  character :: where * 20
  ! are we in the right place?
  integer :: ik, ibnd, ik_, iter
  ! counters
  ! last completed kpoint
  ! last completed iteration
  ! check for bravais lattice
  ! check for number of atoms
  logical :: exst

  real(kind=DP) :: dr2, charge
  call seqopn (iunres, 'restart', 'unformatted', exst)

  if (.not.exst) goto 10
  read (iunres, err = 10, end = 10) where
  !
  ! is this the right place where to restart ?
  !
  if (where.ne.'IONS') then
     close (unit = iunres, status = 'keep')
     write (*,*) where, '.......?'
     call errore ('restart_i', ' we should not be here ...!', 1)
  endif
  !
  !  read recover information
  !
  read (iunres, err=10, end=10) ( (et(ibnd,ik), ibnd=1,nbnd), ik=1,nks)
  read (iunres, err=10, end=10) etot, tr2
  ! vnew = V(in)-V(out) is needed in the scf correction term to forces
  read (iunres, err=10, end=10) vnew
  close (unit = iunres, status = 'keep')
  write (6, '(5x,"Calculation restarted from IONS ",i3)')
  !
  ! store wavefunctions in memory here if there is just one k-point
  ! (otherwise it is never done)
  !
  if (nks.eq.1) call davcio (evc, nwordwfc, iunwfc, 1, -1)
  !
  ! recalculate rho
  !
  call sum_band
  !
  ! recalculate etxc, vtxc, ehart, needed by stress calculation
  !
  call v_of_rho (rho, rho_core, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
       nrxx, nl, ngm, gstart, nspin, g, gg, alat, omega, &
       ehart, etxc, vtxc, charge, psic)
  !
  !  restart procedure completed
  !

  restart = .false.

  return
  !
  ! in case of problems
  !

10 call errore ('restart_i', 'problems in reading recover file', -1)
  return

end subroutine restart_in_ions
