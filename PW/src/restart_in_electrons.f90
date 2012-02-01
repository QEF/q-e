!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine restart_in_electrons (iter, ik_, dr2)
  !-----------------------------------------------------------------------
  USE io_global,  ONLY : stdout
  USE io_files,   ONLY : iunwfc, nwordwfc, iunres, prefix
  USE kinds, ONLY: DP
  USE klist, ONLY: nks
  USE control_flags, ONLY: restart, tr2, ethr
  USE wvfct, ONLY: nbnd, et
  USE noncollin_module, ONLY: noncolin
  USE wavefunctions_module,    ONLY : evc
  USE exx,                     ONLY : exx_restart, &
                                      fock0, fock1, fock2, dexx, x_occupation
  implicit none
  character :: where * 20
  ! are we in the right place?
  integer :: ik, ibnd, ik_, iter_, iter
  ! counters
  ! last completed kpoint
  ! iteration number when program crashed
  ! last completed iteration
  logical :: exst
  logical :: l_exx_was_active

  real(DP) :: dr2
  call seqopn (iunres, 'restart', 'unformatted', exst)

  if (.not.exst) goto 10
  read (iunres, err=10, end=10) where
  !
  ! is  this the right place where to restart ?
  !
  if (where.ne.'ELECTRONS') then
     close (unit = iunres, status = 'keep')
     !
     ! this is a signal for the calling routine saying we are in the wrong place
     !
     ik_ = - 1000
     return
  endif

  read (iunres) ( (et(ibnd,ik), ibnd=1,nbnd), ik=1,nks)
  read (iunres, err=10, end=10) iter_, ik_, dr2, tr2, ethr
  read (iunres, err=10, end=10) l_exx_was_active, fock0, fock1, fock2, dexx
  IF (l_exx_was_active) read (iunres) &
     ( (x_occupation(ibnd,ik), ibnd=1,nbnd), ik=1,nks)
  call exx_restart(l_exx_was_active)

  close (unit = iunres, status = 'keep')
  if (ik_.eq.0) then
     iter = iter_
     WRITE( stdout, '(5x,"Calculation restarted from first kpoint ", &
          &" of iteration #",i3)') iter + 1
  elseif (ik_.ne.nks) then
     iter = iter_ - 1
     WRITE( stdout, '(5x,"Calculation restarted from kpoint #",i6, &
          &" of iteration #",i3)') ik_ + 1, iter + 1
  else
     iter = iter_ - 1
     WRITE( stdout, '(5x,"Calculation restarted from charge/pot", &
          &" of iteration #",i3)') iter + 1
     !
     ! with only one k-point wavefunctions are not read in sum_band
     !
     if (nks.eq.1) call davcio (evc,2*nwordwfc, iunwfc, 1, - 1)
  endif

  WRITE( stdout, '(5x,"tr2 = ",1pe8.2," ethr = ",1pe8.2)') tr2, ethr
  !
  !  restart procedure completed
  !

  restart = .false.

  return
  !
  ! in case of problems
  !

10 call infomsg ('restart_e', 'problems in reading recover file')
  close ( unit=iunres, status='keep')
  restart = .false.
  !
  return
end subroutine restart_in_electrons
