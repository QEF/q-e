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
  USE io_global,  ONLY : stdout
  USE io_files,   ONLY : iunwfc, nwordwfc, iunres, prefix
  USE kinds, ONLY: DP
  USE cell_base, ONLY: omega, alat
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE ener,  ONLY: etot, ehart, etxc, vtxc
  USE gvect, ONLY: gstart, g, gg, nl, ngm, nr1,nr2,nr3, nrx1,nrx2,nrx3, &
       nrxx
  USE klist, ONLY: nks
  USE lsda_mod, ONLY: nspin
  USE scf, ONLY : rho, rho_core, rhog_core
  USE control_flags, ONLY: restart, tr2, ethr
  USE vlocal, ONLY: vnew
  USE wvfct, ONLY: nbnd, et
  USE wavefunctions_module,    ONLY : evc, psic
  implicit none
  character :: where * 20
  ! are we in the right place?
  integer :: ik, is, ibnd, ik_, iter
  ! counters
  ! last completed kpoint
  ! last completed iteration
  ! check for bravais lattice
  ! check for number of atoms
  logical :: exst

  real(DP) :: dr2, charge, etotefield
  call seqopn (iunres, 'restart', 'unformatted', exst)

  if (.not.exst) goto 10
  read (iunres, err = 10, end = 10) where
  !
  ! is this the right place where to restart ?
  !
  if (where.ne.'IONS') then
     close (unit = iunres, status = 'keep')
     WRITE( stdout,*) where, '.......?'
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
  WRITE( stdout, '(5x,"Calculation restarted from IONS ",i3)')
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
  ! ... bring rho to G-space
  !
  DO is = 1, nspin
     !
     psic(:) = rho%of_r(:,is)
     !
     CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
     !
     rho%of_g(:,is) = psic(nl(:))
     !
  END DO
  
  ! recalculate etxc, vtxc, ehart, needed by stress calculation
  !
  CALL v_of_rho( rho%of_r, rho%of_g, rho_core, rhog_core, rho%kin_r, &
                 ehart, etxc, vtxc, etotefield, charge, psic )
  !
  !  restart procedure completed
  !

  restart = .false.

  return
  !
  ! in case of problems
  !

10 call infomsg ('restart_i', 'problems in reading recover file')
  return

end subroutine restart_in_ions
