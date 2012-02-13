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
  USE kinds,      ONLY: DP
  USE io_global,  ONLY : stdout
  USE io_files,   ONLY : iunwfc, nwordwfc, iunres, prefix, seqopn
  USE cell_base,  ONLY: omega, alat
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp
  USE ener,       ONLY: etot, ehart, etxc, vtxc, epaw
  USE fft_base,   ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft
  USE gvect,      ONLY: gstart, g, gg, nl, ngm
  USE klist,      ONLY: nks
  USE lsda_mod,   ONLY: nspin
  USE scf,        ONLY : rho, rho_core, rhog_core, v, vnew
  USE ldaU,       ONLY : eth
  USE control_flags, ONLY: restart, tr2, ethr
  USE wvfct,      ONLY: nbnd, et
  USE wavefunctions_module,    ONLY : evc, psic
  USE uspp,  ONLY: becsum
  USE paw_variables,  ONLY: okpaw, ddd_PAW
  USE paw_onecenter, ONLY : PAW_potential
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
  read (iunres, err=10, end=10) vnew%of_r
  close (unit = iunres, status = 'keep')
  WRITE( stdout, '(5x,"Calculation restarted from IONS ",i3)')
  !
  ! store wavefunctions in memory here if there is just one k-point
  ! (otherwise it is never done)
  !
  if (nks.eq.1) call davcio (evc,2*nwordwfc, iunwfc, 1, -1)
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
     CALL fwfft ('Dense', psic, dfftp)
     !
     rho%of_g(:,is) = psic(nl(:))
     !
  END DO
  
  ! recalculate etxc, vtxc, ehart, needed by stress calculation
  !
  CALL v_of_rho( rho, rho_core, rhog_core, &
                 ehart, etxc, vtxc, eth, etotefield, charge, v )
  IF (okpaw) CALL PAW_potential(rho%bec, ddd_PAW, epaw)
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
