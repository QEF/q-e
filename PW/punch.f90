!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine punch
  !-----------------------------------------------------------------------
  !
  !     This routine is called at the end of the run to save on a file
  !     the information needed to the phonon program.
  !
  !
  USE io_global, ONLY : stdout
  USE klist, ONLY: nks, nkstot
  USE lsda_mod, ONLY: nspin
  USE scf, ONLY: rho
  USE control_flags, ONLY: reduce_io, lscf
  USE wvfct, ONLY: et, wg, nbnd
  USE wavefunctions_module, ONLY : evc, evc_nc
  USE io_files, ONLY: prefix, iunpun, iunwfc, nwordwfc
  USE noncollin_module, ONLY: noncolin
#ifdef __PARA
  use para
#endif
  use restart_module, only: writefile_new
  !
  implicit none
  !
  integer :: ik, i, ibnd, kunittmp
  logical :: exst
  !
  WRITE( stdout, '(/,5x,"Writing file ",a14," for program phonon")') &
     trim(prefix)//'.save'
  !
  kunittmp = 1
  !
  ! if the wavefunction has not been written on file, do it now
  !
  if (noncolin) then
     if (nks.eq.1.and.reduce_io) call davcio (evc_nc, nwordwfc, iunwfc, 1, + 1)
  else
     if (nks.eq.1.and.reduce_io) call davcio (evc, nwordwfc, iunwfc, 1, + 1)
  endif
  !
  ! The following instruction is used  when more k-points are needed
  ! for finite-q phonon calculations (on fine q-grid) then those needed
  ! for self-consistency. In such a case, a self-consistent calculation
  ! with few k-points is followed by a non-self-consistent one with added
  ! k-points, whose weight is set to zero.
  !
  if (.not.lscf) call sum_band
  !
  !  Write: general variables (including dimensions of the arrays),
  !  atomic positions, forces, k-points, eigenvalues
  !
#ifdef __PARA
  !
  ! xk, wk, isk, et, wg are distributed across pools
  ! the first node has a complete copy of xk, wk, isk,
  ! while eigenvalues et and weights wg must be
  ! explicitely collected to the first node
  !
  call poolrecover (et, nbnd, nkstot, nks)
  call poolrecover (wg, nbnd, nkstot, nks)
  !
  ! In parallel execution, only the first node writes this file
  !
  kunittmp = kunit
  !
#endif
  !
  !  Write the charge density on a separate file
  !
  call io_pot ( + 1, trim(prefix)//'.rho', rho, nspin)

  iunpun = 4
  call writefile_new( 'all', iunpun, et, wg, kunittmp )

  return
end subroutine punch
