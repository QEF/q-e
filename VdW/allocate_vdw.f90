
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE allocate_vdw
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays:
  ! local potential for each kind of atom, structure factor
  !
  USE lsda_mod,  ONLY : nspin
  USE grid_dimensions, ONLY : nrxx
  USE wvfct,     ONLY : npwx, nbnd
  USE qpoint,    ONLY : nksq
  USE eff_v
  !
  IMPLICIT NONE
  !
  ALLOCATE (rho_fft (nrxx, nspin))
  ALLOCATE (rho_veff(nrxx, nspin))
  ALLOCATE (veff   (nrxx, nspin))
  !
  ALLOCATE (evc_veff(npwx, nbnd ))
  ALLOCATE (et_c(nbnd, nksq ))
  ALLOCATE (dvext(npwx, 3, nbnd))
  ALLOCATE (dpsi_eff(npwx, 3, nbnd))
  !
  RETURN
END SUBROUTINE allocate_vdw

