
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
subroutine allocate_vdw
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays:
  ! local potential for each kind of atom, structure factor
  !
  USE lsda_mod,  ONLY : nspin
  USE gvect,     ONLY : nrxx
  USE wvfct,     ONLY : npwx, nbnd
  USE qpoint,    ONLY : nksq
  USE eff_v
  !
  implicit none
  !
  allocate (rho_fft (nrxx, nspin))    
  allocate (rho_veff(nrxx, nspin))    
  allocate (veff   (nrxx, nspin))    
  !
  allocate (evc_veff(npwx, nbnd ))    
  allocate (et_c(nbnd, nksq ))    
  allocate (dvext(npwx, 3, nbnd))
  allocate (dpsi_eff(npwx, 3, nbnd))
  !
  return
end subroutine allocate_vdw

