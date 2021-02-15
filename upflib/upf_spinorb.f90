!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE upf_spinorb
  !
  !! Variables needed for calculations with spin-orbit
  !
  USE upf_kinds,   ONLY : DP
  USE upf_params,  ONLY : lmaxx, lqmax 
  !
  !! FIXME: rot_ylm could be dynamically allocated
  !
  IMPLICIT NONE
  SAVE

  LOGICAL :: lspinorb
  !! if .TRUE. this is a spin-orbit calculation
  COMPLEX (DP) :: rot_ylm(lqmax,lqmax)
  !! transform real spherical harmonics into complex ones
  COMPLEX (DP), ALLOCATABLE :: fcoef(:,:,:,:,:)
  !! function needed to account for spinors.

END MODULE upf_spinorb

