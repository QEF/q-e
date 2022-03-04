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
  !!    AF: this module could be merged with uspp, where
  !         other SO variables are
  !
  IMPLICIT NONE
  SAVE

  LOGICAL :: is_spinorbit
  !! if .TRUE. this is a spin-orbit calculation
  !! internal flag, set by uspp_allocate, to be used only in upflib
  COMPLEX (DP) :: rot_ylm(lqmax,lqmax)
  !! transform real spherical harmonics into complex ones
  COMPLEX (DP), ALLOCATABLE :: fcoef(:,:,:,:,:)
  !! function needed to account for spinors.
  !
  ! GPU vars
  COMPLEX(DP), ALLOCATABLE :: fcoef_d(:,:,:,:,:)
#if defined(__CUDA)
  attributes (DEVICE) :: fcoef_d
#endif

END MODULE upf_spinorb

