!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
! ... Common variables for the vdw program
!
MODULE eff_v
  !
  USE kinds,  ONLY : DP
  !
  SAVE
  !
  REAL (kind=DP) :: nelecr
  !
  REAL (kind=DP), ALLOCATABLE ::            &
     !
     rho_fft(:,:),      &  !charge density using cut-off in FFT
     rho_veff(:,:),     &  !charge density from v_eff
     veff(:,:)            !effective pot in TFvW scheme
     !
  COMPLEX (kind=DP), ALLOCATABLE :: &
     !
     evc_veff (:,:),    &     !eigenvectors corespond to v_eff
     et_c (:,:),        &     !eigenvalues plus imag. freq.
     dvext(:,:,:),      &     !r.h.s. of linear equation
     dpsi_eff(:,:,:)             !used to keep dpsi
     !
END MODULE eff_v
!
MODULE control_vdw
  !
  USE kinds,  ONLY : DP
  !
  SAVE
  !
  INTEGER, PARAMETER :: maxter_vdw = 100    ! Max. # of iterations
  !
  INTEGER ::            niter_vdw, nmix_vdw
  ! Max. # of iterations
  ! number of iterations used in potential mixing
  !
  REAL (kind=DP) ::  tr2_vdw, thresh_veff, al_mix_vdw(maxter_vdw)
  ! Threshold for iterative solution of linear system
  ! thresh_hold for iterative optimization of Veff
  ! mixing factor (for each iteration) for updating the scf potential:
  !
END MODULE control_vdw
