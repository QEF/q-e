!
! Copyright (C) 2001-2011 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
! ... Common variables for LR_Modules routines
!
!
MODULE qpoint
  USE kinds, ONLY :  DP
  USE parameters, ONLY : npk
  !
  ! ... The q point
  !
  SAVE
  !
  INTEGER, POINTER :: igkq(:)     ! npwx)
  ! correspondence k+q+G <-> G
  INTEGER :: nksq, npwq, nksqtot
  ! the real number of k points
  ! the number of plane waves for q
  ! the total number of q points 
  INTEGER, ALLOCATABLE :: ikks(:), ikqs(:)
  ! the index of k point in the list of k
  ! the index of k+q point in the list of k
  REAL (DP) :: xq(3)
  ! the coordinates of the q point
  COMPLEX (DP), ALLOCATABLE :: eigqts(:) ! nat)
  ! the phases associated to the q
  REAL (DP), ALLOCATABLE :: xk_col(:,:)
  !
END MODULE qpoint
!
MODULE control_lr
  USE kinds, ONLY :  DP
  USE parameters, ONLY: npk
  !
  ! ... the variable controlling the phonon run
  !
  SAVE
  !
  INTEGER :: nbnd_occ(npk)   ! occupated bands in metals
  !
  LOGICAL :: lgamma          ! if .TRUE. this is a q=0 computation
  !
  REAL(DP) :: alpha_pv       ! the alpha value for shifting the bands
  !
END MODULE control_lr
!
MODULE eqv
  USE kinds, ONLY :  DP
  !
  ! ... The wavefunctions at point k+q
  !
  SAVE
  !
  COMPLEX (DP), POINTER :: evq(:,:)
  !
  ! ... The variable describing the linear response problem
  !
  COMPLEX (DP), ALLOCATABLE :: dvpsi(:,:), dpsi(:,:), drhoscfs (:,:,:)
  ! the product of dV psi
  ! the change of the wavefunctions
  REAL (DP), ALLOCATABLE :: dmuxc(:,:,:)        ! nrxx, nspin, nspin),
  REAL (DP), ALLOCATABLE, TARGET :: vlocq(:,:)  ! ngm, ntyp)
  ! the derivative of the xc potential
  ! the local potential at q+G
  REAL (DP), ALLOCATABLE :: eprec(:,:) ! needed for preconditioning
  !
END MODULE eqv
!
!
MODULE gc_lr
  USE kinds, ONLY :  DP
  !
  ! ... The variables needed for gradient corrected calculations
  !
  SAVE
  !
  REAL (DP), ALLOCATABLE :: &
       grho(:,:,:),              &! 3, nrxx, nspin),
       gmag(:,:,:),              &! 3, nrxx, nspin),
       vsgga(:),                 &! nrxx
       segni(:),                 &! nrxx
       dvxc_rr(:,:,:),           &! nrxx, nspin, nspin), &
       dvxc_sr(:,:,:),           &! nrxx, nspin, nspin),
       dvxc_ss(:,:,:),           &! nrxx, nspin, nspin), &
       dvxc_s(:,:,:)              ! nrxx, nspin, nspin)
  !
  ! in the noncollinear case gmag contains the gradient of the magnetization
  ! grho the gradient of rho+ and of rho-, the eigenvalues of the spin density
  ! vsgga= 0.5* (V_up-V_down) to be used in the calculation of the change
  ! of the exchange and correlation magnetic field.
  ! gradient of the unpert. density
  !
  ! derivatives of the E_xc functiona
  ! r=rho and s=|grad(rho)|
  !
END MODULE gc_lr
!

