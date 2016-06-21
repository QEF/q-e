!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
! ... Common variables for LR_Modules routines
!
MODULE qpoint
  !
  USE kinds,      ONLY : DP
  !
  ! ... The variables needed to specify various indices,
  ! ... number of plane waves and k points and their coordiantes.
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
  !
  USE kinds,      ONLY : DP
  !
  ! ... The variables controlling the run of linear response codes
  !
  SAVE
  !
  INTEGER, ALLOCATABLE :: nbnd_occ(:)  ! occupied bands in metals
  REAL(DP) :: alpha_pv       ! the alpha value for shifting the bands
  LOGICAL  :: lgamma         ! if .TRUE. this is a q=0 computation
  LOGICAL  :: lrpa           ! if .TRUE. uses the Random Phace Approximation
  !
END MODULE control_lr
!
MODULE eqv
  !
  USE kinds,  ONLY : DP
  !
  ! ... The variables describing the linear response problem
  !
  SAVE
  !
  COMPLEX (DP), POINTER :: evq(:,:)
  ! the wavefunctions at point k+q
  COMPLEX (DP), ALLOCATABLE :: dvpsi(:,:), dpsi(:,:), drhoscfs (:,:,:)
  ! the product of dV psi
  ! the change of the wavefunctions
  REAL (DP), ALLOCATABLE :: dmuxc(:,:,:)        ! nrxx, nspin, nspin)
  ! the derivative of the xc potential
  REAL (DP), ALLOCATABLE, TARGET :: vlocq(:,:)  ! ngm, ntyp)
  ! the local potential at q+G
  !
END MODULE eqv
!
MODULE gc_lr
  !
  USE kinds, ONLY : DP
  !
  ! ... The variables needed for gradient corrected calculations
  !
  SAVE
  !
  REAL (DP), ALLOCATABLE :: &
       grho(:,:,:),              &! 3, nrxx, nspin)
       gmag(:,:,:),              &! 3, nrxx, nspin)
       vsgga(:),                 &! nrxx)
       segni(:),                 &! nrxx)
       dvxc_rr(:,:,:),           &! nrxx, nspin, nspin)
       dvxc_sr(:,:,:),           &! nrxx, nspin, nspin)
       dvxc_ss(:,:,:),           &! nrxx, nspin, nspin)
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
MODULE lr_symm_base
  !
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to describe the modes and the small group of q
  !
  SAVE
  !
  INTEGER :: irgq(48), nsymq=0, irotmq
  ! selects the operations of the small group
  ! the number of symmetry of the small group
  ! selects the symmetry sending q <-> -q+G
  REAL (DP), ALLOCATABLE :: rtau(:,:,:) !3, 48, nat)
  ! coordinates of direct translations
  REAL (DP) :: gi(3,48), gimq(3)
  ! the possible G associated to each symmetry
  ! the G associated to the symmetry q<->-q+G
  LOGICAL :: minus_q, & ! if .TRUE. there is the symmetry sending q<->-q
             invsymq    ! if .TRUE. the small group of q has inversion
  !
END MODULE lr_symm_base
!
MODULE lrus
  !
  USE kinds,  ONLY : DP
  USE becmod, ONLY : bec_type
  !
  ! ... These are additional variables needed for the linear response
  ! ... with US pseudopotentials and a generic perturbation Delta Vscf
  !
  SAVE
  !
  COMPLEX (DP), ALLOCATABLE :: &
       int3(:,:,:,:,:),     &! nhm, nhm, nat, nspin, npert)
       int3_paw(:,:,:,:,:), &! nhm, nhm, nat, nspin, npert)
       int3_nc(:,:,:,:,:)    ! nhm, nhm, nat, nspin, npert)
  ! int3 -> \int (Delta V_Hxc) Q d^3r
  ! similarly for int_nc while
  ! int3_paw contains Delta (D^1-\tilde D^1)
  !
  REAL (DP), ALLOCATABLE ::    dpqq(:,:,:,:)       ! nhm, nhm, 3, ntyp)
  COMPLEX (DP), ALLOCATABLE :: dpqq_so(:,:,:,:,:)  ! nhm, nhm, nspin, 3, ntyp)
  ! dpqq and dpqq_so: dipole moment of each Q multiplied by the fcoef factors 
  !
  type (bec_type), ALLOCATABLE, TARGET :: becp1(:) ! nksq)
  ! becp1 contains < beta_n | psi_i >
  !
  REAL (DP),    ALLOCATABLE :: bbg(:,:)      ! nkb, nkb)
  ! for gamma_only     
  COMPLEX (DP), ALLOCATABLE :: bbk(:,:,:)    ! nkb, nkb, nks)
  ! for k points
  COMPLEX (DP), ALLOCATABLE :: bbnc(:,:,:,:) ! nkb, nkb, nspin_mag, nks)
  ! for the noncollinear case
  ! bbg = < beta^N_i | beta^P_j > 
  ! bbg/bbk/bbnc are the scalar products of beta functions 
  ! localized on atoms N and P.
  !
END MODULE lrus
