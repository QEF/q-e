!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! OBM
! 150608  now obsolote nbrx resides here

module lr_variables
  !--------------------------------------------------------------------------
  ! ... sets the dimensions of the variables required in the
  ! ... Lanczos/linear response calculation
  !--------------------------------------------------------------------------
  !
  use kinds,                only : dp
  use control_flags,                only : gamma_only
  !
  implicit none
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Parameters 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER, PARAMETER:: nbrx=14   ! max number of beta functions
  !
  integer :: iund0psi   = 20
  integer :: iunrestart = 20
  integer :: nwordd0psi
  integer :: nwordrestart
  !
  integer :: n_ipol
  !
  !
  integer :: size_evc
  !
  !LOGICAL :: nlcc_any  ! .T. if any atom-type has nlcc
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Global Variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !REAL (DP), ALLOCATABLE :: dmuxc(:,:,:)        ! nrxx, nspin, nspin),
  !
  !REAL (DP), ALLOCATABLE :: &
  !     grho(:,:,:),              &! 3, nrxx, nspin),
  !     dvxc_rr(:,:,:),           &! nrxx, nspin, nspin), &
  !     dvxc_sr(:,:,:),           &! nrxx, nspin, nspin),
  !     dvxc_ss(:,:,:),           &! nrxx, nspin, nspin), &
  !     dvxc_s(:,:,:)              ! nrxx, nspin, nspin)
  !
  !
  real(kind=dp), allocatable :: becp1(:,:)
  complex(kind=dp), allocatable :: becp1_c(:,:,:)
  real(kind=dp), allocatable :: becp1_virt(:,:)
  complex(kind=dp), allocatable :: becp1_c_virt(:,:,:)

  ! 
  complex(kind=dp), allocatable :: & 
       evc0(:,:,:),       &    ! the ground state wavefunctions (plane wave, band, k point)
       evc0_virt(:,:,:),  &    ! unoccupied ground state wavefunctions (plane wave, band, k point)
       sevc0(:,:,:),      &    ! S * ground state wavefunctions
       sevc0_virt(:,:,:), &    ! S * virtual ground state wavefunctions
       evc1_old(:,:,:,:), &    ! response wavefunctions in the pw basis (last
                               ! index 1: q' using rotated SBR 2: p') 
       evc1(:,:,:,:),     &    !  "    "
       evc1_new(:,:,:,:), &    !  "    "
       sevc1_new(:,:,:,:),&    ! S * "    "
       d0psi(:,:,:,:)          ! for saving the original starting vectors
  !
  complex(kind=dp), allocatable :: revc0(:,:,:)
  !
  real(kind=dp), allocatable :: &
       rho_1(:,:),           &         ! response charge density in real space
       !rho_tot(:),        &         ! ground state + resp. charge density in real space (obm: is it used at all?)
       rho_1_tot(:,:)                !response charge density (mode 2) 
  !
  !integer, allocatable :: &
  !     igk_k(:,:),&         
  !     npw_k(:)
  ! 
  integer :: &
        nbnd_total               !Actual number of bands calculated by PWSCF (virtual+ocuppied)
  !
  integer, allocatable :: cube_save(:,:) !used in response charge density mode 1
  !
  !open shell related...
  !
  !real(kind=dp) :: lr_alpha_pv ! Spread in eigenvalues
  !integer, allocatable :: nbnd_occ(:) !number of occupied points for the given k point
  !integer, allocatable :: &
  !     igk_k(:,:),&                   ! The g<->k correspondance for each k point 
  !     npw_k(:)                       ! number of plane waves at each k point
                           ! They are (used many times, it is much better to hold them in memory 
 
  !
  !Lanczos Matrix
  !
  !
  !       | alpha(:,1),  gamma(:,2),                             0|
  !       | beta(:,2) ,  alpha(:,2), gamma(:,3)                   |
  ! T^j = |           ,  beta(:,3) , ...                          |
  !       |                                                       |
  !       |                                       ..., gamma(:,j) |
  !       |0                               beta(:,j) , alpha(:,j) |
  !
  ! Zeta is the \sum_valance (V^T_j * r_i ) where r_i is the density operator acting
  ! on ground state orbitals
  ! 
  ! zeta.w_T gives the polarizability (w_T is the solution of 
  ! (\omega-L)e_1 = w_T , this is handled in a post processing program)
  ! 

  real(kind=dp), allocatable :: &  ! (pol, iter)
       alpha_store(:,:),&      
       beta_store(:,:),&
       gamma_store(:,:)
  complex(kind=dp), allocatable :: zeta_store(:,:,:)  !polarization of external field, polarization of internal field, iteration number.
  !
  !The currently processed polarization direction and Lanczos iteration
  !
  INTEGER           :: LR_iteration, LR_polarization !polarization directions 1:x 2:y 3:z
  !
  ! variables for diagonalising the coefficient matrix
  !
  real(kind=dp) :: lanc_norm
  real(kind=dp), allocatable :: eval1(:),eval2(:)
  real(kind=dp), allocatable :: vl(:,:),vr(:,:)
  !
   REAL(kind=dp) :: norm0(3) 
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! lr_input:
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !
  logical :: restart            ! set True if the calculation is a restart run
  integer :: restart_step       ! the amount of steps to write a restart file
  !
  integer :: lr_verbosity       ! verbosity level for linear response routines
  !
  integer :: test_case_no = 0   ! OBM, this dummy variable performs various tests
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! lr_control:
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  integer :: charge_response    ! A variable for calculating response charge density
  !
  integer :: itermax            ! number of Lanczos vectors to be calculated
  logical :: ltammd             ! Tarn-Darnkhoff approximation
  logical :: no_hxc             ! If .true. no hartree exchange correlation corrections will be considered.
  logical :: project            ! If .true. projections to read virtual states will be calculated
  !
  !  
  !integer :: ipol               ! sets the polarization direction to be calculated. (Used as a variable if 4)
  !                             ! 1=x 2=y 3=z 4=x,y,z
  !
  !integer :: grid_coarsening    ! Coarses the real space grid (by dividing nr1x nr2x and nr3x by this value)
  !
  ! 
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! lr_post:
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !real(kind=dp) :: broadening                         !Broadening 
  integer :: plot_type                                 ! Dumps rho as: 1=xyzd 2=xsf 3=cube
end module lr_variables
!----------------------------------------------------------------------------
