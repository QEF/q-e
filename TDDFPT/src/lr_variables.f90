!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
MODULE lr_variables
  !--------------------------------------------------------------------------
  ! ... sets the dimensions of the variables required in the
  ! ... Lanczos linear-response calculation
  !--------------------------------------------------------------------------
  !
  USE kinds,                ONLY : dp
  USE control_flags,        ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  ! Parameters
  !
  INTEGER, PARAMETER:: nbrx = 14   ! max number of beta functions
  INTEGER :: iund0psi   = 20       ! unit for writing/reding of d0psi
  INTEGER :: iunrestart = 22
  INTEGER :: nwordd0psi, nwordrestart, n_ipol
  CHARACTER (len=10), PARAMETER :: code1 = 'turboTDDFT', code2 = 'turboEELS'
  INTEGER :: size_evc
  CHARACTER (len=24) :: bgz_suffix
  !
  LOGICAL :: lr_exx
  REAL(kind=dp) :: scissor
  !
  !------------------------------------------------------------------------!
  ! Variables for EELS
  !------------------------------------------------------------------------!
  !
  LOGICAL :: eels=.false.      ! If .true. then EELS calculation is activated
  REAL(kind=dp) :: q1,q2,q3    ! Components of the q-vector in units of 2*pi/a
                               ! in Cartesian coordinates
  CHARACTER(len=30) :: approximation ! Level of approximation in TDDFPT  
  CHARACTER(LEN=256) :: tmp_dir_lr   ! Name of a temporary directory 
  !
  !------------------------------------------------------------------------!
  !
  REAL(kind=dp), ALLOCATABLE    :: becp_1(:,:), becp1_virt(:,:)
  COMPLEX(kind=dp), ALLOCATABLE :: becp1_c(:,:,:), becp1_c_virt(:,:,:)
  ! the product of the beta-functions with the unperturbed wavefunctions
  !
  COMPLEX(kind=dp), ALLOCATABLE :: &
       evc0(:,:,:),       &    ! the ground state wavefunctions (plane wave, band, k point)
       evc0_virt(:,:,:),  &    ! unoccupied ground state wavefunctions (plane wave, band, k point)
       sevc0(:,:,:),      &    ! S * ground state wavefunctions
       sevc0_virt(:,:,:), &    ! S * virtual ground state wavefunctions
       evc1_old(:,:,:,:), &    ! response wavefunctions in the pw basis (last
                               ! index 1: q' using rotated SBR 2: p')
       evc1(:,:,:,:),     &    !  "    "
       evc1_new(:,:,:,:), &    !  "    "
       sevc1(:,:,:),      &    ! S * "    "
       sevc1_new(:,:,:),  &    ! S * "    "
       d0psi(:,:,:,:),    &    ! for saving the original starting vectors
       d0psi2(:,:,:,:),   &    ! for saving the original starting vectors (without P^+_c)
       revc0(:,:,:),      &    ! ground state wavefunctions in real space
       tg_revc0(:,:,:)         ! ground state wavefunctions in real space
  REAL(kind=dp), ALLOCATABLE ::    &
       rho_1(:,:)              ! response charge density in real space
  COMPLEX(kind=dp), ALLOCATABLE :: &
       rho_1c(:,:)             ! response charge density in real space
  INTEGER :: nbnd_total        ! actual number of bands calculated by PWSCF (virtual+ocuppied)
  INTEGER, ALLOCATABLE :: cube_save(:,:) !used in response charge density mode 1
  !
  COMPLEX(kind=dp), ALLOCATABLE :: F(:,:,:) ! The intensity of transition from valence state (first index)
                                            ! to conduction  state (second index), for each polarization
                                            ! direction (third index).
  !
  COMPLEX(kind=dp), ALLOCATABLE :: R(:,:,:) ! The oscillator strength from valence state (first index)
                                            ! to conduction  state (second index), for each polarization
                                            ! direction (third index).
  !
  ! Lanczos Matrix
  !
  !
  !       | alpha(:,1),  gamma(:,2),                             0|
  !       | beta(:,2) ,  alpha(:,2), gamma(:,3)                   |
  ! T^j = |           ,  beta(:,3) , ...                          |
  !       |                                                       |
  !       |                                       ..., gamma(:,j) |
  !       |0                               beta(:,j) , alpha(:,j) |
  !
  ! Zeta is the \sum_valance (V^T_j * r_i ), where r_i is the density operator acting
  ! on ground state orbitals
  !
  ! zeta.w_T gives the polarizability (w_T is the solution of
  ! (\omega-L)e_1 = w_T , this is handled in a post processing program)
  !
  REAL(kind=dp), ALLOCATABLE :: &  ! (pol, iter)
               alpha_store(:,:),&
               beta_store(:,:), &
               gamma_store(:,:)
  COMPLEX(kind=dp), ALLOCATABLE :: zeta_store(:,:,:)  ! polarization of external field, 
                                                      ! polarization of internal field, 
                                                      ! iteration number.
  !
  ! The currently processed polarization direction and Lanczos iteration
  !
  INTEGER :: LR_iteration, LR_polarization ! polarization directions 1:x 2:y 3:z
  !
  REAL(kind=dp) :: norm0(3)
  !
  LOGICAL :: davidson = .false.           
  !
#if defined(__ENVIRON)
  CHARACTER(len=80) :: assume_isolated = 'none'
#endif
  !
  ! lr_input
  !
  LOGICAL :: restart            ! .true. if the calculation is a restart run
  INTEGER :: restart_step       ! the amount of steps to write a restart file
  INTEGER :: lr_verbosity       ! verbosity level for linear response routines
  INTEGER :: test_case_no = 0   ! OBM, this dummy variable performs various tests
  INTEGER :: lr_io_level = 1    ! Controls disk io
  !
  ! lr_control
  !
  INTEGER :: charge_response    ! variable for calculating response charge density
  INTEGER :: itermax            ! number of Lanczos vectors to be calculated
  INTEGER :: itermax_int        ! interpolated number of lanczos steps for Ritz vectors
  LOGICAL :: ltammd             ! Tamm-Dancoff approximation
  LOGICAL :: d0psi_rs           ! Calculate d0psi in the real space
  LOGICAL :: lshift_d0psi       ! When calculate d0psi in real space, automatically
                                ! shift the position operator r to the center of the molecule
  LOGICAL :: pseudo_hermitian   ! If psedudo-Hermitian algorithm is used
  LOGICAL :: no_hxc             ! If .true. no Hartree and exchange-correlation corrections will be considered.
  LOGICAL :: project            ! If .true. projections to read virtual states will be calculated
  !
  ! lr_post
  !
  INTEGER :: plot_type          ! dumps rho as: 1=xyzd 2=xsf 3=cube
  INTEGER :: sum_rule           ! currently supported sum rules : -2 for alpha(w->0)
  !
END MODULE lr_variables
