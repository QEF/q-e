!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
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
  !------------------------------------------------------------------------!
  ! Variables for TD-HF using Dario Rocca's BSE routines.
  !------------------------------------------------------------------------!
  !
  LOGICAL :: lr_exx
  REAL(kind=dp), ALLOCATABLE :: xkk_m(:,:), xkk_p(:,:), weight_kk(:)
  REAL(kind=dp) :: exxdiv
  REAL(kind=dp) :: ecutfock_
  !  REAL(kind=dp), ALLOCATABLE :: revc_int(:,:)
  !  COMPLEX(kind=dp), ALLOCATABLE :: revc_int_c(:,:,:)
  !
  INTEGER, ALLOCATABLE :: &
           igkk_kk_m(:,:),&
           npw_kk_m(:),   &
           ind_kk_m(:,:), &
           igkk_kk_p(:,:),&
           npw_kk_p(:),   &
           ind_kk_p(:,:)
  REAL(kind=dp) :: scissor
  CHARACTER(len=200) :: eig_dir
  !
  !------------------------------------------------------------------------!
  ! Variables for EELS
  !------------------------------------------------------------------------!
  !
  LOGICAL :: eels=.false.      ! If .true. then EELS calculation is activated
  REAL(kind=dp) :: q1,q2,q3    ! Components of the q-vector in units of 2*pi/a
                               ! in Cartesian coordinates
  LOGICAL :: lr_periodic       ! If .true. the perturbation is periodic (q=G)
  CHARACTER(len=30) :: approximation ! Level of approximation in TDDFPT  
  LOGICAL :: clfe              ! Crystal Local Field Effects (CLFE) are used 
                               ! if clfe=.true.
  !LOGICAL :: eps               ! If eps=.true. one computes the dielectric
  !                             ! function (not the inverse) for a fixed q
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
       sevc1(:,:,:,:),    &    ! S * "    "
       sevc1_new(:,:,:,:),&    ! S * "    "
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
  ! variables for diagonalising the coefficient matrix
  !real(kind=dp) :: lanc_norm
  !real(kind=dp), ALLOCATABLE :: eval1(:),eval2(:)
  !real(kind=dp), ALLOCATABLE :: vl(:,:),vr(:,:)
  !
  REAL(kind=dp) :: norm0(3)
  !
  LOGICAL :: davidson = .false.           
  !
#ifdef __ENVIRON
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
  LOGICAL :: ltammd             ! Tarn-Darnkhoff approximation
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
  ! I. Timrov's comment: It is not a proper place for debugging subroutines
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Debugging subroutines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
CONTAINS

SUBROUTINE check_vector_gamma (x)
   !----------------------------------------------------------------------------
   ! Checks the inner product for a given vector, and its imaginary and real component
   ! input: evc
   ! output : screen output
   !
   USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
   USE mp,                   ONLY : mp_sum
   USE klist ,               ONLY : npw_k=>ngk
   USE gvect,                ONLY : gstart
   USE io_global,            ONLY : stdout
   !
   IMPLICIT NONE
   !input/output
   COMPLEX(kind=dp),INTENT(in)  :: x(:)
   !
   ! local variables
   !
   REAL(kind=dp) :: temp_gamma
   REAL(kind=dp), EXTERNAL :: DDOT
   !
   temp_gamma = 2.D0*DDOT(2*npw_k(1),x(:),1,x(:),1)
   !    
   IF (gstart==2) temp_gamma = temp_gamma - dble(x(1))*dble(x(1))
   ! 
#ifdef __MPI
   CALL mp_sum(temp_gamma, intra_bgrp_comm)
#endif
   !    
   WRITE(stdout,'("<x> = ",E15.8)') temp_gamma
   !
   RETURN
   !  
END SUBROUTINE check_vector_gamma

SUBROUTINE check_vector_f (x)
   !-----------------------------------------------------------------------
   !
   ! Checks the inner product for a given vector, and its imaginary and real component
   ! input: evc
   ! output: screen output
   ! 
   USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
   USE mp,                   ONLY : mp_sum
   USE klist ,               ONLY : npw_k=>ngk
   USE gvect,                ONLY : gstart
   USE io_global,            ONLY : stdout
   !
   IMPLICIT NONE
   !input/output
   COMPLEX(kind=dp),INTENT(in)  :: x(:)
   !
   ! local variables
   !
   COMPLEX(kind=dp) :: temp_f
   COMPLEX(kind=dp), EXTERNAL :: ZDOTC
   !
   temp_f = ZDOTC(npw_k(1),x(:),1,x(:),1)
   !
#ifdef __MPI
   CALL mp_sum(temp_f, intra_bgrp_comm)
#endif
   !
   WRITE(stdout,'("<x> = ",2E15.8,1X)') temp_f
   !
   RETURN
   ! 
END SUBROUTINE check_vector_f

SUBROUTINE check_all_bands_gamma (x,sx,nbnd1,nbnd2)
  !----------------------------------------------------------------------
  !
  ! Checks all bands of given KS states for orthoganilty
  ! input: evc and sevc
  ! output : screen output
  ! 
  USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE klist ,               ONLY : npw_k=>ngk
  USE io_global,            ONLY : stdout
  USE gvect,                ONLY : gstart
  !
  IMPLICIT NONE
  !input/output
  INTEGER, INTENT(in) :: nbnd1,nbnd2 !Total number of bands for x and sx
  COMPLEX(kind=dp),INTENT(in) :: x(:,:), sx(:,:)
  !
  ! local variables
  !
  INTEGER :: ibnd, jbnd
  REAL(kind=dp) :: temp_gamma
  REAL(kind=dp), EXTERNAL :: DDOT
  !
  DO ibnd=1,nbnd1
     DO jbnd=ibnd,nbnd2
        !
        temp_gamma = 2.D0*DDOT(2*npw_k(1),x(:,ibnd),1,sx(:,jbnd),1)
        !
        IF (gstart==2) temp_gamma = temp_gamma - dble(x(1,ibnd))*dble(sx(1,jbnd))
        !
#ifdef __MPI
        CALL mp_sum(temp_gamma, intra_bgrp_comm)
#endif
        !
        WRITE(stdout,'("<x,",I02,"|S|x,",I02,"> =",E15.8)') ibnd,jbnd,temp_gamma
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE check_all_bands_gamma

SUBROUTINE check_density_gamma (rx,nbnd)
  !---------------------------------------------------------------------------------
  !
  ! Checks the contirbution of a given function transformed into real space
  ! input: revc
  ! output : screen output
  !
  USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE wvfct,                ONLY : wg
  USE fft_base,             ONLY : dfftp
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : omega
  !
  IMPLICIT NONE
  !input/output
  INTEGER, INTENT(in)          :: nbnd !Total number of bands for x and sx
  COMPLEX(kind=dp),INTENT(in)  :: rx(:,:)
  ! 
  ! local variables
  !
  INTEGER :: ibnd
  REAL(kind=dp) :: temp_gamma,w1,w2
  !
  DO ibnd=1,nbnd,2
     w1 = wg(ibnd,1)/omega
     !
     IF (ibnd<nbnd) THEN
        w2 = wg(ibnd+1,1)/omega
     ELSE
        w2 = w1
     ENDIF
     temp_gamma = sum(w1*dble(rx(1:dfftp%nnr,ibnd))*dble(rx(1:dfftp%nnr,ibnd))&
                + w2*aimag(rx(1:dfftp%nnr,ibnd))*aimag(rx(1:dfftp%nnr,ibnd)))
#ifdef __MPI
     CALL mp_sum(temp_gamma, intra_bgrp_comm)
#endif
     WRITE(stdout,'("Contribution of bands ",I02," and ",I02," to total density",E15.8)') ibnd,ibnd+1,temp_gamma
     !    
  ENDDO
  !
  RETURN
  !
END SUBROUTINE check_density_gamma

END MODULE lr_variables
