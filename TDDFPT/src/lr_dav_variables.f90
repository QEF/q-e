!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
MODULE lr_dav_variables
  !--------------------------------------------------------------------------
  ! ... sets the dimensions of the variables required in the
  ! ... lr_davidson algorithm
  !--------------------------------------------------------------------------
  ! Created by Xiaochuan Ge (Oct.2012)
  USE kinds,                ONLY : dp

  IMPLICIT NONE
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Parameters
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real (dp), PARAMETER :: zero=1.0E-10
  real (dp), PARAMETER :: PI = 3.14159265d0

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Namelist
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER :: num_eign, num_init, num_basis_max,max_iter,&
             p_nbnd_occ,p_nbnd_virt ! number of occ and virt bands for projection
  REAL(kind=dp) :: residue_conv_thr, reference, close_pre, broadening,start,&
                   finish,step,turn2planb, vccouple_shift
  logical :: precondition,dav_debug, single_pole,&
             sort_contr,diag_of_h,print_spectrum,if_check_orth,&
             if_random_init,if_check_her,poor_of_ram,poor_of_ram2,&
             conv_assistant,if_dft_spectrum,lplot_drho

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Global Variables
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer , allocatable :: &
        vc_couple(:,:),&  ! for calculating the optimal initial vector
        eign_value_order(:), energy_dif_order(:)

  LOGICAL :: dav_conv
  LOGICAL :: ploted(3)=.false., done_calc_R=.false.
  real(dp) :: max_res
  logical, allocatable :: kill_left(:),kill_right(:) 
  INTEGER :: num_basis_old, &   ! num of basis that has been applied Liouvillian
             num_basis,toadd,       &   ! total number of basis
             lwork,info,&         ! for calling lapack routine 
             dav_iter,num_basis_tot  

  COMPLEX (DP), ALLOCATABLE ::  &
        vec_b(:,:,:,:), svec_b(:,:,:,:),   & ! vectors for the basis of Liouvillian
        left_full(:,:,:,:), &
        right_full(:,:,:,:),&  ! left and right appro. eigen vector of D*C in full space
        left_res(:,:,:,:), &  
        right_res(:,:,:,:), &  ! left and right residual of the appr. eigen-v
        vecwork(:,:,:), & ! Working space for supervectors
        left2(:),right2(:), & ! square of the module of left_res and right_res
        M_C(:,:),M_D(:,:),M(:,:),& ! Matrix represented in the sub space
        ground_state(:,:,:),&
        C_vec_b(:,:,:,:),       &  ! D and C applied to vec_b 
        D_vec_b(:,:,:,:),       &  ! D and C applied to vec_b. see documentation for explaination of D&C
        Fx(:,:),Fy(:,:)
        
  real(dp) :: FxR(3), FyR(3)

  REAL (DP), ALLOCATABLE :: &
        work(:) ,&            ! global work space for linear algorithm routine
        left_M(:,:),        &
        right_M(:,:),       & ! left and right appro. eigen vectors of
        eign_value(:,:),&
        M_shadow_avatar(:,:),&
        inner_matrix(:,:),&
        tr_energy(:),&
        energy_dif(:),contribute(:,:),&
        chi_dav(:,:),total_chi(:),norm_F(:),&
        omegal(:),omegar(:)
  
end module lr_dav_variables
