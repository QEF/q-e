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
  complex(kind=dp), allocatable :: F(:,:,:) !the intensity of transition from valance state (first index)
                                       ! to conduction  state (second index), for each polarization 
                                       !direction (third index)
  complex(kind=dp), allocatable :: R(:,:,:) !the oscillator strength from valanace state (first index)
                                       ! to conduction  state (second index), for each polarization 
                                       !direction (third index)
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
  integer :: itermax_int        ! interpolated number of lanczos steps for Ritz vectors
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Debugging subroutines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  contains
!----------------------------------------------------------------------------
SUBROUTINE check_vector_gamma (x)
! Checks the inner product for a given vector, and its imaginary and real
! component
! input, evc
! output : screen output
   USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
   use mp,                   only : mp_sum
   use realus,               only : npw_k
  use gvect,                only : gstart
  use io_global,            only : stdout

IMPLICIT NONE
  !input/output
  complex(kind=dp),intent(in)  :: x(:)
  !local
  real(kind=dp) :: temp_gamma
  real(kind=dp), external    :: DDOT

      temp_gamma = 2.D0*DDOT(2*npw_k(1),x(:),1,x(:),1)
       if (gstart==2) temp_gamma = temp_gamma - dble(x(1))*dble(x(1))
#ifdef __PARA
       call mp_sum(temp_gamma, intra_pool_comm)
#endif 
      write(stdout,'("<x> = ",E15.8)') temp_gamma
    END SUBROUTINE check_vector_gamma

!----------------------------------------------------------------------------
SUBROUTINE check_vector_f (x)
! Checks the inner product for a given vector, and its imaginary and real
! component
! input, evc
! output : screen output
   USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
   use mp,                   only : mp_sum
   use realus,               only : npw_k
  use gvect,                only : gstart
  use io_global,            only : stdout

IMPLICIT NONE
  !input/output
  complex(kind=dp),intent(in)  :: x(:)
  !local
  complex(kind=dp) :: temp_f
  complex(kind=dp), external    :: ZDOTC

      temp_f = ZDOTC(npw_k(1),x(:),1,x(:),1)
#ifdef __PARA
       call mp_sum(temp_f, intra_pool_comm)
#endif 
      write(stdout,'("<x> = ",2E15.8,1X)') temp_f
 END SUBROUTINE check_vector_f
!----------------------------------------------------------------------------
SUBROUTINE check_all_bands_gamma (x,sx,nbnd1,nbnd2)
! Checks all bands of given KS states for orthoganilty 
! input, evc and sevc
! output : screen output
   USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
   use mp,                   only : mp_sum
   use realus,               only : npw_k
   use io_global,            only : stdout
   use gvect,                only : gstart
IMPLICIT NONE
  !input/output
  integer,intent(in)           :: nbnd1,nbnd2 !Total number of bands for x and sx
  complex(kind=dp),intent(in)  :: x(:,:), sx(:,:)
  !local
  integer          :: ibnd, jbnd
  real(kind=dp) :: temp_gamma
  real(kind=dp), external    :: DDOT

    do ibnd=1,nbnd1
     do jbnd=ibnd,nbnd2      
       !
      temp_gamma = 2.D0*DDOT(2*npw_k(1),x(:,ibnd),1,sx(:,jbnd),1)
       if (gstart==2) temp_gamma = temp_gamma - dble(x(1,ibnd))*dble(sx(1,jbnd))
#ifdef __PARA
       call mp_sum(temp_gamma, intra_pool_comm)
#endif 
      write(stdout,'("<x,",I02,"|S|x,",I02,"> =",E15.8)') ibnd,jbnd,temp_gamma
     enddo
    enddo
    END SUBROUTINE check_all_bands_gamma
!----------------------------------------------------------------------------
SUBROUTINE check_density_gamma (rx,nbnd)
! Checks the contirbution of a given function transformed into real space
! input, revc
! output : screen output
   USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
   use mp,                   only : mp_sum
   use realus,               only : npw_k
   use wvfct,                only : wg
   use gvect,                only : nrxx
   use io_global,            only : stdout
   use cell_base,                    only : omega

IMPLICIT NONE
  !input/output
  integer,intent(in)           :: nbnd !Total number of bands for x and sx
  complex(kind=dp),intent(in)  :: rx(:,:)
  !local
  integer          :: ibnd
  real(kind=dp) :: temp_gamma,w1,w2

    do ibnd=1,nbnd,2
       w1=wg(ibnd,1)/omega
       !
       if(ibnd<nbnd) then
          w2=wg(ibnd+1,1)/omega
       else
          w2=w1
       endif
       temp_gamma=SUM(w1*DBLE(rx(1:nrxx,ibnd))*DBLE(rx(1:nrxx,ibnd))&
               +w2*aimag(rx(1:nrxx,ibnd))*aimag(rx(1:nrxx,ibnd)))
#ifdef __PARA
       call mp_sum(temp_gamma, intra_pool_comm)
#endif 
      write(stdout,'("Contribution of bands ",I02," and ",I02," to total density",E15.8)') ibnd,ibnd+1,temp_gamma
    enddo
    !
    !
END SUBROUTINE check_density_gamma
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
end module lr_variables
!----------------------------------------------------------------------------
