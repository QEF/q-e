!
! Copyright (C) 2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Written by Aadhityan A, Lorenzo Paulatto, Michele Casula, Tommaso Morresi
!
!--------------------------------------------------------------------------
MODULE pimd_variables
  !---------------------------------------------------------------------------
  !
  ! ... This module contains all variables needed by PIOUD
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! ... "general" variables :
  !
  LOGICAL :: &
       verbose,                  &
       averaged_cov,             &
       fixcm,                    &
       yessecond,                &
       yesquantum,               &
       yesturboq,                &
       yesglobal,                &
       nh,                       &
       restart_pimd     
       
  INTEGER :: &
       ndimMD=0,                   &! dimension of the configuration space
       nbeadMD=0                  ! number of images
                                  
  REAL(DP) :: &
       delt,                     &
       tempMD,                   &
       gammaMD,                    &
       delta_force,              &
       delta_harm               
  
  INTEGER :: &
       nblocks,                  &! number of blocks for averages  
       nstep_block,              &! number of steps for block
       iprint,                   &! printing trajectories frequency
       irun                       ! tipe of run: only 0,3 or 4    
                                   ! irun = 0  - Classical MD with Nosé-Hoover thermostat
                                   ! irun = 3  - Classical/Quantum Langevin dyn. with Ceriotti integrator
                                   ! irun = 4  - Quantum Langevin dyn. with PIOUD integrator
       

  INTEGER,PARAMETER :: verlet=0, ceriotti=3, pioud=4

  
  !
  ! ... "general" real space arrays
       
  INTEGER  :: &
       unit_dot_ep,  &
       unit_dot_ep2, &
       unit_dot_epsr, &
       unit_dot_eplr, &
       unit_dot_ek,   &
       unit_dot_positions, &
       unit_dot_positions_cen, &
       unit_dot_forces,    &
       unit_dot_forces_cen,    &
       unit_dot_stress,    &
       unit_dot_stress_cen,    &
       unit_dot_velocities, &
       unit_dot_velocities_cen, &
       unit_dot_xyz,        &
       unit_dot_localtemp,  &
       unit_dot_sigma,      &
       unit_dot_out
       
  integer,parameter :: nmax=7
  integer,parameter :: mav=2000+nmax
       
  real(dp), parameter :: pi=3.14159265358979d0, kbm1=315775.024d0  ! in K/Ha units
  real(dp), parameter :: eps_cov = 1.d-10   ! cut off of lower eigenvalues for the inversion of the covariance matrix 
  real(dp), parameter :: ha2cm1 = 219474.63137d0

! Variables 

  integer :: natMD=0,nspecies,gMD,iblockMD
  integer :: lfnl,lonl,ifl,nav,ikin,nunitcells
  integer, dimension (:), allocatable :: ipot,indx
  
  real(dp) :: alat
  real(8) :: rcut
  real(8) :: SMD,PSMD,QMD,sint,mtot,sigmavar
  real(8) :: ekinq, ekinqp, epot_centroid,tfakeMD,deltahtilde,epot_old
  real(8), dimension (:), allocatable :: amas,av,avp,anorm,anormp
  real(8), dimension (:), allocatable :: dynmat_eig,tmes_bead,dynmatforce_eig
  real(8), dimension (:), allocatable :: friction_mode,omega_mode,cost1,cost2,gamma_eigen
  real(8), dimension (:), allocatable :: alpha_qmc,alphaqmc_eig
  real(8), dimension (:,:), allocatable :: rcm,vcm,dynmat,cmatrix,dynmat_force0
  real(8), dimension (:,:), allocatable :: el,rtilde,rcentroid,forcedyn
  real(8), dimension (:,:,:), allocatable :: rpos,rpos_init,forceMD,vel,pimp,velocity
  real(8), dimension (:,:,:), allocatable :: rtilde_mode, ptilde,forceMD_old,rpos_old

  REAL(dp), ALLOCATABLE :: &
         pes(:),                   &! the potential enrgy along the path
     !   pos(:,:),                 &! reaction path
     !   grad_pes(:,:),            &! gradients acting on the path
       stress_pes_md(:,:)             !Stress

   
  character(len=3), dimension(:), allocatable :: ion_name   
  ! Myreweight and former TurboRVB variables 

  integer :: ieskin,iflagerr,info,nbead,nion
  integer :: iscramax,rankMD
  real(8) :: dt,cost,temp,ener_true,friction,dynmat_cutoff
  real(8) :: normcorr,scalecov,sigma_true,tmes,friction_qmc
  real(8), dimension(:,:), allocatable :: fk,mass_ion,fbead,kdyn,cov_pimd
  real(8), dimension(:), allocatable :: psip,scalpar,kdyn_eig,cov,cov_old
  logical :: write_cov
  
END MODULE pimd_variables
