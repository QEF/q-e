!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! Copyright (C) 2009 Brian Kolb, Timo Thonhauser - Wake Forest University
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------

#define FFTGRADIENT
!#undef FFTGRADIENT

MODULE vdW_DF
 
  !! This module calculates the non-local correlation contribution to the energy
  !! and potential. This method is based on the method of Guillermo Roman-Perez 
  !! and Jose M. Soler described in:
  !!
  !!    G. Roman-Perez and J. M. Soler, PRL 103, 096102 (2009)
  !!
  !! henceforth referred to as SOLER. That method is a new implementation
  !! of the method found in:
  !!
  !!    M. Dion, H. Rydberg, E. Schroeder, D. C. Langreth, and
  !!    B. I. Lundqvist, Phys. Rev. Lett. 92, 246401 (2004).
  !!
  !! henceforth referred to as DION. Further information about the
  !! functional and its corresponding potential can be found in:
  !!
  !!    T. Thonhauser, V.R. Cooper, S. Li, A. Puzder, P. Hyldgaard,
  !!    and D.C. Langreth, Phys. Rev. B 76, 125112 (2007).
  !!
  !! A review article that shows many of the applications vdW-DF has been
  !! applied to so far can be found at:
  !!
  !!    D. C. Langreth et al., J. Phys.: Condens. Matter 21, 084203 (2009).
  !!
  !! There are a number of subroutines in this file. All are used only
  !! by other subroutines here except for the xc_vdW_DF subroutine
  !! which is the driver routine for the vdW-DF calculations and is called
  !! from v_of_rho.  This routine handles setting up the parallel run (if
  !! any) and carries out the calls necessary to calculate the non-local
  !! correlation contributions to the energy and potential.
  
  USE kinds,             ONLY : dp
  USE constants,         ONLY : pi, e2
  USE kernel_table,      ONLY : q_mesh, Nr_points, Nqs, r_max
  USE mp,                ONLY : mp_bcast, mp_sum, mp_barrier
  USE mp_pools,          ONLY : me_pool, nproc_pool, intra_pool_comm, root_pool
  USE io_global,         ONLY : ionode, stdout
  USE fft_base,          ONLY : dfftp
  USE fft_interfaces,    ONLY : fwfft, invfft 
  USE control_flags,     ONLY : iverbosity, gamma_only
 
  IMPLICIT NONE
  
  REAL(DP), PARAMETER :: epsr =1.d-12  ! a small number to cut off points with negative or undefined densities 
  integer :: vdw_type = 1

  private  
  public :: xc_vdW_DF, stress_vdW_DF, interpolate_kernel, vdw_type, numerical_gradient, initialize_spline_interpolation

CONTAINS

!! #################################################################################################
!!                                       |             |
!!                                       |  XC_VDW_DF  |
!!                                       |_____________|

  SUBROUTINE xc_vdW_DF(rho_valence, rho_core, nspin, etxc, vtxc, v)
    
    !! Modules to include
    !! -------------------------------------------------------------------------
    
    use gvect,           ONLY : ngm, nl, g, nlm
    USE fft_base,        ONLY : dfftp
    USE cell_base,       ONLY : omega, tpiba
    !! -------------------------------------------------------------------------
    
    !! Local variables
    !! ----------------------------------------------------------------------------------
    !                                               _
    real(dp), intent(IN) :: rho_valence(:,:)       !
    real(dp), intent(IN) :: rho_core(:)            !  PWSCF input variables 
    INTEGER,  INTENT(IN) :: nspin                  !
    real(dp), intent(inout) :: etxc, vtxc, v(:,:)  !_  
   
    
    integer :: i_grid, theta_i, i_proc, I       !! Indexing variables over grid points,
    !                                           !! theta functions, and processors, and a
    !                                           !! generic index.
                                                
    real(dp) :: grid_cell_volume                !! The volume of the unit cell per G-grid point
                                                
    real(dp), allocatable ::  q0(:)             !! The saturated value of q (equations 11 and 12 of DION)
    !                                           !! This saturation is that of equation 5 in 
    !                                           !! SOLER
    
    real(dp), allocatable :: gradient_rho(:,:)  !! The gradient of the charge density. The 
    !                                           !! format is as follows:
    !                                           !! gradient_rho(grid_point, cartesian_component)
                                                
    real(dp), allocatable :: potential(:)       !! The vdW contribution to the potential
                                                
    real(dp), allocatable :: dq0_drho(:)        !! The derivative of the saturated q0 
    !                                           !! (equation 5 of SOLER) with respect
    !                                           !! to the charge density (sort of.  see
    !                                           !! get_q0_on_grid subroutine below.) 
                                                
    real(dp), allocatable :: dq0_dgradrho(:)    !! The derivative of the saturated q0
    !                                           !! (equation 5 of SOLER) with respect
    !                                           !! to the gradient of the charge density
    !                                           !! (again, see get_q0_on_grid subroutine)
                                                
    complex(dp), allocatable :: thetas(:,:)     !! These are the functions of equation 8 of 
    !                                           !! SOLER.  They will be forward Fourier transformed
    !                                           !! in place to get theta(k) and worked on in
    !                                           !! place to get the u_alpha(r) of equation 11
    !                                           !! in SOLER.  They are formatted as follows:
    !                                           !1 thetas(grid_point, theta_i)
                                                
    real(dp) :: Ec_nl                           !! The non-local vdW contribution to the energy
                                                
    real(dp), allocatable :: total_rho(:)       !! This is the sum of the valence and core 
    !                                           !! charge.  This just holds the piece assigned
    !                                           !! to this processor.
#ifndef FFTGRADIENT
    integer, parameter :: Nneighbors = 4        !! How many neighbors on each side
    !                                           !! to include in numerical derivatives.
    !                                           !! Can be from 1 to 6
                                                
    real(dp), allocatable :: full_rho(:)        !! This is the whole charge density.  It
    !                                           !! is the sum of valence and core density
    !                                           !! over the entire simulation cell.  Each
    !                                           !! processor has a copy of this to do the
    !                                           !! numerical gradients.
    integer, ave :: msy_start_z, my_end_z       !! Starting and ending z-slabs for this processor


    integer, allocatable, save :: procs_Npoints(:) !! The number of grid points assigned to each proc
    integer, allocatable, save :: procs_start(:)   !! The first assigned index into the charge-density array for each proc
    integer, allocatable, save :: procs_end(:)     !! The last assigned index into the charge density array for each proc

#endif                                                
                                                
    logical, save :: first_iteration = .true.      !! Whether this is the first time this
    !                                              !! routine has been called.
    
    !! ---------------------------------------------------------------------------------------------
    !!   Begin calculations
  
    !! Check to make sure we aren't trying to do a spin-polarized run.
    !! Gamma point calculations can be done using the special {gamma} features
    !! stress tensor calcultion and cell relaxation runs are also possible.
    !! --------------------------------------------------------------------------------------------------------

    call errore('xc_vdW_DF','vdW functional not implemented for spin polarized runs', size(rho_valence,2)-1)
    if (nspin>2) call errore('xc_vdW_DF','vdW functional not implemented for nspin > 2', nspin)

    !! --------------------------------------------------------------------------------------------------------

    if (first_iteration) then
       
#ifndef FFTGRADIENT
    !! Here we set up the calculations on the first iteration.  If this is a parallel run, each
    !! processor figures out which element in the charge-density array it should start and stop on.
    !! PWSCF splits the cell up into slabs in the z-direction to distribute over processors.
    !! Thus, each processor figures out what z-planes its region corresponds to.  That is important
    !! for the get_3d_indices and get_potential subroutines below.
    !! --------------------------------------------------------------------------------------------------

       
       allocate( procs_Npoints(0:nproc_pool-1), procs_start(0:nproc_pool-1), procs_end(0:nproc_pool-1) )
       
       procs_Npoints(me_pool) = dfftp%nnr
       procs_start(0) = 1
       
       ! All processors communicate how many points they have been assigned.  Each processor
       ! then calculates for itself what the starting and ending indices should be for every
       ! other processor.
       !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       do i_proc = 0, nproc_pool-1
          
          call mp_bcast(procs_Npoints(i_proc), i_proc, intra_pool_comm)
          call mp_barrier(intra_pool_comm)
          
          procs_end(i_proc) = procs_start(i_proc) + procs_Npoints(i_proc) - 1
          
          if (i_proc .ne. nproc_pool-1) then
             
             procs_start(i_proc+1) = procs_end(i_proc)+1
             
          end if
          
       end do
       
       !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       ! Each processor finds the starting and ending z-planes assined to them.  Since
       ! PWSCF splits the cell into slabs in the z-direction, the beginning (ending)
       ! z slabs can be found by dividing the starting (ending) index into the charge density
       ! array by the number of points in a slab of thickness 1.  We add 1 to the starting
       ! z plane because of the integer division and the fact that arrays in Fortran start at 1.
       ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       my_start_z = procs_start(me_pool)/(dfftp%nr1x*dfftp%nr2x)+1
       my_end_z = procs_end(me_pool)/(dfftp%nr1x*dfftp%nr2x)

       !write(*,'(A,3I5)') "Parall en [proc, my_start_z, my_end_z]", me_pool, my_start_z, my_end_z
       ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#endif
       first_iteration = .false.
       
       !! Here we output some of the parameters being used in the run.  This is important because
       !! these parameters are read from the vdW_kernel_table file.  The user should ensure that
       !! these are the parameters they were intending to use on each run.
       
       !! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       
       if (ionode ) then
          
          WRITE( stdout, '(//,5x,"************************************************************************")')
          WRITE( stdout, '(5x,"*")')
          WRITE( stdout, '(5x,"* You are using vdW-DF for exchange-correlation in this calculation.")')
          WRITE( stdout, '(5x,"* Please cite the following three references that made this development")')
          WRITE( stdout, '(5x,"* possible:")')
          WRITE( stdout, '(5x,"*")')
          WRITE( stdout, '(5x,"*   M. Dion, H. Rydberg, E. Schroder, D. C. Langreth, and")')
          WRITE( stdout, '(5x,"*   B. I. Lundqvist, Phys. Rev. Lett. 92, 246401 (2004).")')
          WRITE( stdout, '(5x,"*")')
          WRITE( stdout, '(5x,"*   T. Thonhauser, V. R. Cooper, S. Li, A. Puzder, P. Hyldgaard, and")')
          WRITE( stdout, '(5x,"*   D. C. Langreth, Phys. Rev. B 76, 125112 (2007).")')
          WRITE( stdout, '(5x,"*")')
          WRITE( stdout, '(5x,"*   G. Roman-Perez and J. M. Soler, Phys. Rev. Lett. 103, 096102 (2009).")')
          WRITE( stdout, '(5x,"*")')
          WRITE( stdout, '(5x,"*")')
          WRITE( stdout, '(5x,"* If you are using vdW-DF2, please also cite:")')
          WRITE( stdout, '(5x,"*")')
          WRITE( stdout, '(5x,"*   K. Lee, E. D. Murray, L. Kong, B. I. Lundqvist, and")')
          WRITE( stdout, '(5x,"*   D. C. Langreth, Phys. Rev. B 82, 081101(R) (2010).")')
          WRITE( stdout, '(5x,"*")')
          WRITE( stdout, '(5x,"*")')
          WRITE( stdout, '(5x,"* If you are calculating the stress with vdW-DF, please also cite:")')
          WRITE( stdout, '(5x,"*")')
          WRITE( stdout, '(5x,"*   R. Sabatini, E. Kucukbenli, B. Kolb, T. Thonhauser, and")')
          WRITE( stdout, '(5x,"*   S. de Gironcoli, J. Phys.: Condens. Matter 24, 424209 (2012).")')
          WRITE( stdout, '(5x,"*")')
          WRITE( stdout, '(5x,"************************************************************************",/)')
          WRITE( stdout, '(5x,"Carrying out vdW-DF run using the following parameters:",/)')
          WRITE( stdout, '(5X,A,I6,A,I6,A,F8.3)' ) "Nqs =  ", Nqs, "    Nr_points =  ", Nr_points, "   r_max =  ", r_max
          WRITE( stdout, '(5X,"q_mesh =")')
          WRITE( stdout, '(10X,4F15.8)' ) (q_mesh(I), I=1, Nqs)
                 
#ifdef FFTGRADIENT
          WRITE( stdout, '(/,5x,"Gradients computed in reciprocal space",/)')
#else
          WRITE( stdout, '(/,5x,"Gradients computed in real space",/)')
#endif
          
          WRITE( stdout, '(5x,"************************************************************************",//)')

       end if
       
       !! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       
    end if

    !! --------------------------------------------------------------------------------------------------
    
    !! Allocate arrays.  nnr is a PWSCF variable that holds the number of points assigned to 
    !! a given processor.  
    !! ---------------------------------------------------------------------------------------

    allocate( q0(dfftp%nnr) )
    allocate( gradient_rho(dfftp%nnr, 3) )
    allocate( dq0_drho(dfftp%nnr), dq0_dgradrho(dfftp%nnr) )
    allocate( total_rho(dfftp%nnr) )
    
    !! ---------------------------------------------------------------------------------------
    !! Add together the valence and core charge densities to get the total charge density    
    total_rho = rho_valence(:,1) + rho_core(:)
    if (nspin == 2) then
      total_rho = rho_valence(:,1) + rho_valence(:,2) + rho_core(:)
    else
      total_rho = rho_valence(:,1) + rho_core(:)
    endif

#ifdef FFTGRADIENT
    !! -------------------------------------------------------------------------
    !! Here we calculate the gradient in reciprocal space using FFT
    !! -------------------------------------------------------------------------
    call numerical_gradient(total_rho,gradient_rho)

#else
    !! -------------------------------------------------------------------------
    !! Here we calculate the gradient numerically in real space
    !! The neighbors variable is set above and gives the number of points in 
    !! each direction to consider when taking the numerical derivatives.
    !! -------------------------------------------------------------------------
    !! If there is only 1 processor the needed information is held by the 
    !! total_rho array, otherwise we need to allocate the full_rho array that 
    !! will be deallocated the call since it is no longer needed.  
    !!
    !! The full_rho array holds the charge density at every point in the 
    !! simulation cell. Each processor needs this because the numerical
    !! gradients require knowledge of the charge density on points outside the
    !! slab one has been given.  We don't allocate this in the case of using a single 
    !! processor since total_rho would already hold this information.
    !! nr1x, nr2x, and nr3x are PWSCF variables that hold the TOTAL number of 
    !! divisions along each lattice vector. Thus, their product is the total 
    !! number of points in the cell (not just those assigned to a particular 
    !! processor).
    !! ------------------------------------------------------------------------
    
    if (nproc_pool > 1) then
       
       allocate( full_rho(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x) )

       full_rho(procs_start(me_pool):procs_end(me_pool)) = total_rho
       
       ! All the processors broadcast their piece of the charge density to fill in the full_rho
       ! arrays of all processors
       ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       do i_proc = 0, nproc_pool - 1
          
          call mp_barrier(intra_pool_comm) 
          call mp_bcast(full_rho(procs_start(i_proc):procs_end(i_proc)), i_proc, intra_pool_comm)
                   
       end do
       
       ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       
       ! Here we calculate the gradient numerically in real space

       call numerical_gradient(full_rho, Nneighbors, gradient_rho, my_start_z, my_end_z)
       deallocate(full_rho)
       
    else 
       
       ! Here we calculate the gradient numerically in real space
       call numerical_gradient(total_rho, Nneighbors, gradient_rho, my_start_z, my_end_z)
       
    end if

#endif
    !! -------------------------------------------------------------------------------------------------------------

    !! Find the value of q0 for all assigned grid points.  q is defined in equations
    !! 11 and 12 of DION and q0 is the saturated version of q defined in equation
    !! 5 of SOLER.  This routine also returns the derivatives of the q0s with respect
    !! to the charge-density and the gradient of the charge-density.  These are needed
    !! for the potential calculated below.
    !! ---------------------------------------------------------------------------------

    CALL get_q0_on_grid(total_rho, gradient_rho, q0, dq0_drho, dq0_dgradrho)

    !! ---------------------------------------------------------------------------------
    
    !! Here we allocate and calculate the theta functions of SOLER equation 8.  These are defined as 
    !! rho * P_i(q0(rho, gradient_rho)) where P_i is a polynomial that interpolates a Kroneker delta 
    !! function at the point q_i (taken from the q_mesh) and q0 is the saturated version of q.  q is
    !! defined in equations 11 and 12 of DION and the saturation proceedure is defined in equation 5 
    !! of SOLER.  This is the biggest memory consumer in the method since the thetas array is 
    !! (total # of FFT points)*Nqs complex numbers.  In a parallel run, each processor will hold the 
    !! values of all the theta functions on just the points assigned to it.
    !! --------------------------------------------------------------------------------------------------
    !! thetas are stored in reciprocal space as  theta_i(k) because this is the way they are used later
    !! for the convolution (equation 8 of SOLER).  The ffts used here are timed.
    !! --------------------------------------------------------------------------------------------------

    allocate( thetas(dfftp%nnr, Nqs) )
    CALL get_thetas_on_grid(total_rho, q0, thetas)
    !! ---------------------------------------------------------------------------------------------

    !! Carry out the integration in equation 8 of SOLER.  This also turns the thetas array into the 
    !! precursor to the u_i(k) array which is inverse fourier transformed to get the u_i(r) functions
    !! of SOLER equation 11.  Add the energy we find to the output variable etxc.  This process is timed.
    !! --------------------------------------------------------------------------------------------------

    call start_clock( 'vdW_energy')

    call vdW_energy(thetas, Ec_nl)
    
    etxc = etxc + Ec_nl

    call stop_clock( 'vdW_energy')

    !! --------------------------------------------------------------------------------------------------

    !! If verbosity is set to high we output the total non-local correlation energy found
    !! ---------------------------------------------------------------------------------------

    if (iverbosity > 0) then

       call mp_sum(Ec_nl,intra_pool_comm)
    
       if (ionode) write(*,'(/ / A /)') "     ----------------------------------------------------------------"

       if (ionode) write(*,'(A, F22.15 /)') "     Non-local correlation energy =         ", Ec_nl

       if (ionode) write(*,'(A /)') "     ----------------------------------------------------------------"
       
    end if

    !! ----------------------------------------------------------------------------------------

    !! Inverse Fourier transform the u_i(k) to get the u_i(r) of SOLER equation 11.  These FFTs
    !! are also timed and added to the timing of the forward FFTs done earlier.
    !!---------------------------------------------------------------------------------------

    call start_clock( 'vdW_ffts')
    
    do theta_i = 1, Nqs
       
       !call cft3(thetas(:,theta_i), dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, 1)
       CALL invfft('Dense', thetas(:,theta_i), dfftp) 
    end do

    call stop_clock( 'vdW_ffts')

    !! -------------------------------------------------------------------------

    !! Here we allocate the array to hold the potential. This is calculated via 
    !! equation 10 of SOLER, using the u_i(r) calculated from quations 11 and 
    !! 12 of SOLER.  Each processor allocates the array to be the size of the 
    !! full grid  because, as can be seen in SOLER equation 10, processors need
    !! to access grid points outside their allocated regions.
    !! This process is timed.  The timer is stopped below after the v output 
    !! variable has been updated with the non-local corelation potential.  
    !! That is, the timer includes the communication time necessary in a 
    !! parallel run.
    !! -------------------------------------------------------------------------

#ifdef FFTGRADIENT

    call start_clock( 'vdW_v' )

    allocate( potential(dfftp%nnr) )
    
    call get_potential(q0, dq0_drho, dq0_dgradrho, gradient_rho, thetas, potential)
    
    !! -------------------------------------------------------------------------

    
    v(:,1) = v(:,1) + e2*potential(:)
    if (nspin == 2) v(:,2) = v(:,2) + e2*potential(:)

    call stop_clock( 'vdW_v' )

    !! -----------------------------------------------------------------------
    
    !! The integral of rho(r)*potential(r) for the vtxc output variable
    !! --------------------------------------------------------------------

    grid_cell_volume = omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)  
 
    do i_grid = 1, dfftp%nnr
       
       vtxc = vtxc + e2*grid_cell_volume*rho_valence(i_grid,1)*potential(i_grid)
    
    end do
    if (nspin==2) then
       do i_grid = 1, dfftp%nnr
       
          vtxc = vtxc + e2*grid_cell_volume*rho_valence(i_grid,2)*potential(i_grid)
    
       end do
    endif 
       

    deallocate(potential)  

#else

    call start_clock( 'vdW_v' )

    allocate( potential(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x) )
    
    call get_potential(q0, dq0_drho, dq0_dgradrho, Nneighbors, gradient_rho, thetas, potential, my_start_z, my_end_z)
    
    !! -------------------------------------------------------------------------


    !! Reduction process to sum all the potentials of all the processors.  
    !! ----------------------------------------------------------------------

    !   call mp_barrier( intra_pool_comm )
    call mp_sum(potential, intra_pool_comm)

    !! ----------------------------------------------------------------------

    !! Here, the potential is rebroadcast. Since each processor has part of the output v array it is easier if
    !! each processor adds only its assigned points to the v array.  After this step, however, all
    !! processors hold the vdW potential over the entire grid.
    !! ------------------------------------------------------------------------------------------------------

    !    call mp_barrier( intra_pool_comm )
    call mp_bcast(potential, root_pool, intra_pool_comm)

    !! ------------------------------------------------------------------------------------------------------
    
    !! Each processor adds its piece of the potential to the output v array.    
    !! Stop the timer for the potential.
    !! -----------------------------------------------------------------------
    
    v(:,1) = v(:,1) + e2*potential(procs_start(me_pool):procs_end(me_pool))
    
    call stop_clock( 'vdW_v' )

    !! -----------------------------------------------------------------------
    
    !! The integral of rho(r)*potential(r) for the vtxc output variable
    !! --------------------------------------------------------------------

    grid_cell_volume = omega/(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x)  
 
    do i_grid = 1, dfftp%nnr
       
       vtxc = vtxc + e2*grid_cell_volume * rho_valence(i_grid,1)*potential(procs_start(me_pool)+i_grid-1)
    
    end do
    if (nspin==2) then
       do i_grid = 1, dfftp%nnr
       
          vtxc = vtxc + e2*grid_cell_volume * rho_valence(i_grid,2)*potential(procs_start(me_pool)+i_grid-1)
    
       end do
    endif
    deallocate(potential)  

#endif
    !! ----------------------------------------------------------------------
    
    !! Deallocate all arrays.
    deallocate(q0, gradient_rho, dq0_drho, dq0_dgradrho, total_rho, thetas)  
    
  END SUBROUTINE xc_vdW_DF

!! #################################################################################################
!!                                       |                 |
!!                                       |  STRESS_VDW_DF  |
!!                                       |_________________|

  SUBROUTINE stress_vdW_DF(rho_valence, rho_core, nspin, sigma)

      USE fft_base,        ONLY : dfftp
      use gvect,           ONLY : ngm, nl, g, nlm
      USE cell_base,       ONLY : tpiba

      implicit none

      real(dp), intent(IN) :: rho_valence(:,:)           !
      real(dp), intent(IN) :: rho_core(:)                ! Input variables 
      INTEGER,  INTENT(IN) :: nspin                      !
      real(dp), intent(inout) :: sigma(3,3)              !  

      real(dp), allocatable :: gradient_rho(:,:)         !
      real(dp), allocatable :: total_rho(:)              ! Rho values

      real(dp), allocatable :: q0(:)                     !
      real(dp), allocatable :: dq0_drho(:)               ! Q-values
      real(dp), allocatable :: dq0_dgradrho(:)           !

      complex(dp), allocatable :: thetas(:,:)            ! Thetas
#ifndef FFTGRADIENT
      real(dp), allocatable :: full_rho(:)               ! additional Rho values on the full grid

      integer, save :: my_start_z, my_end_z              ! 
      integer, allocatable, save :: procs_Npoints(:)     ! 
      integer, allocatable, save :: procs_start(:)       !
      integer, allocatable, save :: procs_end(:)         !

      logical,  save :: first_stress_iteration = .true.  !
      integer  :: Nneighbors = 4
#endif
      integer :: i_proc, theta_i, l, m

      real(dp)  :: sigma_grad(3,3)
      real(dp)  :: sigma_ker(3,3)

      !! ---------------------------------------------------------------------------------------------
      !!   Tests
      !! --------------------------------------------------------------------------------------------------------

      !call errore('xc_vdW_DF','vdW functional not implemented for spin polarized runs', size(rho_valence,2)-1)
      if (nspin>2) call errore('xc_vdW_DF','vdW functional not implemented for nspin > 2', nspin)

      !IF ( gamma_only) CALL errore ('xc_vdW_DF', &
      !  & 'vdW functional not implemented for gamma point calculations. &
      !  & Use kpoints automatic and specify the gamma point explicitly', 2)
      sigma(:,:) = 0.0_DP
      sigma_grad(:,:) = 0.0_DP
      sigma_ker(:,:) = 0.0_DP

#ifndef FFTGRADIENT
      !! ---------------------------------------------------------------------------------------------
      !!   Parallel setup
      !! ---------------------------------------------------------------------------
      if (first_stress_iteration) then
         
         allocate( procs_Npoints(0:nproc_pool-1), procs_start(0:nproc_pool-1), procs_end(0:nproc_pool-1) )
         
         procs_Npoints(me_pool) = dfftp%nnr
         procs_start(0) = 1
         
         do i_proc = 0, nproc_pool-1
            
            call mp_bcast(procs_Npoints(i_proc), i_proc, intra_pool_comm)
            call mp_barrier(intra_pool_comm)
            
            procs_end(i_proc) = procs_start(i_proc) + procs_Npoints(i_proc) - 1
            
            if (i_proc .ne. nproc_pool-1) then
               procs_start(i_proc+1) = procs_end(i_proc)+1
            end if
            
         end do
         
         my_start_z = procs_start(me_pool)/(dfftp%nr1x*dfftp%nr2x)+1
         my_end_z = procs_end(me_pool)/(dfftp%nr1x*dfftp%nr2x)

         !write(*,'(A,3I5)') "Parall stress [proc, my_start_z, my_end_z]", me_pool, my_start_z, my_end_z

         first_stress_iteration = .false.

      end if
#endif
      !! ---------------------------------------------------------------------------------------
      !! Allocations
      !! ---------------------------------------------------------------------------------------

      allocate( gradient_rho(dfftp%nnr, 3) )
      allocate( total_rho(dfftp%nnr) )
      allocate( q0(dfftp%nnr) )
      allocate( dq0_drho(dfftp%nnr), dq0_dgradrho(dfftp%nnr) )
      allocate( thetas(dfftp%nnr, Nqs) )
 
      !! ---------------------------------------------------------------------------------------
      !! Charge
      !! ---------------------------------------------------------------------------------------

      if (nspin == 2) then
        total_rho = rho_valence(:,1) + rho_valence(:,2) + rho_core(:)
      else
        total_rho = rho_valence(:,1) + rho_core(:)
      endif

#ifdef FFTGRADIENT
      !! -------------------------------------------------------------------------
      !! Here we calculate the gradient in reciprocal space using FFT
      !! -------------------------------------------------------------------------
      call numerical_gradient(total_rho,gradient_rho)
#else
      !! ---------------------------------------------------------------------------------------
      !! Here we calculate the gradient in Real space 
      !! ---------------------------------------------------------------------------------------

      if (nproc_pool > 1) then
         
         allocate( full_rho(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x) )

         full_rho(procs_start(me_pool):procs_end(me_pool)) = total_rho
         
         do i_proc = 0, nproc_pool - 1
            
            call mp_barrier(intra_pool_comm) 
            call mp_bcast(full_rho(procs_start(i_proc):procs_end(i_proc)), i_proc, intra_pool_comm)
                     
         end do
         
         call numerical_gradient(full_rho, Nneighbors, gradient_rho, my_start_z, my_end_z)

         deallocate(full_rho)
         
      else 
         
         call numerical_gradient(total_rho, Nneighbors, gradient_rho, my_start_z, my_end_z)
         
      end if
#endif
      !! -------------------------------------------------------------------------------------------------------------
      !! Get q0.
      !! ---------------------------------------------------------------------------------

      CALL get_q0_on_grid(total_rho, gradient_rho, q0, dq0_drho, dq0_dgradrho)

      !! ---------------------------------------------------------------------------------
      !! Get thetas in reciprocal space.
      !! ---------------------------------------------------------------------------------

      CALL get_thetas_on_grid(total_rho, q0, thetas)

      !! ---------------------------------------------------------------------------------------
      !! Stress
      !! ---------------------------------------------------------------------------------------
      CALL stress_vdW_DF_gradient(total_rho, gradient_rho, q0, dq0_drho, &
                                  dq0_dgradrho, thetas, sigma_grad)
      
      CALL stress_vdW_DF_kernel(total_rho, q0, thetas, sigma_ker)
      
      sigma = - (sigma_grad + sigma_ker) 

      do l = 1, 3
         do m = 1, l - 1
            sigma (m, l) = sigma (l, m)
         enddo
      enddo

      deallocate( gradient_rho, total_rho, q0, dq0_drho, dq0_dgradrho, thetas )
 
   END SUBROUTINE stress_vdW_DF

   !! ###############################################################################################################
   !!                                       |                          |
   !!                                       |  STRESS_VDW_DF_GRADIENT  |
   !!                                       |                          |

   SUBROUTINE stress_vdW_DF_gradient (total_rho, gradient_rho, q0, dq0_drho, &
                                      dq0_dgradrho, thetas, sigma)

      !!-----------------------------------------------------------------------------------
      !! Modules to include
      !! ----------------------------------------------------------------------------------
      use gvect,                 ONLY : ngm, nl, g, nlm, nl, gg, igtongl, &
                                        gl, ngl, gstart
      USE fft_base,              ONLY : dfftp
      USE cell_base,             ONLY : omega, tpiba, alat, at, tpiba2

      !! ----------------------------------------------------------------------------------

      implicit none

      real(dp), intent(IN) :: total_rho(:)               !
      real(dp), intent(IN) :: gradient_rho(:, :)         ! Input variables
      real(dp), intent(inout) :: sigma(:,:)              !  
      real(dp), intent(IN) :: q0(:)                      !
      real(dp), intent(IN) :: dq0_drho(:)                ! 
      real(dp), intent(IN) :: dq0_dgradrho(:)            !
      complex(dp), intent(IN) :: thetas(:,:)             !

      complex(dp), allocatable :: u_vdW(:,:)             !

      real(dp), allocatable    :: d2y_dx2(:,:)           !
      real(dp) :: y(Nqs), dP_dq0, P, a, b, c, d, e, f    ! Interpolation
      real(dp) :: dq                                     !

      integer  :: q_low, q_hi, q, q1_i, q2_i , g_i       ! Loop and q-points

      integer  :: l, m
      real(dp) :: prefactor                              ! Final summation of sigma

      integer  :: i_proc, theta_i, i_grid, q_i, &        !
                  ix, iy, iz                             ! Iterators
      
      character(LEN=1) :: intvar

      !real(dp)       :: at_inverse(3,3)

      allocate( d2y_dx2(Nqs, Nqs) ) 
      allocate( u_vdW(dfftp%nnr, Nqs) )

      sigma(:,:) = 0.0_DP
      prefactor = 0.0_DP
      
      !! --------------------------------------------------------------------------------------------------
      !! Get u in k-space.
      !! ---------------------------------------------------------------------------------------------------

      call thetas_to_uk(thetas, u_vdW)

      !! --------------------------------------------------------------------------------------------------
      !! Get u in real space.
      !! ---------------------------------------------------------------------------------------------------

      call start_clock( 'vdW_ffts')

      do theta_i = 1, Nqs
         CALL invfft('Dense', u_vdW(:,theta_i), dfftp) 
      end do

      call stop_clock( 'vdW_ffts')

      !! --------------------------------------------------------------------------------------------------
      !! Get the second derivatives for interpolating the P_i
      !! ---------------------------------------------------------------------------------------------------

      call initialize_spline_interpolation(q_mesh, d2y_dx2(:,:))

      !! ---------------------------------------------------------------------------------------------

      i_grid = 0

      !! ----------------------------------------------------------------------------------------------------	
      !! Do the real space integration to obtain the stress component
      !! ----------------------------------------------------------------------------------------------------


      do i_grid = 1, dfftp%nnr

                  q_low = 1
                  q_hi = Nqs 

                  !
                  ! Figure out which bin our value of q0 is in the q_mesh
                  !
                  do while ( (q_hi - q_low) > 1)

                      q = int((q_hi + q_low)/2)

                      if (q_mesh(q) > q0(i_grid)) then
                          q_hi = q
                      else 
                          q_low = q
                      end if

                  end do

                  if (q_hi == q_low) call errore('stress_vdW_gradient','qhi == qlow',1)

                  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                  dq = q_mesh(q_hi) - q_mesh(q_low)

                  a = (q_mesh(q_hi) - q0(i_grid))/dq
                  b = (q0(i_grid) - q_mesh(q_low))/dq
                  c = (a**3 - a)*dq**2/6.0D0
                  d = (b**3 - b)*dq**2/6.0D0
                  e = (3.0D0*a**2 - 1.0D0)*dq/6.0D0
                  f = (3.0D0*b**2 - 1.0D0)*dq/6.0D0
                 
                  do q_i = 1, Nqs

                      y(:) = 0.0D0
                      y(q_i) = 1.0D0

                      dP_dq0 = (y(q_hi) - y(q_low))/dq - e*d2y_dx2(q_i,q_low) + f*d2y_dx2(q_i,q_hi)

                      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                      prefactor = u_vdW(i_grid,q_i) * dP_dq0 * dq0_dgradrho(i_grid)

                      do l = 1, 3
                          do m = 1, l
                                        
                              sigma (l, m) = sigma (l, m) -  e2 * prefactor * &
                                             (gradient_rho(i_grid,l) * gradient_rho(i_grid,m))
                          enddo
                     enddo
              
                     !! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                 end do
           
      end do

#ifdef __MPI
      call mp_sum(  sigma, intra_pool_comm )
#endif

      call dscal (9, 1.d0 / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3), sigma, 1)

      deallocate( d2y_dx2, u_vdW )

   END SUBROUTINE stress_vdW_DF_gradient



   !! ###############################################################################################################
   !!                                       |                          |
   !!                                       |  STRESS_VDW_DF_KERNEL    |
   !!                                       |                          |

   SUBROUTINE stress_vdW_DF_kernel (total_rho, q0, thetas, sigma)

      !! Modules to include
      !! ----------------------------------------------------------------------------------
      use gvect,                 ONLY : ngm, nl, g, nl, gg, igtongl, gl, ngl, gstart 
      USE fft_base,              ONLY : dfftp
      USE cell_base,             ONLY : omega, tpiba, tpiba2
      USE constants, ONLY: pi

      implicit none
     
      real(dp), intent(IN) :: q0(:) 
      real(dp), intent(IN) :: total_rho(:)
      real(dp), intent(inout) :: sigma(3,3)                     !  
      complex(dp), intent(IN) :: thetas(:,:) 

      real(dp), allocatable :: dkernel_of_dk(:,:)               !
      
      integer               :: l, m, q1_i, q2_i , g_i           !
      real(dp)              :: g2, ngmod2, g_kernel, G_multiplier             ! 
      integer               :: last_g, theta_i

      allocate( dkernel_of_dk(Nqs, Nqs) )

      sigma(:,:) = 0.0_DP

      !! --------------------------------------------------------------------------------------------------
      !! Integration in g-space
      !! ---------------------------------------------------------------------------------------------------

      last_g = -1

      G_multiplier = 1.0D0
      if (gamma_only) G_multiplier = 2.0D0

      do g_i = gstart, ngm

          g2 = gg (g_i) * tpiba2
          g_kernel = sqrt(g2)

          if ( igtongl(g_i) .ne. last_g) then

             call interpolate_Dkernel_Dk(g_kernel, dkernel_of_dk)  ! Gets the derivatives
             last_g = igtongl(g_i)

          end if
          
          do q2_i = 1, Nqs
             do q1_i = 1, Nqs
                 do l = 1, 3
                     do m = 1, l

                     sigma (l, m) = sigma (l, m) - G_multiplier * 0.5 * e2 *&
                                     thetas(nl(g_i),q1_i)*dkernel_of_dk(q1_i,q2_i)*conjg(thetas(nl(g_i),q2_i))* &
                                     (g (l, g_i) * g (m, g_i) * tpiba2) / g_kernel 
                     end do
                 end do 
             enddo
         end do      

         if (g_i < gstart ) sigma(:,:) = sigma(:,:) / G_multiplier
         
      enddo

#ifdef __MPI
      call mp_sum(  sigma, intra_pool_comm )
#endif
      
      deallocate( dkernel_of_dk )
      
   END SUBROUTINE stress_vdW_DF_kernel

  !! ###############################################################################################################
  !!                                    |                  |
  !!                                    |  GET_Q0_ON_GRID  |
  !!                                    |__________________|

  !! This routine first calculates the q value defined in (DION equations 11 and 12), then
  !! saturates it according to (SOLER equation 5).  
  
  SUBROUTINE get_q0_on_grid (total_rho, gradient_rho, q0, dq0_drho, dq0_dgradrho)
  !!
  !! more specifically it calculates the following
  !!
  !!     q0(ir) = q0 as defined above
  !!     dq0_drho(ir) = total_rho * d q0 /d rho 
  !!     dq0_dgradrho = total_rho / |gradient_rho| * d q0 / d |gradient_rho|
  !!
    
    USE fft_base,        ONLY : dfftp
    USE kernel_table,    ONLY : q_cut, q_min
    
    real(dp),  intent(IN)    :: total_rho(:), gradient_rho(:,:)      !! Input variables needed
    
    real(dp),  intent(OUT) :: q0(:), dq0_drho(:), dq0_dgradrho(:)    !! Output variables that have been allocated
  !                                                                  !! outside this routine but will be set here.
  !                                                                        _
  real(dp),   parameter      :: LDA_A  = 0.031091D0, LDA_a1 = 0.2137D0    !
  real(dp),   parameter      :: LDA_b1 = 7.5957D0  , LDA_b2 = 3.5876D0    ! see J.P. Perdew and Yue Wang, Phys. Rev. B 45, 13244 (1992). 
  real(dp),   parameter      :: LDA_b3 = 1.6382D0  , LDA_b4 = 0.49294D0   !_ 
  real(dp)                   :: Z_ab = -0.8491D0                          !! see DION

  integer,    parameter      :: m_cut = 12                                !! How many terms to include in the sum
  !                                                                       !! of SOLER equation 5
  
  
  real(dp)                   :: kF, r_s, sqrt_r_s, gradient_correction    !! Intermediate variables needed to get q and q0
  real(dp)                   :: LDA_1, LDA_2, q, exponent                 !!
 
  real(dp)                   :: dq0_dq                                    !! The derivative of the saturated q0 with respect to q.
  !                                                                       !! Needed by dq0_drho and dq0_dgradrho by the chain rule.
  
  integer                    :: i_grid, idx, count=0                    !! Indexing variables
 
  if (vdw_type==1) Z_ab = -0.8491D0
  if (vdw_type==2) Z_ab = -1.887D0
  
 
  ! initialize q0-related arrays ... 
  q0(:) = q_cut
  dq0_drho(:) = 0.0_DP
  dq0_dgradrho(:) = 0.0_DP

  do i_grid = 1, dfftp%nnr
          
     !! This prevents numerical problems.  If the charge density is negative (an
     !! unphysical situation), we simply treat it as very small.  In that case,
     !! q0 will be very large and will be saturated.  For a saturated q0 the derivative
     !! dq0_dq will be 0 so we set q0 = q_cut and dq0_drho = dq0_dgradrho = 0 and go on
     !! to the next point.
     !! ------------------------------------------------------------------------------------
     
     if (total_rho(i_grid) < epsr) cycle

     !! ------------------------------------------------------------------------------------
     
     !! Calculate some intermediate values needed to find q
     !! ------------------------------------------------------------------------------------

     kF = (3.0D0*pi*pi*total_rho(i_grid))**(1.0D0/3.0D0)
     r_s = (3.0D0/(4.0D0*pi*total_rho(i_grid)))**(1.0D0/3.0D0)
     sqrt_r_s = sqrt(r_s)
     
     gradient_correction = -Z_ab/(36.0D0*kF*total_rho(i_grid)**2) &
          * (gradient_rho(i_grid,1)**2+gradient_rho(i_grid,2)**2+gradient_rho(i_grid,3)**2)
     
     LDA_1 =  8.0D0*pi/3.0D0*(LDA_A*(1.0D0+LDA_a1*r_s))
     LDA_2 =  2.0D0*LDA_A * (LDA_b1*sqrt_r_s + LDA_b2*r_s + LDA_b3*r_s*sqrt_r_s + LDA_b4*r_s*r_s)
     
     !! ------------------------------------------------------------------------------------

     !! This is the q value defined in equations 11 and 12 of DION
     !! ---------------------------------------------------------------
     
     q = kF + LDA_1 * log(1.0D0+1.0D0/LDA_2) + gradient_correction
     
     !! ---------------------------------------------------------------
     
     !! Here, we calculate q0 by saturating q according to equation 5 of SOLER.  Also, we find 
     !! the derivative dq0_dq needed for the derivatives dq0_drho and dq0_dgradrh0 discussed below.
     !! ---------------------------------------------------------------------------------------

     exponent = 0.0D0
     dq0_dq = 0.0D0
     
     do idx = 1, m_cut
        
        exponent = exponent + ( (q/q_cut)**idx)/idx
        dq0_dq = dq0_dq + ( (q/q_cut)**(idx-1))
        
     end do
     
     q0(i_grid) = q_cut*(1.0D0 - exp(-exponent))
     dq0_dq = dq0_dq * exp(-exponent)
     
     !! ---------------------------------------------------------------------------------------
     
     !! This is to handle a case with q0 too small.  We simply set it to the smallest q value in
     !! out q_mesh.  Hopefully this doesn't get used often (ever)
     !! ---------------------------------------------------------------------------------------

     if (q0(i_grid) < q_min) then
        
        q0(i_grid) = q_min
        
     end if

     !! ---------------------------------------------------------------------------------------
     
     !! Here we find derivatives.  These are actually the density times the derivative of q0 with respect
     !! to rho and gradient_rho.  The density factor comes in since we are really differentiating
     !! theta = (rho)*P(q0) with respect to density (or its gradient) which will be
     !! dtheta_drho = P(q0) + dP_dq0 * [rho * dq0_dq * dq_drho]   and
     !! dtheta_dgradient_rho =  dP_dq0  * [rho * dq0_dq * dq_dgradient_rho]
     !! The parts in square brackets are what is calculated here.  The dP_dq0 term will be interpolated
     !! later.  There should actually be a factor of the magnitude of the gradient in the gradient_rho derivative
     !! but that cancels out when we differentiate the magnitude of the gradient with respect to a particular
     !! component.
     !! -------------------------------------------------------------------------------------------------------------------------

     dq0_drho(i_grid) = dq0_dq * (kF/3.0D0 - 7.0D0/3.0D0*gradient_correction &
          - 8.0D0*pi/9.0D0 * LDA_A*LDA_a1*r_s*log(1.0D0+1.0D0/LDA_2) &
          + LDA_1/(LDA_2*(1.0D0 + LDA_2)) & 
          * (2.0D0*LDA_A*(LDA_b1/6.0D0*sqrt_r_s + LDA_b2/3.0D0*r_s + LDA_b3/2.0D0*r_s*sqrt_r_s + 2.0D0*LDA_b4/3.0D0*r_s**2)))
     
     dq0_dgradrho(i_grid) = total_rho(i_grid) * dq0_dq * 2.0D0 * (-Z_ab)/(36.0D0*kF*total_rho(i_grid)**2)
     
     !! --------------------------------------------------------------------------------------------------------------------------

  end do

end SUBROUTINE get_q0_on_grid


!! ###############################################################################################################


!! ###############################################################################################################
!!                                      |                      |
!!                                      |  GET_THETAS_ON_GRID  |
!!                                      |______________________|


SUBROUTINE get_thetas_on_grid (total_rho, q0_on_grid, thetas)

  real(dp), intent(in) :: total_rho(:), q0_on_grid(:)  !! Input arrays

  complex(dp), intent(inout):: thetas(:,:)             !! value of thetas for the grid points
  !                                                    !! assigned to this processor. The format
  !                                                    !! is thetas(grid_point, theta_i)
  !     NB: thetas are returned in reciprocal space
  
  integer :: i_grid, Ngrid_points                      !! An index for the point on the grid and the total
  !                                                    !! number of grid points
  integer :: theta_i                                   !! an index
  
  Ngrid_points = size(q0_on_grid)
  
  
  !! Interpolate the P_i polynomials defined in equation 3 in SOLER for the particular
  !! q0 values we have.
  CALL spline_interpolation(q_mesh, q0_on_grid, thetas)
  
  !! Form the thetas where theta is defined as rho*p_i(q0)
  !! ------------------------------------------------------------------------------------

  do i_grid = 1, Ngrid_points
     
     thetas(i_grid,:) = thetas(i_grid,:) * total_rho(i_grid)
     
  end do

  !! ------------------------------------------------------------------------------------
  !! Get thetas in reciprocal space.
  
  call start_clock( 'vdW_ffts')

  do theta_i = 1, Nqs

     CALL fwfft ('Dense', thetas(:,theta_i), dfftp)
  end do

  call stop_clock( 'vdW_ffts')
  
END SUBROUTINE get_thetas_on_grid


!! ###############################################################################################################

  
  
  




!! ###############################################################################################################
!!                                     |                        | 
!!                                     |  SPLINE_INTERPOLATION  |
!!                                     |________________________|


!! This routine is modeled after an algorithm from "Numerical Recipes in C" by Cambridge University
!! press, page 97.  It was adapted for Fortran, of course and for the problem at hand, in that it finds
!! the bin a particular x value is in and then loops over all the P_i functions so we only have to find
!! the bin once.


SUBROUTINE spline_interpolation (x, evaluation_points, values)
  
  
  real(dp), intent(in) :: x(:), evaluation_points(:)          !! Input variables.  The x values used to form the interpolation
  !                                                           !! (q_mesh in this case) and the values of q0 for which we are 
  !                                                           !! interpolating the function
  
  complex(dp), intent(inout) :: values(:,:)                   !! An output array (allocated outside this routine) that stores the
  !                                                           !! interpolated values of the P_i (SOLER equation 3) polynomials.  The
  !                                                           !! format is values(grid_point, P_i)
  
  integer :: Ngrid_points, Nx                                 !! Total number of grid points to evaluate and input x points
  
  real(dp), allocatable, save :: d2y_dx2(:,:)                 !! The second derivatives required to do the interpolation
  
  integer :: i_grid, lower_bound, upper_bound, idx, P_i     !! Some indexing variables
  
  real(dp), allocatable :: y(:)                               !! Temporary variables needed for the interpolation
  real(dp) :: a, b, c, d, dx                                  !!
  
  Nx = size(x)
  Ngrid_points = size(evaluation_points)
  

  !! Allocate the temporary array
  allocate( y(Nx) )

  !! If this is the first time this routine has been called we need to get the second
  !! derivatives (d2y_dx2) required to perform the interpolations.  So we allocate the
  !! array and call initialize_spline_interpolation to get d2y_dx2.
  !! ------------------------------------------------------------------------------------

  if (.not. allocated(d2y_dx2) ) then

     allocate( d2y_dx2(Nx,Nx) )
     call initialize_spline_interpolation(x, d2y_dx2)
     
  end if

  !! ------------------------------------------------------------------------------------
  
  
  do i_grid=1, Ngrid_points
     
     lower_bound = 1
     upper_bound = Nx
     
     do while ( (upper_bound - lower_bound) > 1 )
        
        idx = (upper_bound+lower_bound)/2
        
        if ( evaluation_points(i_grid) > x(idx) ) then
           lower_bound = idx 
        else
           upper_bound = idx
        end if
        
     end do
     
     dx = x(upper_bound)-x(lower_bound)
     
     a = (x(upper_bound) - evaluation_points(i_grid))/dx
     b = (evaluation_points(i_grid) - x(lower_bound))/dx
     c = ((a**3-a)*dx**2)/6.0D0
     d = ((b**3-b)*dx**2)/6.0D0
     
     do P_i = 1, Nx
        
        y = 0
        y(P_i) = 1
        
        values(i_grid, P_i) = a*y(lower_bound) + b*y(upper_bound) &
             + (c*d2y_dx2(P_i,lower_bound) + d*d2y_dx2(P_i, upper_bound))
     end do
     
  end do

  deallocate( y )

END SUBROUTINE spline_interpolation

  
!! ###############################################################################################################






!! ###############################################################################################################
!!                                |                                   |
!!                                |  INITIALIZE_SPLINE_INTERPOLATION  |
!!                                |___________________________________|


!! This routine is modeled after an algorithm from "Numerical Recipes in C" by Cambridge
!! University Press, pages 96-97.  It was adapted for Fortran and for the problem at hand.

SUBROUTINE initialize_spline_interpolation (x, d2y_dx2)
  
  real(dp), intent(in)  :: x(:)                    !! The input abscissa values 
  real(dp), intent(inout) :: d2y_dx2(:,:)          !! The output array (allocated outside this routine)
  !                                                !! that holds the second derivatives required for 
  !                                                !! interpolating the function

  integer :: Nx, P_i, idx                        !! The total number of x points and some indexing
  !                                                !! variables

  real(dp), allocatable :: temp_array(:), y(:)     !! Some temporary arrays required.  y is the array
  !                                                !! that holds the funcion values (all either 0 or 1 here).

  real(dp) :: temp1, temp2                         !! Some temporary variables required


  
  Nx = size(x)
  
  allocate( temp_array(Nx), y(Nx) )

  do P_i=1, Nx
     

     !! In the Soler method, the polynomicals that are interpolated are Kroneker delta funcions
     !! at a particular q point.  So, we set all y values to 0 except the one corresponding to 
     !! the particular function P_i.
     !! ----------------------------------------------------------------------------------------

     y = 0.0D0
     y(P_i) = 1.0D0

     !! ----------------------------------------------------------------------------------------
     
     d2y_dx2(P_i,1) = 0.0D0
     temp_array(1) = 0.0D0
     
     do idx = 2, Nx-1
        
        temp1 = (x(idx)-x(idx-1))/(x(idx+1)-x(idx-1))
        temp2 = temp1 * d2y_dx2(P_i,idx-1) + 2.0D0
        d2y_dx2(P_i,idx) = (temp1-1.0D0)/temp2
        temp_array(idx) = (y(idx+1)-y(idx))/(x(idx+1)-x(idx)) &
             - (y(idx)-y(idx-1))/(x(idx)-x(idx-1))
        temp_array(idx) = (6.0D0*temp_array(idx)/(x(idx+1)-x(idx-1)) &
             - temp1*temp_array(idx-1))/temp2
        
     end do
     
     d2y_dx2(P_i,Nx) = 0.0D0
     
     do idx=Nx-1, 1, -1
        
        d2y_dx2(P_i,idx) = d2y_dx2(P_i,idx) * d2y_dx2(P_i,idx+1) + temp_array(idx)
        
     end do

  end do

  deallocate( temp_array, y)

end SUBROUTINE initialize_spline_interpolation


!! ###############################################################################################################






!! ###############################################################################################################
!!                                         |                    |
!!                                         | INTERPOLATE_KERNEL |
!!                                         |____________________|


!! This routine is modeled after an algorithm from "Numerical Recipes in C" by Cambridge
!! University Press, page 97.  Adapted for Fortran and the problem at hand.  This function is used to 
!! find the Phi_alpha_beta needed for equations 8 and 11 of SOLER.


subroutine interpolate_kernel(k, kernel_of_k)
  
  USE kernel_table,             ONLY : r_max, Nr_points, kernel, d2phi_dk2, dk

  real(dp), intent(in) :: k                     !! Input value, the magnitude of the g-vector for the 
  !                                             !! current point.
  
  real(dp), intent(inout) :: kernel_of_k(:,:)   !! An output array (allocated outside this routine)
  !                                             !! that holds the interpolated value of the kernel
  !                                             !! for each pair of q points (i.e. the phi_alpha_beta 
  !                                             !! of the Soler method.

  integer :: q1_i, q2_i, k_i                    !! Indexing variables
 
  real(dp) :: A, B, C, D                        !! Intermediate values for the interpolation
  

  !! Check to make sure that the kernel table we have is capable of dealing with this
  !! value of k.  If k is larger than Nr_points*2*pi/r_max then we can't perform the 
  !! interpolation.  In that case, a kernel file should be generated with a larger number
  !! of radial points.
  !! -------------------------------------------------------------------------------------

  if ( k >= Nr_points*dk ) then
     
     write(*,'(A,F10.5,A,F10.5)') "k =  ", k, "     k_max =  ",Nr_points*dk
     call errore('interpolate kernel', 'k value requested is out of range',1)
     
  end if

  !! -------------------------------------------------------------------------------------
  
  kernel_of_k = 0.0D0
  
  !! This integer division figures out which bin k is in since the kernel
  !! is set on a uniform grid.
  k_i = int(k/dk)
  
  !! Test to see if we are trying to interpolate a k that is one of the actual
  !! function points we have.  The value is just the value of the function in that
  !! case.
  !! ----------------------------------------------------------------------------------------

  if (mod(k,dk) == 0) then
     
     do q1_i = 1, Nqs
        do q2_i = 1, q1_i
           
           kernel_of_k(q1_i, q2_i) = kernel(k_i,q1_i, q2_i)
           kernel_of_k(q2_i, q1_i) = kernel(k_i,q2_i, q1_i)
           
        end do
     end do
     
     return
     
  end if

  !! ----------------------------------------------------------------------------------------


  !! If we are not on a function point then we carry out the interpolation
  !! ----------------------------------------------------------------------------------------
  
  A = (dk*(k_i+1.0D0) - k)/dk
  B = (k - dk*k_i)/dk
  C = (A**3-A)*dk**2/6.0D0
  D = (B**3-B)*dk**2/6.0D0
  
  do q1_i = 1, Nqs
     do q2_i = 1, q1_i
        
        kernel_of_k(q1_i, q2_i) = A*kernel(k_i, q1_i, q2_i) + B*kernel(k_i+1, q1_i, q2_i) &
             +(C*d2phi_dk2(k_i, q1_i, q2_i) + D*d2phi_dk2(k_i+1, q1_i, q2_i))
        
        kernel_of_k(q2_i, q1_i) = kernel_of_k(q1_i, q2_i)
        
     end do
  end do

  !! ----------------------------------------------------------------------------------------

  
end subroutine interpolate_kernel


!! ###############################################################################################################




!! ###############################################################################################################
!!                                         |                        |
!!                                         | INTERPOLATE_DKERNEL_DK |
!!                                         |________________________|



subroutine interpolate_Dkernel_Dk(k, dkernel_of_dk)
  
  USE kernel_table,             ONLY : r_max, Nr_points, kernel, d2phi_dk2, dk

  implicit none 

  real(dp), intent(in) :: k                     !! Input value, the magnitude of the g-vector for the 
  !                                             !! current point.
  
  real(dp), intent(inout) :: dkernel_of_dk(Nqs,Nqs)   !! An output array (allocated outside this routine)
  !                                             !! that holds the interpolated value of the kernel
  !                                             !! for each pair of q points (i.e. the phi_alpha_beta 
  !                                             !! of the Soler method.

  integer :: q1_i, q2_i, k_i                    !! Indexing variables
 
  real(dp) :: A, B, dAdk, dBdk, dCdk, dDdk      !! Intermediate values for the interpolation
  

  !! -------------------------------------------------------------------------------------

  if ( k >= Nr_points*dk ) then
     
     write(*,'(A,F10.5,A,F10.5)') "k =  ", k, "     k_max =  ",Nr_points*dk
     call errore('interpolate kernel', 'k value requested is out of range',1)
     
  end if

  !! -------------------------------------------------------------------------------------

  dkernel_of_dk = 0.0D0

  k_i = int(k/dk)


  !! ----------------------------------------------------------------------------------------

  A = (dk*(k_i+1.0D0) - k)/dk
  B = (k - dk*k_i)/dk

  dAdk = -1.0D0/dk
  dBdk = 1.0D0/dk
  dCdk = -((3*A**2 -1.0D0)/6.0D0)*dk
  dDdk = ((3*B**2 -1.0D0)/6.0D0)*dk

  do q1_i = 1, Nqs
     do q2_i = 1, q1_i

        dkernel_of_dk(q1_i, q2_i) = dAdk*kernel(k_i, q1_i, q2_i) + dBdk*kernel(k_i+1, q1_i, q2_i) &
             + dCdk*d2phi_dk2(k_i, q1_i, q2_i) + dDdk*d2phi_dk2(k_i+1, q1_i, q2_i)

        dkernel_of_dk(q2_i, q1_i) = dkernel_of_dk(q1_i, q2_i)

     end do
  end do

  !! ----------------------------------------------------------------------------------------

  
end subroutine interpolate_Dkernel_Dk 



!! ###############################################################################################################
!!                                       |                       |
!!                                       |   NUMERICAL_GRADIENT  |
!!                                       |_______________________|


#ifdef FFTGRADIENT
!! Calculates the gradient of the charge density numerically on the grid.  We use
!! the PWSCF gradient style.

subroutine numerical_gradient(total_rho, gradient_rho)

   use gvect,             ONLY : ngm, nl, g, nlm
   USE cell_base,         ONLY : tpiba
   USE fft_base,          ONLY : dfftp
   USE fft_interfaces,    ONLY : fwfft, invfft 
   !
   ! I/O variables
   !
   real(dp), intent(in) :: total_rho(:)       !! Input array holding total charge density.
 
   real(dp), intent(out) :: gradient_rho(:,:) !! Output array that will holds the gradient
   !                                          !! of the charge density.
   ! local variables
   !
   integer :: icar                            !! counter on cartesian components
   complex(dp), allocatable :: c_rho(:)       !! auxiliary complex array for rho
   complex(dp), allocatable :: c_grho(:)      !! auxiliary complex array for grad rho
 
   ! rho in G space
   allocate ( c_rho(dfftp%nnr), c_grho(dfftp%nnr) )
   c_rho(1:dfftp%nnr) = CMPLX(total_rho(1:dfftp%nnr),0.0_DP)
   CALL fwfft ('Dense', c_rho, dfftp) 
 
   do icar=1,3
      ! compute gradient in G space
      c_grho(:) =CMPLX(0.0_DP,0.0_DP)
      c_grho(nl(:)) = CMPLX (0.0_DP,1.0_DP) * tpiba * g(icar,:) * c_rho(nl(:))
      if (gamma_only) c_grho( nlm(:) ) = CONJG( c_grho( nl(:) ) )
      ! back in real space
      CALL invfft ('Dense', c_grho, dfftp) 
      gradient_rho(:,icar) = REAL( c_grho(:) )
   end do
   deallocate ( c_rho, c_grho )

   return

end subroutine numerical_gradient


#else
!! Calculates the gradient of the charge density numerically on the grid.  We could simply
!! use the PWSCF gradient routine but we need the derivative of the gradient at point j 
!! with respect to the density at point i for the potential (SOLER equation 10).  This is 
!! difficult to do with the standard means of calculating the density gradient but trivial
!! in the case of the numerical formula because the derivative of the gradient at point j
!! with respect to the density at point i is just whatever the coefficient is in the numerical
!! derivative formula.


subroutine numerical_gradient(full_rho, Nneighbors, gradient_rho, my_start_z, my_end_z)

  USE fft_base,          ONLY : dfftp
  USE cell_base,         ONLY : alat, at
  
  real(dp), intent(in) :: full_rho(:)                        !! Input array holding the value of the total charge density
  !                                                          !! on all grid points of the simulation cell
  
  integer, intent(in) :: Nneighbors, my_start_z, my_end_z    !! Input variables giving the order of the numerical derivative,
  !                                                          !! and the starting and ending z-slabs for the given processor.

  real(dp), intent(inout) :: gradient_rho(:,:)               !! Output array (allocated outside the routine) that holds the 
  !                                                          !! gradient of the charge density only in the region assigned to
  !                                                          !! the given processor in the format: 
  !                                                          !! gradient_rho(grid_point, cartesian_component)

  real(dp), pointer, save ::  coefficients(:)                !! A pointer to an array of coefficients used for the numerical
  !                                                          !! differentiation.  See gradient_coefficients function for more
  !                                                          !! detail.

  integer, pointer, save :: indices3d(:,:,:)                 !! A pointer to a rank 3 array that gives the relation between the
  !                                                          !! x, y, and z indices of a point and its index in the charge density
  !                                                          !! array.  Used to easily find neighbors in the x, y, and z directions.

  integer :: i_grid, ix1, ix2, ix3, nx                       !! Indexing variables
  
  real(dp) :: temp(3)                                        !! A temporary array for the gradient at a point 

  real(dp), save :: at_inverse(3,3)                          !! The inverse of the matrix of unit cell basis vectors
  
  logical, save :: have_at_inverse = .false.                 !! Flag to determine if we have found the inverse matrix yet


  gradient_rho = 0.0D0

  !! Get pointers to the gradient coefficients and the 3d index array needed to find
  !! the gradient if we don't have them already.
  !! ----------------------------------------------------------------------------------
  
  if (.not. associated(indices3d) ) then

     indices3d => get_3d_indices(Nneighbors)
     coefficients => gradient_coefficients(Nneighbors)
  
  end if

  !! ----------------------------------------------------------------------------------

  !! Here we need to get the transformation matrix that takes our calculated "gradient"
  !! , gradient_rho()!! to the real thing.  It is just the (normalized) inverse of the matrix of unit cell 
  !! basis vectors.  If the unit cell has orthogonal basis vectors then this will be a 
  !! diagonal matrix with the diagonal elements bein 1/(basis vector length).  In the 
  !! general case this will not be diagonal (e.g. for hexagonal unit cells). 
  !! ----------------------------------------------------------------------------------
  
  if (.not. have_at_inverse) then

     at_inverse = alat*at
     call invert_3x3_matrix(at_inverse)

     ! Normalize by the number of grid points in each direction
     ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     at_inverse(1,:) = at_inverse(1,:) * dble(dfftp%nr1x)
     at_inverse(2,:) = at_inverse(2,:) * dble(dfftp%nr2x)
     at_inverse(3,:) = at_inverse(3,:) * dble(dfftp%nr3x)

     ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     
     ! Take the transpose because of the way Fortran does matmul() (used below)
     at_inverse = transpose(at_inverse)

     ! Mark that we have gotten the transformation matrix so we don't have to find it again
     have_at_inverse = .true.

  end if

  !! ----------------------------------------------------------------------------------



  i_grid = 0

  !! Here we loop over all of the points assigned to a given processor.  For each point we loop
  !! over all relavant neighbors (determined by the variable Nneighbors) and multiply the value
  !! of the density of each by the corresponding coefficient.  We then tranform the vector by 
  !! multiplying it by the inverse of the unit cell matrix found above.  This takes care of cases
  !! where the basis vectors are not the same length or are not even orthogonal.
  !! -----------------------------------------------------------------------------------------------

  do ix3 = my_start_z, my_end_z
     do ix2 = 1, dfftp%nr2x
        do ix1 = 1, dfftp%nr1x
           
           i_grid = i_grid + 1
          
           temp = 0.0D0

           do nx = -Nneighbors, Nneighbors
              
              temp(1) = temp(1) + coefficients(nx) * full_rho(indices3d(ix1+nx,ix2,ix3))
              temp(2) = temp(2) + coefficients(nx) * full_rho(indices3d(ix1,ix2+nx,ix3))
              temp(3) = temp(3) + coefficients(nx) * full_rho(indices3d(ix1,ix2,ix3+nx))
              
           end do
           
           gradient_rho(i_grid,:) = matmul(at_inverse,temp)
           
        end do
     end do
  end do
  
  !! -----------------------------------------------------------------------------------------------
  !! FAKE PATCH !!
  !gradient_rho = 0.0D0

end subroutine numerical_gradient

#endif

!! #################################################################################################
!!                                          |              |
!!                                          | thetas_to_uk |
!!                                          |______________|


subroutine thetas_to_uk(thetas, u_vdW)
  
  USE gvect,           ONLY : nl, nlm, gg, ngm, igtongl, gl, ngl, gstart
  USE fft_base,        ONLY : dfftp
  USE cell_base,       ONLY : tpiba, omega

  complex(dp), intent(in) :: thetas(:,:)    !! On input this variable holds the theta functions (equation 8, SOLER)
  !                                         !! in the format thetas(grid_point, theta_i).  
  complex(dp), intent(out) :: u_vdW(:,:)
  !                                            !! On output this array holds u_alpha(k) = Sum_j[theta_beta(k)phi_alpha_beta(k)]

  real(dp), allocatable :: kernel_of_k(:,:)    !! This array will hold the interpolated kernel values for each pair of q values
  !                                            !! in the q_mesh.

  real(dp) :: g
  integer :: last_g, g_i, q1_i, q2_i, count, i_grid    !! Index variables

  complex(dp) :: theta(Nqs)                    !! Temporary storage vector used since we are overwriting the thetas array here.
  
  !! -------------------------------------------------------------------------------------------------
  
  allocate( kernel_of_k(Nqs, Nqs) )

  u_vdW(:,:) = CMPLX(0.0_DP,0.0_DP)
  
  last_g = -1 

  do g_i = 1, ngm
    
     if ( igtongl(g_i) .ne. last_g) then
        
        g = sqrt(gl(igtongl(g_i))) * tpiba
        call interpolate_kernel(g, kernel_of_k)
        last_g = igtongl(g_i)
        
     end if
     
     theta = thetas(nl(g_i),:)
     
     do q2_i = 1, Nqs
        do q1_i = 1, Nqs
           u_vdW(nl(g_i),q2_i) = u_vdW(nl(g_i),q2_i) + kernel_of_k(q2_i,q1_i)*theta(q1_i)
        end do
     end do

  end do

  if (gamma_only) u_vdW(nlm(:),:) = CONJG(u_vdW(nl(:),:))
  
  deallocate( kernel_of_k )
     
  !! -----------------------------------------------------------------------------------------------
  
end subroutine thetas_to_uk

!! #################################################################################################
!!                                              |             |
!!                                              | VDW_ENERGY  |
!!                                              |_____________|

!! This routine carries out the integration of equation 8 of SOLER.  It returns the non-local 
!! exchange-correlation energy and the u_alpha(k) arrays used to find the u_alpha(r) arrays via 
!! equations 11 and 12 in SOLER.

subroutine vdW_energy(thetas, vdW_xc_energy)
  
  USE gvect,           ONLY : nl, nlm, gg, ngm, igtongl, gl, ngl, gstart
  USE fft_base,        ONLY : dfftp
  USE cell_base,       ONLY : tpiba, omega

  complex(dp), intent(inout) :: thetas(:,:)    !! On input this variable holds the theta functions 
  !                                            !! (equation 8, SOLER) in the format thetas(grid_point, theta_i).  
  !                                            !! On output this array holds 
  !                                            !! u_alpha(k) = Sum_j[theta_beta(k)phi_alpha_beta(k)]

  real(dp), intent(out) :: vdW_xc_energy       !! The non-local correlation energy.  An output variable.

  real(dp), allocatable :: kernel_of_k(:,:)    !! This array will hold the interpolated kernel values 
  !                                            !! for each pair of q values in the q_mesh.

  real(dp) :: g                                !! The magnitude of the current g vector 
  integer  :: last_g                           !! The shell number of the last g vector
  !                                          
  
  integer :: g_i, q1_i, q2_i, count, i_grid    !! Index variables

  complex(dp) :: theta(Nqs), thetam(Nqs), theta_g(Nqs)      !! Temporary storage vector used since we are overwriting the thetas array here.
  real(dp)    :: G0_term, G_multiplier 

  complex(dp), allocatable :: u_vdw(:,:)       !! temporary array holding u_alpha(k)

  vdW_xc_energy = 0.0D0
  allocate (u_vdW(dfftp%nnr,Nqs))
  u_vdW(:,:) = CMPLX(0.0_DP,0.0_DP)

  allocate( kernel_of_k(Nqs, Nqs) )
  
  !! Loop over PWSCF's array of magnitude-sorted g-vector shells.  For each shell, interpolate
  !! the kernel at this magnitude of g, then find all points on the shell and carry out the 
  !! integration over those points.  The PWSCF variables used here are 
  !! ngm = number of g-vectors on this processor, nl = an array that gives the indices into the 
  !! FFT grid for a particular g vector, igtongl = an array that gives the index of which shell a 
  !! particular g vector is in, gl = an array that gives the magnitude of the g vectors for each shell.
  !! In essence, we are forming the reciprocal-space u(k) functions of SOLER equation 11.  These are
  !! kept in thetas array.
  !! -------------------------------------------------------------------------------------------------

  !!
  !! Here we should use gstart,ngm but all the cases are handled by conditionals inside the loop
  !!
  G_multiplier = 1.0D0
  if (gamma_only) G_multiplier = 2.0D0

  last_g = -1 

  do g_i = 1, ngm
     
     if ( igtongl(g_i) .ne. last_g) then
        
        g = sqrt(gl(igtongl(g_i))) * tpiba
        call interpolate_kernel(g, kernel_of_k)
        last_g = igtongl(g_i)
        
     end if
     
     theta = thetas(nl(g_i),:)

     do q2_i = 1, Nqs
        do q1_i = 1, Nqs
           u_vdW(nl(g_i),q2_i)  = u_vdW(nl(g_i),q2_i) + kernel_of_k(q2_i,q1_i)*theta(q1_i)
        end do
        vdW_xc_energy = vdW_xc_energy + G_multiplier * (u_vdW(nl(g_i),q2_i)*conjg(theta(q2_i)))
     end do
     
     if (g_i < gstart ) vdW_xc_energy = vdW_xc_energy / G_multiplier

  end do

  if (gamma_only) u_vdW(nlm(:),:) = CONJG(u_vdW(nl(:),:))
 
  !! ---------------------------------------------------------------------------------------------------
 
  !! Apply scaling factors.  The e2 comes from PWSCF's choice of units.  This should be 
  !! 0.5 * e2 * vdW_xc_energy * (2pi)^3/omega * (omega)^2, with the (2pi)^3/omega being
  !! the volume element for the integral (the volume of the reciprocal unit cell) and the 
  !! 2 factors of omega being used to cancel the factor of 1/omega PWSCF puts on forward  
  !! FFTs of the 2 theta factors.  1 omega cancels and the (2pi)^3 cancels because there should
  !! be a factor of 1/(2pi)^3 on the radial Fourier transform of phi that was left out to cancel
  !! with this factor.
  !! ---------------------------------------------------------------------------------------------------
 
  vdW_xc_energy = 0.5D0 * e2 * omega * vdW_xc_energy 
  
  deallocate( kernel_of_k )
  thetas(:,:) = u_vdW(:,:)
  deallocate (u_vdW)
  !! ---------------------------------------------------------------------------------------------------
  
end subroutine vdW_energy


!! ###############################################################################################################

!! ###############################################################################################################
!!                                             |                 |
!!                                             |  dv_drho_vdw    |
!!                                             |_________________|

#ifdef FFTGRADIENT

subroutine dv_drho_vdw(rho_valence, rho_core, drho, nspin, dv_drho)

    USE gvect,               ONLY : nl, g, nlm, ngm
    USE fft_base,            ONLY : dfftp
    USE cell_base,           ONLY : alat, tpiba, omega

    integer :: nspin
    real(dp), intent(IN) :: rho_valence(:,:)       !
    real(dp), intent(IN) :: rho_core(:)    
    complex(DP), intent(IN) :: drho (dfftp%nnr, nspin)
    complex(DP), intent(INOUT) :: dv_drho(dfftp%nnr, nspin)

    !! -------------------------------------------------------------------------
    !! For the potential
    !! -------------------------------------------------------------------------
    integer :: i_grid, theta_i, i_proc, I
    real(dp) :: grid_cell_volume
    real(dp), allocatable ::  q0(:)
    real(dp), allocatable :: gradient_rho(:,:)
    real(dp), allocatable :: potential(:)
    real(dp), allocatable :: dq0_drho(:)
    real(dp), allocatable :: dq0_dgradrho(:)
    complex(dp), allocatable :: thetas(:,:)
    real(dp), allocatable :: total_rho(:)
    complex(dp), allocatable :: u_vdW(:,:)

    !! -------------------------------------------------------------------------
    !! For the derivative
    !! -------------------------------------------------------------------------
    real(dp), allocatable :: potential_plus(:), potential_minus(:)
    real(dp) :: lambda
    real(DP), allocatable :: drho_real(:)

 
    allocate( q0(dfftp%nnr) )
    allocate( gradient_rho(dfftp%nnr, 3) )
    allocate( dq0_drho(dfftp%nnr), dq0_dgradrho(dfftp%nnr) )
    allocate( total_rho(dfftp%nnr) )
    allocate( drho_real(dfftp%nnr) )
    allocate( thetas(dfftp%nnr, Nqs) )
    allocate( u_vdW(dfftp%nnr, Nqs) )
    allocate( potential_plus(dfftp%nnr), potential_minus(dfftp%nnr) )
   

    !! Derivative parameter 
    lambda = 0.01D0


    !! Delta rho in real space 
    CALL invfft ('Dense', drho(:,1), dfftp)
    drho_real(:) = REAL( drho(:,1) )


    !! -------------------------------------------------------------------------
    !! Potential plus
    !! -------------------------------------------------------------------------
    
    total_rho = rho_valence(:,1) + rho_core(:) + lambda*drho_real(:)

    call numerical_gradient(total_rho,gradient_rho)
    CALL get_q0_on_grid(total_rho, gradient_rho, q0, dq0_drho, dq0_dgradrho)
    CALL get_thetas_on_grid(total_rho, q0, thetas)
    
    !!call vdW_energy(thetas, Ec_nl)
    call thetas_to_uk(thetas, u_vdW)
    
    call start_clock( 'vdW_ffts')
    do theta_i = 1, Nqs

       CALL invfft('Dense', u_vdW(:,theta_i), dfftp)
    end do
    call stop_clock( 'vdW_ffts')
    
    !!call get_potential(q0, dq0_drho, dq0_dgradrho, gradient_rho, thetas, potential)
    call get_potential(q0, dq0_drho, dq0_dgradrho, gradient_rho, u_vdW, potential_plus)
  

    !! -------------------------------------------------------------------------
    !! Potential minus
    !! -------------------------------------------------------------------------

    total_rho = rho_valence(:,1) + rho_core(:) - lambda*drho_real(:)

    call numerical_gradient(total_rho,gradient_rho)
    CALL get_q0_on_grid(total_rho, gradient_rho, q0, dq0_drho, dq0_dgradrho)
    CALL get_thetas_on_grid(total_rho, q0, thetas)

    !!call vdW_energy(thetas, Ec_nl)
    call thetas_to_uk(thetas, u_vdW)

    call start_clock( 'vdW_ffts')
    do theta_i = 1, Nqs

       CALL invfft('Dense', u_vdW(:,theta_i), dfftp)
    end do
    call stop_clock( 'vdW_ffts')

    !!call get_potential(q0, dq0_drho, dq0_dgradrho, gradient_rho, thetas, potential)
    call get_potential(q0, dq0_drho, dq0_dgradrho, gradient_rho, u_vdW, potential_minus)


    !! -------------------------------------------------------------------------
    !! Derivative
    !! -------------------------------------------------------------------------
  
    dv_drho(:,1) = (potential_plus(:) - potential_minus(:))/(2*lambda)


    !! -------------------------------------------------------------------------
    !! Deallocate
    !! -------------------------------------------------------------------------

    CALL fwfft ('Dense', drho(:,1), dfftp)

    deallocate( q0, gradient_rho, dq0_drho, dq0_dgradrho, total_rho)
    deallocate( drho_real,thetas, u_vdW)
    deallocate( potential_plus, potential_minus )
    

end subroutine dv_drho_vdw


#endif

!! ###############################################################################################################
!!                                             |                 |
!!                                             |  GET_POTENTIAL  |
!!                                             |_________________|

!! This routine finds the non-local correlation contribution to the potential (i.e. the derivative of the non-local
!! piece of the energy with respect to density) given in SOLER equation 10.  The u_alpha(k) functions were found
!! while calculating the energy.  They are passed in as the matrix u_vdW.  Most of the required derivatives were
!! calculated in the "get_q0_on_grid" routine, but the derivative of the interpolation polynomials, P_alpha(q),
!! (SOLER equation 3) with respect to q is interpolated here, along with the polynomials themselves.

#ifdef FFTGRADIENT

subroutine get_potential(q0, dq0_drho, dq0_dgradrho, gradient_rho, u_vdW, potential)

  use gvect,               ONLY : nl, g, nlm
  USE fft_base,            ONLY : dfftp
  USE cell_base,           ONLY : alat, tpiba
  
  real(dp), intent(in) ::  q0(:), gradient_rho(:,:)   !! Input arrays holding the value of q0 for all points assigned
  !                                                   !! to this processor and the gradient of the charge density for
  !                                                   !! points assigned to this processor.

  real(dp), intent(in) :: dq0_drho(:), dq0_dgradrho(:)!! The derivative of q0 with respect to the charge density and
  !                                                   !! gradient of the charge density (almost).  See comments in
  !                                                   !! the get_q0_on_grid subroutine above.

  complex(dp), intent(in) :: u_vdW(:,:)               !! The functions u_alpha(r) obtained by inverse transforming the 
  !                                                   !! functions u_alph(k).  See equations 11 and 12 in SOLER

  real(dp), intent(inout) :: potential(:)             !! The non-local correlation potential for points on the grid over
  !                                                   !! the whole cell (not just those assigned to this processor).

  real(dp), allocatable, save :: d2y_dx2(:,:)         !! Second derivatives of P_alpha polynomials for interpolation

  integer :: i_grid, P_i,icar                         !! Index variables

  integer :: q_low, q_hi, q                           !! Variables to find the bin in the q_mesh that a particular q0
  !                                                   !! belongs to (for interpolation).
  real(dp) :: dq, a, b, c, d, e, f                    !! Intermediate variables used in the interpolation of the polynomials
  
  real(dp) :: y(Nqs), dP_dq0, P                       !! The y values for a given polynomial (all 0 exept for element i of P_i)
  !                                                   !! The derivative of P at a given q0 and the value of P at a given q0.  Both
  !                                                   !! of these are interpolated below

  real(dp), allocatable ::h_prefactor(:)
  complex(dp), allocatable ::h(:)

  allocate (h_prefactor(dfftp%nnr),h(dfftp%nnr))

  potential = 0.0D0
  h_prefactor   = 0.0D0

  !! -------------------------------------------------------------------------------------------
  
  !! Get the second derivatives of the P_i functions for interpolation.  We have already calculated
  !! this once but it is very fast and it's just as easy to calculate it again.
  !! ---------------------------------------------------------------------------------------------

  if (.not. allocated( d2y_dx2) ) then
     
     allocate( d2y_dx2(Nqs, Nqs) )
     
     call initialize_spline_interpolation(q_mesh, d2y_dx2(:,:))
     
  end if
  
  !! ---------------------------------------------------------------------------------------------
  

  do i_grid = 1,dfftp%nnr
           
     q_low = 1
     q_hi = Nqs 

     ! Figure out which bin our value of q0 is in in the q_mesh
     ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     do while ( (q_hi - q_low) > 1)
              
        q = int((q_hi + q_low)/2)
              
        if (q_mesh(q) > q0(i_grid)) then
           q_hi = q
        else 
           q_low = q
        end if
              
     end do
           
     if (q_hi == q_low) call errore('get_potential','qhi == qlow',1)
           
     ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     dq = q_mesh(q_hi) - q_mesh(q_low)
           
     a = (q_mesh(q_hi) - q0(i_grid))/dq
     b = (q0(i_grid) - q_mesh(q_low))/dq
     c = (a**3 - a)*dq**2/6.0D0
     d = (b**3 - b)*dq**2/6.0D0
     e = (3.0D0*a**2 - 1.0D0)*dq/6.0D0
     f = (3.0D0*b**2 - 1.0D0)*dq/6.0D0
           
     do P_i = 1, Nqs
        y = 0.0D0
        y(P_i) = 1.0D0
              
        dP_dq0 = (y(q_hi) - y(q_low))/dq - e*d2y_dx2(P_i,q_low) + f*d2y_dx2(P_i,q_hi)
             
        P = a*y(q_low) + b*y(q_hi) + c*d2y_dx2(P_i,q_low) + d*d2y_dx2(P_i,q_hi)
              
        !! The first term in equation 10 of SOLER
        potential(i_grid) = potential(i_grid) + u_vdW(i_grid,P_i)* (P + dP_dq0 * dq0_drho(i_grid))
        if (q0(i_grid) .ne. q_mesh(Nqs)) then
           h_prefactor(i_grid) = h_prefactor(i_grid) +  u_vdW(i_grid,P_i)* dP_dq0 * dq0_dgradrho(i_grid)
        end if
     end do
  end do
  do icar = 1,3
     h(:) = CMPLX(h_prefactor(:) * gradient_rho(:,icar),0.0_DP)
     CALL fwfft ('Dense', h, dfftp) 
     h(nl(:)) = CMPLX(0.0_DP,1.0_DP) * tpiba * g(icar,:) * h(nl(:))
     if (gamma_only) h(nlm(:)) = CONJG(h(nl(:)))
     CALL invfft ('Dense', h, dfftp) 
     potential(:) = potential(:) - REAL(h(:))
  end do

  !! ------------------------------------------------------------------------------------------------------------------------
  deallocate (h_prefactor,h)

end subroutine get_potential

#else

subroutine get_potential(q0, dq0_drho, dq0_dgradrho, N, gradient_rho, u_vdW, potential, my_start_z, my_end_z)

  USE fft_base,            ONLY : dfftp
  USE cell_base,           ONLY : alat, at
  
  real(dp), intent(in) ::  q0(:), gradient_rho(:,:)       !! Input arrays holding the value of q0 for all points assigned
  !                                                       !! to this processor and the gradient of the charge density for
  !                                                       !! points assigned to this processor.

  real(dp), intent(in) :: dq0_drho(:), dq0_dgradrho(:)    !! The derivative of q0 with respect to the charge density and
  !                                                       !! gradient of the charge density (almost).  See comments in
  !                                                       !! the get_q0_on_grid subroutine above.

  real(dp), intent(inout) :: potential(:)                 !! The non-local correlation potential for points on the grid over
  !                                                       !! the whole cell (not just those assigned to this processor).

  integer, intent(in) :: N, my_start_z, my_end_z          !! The number of neighbors used in the numerical gradient formula
  !                                                       !! and the starting and ending z planes for this processor

  complex(dp), intent(in) :: u_vdW(:,:)                   !! The functions u_alpha(r) obtained by inverse transforming the 
  !                                                       !! functions u_alph(k).  See equations 11 and 12 in SOLER

  real(dp), allocatable, save :: d2y_dx2(:,:)             !! Second derivatives of P_alpha polynomials for interpolation

  integer :: i_grid, ix1, ix2, ix3, P_i, nx               !! Index variables

  integer :: q_low, q_hi, q                               !! Variables to find the bin in the q_mesh that a particular q0
  !                                                       !! belongs to (for interpolation).
  
  real(dp) :: prefactor                                   !! Intermediate variable used to minimize calculations

  real(dp), pointer, save :: coefficients(:)              !! Pointer to the gradient coefficients.  Used to find the derivative
  !                                                       !! of the magnitude of the gradient of the charge density with 
  !                                                       !! respect to the charge density at another point.  Equation 10 in SOLER

  integer, pointer, save :: indices3d(:,:,:)              !! A pointer to a rank 3 array that gives the relation between the
  !                                                       !! x, y, and z indices of a point and its index in the charge density
  !                                                       !! array.  Used to easily find neighbors in the x, y, and z directions.
  
  real(dp) :: dq, a, b, c, d, e, f                        !! Inermediate variables used in the interpolation of the polynomials
  
  real(dp) :: y(Nqs), dP_dq0, P                           !! The y values for a given polynomial (all 0 exept for element i of P_i)
  !                                                       !! The derivative of P at a given q0 and the value of P at a given q0.  Both
  !                                                       !! of these are interpolated below

  real(dp), save :: at_inverse(3,3)     

  logical, save :: have_at_inverse = .false.


  if (.not. have_at_inverse) then
     
     at_inverse = alat * at
     call invert_3x3_matrix(at_inverse)
     
     at_inverse(1,:) = at_inverse(1,:) * dble(dfftp%nr1x)
     at_inverse(2,:) = at_inverse(2,:) * dble(dfftp%nr2x)
     at_inverse(3,:) = at_inverse(3,:) * dble(dfftp%nr3x)

     at_inverse = transpose(at_inverse)
     
     have_at_inverse = .true.

  end if

  
  potential = 0.0D0

  !! Find the gradient coefficients and the 3d index mapping array if we don't already have it.
  !! -------------------------------------------------------------------------------------------

  if (.not. associated(indices3d) ) then
     
     indices3d => get_3d_indices()
     coefficients => gradient_coefficients()
  
  end if

  !! -------------------------------------------------------------------------------------------

  
  !! Get the second derivatives of the P_i functions for interpolation.  We have already calculated
  !! this once but it is very fast and it's just as easy to calculate it again.
  !! ---------------------------------------------------------------------------------------------

  if (.not. allocated( d2y_dx2) ) then
     
     allocate( d2y_dx2(Nqs, Nqs) )
     
     call initialize_spline_interpolation(q_mesh, d2y_dx2(:,:))
     
  end if
  
  !! ---------------------------------------------------------------------------------------------
  
i_grid = 0
  
  !! Loop over all the points assigned to this processor.  For each point and each q value in the q_mesh,
  !! interpolate P_i and dP_dq0.  
  !! --------------------------------------------------------------------------------------------------------------------

  do ix3 = my_start_z, my_end_z
     do ix2 = 1, dfftp%nr2x
        do ix1 = 1, dfftp%nr1x
           
           i_grid = i_grid + 1
           
           q_low = 1
           q_hi = Nqs 


           ! Figure out which bin our value of q0 is in in the q_mesh
           ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

           do while ( (q_hi - q_low) > 1)
              
              q = int((q_hi + q_low)/2)
              
              if (q_mesh(q) > q0(i_grid)) then
                 q_hi = q
              else 
                 q_low = q
              end if
              
           end do
           
           if (q_hi == q_low) call errore('get_potential','qhi == qlow',1)
           
           ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

           dq = q_mesh(q_hi) - q_mesh(q_low)
           
           a = (q_mesh(q_hi) - q0(i_grid))/dq
           b = (q0(i_grid) - q_mesh(q_low))/dq
           c = (a**3 - a)*dq**2/6.0D0
           d = (b**3 - b)*dq**2/6.0D0
           e = (3.0D0*a**2 - 1.0D0)*dq/6.0D0
           f = (3.0D0*b**2 - 1.0D0)*dq/6.0D0
           
           do P_i = 1, Nqs
              
              y = 0.0D0
              y(P_i) = 1.0D0
              
              dP_dq0 = (y(q_hi) - y(q_low))/dq - e*d2y_dx2(P_i,q_low) + f*d2y_dx2(P_i,q_hi)
             
              P = a*y(q_low) + b*y(q_hi) + c*d2y_dx2(P_i,q_low) + d*d2y_dx2(P_i,q_hi)
              
              !! The first term in equation 10 of SOLER
              potential(indices3d(ix1,ix2,ix3)) = potential(indices3d(ix1,ix2,ix3)) + &
                                                  u_vdW(i_grid,P_i)* (P + dP_dq0 * dq0_drho(i_grid))


              ! Now, loop over all relevant neighbors and calculate the second term in equation 10 of SOLER.  Note,
              ! that we are using our value of u_vdW and gradients and adding the piece of the potential point i_grid
              ! contributes to the neighbor's potential.  If the value of q0 at point i_grid is equal to q_cut, the
              ! derivative dq0_dq will be 0 so both of dq0_drho and dq0_dgradrho will be 0.  Thus, we can safely
              ! skip these points.
              ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

              if (q0(i_grid) .ne. q_mesh(Nqs)) then

                 prefactor = u_vdW(i_grid,P_i) * dP_dq0 * dq0_dgradrho(i_grid)
                 
                 do nx = -N,N
                    
                    potential(indices3d(ix1+nx,ix2,ix3)) = potential(indices3d(ix1+nx,ix2,ix3)) &
                         + prefactor * coefficients(nx) &
                         * (gradient_rho(i_grid,1)*at_inverse(1,1) + gradient_rho(i_grid,2)*at_inverse(2,1) &
                         + gradient_rho(i_grid,3)*at_inverse(3,1))
                    
                    potential(indices3d(ix1,ix2+nx,ix3)) = potential(indices3d(ix1,ix2+nx,ix3)) &
                         + prefactor * coefficients(nx) &
                         * (gradient_rho(i_grid,1)*at_inverse(1,2) + gradient_rho(i_grid,2)*at_inverse(2,2) &
                         + gradient_rho(i_grid,3)*at_inverse(3,2))
                         
                    potential(indices3d(ix1,ix2,ix3+nx)) = potential(indices3d(ix1,ix2,ix3+nx)) &
                         + prefactor * coefficients(nx) &
                         * (gradient_rho(i_grid,1)*at_inverse(1,3) + gradient_rho(i_grid,2)*at_inverse(2,3) &
                         + gradient_rho(i_grid,3)*at_inverse(3,3))
                    
                 end do
                 
              end if
              
              !! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

           end do
           
        end do
     end do
  end do

  !! ------------------------------------------------------------------------------------------------------------------------

end subroutine get_potential
#endif

!! ###############################################################################################################


!! ###############################################################################################################
!!                                            |                         |
!!                                            |  GRADIENT_COEFFICIENTS  |
!!                                            |_________________________|


!!  This routine returns a pointer to an array holding the coefficients for a derivative expansion to some order.
!!  The derivative is found by multiplying the value of the function at a point + or - n away from the sample point by
!!  the coefficient gradient_coefficients(+ or - n) and dividing by the appropriate dx for that direction.


function gradient_coefficients(N)
  
  real(dp), allocatable, target, save:: coefficients(:)  !! The local array that will hold the coefficients.  A pointer to this
  !                                                      !! array will be returned by the function

  integer, intent(in), optional :: N                     !! The number of neighbors to use on each side for the gradient
  !                                                      !! calculation.  Can be between 1 (i.e. 3 point derivative formula) 
  !                                                      !! and 6 (i.e. 13 point derivative formula).

  real(dp), pointer :: gradient_coefficients(:)          !! Pointer to the coefficients array that will be returned
  

  if (.not. allocated(coefficients) ) then

     if (.not. present(N) ) call errore('gradient_coefficients', 'Number of neighbors for gradient must be specified',2)
  
     allocate( coefficients(-N:N) )
        
     select case (N)
        
     case (1) 
        coefficients(-1:1) = (/-0.5D0, 0.0D0, 0.5D0/)
     case (2)
        coefficients(-2:2) = (/0.0833333333333333D0, -0.6666666666666666D0, 0.0D0, &
                               0.6666666666666666D0, -0.0833333333333333D0/)
     case (3) 
        coefficients(-3:3) = (/-0.0166666666666666D0, 0.15D0, -0.75D0, 0.0D0, 0.75D0, &
                               -0.15D0, 0.016666666666666666D0/)
     case (4)
        coefficients(-4:4) = (/0.00357142857143D0, -0.03809523809524D0, 0.2D0, -0.8D0, 0.0D0, &
             0.8D0, -0.2D0, 0.03809523809524D0, -0.00357142857143D0/)
     case (5)
        coefficients(-5:5) = (/-0.00079365079365D0, 0.00992063492063D0, -0.05952380952381D0, &
                                0.23809523809524D0, -0.8333333333333333D0, 0.0D0, 0.8333333333333333D0, &
                               -0.23809523809524D0, 0.05952380952381D0, -0.00992063492063D0, 0.00079365079365D0/)
     case (6) 
        coefficients(-6:6) = (/0.00018037518038D0, -0.00259740259740D0, 0.01785714285714D0, &
                              -0.07936507936508D0, 0.26785714285714D0, -0.85714285714286D0, 0.0D0, &
                               0.85714285714286D0, -0.26785714285714D0, 0.07936507936508D0, &
                              -0.01785714285714D0, 0.00259740259740D0, -0.00018037518038D0/)
     case default
        
        call errore('xc_vdW_DF', 'Order of numerical gradient not implemented', 2)
        
     end select
     
  end if
     
  gradient_coefficients => coefficients  

  
end function gradient_coefficients


!! ###############################################################################################################


!! ###############################################################################################################
!!                                                 |                  |
!!                                                 |  GET_3D_INDICES  |
!!                                                 |__________________|

!! This routine builds a rank 3 array that holds the indices into the FFT grid for a point with a given
!! set of x, y, and z indices.  The array holds an extra 2N points in each dimension (N to the left and N
!! to the right) so the code can find the neighbors of edge points easily.  This is done by just copying the
!! first N points in each dimension to the end of that dimension and the end N points to the beginning.


function get_3d_indices(N)
  
  USE fft_base,            ONLY : dfftp
   

  integer, intent(in), optional :: N                     !! The number of neighbors in each direction that will
  !                                                      !! be used for the gradient formula.  If not supplied,
  !                                                      !! the code just returns the pointer to the already
  !                                                      !! allocated rho_3d array.

  real(dp) :: dx, dy, dz                                 !! 
  integer :: ix1, ix2, ix3, i_grid                       !! Index variables

  integer, allocatable, target, save :: rho_3d(:,:,:)    !! The local array that will store the indices.  Only a pointer
  !                                                      !! to this array will be returned.

  integer, pointer :: get_3d_indices(:,:,:)              !! The returned pointer to the rho_3d array of indices.
  

  
  !! If the routine has not already been run we set up the rho_3d array by looping over it
  !! and assigning indices to its elements.  If this routine has already been run we simply
  !! return a pointer to the existing array.
  !! --------------------------------------------------------------------------------
  
  if (.not. allocated(rho_3d)) then
     
     ! Check to make sure we have been given the number of neighbors since the routine has
     ! not been run yet.
     ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     
     if (.not. present(N)) then
        
        call errore('get_3d_rho','Number of neighbors for numerical derivatives &
&             must be specified',2)
        
     end if

     ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     
     allocate( rho_3d(-N+1:dfftp%nr1x+N, -N+1:dfftp%nr2x+N, -N+1:dfftp%nr3x+N) )
     
     i_grid = 0
     
     do ix3 = 1, dfftp%nr3x
        do ix2 = 1, dfftp%nr2x
           do ix1 = 1, dfftp%nr1x
              
              i_grid = i_grid + 1
              
              rho_3d(ix1, ix2, ix3) = i_grid
              
           end do
        end do
     end do
     
     
     
     ! Apply periodic boundary conditions to extend the array by N places in each
     ! direction
     ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

     rho_3d(-N+1:0,:,:) = rho_3d(dfftp%nr1x-N+1:dfftp%nr1x, :, :)
     rho_3d(:,-N+1:0,:) = rho_3d(:, dfftp%nr2x-N+1:dfftp%nr2x, :)
     rho_3d(:,:,-N+1:0) = rho_3d(:, :, dfftp%nr3x-N+1:dfftp%nr3x)
     
     rho_3d(dfftp%nr1x+1:dfftp%nr1x+N, :, :) = rho_3d(1:N, :, :)
     rho_3d(:, dfftp%nr2x+1:dfftp%nr2x+N, :) = rho_3d(:, 1:N, :)
     rho_3d(:, :, dfftp%nr3x+1:dfftp%nr3x+N) = rho_3d(:, :, 1:N)
     
     ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end if

  !! ------------------------------------------------------------------------------------------

  
  !! Return the point to rho_3d
  get_3d_indices => rho_3d
  
  
end function  get_3d_indices


!! ###############################################################################################################







!! ###############################################################################################################
!!                                           |                     |
!!                                           |  INVERT_3X3_MATRIX  |
!!                                           |_____________________|

!! This routine is just a hard-wired subroutine to invert a 3x3 matrix.  It is used to invert the matrix of
!! unit cell basis vectors to find the gradient and the derivative of the gradient with respect to the
!! density.

subroutine invert_3x3_matrix(M) 
  
  real(dp), intent(inout) :: M(3,3)     !! On input, the 3x3 matrix to be inverted
  !                                     !! On output, the inverse of the 3x3 matrix given
  
  real(dp) :: temp(3,3)                 !! Temporary storage

  real(dp) :: determinant_M             !! The determinant of the input 3x3 matrix


  temp = 0.0D0

  temp(1,1) = M(2,2)*M(3,3) - M(2,3)*M(3,2)
  temp(1,2) = M(1,3)*M(3,2) - M(1,2)*M(3,3)
  temp(1,3) = M(1,2)*M(2,3) - M(1,3)*M(2,2)
  temp(2,1) = M(2,3)*M(3,1) - M(2,1)*M(3,3)
  temp(2,2) = M(1,1)*M(3,3) - M(1,3)*M(3,1)
  temp(2,3) = M(1,3)*M(2,1) - M(1,1)*M(2,3)
  temp(3,1) = M(2,1)*M(3,2) - M(2,2)*M(3,1)
  temp(3,2) = M(1,2)*M(3,1) - M(1,1)*M(3,2)
  temp(3,3) = M(1,1)*M(2,2) - M(1,2)*M(2,1)

  determinant_M = M(1,1) * (M(2,2)*M(3,3) - M(2,3)*M(3,2)) &
       - M(1,2) * (M(2,1)*M(3,3) - M(2,3)*M(3,1)) &
       + M(1,3) * (M(2,1)*M(3,2) - M(2,2)*M(3,1))

  if (abs(determinant_M) > 1e-6) then
     
     M = 1.0D0/determinant_M*temp

  else

     call errore('invert_3x3_matrix','Matrix is close to singular',1)

  end if

end subroutine invert_3x3_matrix

SUBROUTINE print_sigma(sigma, title)
  
  real(dp), intent(in) :: sigma(:,:)
  character(len=*), intent(in) :: title
  integer :: l

  WRITE( stdout, '(10x,A)') TRIM(title)//" stress"
  WRITE( stdout, '(10x,3F13.8)')  sigma(1,1), sigma(1,2), sigma(1,3)
  WRITE( stdout, '(10x,3F13.8)')  sigma(2,1), sigma(2,2), sigma(2,3)
  WRITE( stdout, '(10x,3F13.8)')  sigma(3,1), sigma(3,2), sigma(3,3)
  WRITE( stdout, '(10x)')

END SUBROUTINE print_sigma  

!! ###############################################################################################################





END MODULE vdW_DF

