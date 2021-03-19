!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE rVV10 
  !--------------------------------------------------------------------------
  !! This module is modeled after the vdW-DF implementation in
  !! 'Modules/xc_vdW_DF.f90'. See that file for references, explanations, and
  !! many useful comments.

  USE kinds,             ONLY : dp
  USE constants,         ONLY : pi
  USE mp,                ONLY : mp_sum
  USE mp_bands,          ONLY : intra_bgrp_comm
  USE io_global,         ONLY : ionode, stdout
  USE fft_base,          ONLY : dfftp
  USE fft_interfaces,    ONLY : fwfft, invfft 
  USE control_flags,     ONLY : gamma_only, iverbosity

  IMPLICIT NONE
  SAVE

  real(dp), parameter :: epsr      = 1.d-12
  real(dp), parameter :: epsg      = 1.D-10
  integer,  parameter :: Nr_points = 1024
  real(dp), parameter :: r_max     = 100.0D0
  real(dp), parameter :: dr        = r_max/Nr_points
  real(dp), parameter :: dk        = 2.0D0*pi/r_max
  real(dp), parameter :: q_min     = 1.0D-4
  real(dp), parameter :: q_cut     = 0.5D0
  integer,  parameter :: Nqs       = 20
  real(dp), parameter, dimension(Nqs):: q_mesh= (/ q_min, 3.0D-4, 5.893850845618885D-4, 1.008103720396345D-3, &
            1.613958359589310D-3, 2.490584839564653D-3, 3.758997979748929D-3, 5.594297198907115D-3, &
            8.249838297569416D-3, 1.209220822453922D-2, 1.765183095571029D-2, 2.569619042667097D-2, &
            3.733577865542191D-2, 5.417739477463518D-2, 7.854595729872216D-2, 0.113805449932145D0,  &
            0.164823306218807D0 , 0.238642339497217D0 , 0.345452975434964D0 , q_cut /)

  real(dp) :: kernel( 0:Nr_points, Nqs, Nqs ), d2phi_dk2( 0:Nr_points, Nqs, Nqs )

  real(dp) :: b_value = 6.3_DP
  real(dp) :: C_value = 0.0093 

  private  
  public :: xc_rVV10,  &
            interpolate_kernel, &
            initialize_spline_interpolation, &
            rVV10_stress, b_value, &
            q_mesh, Nr_points, r_max, q_min, q_cut, Nqs

CONTAINS

! #################################################################################################
!                                       |             |
!                                       |  xc_rVV10   |
!                                       |_____________|

  SUBROUTINE xc_rVV10(rho_valence, rho_core, nspin, etxc, vtxc, v, b_value_)
  
    !! Calculate exchange-correlation energy and potential for rVV10.

    ! Modules to include
    ! -------------------------------------------------------------------------
    
    use gvect,           ONLY : ngm, g
    USE fft_base,        ONLY : dfftp
    USE cell_base,       ONLY : omega, tpiba
    ! -------------------------------------------------------------------------
    
    real(dp), intent(IN) :: rho_valence(:)
    !! valence charge density
    real(dp), intent(IN) :: rho_core(:)
    !! core charge density
    INTEGER,  INTENT(IN) :: nspin
    !! number of spin components
    real(dp), intent(inout) :: etxc
    !! total XC energy
    real(dp), intent(inout) :: vtxc
    !! total XC potential
    real(dp), intent(inout) :: v(:,:)
    !! XC potential on rho grid
    real(DP),optional,intent(in) :: b_value_
    
    !
    ! Local variables
    ! ----------------------------------------------------------------------------------
    !   
    
    integer :: i_grid, theta_i, i_proc, I      
    real(dp) :: grid_cell_volume                
    
    real(dp), allocatable :: total_rho(:)
    real(dp), allocatable :: gradient_rho(:,:)  

    real(dp), allocatable :: q0(:)                                        
    real(dp), allocatable :: dq0_drho(:) 
    real(dp), allocatable :: dq0_dgradrho(:) 
    complex(dp), allocatable :: thetas(:,:)    
                                                
    real(dp) :: Ec_nl                          
    real(dp), allocatable :: potential(:)                                                 

    logical, save :: first_iteration = .true.  

    real(dp) ::  beta
 
 
    ! ---------------------------------------------------------------------------------------------
    !   Begin calculations
  
    !call errore('xc_rVV10','rVV10 functional not implemented for spin polarized runs', size(rho_valence,2)-1)
    if (nspin>2) call errore('xc_vdW_DF','vdW functional not implemented for nspin > 2', nspin)

    if(present(b_value_)) b_value = b_value_
    
    ! --------------------------------------------------------------------------------------------------------

    call start_clock( 'rVV10' )

    beta = 0.0625d0 * (3.0d0 / (b_value**2.0D0) )**(0.75d0)


    ! Write parameters during the first iteratio
    !
    if (first_iteration) then

       first_iteration = .false.
      
       CALL generate_kernel

       if (ionode .and. iverbosity > -1 ) then

          WRITE(stdout,'(/ /A )') "---------------------------------------------------------------------------------"
          WRITE(stdout,'(A)') "Carrying out rVV10 run using the following parameters:"
          WRITE(stdout,'(A,I6,A,I6,A,F8.3)') "Nqs =  ",Nqs, "    Nr_points =  ", Nr_points,"   r_max =  ",r_max
          WRITE(stdout, '(A, F8.5, A, F8.5 )') "b_value = ", b_value, "    beta = ", beta        
          WRITE(stdout,'(5X,"q_mesh =",4F12.8)') (q_mesh(I), I=1, 4)
          WRITE(stdout,'(13X,4F12.8)') (q_mesh(I), I=5, Nqs)
                 
          WRITE(stdout,'(/ A )') "Gradients computed in Reciprocal space"
          WRITE(stdout,'(/ A / /)') "---------------------------------------------------------------------------------"

          
       end if
       
    end if

    ! --------------------------------------------------------------------------------------------------
    ! Allocate arrays.   
    ! ---------------------------------------------------------------------------------------

    allocate( q0(dfftp%nnr) )
    allocate( gradient_rho(3,dfftp%nnr) )
    allocate( dq0_drho(dfftp%nnr), dq0_dgradrho(dfftp%nnr) )
    allocate( total_rho(dfftp%nnr) )
   
 
    ! ---------------------------------------------------------------------------------------
    ! Add together the valence and core charge densities to get the total charge density    
    !
    total_rho = rho_valence(:) + rho_core(:)

    ! -------------------------------------------------------------------------
    ! Here we calculate the gradient in reciprocal space using FFT
    ! -------------------------------------------------------------------------
    call fft_gradient_r2r( dfftp, total_rho, g, gradient_rho)

    ! -------------------------------------------------------------------------
    ! Get Q and all the derivatives
    ! -------------------------------------------------------------------------
    CALL get_q0_on_grid(total_rho, gradient_rho, q0, dq0_drho, dq0_dgradrho)

    ! ---------------------------------------------------------------------------------    

    allocate( thetas(dfftp%nnr, Nqs) )
    CALL get_thetas_on_grid(total_rho, q0, thetas)
    
    call start_clock( 'rVV10_energy')

    call vdW_energy(thetas, Ec_nl)
   
    Ec_nl = Ec_nl + beta * SUM(total_rho) * (omega/(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
    etxc = etxc + Ec_nl

    call stop_clock( 'rVV10_energy')

    ! Print stuff if verbose run
    !
    if (iverbosity > 0) then

       call mp_sum(Ec_nl,intra_bgrp_comm)
       if (ionode) write(*,'(/ / A /)') "     ----------------------------------------------------------------"
       if (ionode) write(*,'(A, F22.15 /)') "     Non-local correlation energy =         ", Ec_nl
       if (ionode) write(*,'(A /)') "     ----------------------------------------------------------------"
       
    end if

    ! ----------------------------------------------------------------------------------------
    ! Inverse Fourier transform the u_i(k) to get the u_i(r) 
    !---------------------------------------------------------------------------------------

    call start_clock( 'rVV10_ffts')
    
    do theta_i = 1, Nqs
       CALL invfft('Rho', thetas(:,theta_i), dfftp) 
    end do

    call stop_clock( 'rVV10_ffts')

    ! -------------------------------------------------------------------------

    call start_clock( 'rVV10_v' )

    allocate( potential(dfftp%nnr) )
    call get_potential(q0, dq0_drho, dq0_dgradrho, total_rho, gradient_rho, thetas, potential)
    
    ! -------------------------------------------------------------------------
    ! Add beta
    ! ------------------------------------------------------------------------- 
    potential = potential + beta
    
    v(:,1) = v(:,1) + potential(:)
    if (nspin==2) v(:,2) = v(:,2) + potential(:) 

    call stop_clock( 'rVV10_v' )

    ! -----------------------------------------------------------------------
    ! The integral of rho(r)*potential(r) for the vtxc output variable
    ! --------------------------------------------------------------------

    grid_cell_volume = omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)  
 
    do i_grid = 1, dfftp%nnr
       vtxc = vtxc + grid_cell_volume*rho_valence(i_grid)*potential(i_grid)
    end do

    deallocate(potential)  

    ! ----------------------------------------------------------------------
    
    ! Deallocate all arrays.
    deallocate(q0, gradient_rho, dq0_drho, dq0_dgradrho, total_rho, thetas)  

    call stop_clock('rVV10')

  END SUBROUTINE xc_rVV10 


  ! #################################################################################################
  !                   |                 |
  !                   |  rVV10_STRESS   |
  !                   |_________________|

  SUBROUTINE rVV10_stress (rho_valence, rho_core, nspin, sigma)

      !! Calculate the stress tensor for rVV10.

      USE fft_base,        ONLY : dfftp
      use gvect,           ONLY : ngm, g
      USE cell_base,       ONLY : tpiba

      implicit none

      real(dp), intent(IN) :: rho_valence(:)
      !! valence charge density
      real(dp), intent(IN) :: rho_core(:)
      !! core charge density
      INTEGER,  INTENT(IN) :: nspin
      !! number of spin components
      real(dp), intent(inout) :: sigma(3,3)
      !! stress tensor

      real(dp), allocatable :: gradient_rho(:,:)         !
      real(dp), allocatable :: total_rho(:)              ! Rho values

      real(dp), allocatable :: q0(:)                     !
      real(dp), allocatable :: dq0_drho(:)               ! Q-values
      real(dp), allocatable :: dq0_dgradrho(:)           !

      complex(dp), allocatable :: thetas(:,:)            ! Thetas
      integer :: i_proc, theta_i, l, m

      real(dp)  :: sigma_grad(3,3)
      real(dp)  :: sigma_ker(3,3)

      ! ---------------------------------------------------------------------------------------------
      !   Tests
      ! --------------------------------------------------------------------------------------------------------

      if (nspin>2) call errore('rV10_stress',' rVV10 stress not implemented for nspin > 2', nspin)

      sigma(:,:) = 0.0_DP
      sigma_grad(:,:) = 0.0_DP
      sigma_ker(:,:) = 0.0_DP

      ! ---------------------------------------------------------------------------------------
      ! Allocations
      ! ---------------------------------------------------------------------------------------

      allocate( gradient_rho(3,dfftp%nnr) )
      allocate( total_rho(dfftp%nnr) )
      allocate( q0(dfftp%nnr) )
      allocate( dq0_drho(dfftp%nnr), dq0_dgradrho(dfftp%nnr) )
      allocate( thetas(dfftp%nnr, Nqs) )
 
      ! ---------------------------------------------------------------------------------------
      ! Charge
      ! ---------------------------------------------------------------------------------------
      total_rho = rho_valence(:) + rho_core(:)

      ! -------------------------------------------------------------------------
      ! Here we calculate the gradient in reciprocal space using FFT
      ! -------------------------------------------------------------------------
      call fft_gradient_r2r( dfftp, total_rho, g, gradient_rho)
      
      ! -------------------------------------------------------------------------------------------------------------
      ! Get q0.
      ! ---------------------------------------------------------------------------------

      CALL get_q0_on_grid(total_rho, gradient_rho, q0, dq0_drho, dq0_dgradrho)

      ! ---------------------------------------------------------------------------------
      ! Get thetas in reciprocal space.
      ! ---------------------------------------------------------------------------------

      CALL get_thetas_on_grid(total_rho, q0, thetas)

      ! ---------------------------------------------------------------------------------------
      ! Stress
      ! ---------------------------------------------------------------------------------------
      CALL rVV10_stress_gradient(total_rho, gradient_rho, q0, dq0_drho, &
                                  dq0_dgradrho, thetas, sigma_grad)

      CALL rVV10_stress_kernel(total_rho, q0, thetas, sigma_ker)

      sigma = - (sigma_grad + sigma_ker) 

      do l = 1, 3
         do m = 1, l - 1
            sigma (m, l) = sigma (l, m)
         enddo
      enddo

      deallocate( gradient_rho, total_rho, q0, dq0_drho, dq0_dgradrho, thetas )
 
   END SUBROUTINE rVV10_stress

   ! ###############################################################################################################
   !                             |                          |
   !                             |  rVV10_stress_gradient   |

   SUBROUTINE rVV10_stress_gradient (total_rho, gradient_rho, q0, dq0_drho, &
                                      dq0_dgradrho, thetas, sigma)

      !! Calculate rVV10 stress with gradient correction.
                                      
      !-----------------------------------------------------------------------------------
      ! Modules to include
      ! ----------------------------------------------------------------------------------
      use gvect,                 ONLY : ngm, g, gg, igtongl, &
                                        gl, ngl, gstart
      USE fft_base,              ONLY : dfftp
      USE cell_base,             ONLY : omega, tpiba, alat, at, tpiba2

      ! ----------------------------------------------------------------------------------

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
      real(dp) :: const

      !real(dp)       :: at_inverse(3,3)

      allocate( d2y_dx2(Nqs, Nqs) ) 
      allocate( u_vdW(dfftp%nnr, Nqs) )

      const = 1.0D0 / (3.0D0 * b_value**(3.0D0/2.0D0) * pi**(5.0D0/4.0D0) )
      sigma(:,:) = 0.0_DP
      prefactor = 0.0_DP
      
      ! --------------------------------------------------------------------------------------------------
      ! Get u in k-space.
      ! ---------------------------------------------------------------------------------------------------

      call thetas_to_uk(thetas, u_vdW)

      ! --------------------------------------------------------------------------------------------------
      ! Get u in real space.
      ! ---------------------------------------------------------------------------------------------------

      call start_clock( 'rVV10_ffts')

      do theta_i = 1, Nqs
         CALL invfft('Rho', u_vdW(:,theta_i), dfftp) 
      end do

      call stop_clock( 'rVV10_ffts')

      ! --------------------------------------------------------------------------------------------------
      ! Get the second derivatives for interpolating the P_i
      ! ---------------------------------------------------------------------------------------------------

      call initialize_spline_interpolation(q_mesh, d2y_dx2(:,:))

      ! ---------------------------------------------------------------------------------------------

      i_grid = 0

      ! ----------------------------------------------------------------------------------------------------
      ! Do the real space integration to obtain the stress component
      ! ----------------------------------------------------------------------------------------------------

      do i_grid = 1, dfftp%nnr

                  q_low = 1
                  q_hi = Nqs 

                  !
                  ! Figure out which bin our value of q0 is in in the q_mesh
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
                      
                      if (total_rho(i_grid) > epsr) then

                        prefactor = u_vdW(i_grid,q_i) * const * total_rho(i_grid)**(3.0D0/4.0D0) * dP_dq0 * dq0_dgradrho(i_grid) 
                      
                        do l = 1, 3
                          do m = 1, l
                                        
                              sigma (l, m) = sigma (l, m) -  prefactor * &
                                             (gradient_rho(l,i_grid) * gradient_rho(m,i_grid))
                           enddo
                        enddo
                     endif

                     ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

                 end do
           
      end do

      call mp_sum(  sigma, intra_bgrp_comm )

      call dscal (9, 1.d0 / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3), sigma, 1)

      deallocate( d2y_dx2, u_vdW )

   END SUBROUTINE rVV10_stress_gradient



   ! ###############################################################################################################
   !                      |                        |
   !                      |  rVV10_stress_kernel   |
   !                      |                        |

   SUBROUTINE rVV10_stress_kernel (total_rho, q0, thetas, sigma)

      ! Modules to include
      ! ----------------------------------------------------------------------------------
      use gvect,                 ONLY : ngm, g, gg, igtongl, gl, ngl, gstart 
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

      ! --------------------------------------------------------------------------------------------------
      ! Integration in g-space
      ! ---------------------------------------------------------------------------------------------------

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

                     sigma (l, m) = sigma (l, m) - G_multiplier * 0.5 * &
                                     thetas(dfftp%nl(g_i),q1_i)*dkernel_of_dk(q1_i,q2_i)*conjg(thetas(dfftp%nl(g_i),q2_i))* &
                                     (g (l, g_i) * g (m, g_i) * tpiba2) / g_kernel 
                     end do
                 end do 
             enddo
         end do      

         if (g_i < gstart ) sigma(:,:) = sigma(:,:) / G_multiplier
         
      enddo

      call mp_sum(  sigma, intra_bgrp_comm )
      
      deallocate( dkernel_of_dk )
      
   END SUBROUTINE rVV10_stress_kernel


  ! ###############################################################################################################
  !                                    |                  |
  !                                    |  GET_Q0_ON_GRID  |
  !                                    |__________________|

  SUBROUTINE get_q0_on_grid (total_rho, gradient_rho, q0, dq0_drho, dq0_dgradrho)
    
    USE fft_base,        ONLY : dfftp
    
    real(dp),  intent(IN)    :: total_rho(:), gradient_rho(:,:)   
    real(dp),  intent(OUT) :: q0(:), dq0_drho(:), dq0_dgradrho(:)

    integer,    parameter      :: m_cut = 12

    real(dp) :: dw0_dn, dk_dn, gmod2
    real(dp) :: mod_grad, wp2, wg2, w0, k   
    real(dp) :: q, exponent, dq0_dq                
    integer  :: i_grid, index, count=0 
 
    ! initialize q0-related arrays ... 
    q0(:) = q_cut
    dq0_drho(:) = 0.0_DP
    dq0_dgradrho(:) = 0.0_DP
  
    do i_grid = 1, dfftp%nnr
          
     if (total_rho(i_grid) > epsr) then

      gmod2 = gradient_rho(1,i_grid)**2 + &
              gradient_rho(2,i_grid)**2 + &
              gradient_rho(3,i_grid)**2
 
       ! Calculate some intermediate values needed to find q
       ! ------------------------------------------------------------------------------------
       mod_grad = sqrt(gmod2)

       wp2= 16.0_dp*pi*total_rho(i_grid)
       wg2 = 4.0_dp*C_value * (mod_grad/total_rho(i_grid))**4

       k = b_value*3.0_dp*pi* ((total_rho(i_grid)/(9.0_dp*pi))**(1.0_dp/6.0_dp))
       w0 = sqrt( wg2 + wp2/3.0_dp  )

       q = w0 / k 
     
       ! Here, we calculate q0 by saturating q according 
       ! ---------------------------------------------------------------------------------------

       exponent = 0.0_dp
       dq0_dq = 0.0_dp
     
       do index = 1, m_cut
        
          exponent = exponent + ( (q/q_cut)**index)/index
          dq0_dq = dq0_dq + ( (q/q_cut)**(index-1))
        
       end do

       q0(i_grid) = q_cut*(1.0_dp - exp(-exponent))
       dq0_dq = dq0_dq * exp(-exponent)

       ! ---------------------------------------------------------------------------------------

       if (q0(i_grid) < q_min) then
         q0(i_grid) = q_min
       end if

       !---------------------------------Final values---------------------------------  

       dw0_dn = 1.0_dp/(2.0_dp*w0) * (16.0_dp/3.0_dp*pi - 4.0_dp*wg2 / total_rho(i_grid) )                
       dk_dn = k / ( 6.0_dp * total_rho(i_grid) )

       dq0_drho(i_grid) = dq0_dq / (k**2) * (dw0_dn * k - dk_dn * w0 )   
       IF ( gmod2 > epsr) THEN
          dq0_dgradrho(i_grid) = dq0_dq  / ( 2.0_dp*k*w0 ) * 4.0_dp*wg2 / (mod_grad**2)
       ELSE
          dq0_dgradrho(i_grid) = 0.0_dp
       ENDIF

     endif  
 
    end do

  end SUBROUTINE get_q0_on_grid


! ###############################################################################################################
!                                      |                      |
!                                      |  GET_THETAS_ON_GRID  |


  SUBROUTINE get_thetas_on_grid (total_rho, q0_on_grid, thetas)

    real(dp), intent(in) :: total_rho(:), q0_on_grid(:)  

    complex(dp), intent(inout):: thetas(:,:)         

    integer :: i_grid, Ngrid_points                
    integer :: theta_i                        
  
    Ngrid_points = size(q0_on_grid)
  
    ! Interpolate the P_i polynomials 
    CALL spline_interpolation(q_mesh, q0_on_grid, thetas)
  
    ! Form the thetas where theta is defined as rho*p_i(q0)
    ! ------------------------------------------------------------------------------------

    do i_grid = 1, Ngrid_points
  
      if (total_rho(i_grid) > epsr ) then 
       thetas(i_grid,:) = thetas(i_grid,:) * (1.0 / (3.0 * sqrt(pi) &
         * ( b_value**(3.0/2.0) ) ) ) * (total_rho(i_grid) / pi)**(3.0/4.0)
      else
       thetas(i_grid,:) = 0.0d0
      endif
     
    end do

    ! ------------------------------------------------------------------------------------
  
    call start_clock( 'rVV10_ffts')

    do theta_i = 1, Nqs

     CALL fwfft ('Rho', thetas(:,theta_i), dfftp)
    end do

    call stop_clock( 'rVV10_ffts')
  
  END SUBROUTINE get_thetas_on_grid


! ###############################################################################################################
!                                     |                        | 
!                                     |  SPLINE_INTERPOLATION  |
!                                     |________________________|


SUBROUTINE spline_interpolation (x, evaluation_points, values)
  
  real(dp), intent(in) :: x(:)
  !! The x values used to form the interpolation
  real(dp), intent(in) :: evaluation_points(:)
  !! (q_mesh in this case) and the values of q0 for which we are 
  !! interpolating the function.
  complex(dp), intent(inout) :: values(:,:)
  !! An output array (allocated outside this routine) that stores the
  !! interpolated values of the P_i (SOLER equation 3) polynomials. The
  !! format is values(grid_point, P_i).
  !
  integer :: Ngrid_points, Nx                                 ! Total number of grid points to evaluate and input x points
  
  real(dp), allocatable, save :: d2y_dx2(:,:)                 ! The second derivatives required to do the interpolation
  
  integer :: i_grid, lower_bound, upper_bound, index, P_i     ! Some indexing variables
  
  real(dp), allocatable :: y(:)                               ! Temporary variables needed for the interpolation
  real(dp) :: a, b, c, d, dx                                  !
 
 
  Nx = size(x)
  Ngrid_points = size(evaluation_points)
  

  ! Allocate the temporary array
  allocate( y(Nx) )

  ! If this is the first time this routine has been called we need to get the second
  ! derivatives (d2y_dx2) required to perform the interpolations.  So we allocate the
  ! array and call initialize_spline_interpolation to get d2y_dx2.
  ! ------------------------------------------------------------------------------------
  if (.not. allocated(d2y_dx2) ) then

     allocate( d2y_dx2(Nx,Nx) )
     call initialize_spline_interpolation(x, d2y_dx2)
     
  end if

  ! ------------------------------------------------------------------------------------
  
  
  do i_grid=1, Ngrid_points
     
     lower_bound = 1
     upper_bound = Nx
     
     do while ( (upper_bound - lower_bound) > 1 )
        
        index = (upper_bound+lower_bound)/2
        
        if ( evaluation_points(i_grid) > x(index) ) then
           lower_bound = index 
        else
           upper_bound = index
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

  
! ###############################################################################################################
!                                |                                   |
!                                |  INITIALIZE_SPLINE_INTERPOLATION  |
!                                |___________________________________|


SUBROUTINE initialize_spline_interpolation (x, d2y_dx2)

  !! This routine is modeled after an algorithm from "Numerical Recipes in C" by Cambridge
  !! University Press, pages 96-97.  It was adapted for Fortran and for the problem at hand.

  real(dp), intent(in)  :: x(:)
  !! The input abscissa values 
  real(dp), intent(inout) :: d2y_dx2(:,:)
  !! The output array (allocated outside this routine) that holds the second derivatives
  !! required for interpolating the function.

  integer :: Nx, P_i, index                        ! The total number of x points and some indexing
  !                                                ! variables

  real(dp), allocatable :: temp_array(:), y(:)     ! Some temporary arrays required.  y is the array
  !                                                ! that holds the funcion values (all either 0 or 1 here).

  real(dp) :: temp1, temp2                         ! Some temporary variables required


  
  Nx = size(x)
  
  allocate( temp_array(Nx), y(Nx) )

  do P_i=1, Nx
     

     ! In the Soler method, the polynomicals that are interpolated are Kroneker delta funcions
     ! at a particular q point.  So, we set all y values to 0 except the one corresponding to 
     ! the particular function P_i.
     ! ----------------------------------------------------------------------------------------

     y = 0.0D0
     y(P_i) = 1.0D0

     ! ----------------------------------------------------------------------------------------
     
     d2y_dx2(P_i,1) = 0.0D0
     temp_array(1) = 0.0D0
     
     do index = 2, Nx-1
        
        temp1 = (x(index)-x(index-1))/(x(index+1)-x(index-1))
        temp2 = temp1 * d2y_dx2(P_i,index-1) + 2.0D0
        d2y_dx2(P_i,index) = (temp1-1.0D0)/temp2
        temp_array(index) = (y(index+1)-y(index))/(x(index+1)-x(index)) &
             - (y(index)-y(index-1))/(x(index)-x(index-1))
        temp_array(index) = (6.0D0*temp_array(index)/(x(index+1)-x(index-1)) &
             - temp1*temp_array(index-1))/temp2
        
     end do
     
     d2y_dx2(P_i,Nx) = 0.0D0
     
     do index=Nx-1, 1, -1
        
        d2y_dx2(P_i,index) = d2y_dx2(P_i,index) * d2y_dx2(P_i,index+1) + temp_array(index)
        
     end do

  end do

  deallocate( temp_array, y)

end SUBROUTINE initialize_spline_interpolation


! ###############################################################################################################
!                                         |                    |
!                                         | INTERPOLATE_KERNEL |
!                                         |____________________|


subroutine interpolate_kernel(k, kernel_of_k)

  !! This routine is modeled after an algorithm from "Numerical Recipes in C" by Cambridge
  !! University Press, page 97.  Adapted for Fortran and the problem at hand.  This function is used to 
  !! find the Phi_alpha_beta needed for equations 11 and 14 of SOLER.

  real(dp), intent(in) :: k
  !! Input value, the magnitude of the g-vector for the current point.
  
  real(dp), intent(inout) :: kernel_of_k(:,:)
  !! An output array (allocated outside this routine) that holds the interpolated value of
  !! the kernel for each pair of q points (i.e. the phi_alpha_beta of the Soler method.

  integer :: q1_i, q2_i, k_i                    ! Indexing variables
 
  real(dp) :: A, B, C, D                        ! Intermediate values for the interpolation
  

  ! Check to make sure that the kernel table we have is capable of dealing with this
  ! value of k.  If k is larger than Nr_points*2*pi/r_max then we can't perform the 
  ! interpolation.  In that case, a kernel file should be generated with a larger number
  ! of radial points.
  ! -------------------------------------------------------------------------------------

  if ( k >= Nr_points*dk ) then
     
     write(*,'(A,F10.5,A,F10.5)') "k =  ", k, "     k_max =  ",Nr_points*dk
     call errore('interpolate kernel', 'k value requested is out of range',1)
     
  end if

  ! -------------------------------------------------------------------------------------
  
  kernel_of_k = 0.0D0
  
  ! This integer division figures out which bin k is in since the kernel
  ! is set on a uniform grid.
  k_i = int(k/dk)
  
  ! Test to see if we are trying to interpolate a k that is one of the actual
  ! function points we have.  The value is just the value of the function in that
  ! case.
  ! ----------------------------------------------------------------------------------------

  if (mod(k,dk) == 0) then
     
     do q1_i = 1, Nqs
        do q2_i = 1, q1_i
           
           kernel_of_k(q1_i, q2_i) = kernel(k_i,q1_i, q2_i)
           kernel_of_k(q2_i, q1_i) = kernel(k_i,q2_i, q1_i)
           
        end do
     end do
     
     return
     
  end if

  ! ----------------------------------------------------------------------------------------


  ! If we are not on a function point then we carry out the interpolation
  ! ----------------------------------------------------------------------------------------
  
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

  ! ----------------------------------------------------------------------------------------

  
end subroutine interpolate_kernel


! ###############################################################################################################
!                                         |                        |
!                                         | INTERPOLATE_DKERNEL_DK |
!                                         |________________________|



subroutine interpolate_Dkernel_Dk(k, dkernel_of_dk)
  
  implicit none 

  real(dp), intent(in) :: k        
  real(dp), intent(inout) :: dkernel_of_dk(Nqs,Nqs)

  integer :: q1_i, q2_i, k_i  
  real(dp) :: A, B, dAdk, dBdk, dCdk, dDdk    
  

  ! -------------------------------------------------------------------------------------

  if ( k >= Nr_points*dk ) then
     
     write(*,'(A,F10.5,A,F10.5)') "k =  ", k, "     k_max =  ",Nr_points*dk
     call errore('interpolate kernel', 'k value requested is out of range',1)
     
  end if

  ! -------------------------------------------------------------------------------------

  dkernel_of_dk = 0.0D0

  k_i = int(k/dk)

  ! ----------------------------------------------------------------------------------------

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

  ! ----------------------------------------------------------------------------------------

  
end subroutine interpolate_Dkernel_Dk 

! #################################################################################################
!                                          |              |
!                                          | thetas_to_uk |
!                                          |______________|


subroutine thetas_to_uk(thetas, u_vdW)
  
  USE gvect,           ONLY : gg, ngm, igtongl, gl, ngl, gstart
  USE fft_base,        ONLY : dfftp
  USE cell_base,       ONLY : tpiba, omega

  complex(dp), intent(in) :: thetas(:,:)  
  complex(dp), intent(out) :: u_vdW(:,:)

  real(dp), allocatable :: kernel_of_k(:,:)    

  real(dp) :: g
  integer :: last_g, g_i, q1_i, q2_i, count, i_grid   

  complex(dp) :: theta(Nqs)        
  
  ! -------------------------------------------------------------------------------------------------
  
  allocate( kernel_of_k(Nqs, Nqs) )

  u_vdW(:,:) = CMPLX(0.0_DP,0.0_DP, kind=dp)
  
  last_g = -1 

  do g_i = 1, ngm
    
     if ( igtongl(g_i) .ne. last_g) then
        
        g = sqrt(gl(igtongl(g_i))) * tpiba
        call interpolate_kernel(g, kernel_of_k)
        last_g = igtongl(g_i)
        
     end if
     
     theta = thetas(dfftp%nl(g_i),:)
     
     do q2_i = 1, Nqs
        do q1_i = 1, Nqs
           u_vdW(dfftp%nl(g_i),q2_i) = u_vdW(dfftp%nl(g_i),q2_i) + kernel_of_k(q2_i,q1_i)*theta(q1_i)
        end do
     end do

  end do

  if (gamma_only) u_vdW(dfftp%nlm(:),:) = CONJG(u_vdW(dfftp%nl(:),:))
  
  deallocate( kernel_of_k )
     
  ! -----------------------------------------------------------------------------------------------
  
end subroutine thetas_to_uk

! #################################################################################################
!                                              |             |
!                                              | VDW_ENERGY  |
!                                              |_____________|

subroutine vdW_energy(thetas, vdW_xc_energy)
  
  USE gvect,           ONLY : gg, ngm, igtongl, gl, ngl, gstart
  USE fft_base,        ONLY : dfftp
  USE cell_base,       ONLY : tpiba, omega

  complex(dp), intent(inout) :: thetas(:,:)  
  real(dp), intent(out) :: vdW_xc_energy     

  real(dp), allocatable :: kernel_of_k(:,:)   

  real(dp) :: g                      
  integer  :: last_g               
  
  integer :: g_i, q1_i, q2_i, count, i_grid  

  complex(dp) :: theta(Nqs), thetam(Nqs), theta_g(Nqs)  
  real(dp)    :: G0_term, G_multiplier 

  complex(dp), allocatable :: u_vdw(:,:)   

  vdW_xc_energy = 0.0D0
 
  allocate (u_vdW(dfftp%nnr,Nqs))
  u_vdW(:,:) = CMPLX(0.0_DP,0.0_DP, kind=dp)

  allocate( kernel_of_k(Nqs, Nqs) )
  
  !
  ! Here we should use gstart,ngm but all the cases are handeld by conditionals inside the loop
  !
  G_multiplier = 1.0D0
  if (gamma_only) G_multiplier = 2.0D0

  last_g = -1 

  do g_i = 1, ngm
     
     if ( igtongl(g_i) .ne. last_g) then
        
        g = sqrt(gl(igtongl(g_i))) * tpiba
        call interpolate_kernel(g, kernel_of_k)
        last_g = igtongl(g_i)
        
     end if
     
     theta = thetas(dfftp%nl(g_i),:)

     do q2_i = 1, Nqs
        do q1_i = 1, Nqs
           u_vdW(dfftp%nl(g_i),q2_i)  = u_vdW(dfftp%nl(g_i),q2_i) + kernel_of_k(q2_i,q1_i)*theta(q1_i)
        end do
        vdW_xc_energy = vdW_xc_energy + G_multiplier * (u_vdW(dfftp%nl(g_i),q2_i)*conjg(theta(q2_i)))
     end do
     
     if (g_i < gstart ) vdW_xc_energy = vdW_xc_energy / G_multiplier

  end do

  if (gamma_only) u_vdW(dfftp%nlm(:),:) = CONJG(u_vdW(dfftp%nl(:),:))

  ! Final value 
  vdW_xc_energy = 0.5D0 * omega * vdW_xc_energy   
  
  deallocate( kernel_of_k )
  thetas(:,:) = u_vdW(:,:)
  deallocate (u_vdW)
  ! ---------------------------------------------------------------------------------------------------
  
end subroutine vdW_energy


! ###############################################################################################################
!                                             |                 |
!                                             |  GET_POTENTIAL  |
!                                             |_________________|

  subroutine get_potential(q0, dq0_drho, dq0_dgradrho, total_rho, gradient_rho, u_vdW, potential)

    use gvect,               ONLY : g
    USE fft_base,            ONLY : dfftp
    USE cell_base,           ONLY : alat, tpiba

    real(dp), intent(in) ::  q0(:), gradient_rho(:,:)   
    real(dp), intent(in) :: dq0_drho(:), dq0_dgradrho(:)
    real(dp), intent(in) :: total_rho(:)
    complex(dp), intent(in) :: u_vdW(:,:)    
    real(dp), intent(inout) :: potential(:)
  
    real(dp), allocatable, save :: d2y_dx2(:,:)         

    integer :: i_grid, P_i,icar                         
    integer :: q_low, q_hi, q                           
    real(dp) :: dq, a, b, c, d, e, f                    
    real(dp) :: y(Nqs), dP_dq0, P                       
    !                                                   

    real(dp), allocatable ::h_prefactor(:)
    complex(dp), allocatable ::h(:)
    real(dp) :: dtheta_dn, dtheta_dgradn

    real(dp) :: const 


    allocate (h_prefactor(dfftp%nnr),h(dfftp%nnr))

    const = 1.0D0 / (3.0D0 * b_value**(3.0D0/2.0D0) * pi**(5.0D0/4.0D0) )
    potential = 0.0D0
    h_prefactor   = 0.0D0

    ! -------------------------------------------------------------------------------------------
    ! Get the second derivatives of the P_i functions for interpolation
    ! ---------------------------------------------------------------------------------------------

    if (.not. allocated( d2y_dx2) ) then
     
     allocate( d2y_dx2(Nqs, Nqs) )
     call initialize_spline_interpolation(q_mesh, d2y_dx2(:,:))
     
    end if
  
    ! ---------------------------------------------------------------------------------------------
  

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
 

        ! IF THE CHARGE DENSITY IS NEGATIVE WE PUT POTENTIAL = 0, OUTSIDE THE SUBROUTINE WE ADD BETA. 
        if (total_rho(i_grid) > epsr) then

          dtheta_dn = const * (3.0D0/4.0D0) / (total_rho(i_grid)**(1.0D0/4.0D0))  * P + &
                      const * total_rho(i_grid)**(3.0D0/4.0D0) *  dP_dq0 * dq0_drho(i_grid)
          dtheta_dgradn = const * total_rho(i_grid)**(3.0D0/4.0D0) * dP_dq0 * dq0_dgradrho(i_grid)

          potential(i_grid) = potential(i_grid) + u_vdW(i_grid,P_i)* dtheta_dn 

          if (q0(i_grid) .ne. q_mesh(Nqs)) then
            h_prefactor(i_grid) = h_prefactor(i_grid) +  u_vdW(i_grid,P_i)* dtheta_dgradn
          end if

        end if

     end do
    end do

    do icar = 1,3
      h(:) = CMPLX( h_prefactor(:)*gradient_rho(icar,:), 0.0_DP, kind=dp)
      CALL fwfft ('Rho', h, dfftp) 
      h(dfftp%nl(:)) = CMPLX(0.0_DP,1.0_DP,kind=dp)*tpiba*g(icar,:)*h(dfftp%nl(:))
      if (gamma_only) h(dfftp%nlm(:)) = CONJG(h(dfftp%nl(:)))
      CALL invfft ('Rho', h, dfftp) 
      potential(:) = potential(:) - REAL(h(:))
    end do

    ! ------------------------------------------------------------------------------------------------------------------------
    deallocate (h_prefactor,h)

  end subroutine get_potential


! ###############################################################################################################
!                                                 |                   |
!                                                 |  generate_kernel  |
!                                                 |___________________|
SUBROUTINE generate_kernel

  implicit none

  integer  :: q1_i, q2_i, r_i                   ! Indexing variables
  real(dp) :: d1, d2                            ! Intermediate values


  kernel    = 0.0D0
  d2phi_dk2 = 0.0D0

  do q1_i = 1, Nqs
     do q2_i = 1, q1_i

        do r_i = 1, Nr_points
           d1 = q_mesh(q1_i) * (dr * r_i)**2    ! Different definition of d1 and d2 for vv10
           d2 = q_mesh(q2_i) * (dr * r_i)**2    ! Different definition of d1 and d2 for vv10
           kernel(r_i, q1_i, q2_i) = -24.0D0 / ( ( d1+1.0 ) * ( d2+1.0 ) * ( d1+d2+2.0 ) )
        end do

        call radial_fft( kernel(:, q1_i, q2_i) )
        call set_up_splines( kernel(:, q1_i, q2_i), d2phi_dk2(:, q1_i, q2_i) )

        kernel    (:, q2_i, q1_i) = kernel   (:, q1_i, q2_i)
        d2phi_dk2 (:, q2_i, q1_i) = d2phi_dk2(:, q1_i, q2_i)

     end do
  end do

END SUBROUTINE generate_kernel


! ###############################################################################################################
!                                                 |              |
!                                                 |  radial_fft  |
!                                                 |______________|
SUBROUTINE radial_fft(phi)

  REAL(DP), INTENT(INOUT) :: phi(0:Nr_points)
  REAL(DP)                :: phi_k(0:Nr_points)
  INTEGER                 :: k_i, r_i
  REAL(DP)                :: r, k
  
  
  phi_k = 0.0D0
  
  DO r_i = 1, Nr_points
     r        = r_i * dr
     phi_k(0) = phi_k(0) + phi(r_i)*r**2
  END DO
  
  phi_k(0) = phi_k(0) - 0.5D0 * (Nr_points*dr)**2 * phi(Nr_points)
  
  DO k_i = 1, Nr_points
     k = k_i * dk
     DO r_i = 1, Nr_points
        r          = r_i * dr
        phi_k(k_i) = phi_k(k_i) + phi(r_i) * r * SIN(k*r) / k
     END DO
     phi_k(k_i) = phi_k(k_i) - 0.5D0 * phi(Nr_points) * r * SIN(k*r) / k
  END DO
  
  phi = 4.0D0 * pi * phi_k * dr

END SUBROUTINE radial_fft


! ###############################################################################################################
!                                              |                  |
!                                              |  set_up_splines  |
!                                              |__________________|
SUBROUTINE set_up_splines(phi, D2)

  REAL(DP), INTENT(IN)    :: phi(0:Nr_points)
  REAL(DP), INTENT(INOUT) :: D2(0:Nr_points)
  REAL(DP), ALLOCATABLE   :: temp_array(:)
  REAL(DP)                :: temp_1, temp_2
  INTEGER                 :: r_i
  
  
  ALLOCATE( temp_array(0:Nr_points) )
  
  D2         = 0
  temp_array = 0
  
  DO r_i = 1, Nr_points - 1
     temp_1  = DBLE(r_i - (r_i - 1))/DBLE( (r_i + 1) - (r_i - 1) )
     temp_2  = temp_1 * D2(r_i-1) + 2.0D0
     D2(r_i) = (temp_1 - 1.0D0)/temp_2
     temp_array(r_i) = ( phi(r_i+1) - phi(r_i))/DBLE( dk*((r_i+1) - r_i) ) - &
          ( phi(r_i) - phi(r_i-1))/DBLE( dk*(r_i - (r_i-1)) )
     temp_array(r_i) = (6.0D0*temp_array(r_i)/DBLE( dk*((r_i+1) - (r_i-1)) )-&
          temp_1*temp_array(r_i-1))/temp_2
  END DO
  
  D2(Nr_points) = 0.0D0
  DO  r_i = Nr_points-1, 0, -1
     D2(r_i) = D2(r_i)*D2(r_i+1) + temp_array(r_i)
  END DO
  
  DEALLOCATE( temp_array )

END SUBROUTINE set_up_splines


END MODULE rVV10
