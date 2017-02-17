!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------

MODULE ph_rVV10
  
  USE kinds,             ONLY : dp
  USE constants,         ONLY : pi, e2
  USE kernel_table,      ONLY : q_mesh, Nr_points, Nqs, r_max
  USE mp,                ONLY : mp_bcast, mp_sum, mp_barrier
  USE mp_global,         ONLY : me_pool, nproc_pool, intra_pool_comm, root_pool
  USE io_global,         ONLY : ionode
  USE fft_base,          ONLY : dfftp
  USE fft_interfaces,    ONLY : fwfft, invfft 
  USE control_flags,     ONLY : iverbosity, gamma_only
  USE io_global,         ONLY : stdout
  USE rVV10,             ONLY : b_value, initialize_spline_interpolation, &
                                numerical_gradient, interpolate_kernel
  USE gc_lr,             ONLY : grho
  

  IMPLICIT NONE
  
  !! -------------------------------------------------------------------------
  !! Rho and gradient rhos
  !! -------------------------------------------------------------------------
  real(dp), allocatable :: total_rho(:)
  real(dp), allocatable :: gradient_rho(:,:)
  complex(dp), allocatable :: gradient_drho(:,:)

  !! ------------------------------------------------------------------------- 
  real(dp), allocatable       :: q0(:), q(:)
  real(dp), allocatable       :: dq0_dq(:), d2q0_dq2(:)
  real(dp), allocatable       :: dq_dn_n(:), dn_dq_dn_n_n(:), dq_dgradn_n_gmod(:) 
  real(dp), allocatable       :: dn_dq_dgradn_n_gmod_n(:), dgradn_dq_dgradn_n_gmod_n_gmod(:)
  real(dp), allocatable, save :: d2y_dx2(:,:)
 
  real(DP), parameter :: epsr = 1.0d-6

  private  
  public :: dv_drho_rvv10

CONTAINS

!! #####################################################################################################
!!                                        |                      |
!!                                        |  dv_drho_vdw_test    |
!!                                        |______________________|


subroutine dv_drho_rvv10(rho, drho, nspin, q_point, dv_drho)


    USE gvect,               ONLY : nl, g, nlm, ngm
    USE cell_base,           ONLY : alat, tpiba, omega

    integer, intent(IN)        :: nspin
    real(dp), intent(IN)       :: rho(:,:), q_point(3)
    complex(DP), intent(IN)    :: drho (dfftp%nnr, nspin)
    complex(DP), intent(INOUT) :: dv_drho(dfftp%nnr, nspin)
    
    !!
    complex(dp), allocatable :: delta_v(:)
    integer  :: i_grid
    character(len=70) :: fn
 
    !!
    allocate(delta_v(dfftp%nnr))

    CALL get_delta_v(rho, drho, nspin, q_point, delta_v) 
    
    dv_drho(:,1) = delta_v(:)

    deallocate(delta_v)

end subroutine dv_drho_rvv10

!! ###############################################################################################################
!!                            |                          |
!!                            |  get_thetas_derivatives  |
!!                            |__________________________|

subroutine get_delta_v(rho, drho, nspin, q_point, delta_v)

    USE gvect,               ONLY : nl, g, nlm, ngm
    USE cell_base,           ONLY : alat, tpiba, omega

    integer, intent(IN) :: nspin
    real(dp),    intent(IN) :: rho(:,:), q_point(3)       !
    complex(DP), intent(IN) :: drho (dfftp%nnr, nspin)

    complex(DP), intent(OUT) :: delta_v(dfftp%nnr)

    !! Varables needed for calcualtions

    real(dp)    :: gmod, gmod2
    real(dp)    :: theta, dtheta_dn, dtheta_dgradn, d2theta_dn2, dn_dtheta_dgradn, dgradn_dtheta_dgradn 
    complex(dp) :: gradn_graddeltan

    !! -------------------------------------------------------------------------
    !! Terms for the delta_b part
    !! -------------------------------------------------------------------------        
    real(dp), allocatable :: b1(:,:)
    complex(dp), allocatable :: b2(:,:)    

    !! -------------------------------------------------------------------------
    !! Terms for the delta_h part
    !! -------------------------------------------------------------------------        
    complex(dp) :: h1, h1part2
    complex(dp), allocatable :: h1t(:), h2t(:)

    !! -------------------------------------------------------------------------
    !! For the interpolation
    !! -------------------------------------------------------------------------        
    integer :: q_low, q_hi, qbin
    real(dp) :: dq, a, b, c, d, e, f, temp
    
    !! -------------------------------------------------------------------------
    !! Indexes and 
    !! -------------------------------------------------------------------------        
    integer :: P_i, icar, i_grid, theta_i, i_proc, I

    !! -------------------------------------------------------------------------
    !! Delta u and delta_h
    !! -------------------------------------------------------------------------        

    complex(dp), allocatable :: u(:,:), delta_u(:,:)
    complex(dp), allocatable :: delta_h(:), delta_h1(:), delta_h2(:) , delta_h_aux(:), delta_h1_aux(:), delta_h2_aux(:)

    !! -------------------------------------------------------------------------
    !! Delta u and delta_h
    !! -------------------------------------------------------------------------        
    
    character(len=70) :: fn
    integer :: temp_unit
    
    !! -------------------------------------------------------------------------
    !! Allocations
    !! -------------------------------------------------------------------------

    !! Global variables
    allocate(total_rho(dfftp%nnr) )
    allocate(gradient_rho(dfftp%nnr, 3))
    allocate(gradient_drho(dfftp%nnr, 3))
    allocate(q0(dfftp%nnr), q(dfftp%nnr))    
    allocate(dq0_dq(dfftp%nnr), d2q0_dq2(dfftp%nnr))
    allocate(dq_dn_n(dfftp%nnr), dn_dq_dn_n_n(dfftp%nnr), dq_dgradn_n_gmod(dfftp%nnr))
    allocate(dn_dq_dgradn_n_gmod_n(dfftp%nnr), dgradn_dq_dgradn_n_gmod_n_gmod(dfftp%nnr))

    !! Local variables 
    allocate(b1(dfftp%nnr, Nqs), b2(dfftp%nnr, Nqs))
    allocate(u(dfftp%nnr, Nqs), delta_u(dfftp%nnr, Nqs))

    !! -------------------------------------------------------------------------
    !! Zero all values
    !! -------------------------------------------------------------------------
    theta = 0.0D0
    dtheta_dn = 0.0D0
    dtheta_dgradn = 0.0D0
    d2theta_dn2 = 0.0D0
    dn_dtheta_dgradn = 0.0D0
    dgradn_dtheta_dgradn = 0.0D0
    
    b1(:,:) = 0.0D0
    b2(:,:) = 0.0D0
    u(:,:) = (0.0D0, 0.0D0)
    delta_u(:,:) = (0.0D0, 0.0D0)

    ! Empty the output vector    
    delta_v(:) = (0.0D0, 0.0_DP)
    gradient_drho(:,:) = (0.0D0, 0.0D0)
    
    !! -------------------------------------------------------------------------
    !! Gradients
    !! -------------------------------------------------------------------------     

    total_rho(:) = rho(:,1)
    call numerical_gradient(total_rho,gradient_rho)

    CALL qgradient (q_point, dfftp%nnr, drho(:,1), ngm, g, nl, alat, gradient_drho)

    call fill_q0_extended_on_grid ()

    call mp_barrier(intra_pool_comm)
    
    !! ---------------------------------------------------------------------------------------------
    !!  Initialize spline
    !! ---------------------------------------------------------------------------------------------
    
    if (.not. allocated( d2y_dx2) ) then
        allocate( d2y_dx2(Nqs, Nqs) )
        call initialize_spline_interpolation(q_mesh, d2y_dx2(:,:))
    end if

    !! --------------------------------------------------------------------------------------------- 
    !! Begin integral for the delta_b part
    !!---------------------------------------------------------------------------------------------    	
 
    do i_grid = 1,dfftp%nnr

        gmod2 = gradient_rho(i_grid,1)**2+gradient_rho(i_grid,2)**2+gradient_rho(i_grid,3)**2 
        if (total_rho(i_grid) <= epsr ) cycle
        
        CALL get_abcdef (q0, i_grid, q_hi, q_low, dq, a,b,c,d,e,f ) 
 
        do P_i = 1, Nqs
        
          CALL get_thetas_exentended( q_hi, q_low, dq, a,b,c,d,e,f, P_i, i_grid, &   ! Input
                                      gmod, gradn_graddeltan,                    &   ! Output
                                      theta, dtheta_dn, dtheta_dgradn,           &   ! Output - first derivatives
                                      d2theta_dn2, dn_dtheta_dgradn, dgradn_dtheta_dgradn, .true., total_rho) ! Output - second derivatives
          !!
          !! Terms needed later
          !!         
          b1(i_grid, P_i) = dtheta_dn 
          b2(i_grid, P_i) = d2theta_dn2*(drho(i_grid,1)/total_rho(i_grid)) + &
                          dn_dtheta_dgradn*(gradn_graddeltan/total_rho(i_grid))          
          
          !! I need complex variable
          u(i_grid, P_i) =  CMPLX(theta, 0.0D0, KIND=dp)  

          !! Here gradn_graddeltan IS complex, the cast is automatic
          delta_u(i_grid, P_i) =  dtheta_dn*drho(i_grid,1) +  dtheta_dgradn*gradn_graddeltan

        end do
        
    end do
  
    !! -------------------------------------------------------------------------
    !! Delta u part
    !! -------------------------------------------------------------------------

    CALL get_u_delta_u(u, delta_u, q_point)

    do i_grid = 1,dfftp%nnr
        do P_i = 1, Nqs
            delta_v(i_grid) = delta_v(i_grid) + &
                              delta_u(i_grid, P_i) * b1(i_grid, P_i) + &
                              u(i_grid, P_i) * b2(i_grid, P_i)
        enddo
    enddo

    call mp_barrier(intra_pool_comm)
   
    !! -------------------------------------------------------------------------
    !! Deallocate something
    !! -------------------------------------------------------------------------
 
    deallocate(b1,b2)
    
    allocate(h1t(dfftp%nnr),h2t(dfftp%nnr))
    allocate(delta_h(dfftp%nnr))
    allocate(delta_h1(dfftp%nnr))
    allocate(delta_h2(dfftp%nnr))

    !! --------------------------------------------------------------------------------------------- 
    !! Begin h
    !!---------------------------------------------------------------------------------------------    	
   
    delta_h(:) = 0.0_DP
    delta_h1(:) = 0.0_DP
    delta_h2(:) = 0.0_DP

    h1t(:) = (0.0D0, 0.0D0)
    h2t(:) = (0.0D0, 0.0D0)
 
    do i_grid = 1,dfftp%nnr

        gmod2 = gradient_rho(i_grid,1)**2+gradient_rho(i_grid,2)**2+gradient_rho(i_grid,3)**2 
        if (total_rho(i_grid) <= epsr) cycle
    
        CALL get_abcdef (q0, i_grid, q_hi, q_low, dq, a,b,c,d,e,f )
 
        do P_i = 1, Nqs
          
          CALL get_thetas_exentended( q_hi, q_low, dq, a,b,c,d,e,f, P_i, i_grid, &   ! Input
                                      gmod, gradn_graddeltan,                    &   ! Output
                                      theta, dtheta_dn, dtheta_dgradn,           &   ! Output - first derivatives
                                      d2theta_dn2, dn_dtheta_dgradn, dgradn_dtheta_dgradn, .false., total_rho) ! Output - second derivatives
          !!
          !! Terms nedded later
          !!
          h1part2 = dn_dtheta_dgradn*(drho(i_grid,1)/total_rho(i_grid)) + dgradn_dtheta_dgradn*(gradn_graddeltan/total_rho(i_grid))
          
          h1t(i_grid) = h1t(i_grid) + delta_u(i_grid,P_i)*dtheta_dgradn + u(i_grid,P_i)*h1part2
          h2t(i_grid) = h2t(i_grid) + u(i_grid,P_i)*dtheta_dgradn
        
        end do

    end do

    allocate (delta_h_aux(dfftp%nnr))
    allocate (delta_h1_aux(dfftp%nnr))
    allocate (delta_h2_aux(dfftp%nnr))

    do icar = 1,3
       delta_h(:) = (h1t(:) * gradient_rho(:,icar)+ h2t(:) * gradient_drho(:,icar))

       CALL fwfft ('Dense', delta_h, dfftp) 

       delta_h_aux(:) = (0.0_DP, 0.0_DP)
       delta_h_aux(nl(:)) = CMPLX(0.0_DP,(g(icar,:)+q_point(icar)),kind=DP ) * delta_h(nl(:))
       
       if (gamma_only) delta_h_aux(nlm(:)) = CONJG(delta_h_aux(nl(:)))

       CALL invfft ('Dense', delta_h_aux, dfftp) 

       delta_h_aux(:) = delta_h_aux(:)*tpiba

       delta_v(:) = delta_v(:) - delta_h_aux(:)


    end do

    !! -------------------------------------------------------------------------
    !! Deallocate everything
    !! -------------------------------------------------------------------------

    call mp_barrier(intra_pool_comm)

    call mp_barrier(intra_pool_comm)

    deallocate(total_rho, gradient_drho)
    deallocate(gradient_rho)
    deallocate(q0, q, dq0_dq, d2q0_dq2)
    deallocate(dq_dn_n, dn_dq_dn_n_n, dq_dgradn_n_gmod)
    deallocate(dn_dq_dgradn_n_gmod_n, dgradn_dq_dgradn_n_gmod_n_gmod)
    deallocate(h1t, h2t)
    deallocate(delta_h_aux, delta_h)
    deallocate(delta_h1_aux, delta_h1)
    deallocate(delta_h2_aux, delta_h2)
    deallocate(u, delta_u)
 
end subroutine get_delta_v

  !! ###############################################################################################################
  !!                                    |                  |
  !!                                    |    get_abcdef    |
  !!                                    |                  |
  !! ###############################################################################################################

  SUBROUTINE get_abcdef (q0, i_grid, q_hi, q_low, dq, a,b,c,d,e,f )

    USE kernel_table,    ONLY : q_cut, q_min

    real(dp),  intent(IN)    :: q0(:)
    integer, INTENT(IN)      :: i_grid
    integer, INTENT(OUT)     :: q_hi, q_low
    real(dp),  intent(OUT)   :: a,b,c,d,e,f, dq

    integer  :: qbin

    q_low = 1
    q_hi = Nqs 

    do while ( (q_hi - q_low) > 1)

       qbin = int((q_hi + q_low)/2)
       if (q_mesh(qbin) > q0(i_grid)) then
           q_hi = qbin
       else
           q_low = qbin
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

  END SUBROUTINE get_abcdef

  SUBROUTINE get_thetas_exentended (q_hi, q_low, dq, a,b,c,d,e,f, P_i, i_grid, &
                                    gmod, gradn_graddeltan, &
                                    theta, dtheta_dn, dtheta_dgradn, d2theta_dn2, &
                                    dn_dtheta_dgradn, dgradn_dtheta_dgradn, do_write, total_rho)
  
    integer,  intent(IN)  :: q_low, q_hi, P_i, i_grid
    real(dp),  intent(IN) :: dq, a, b, c, d, e, f
    real(dp),  intent(IN) :: total_rho(dfftp%nnr) 
    logical :: do_write

    real(dp), INTENT(OUT)    :: gmod, theta, dtheta_dn, dtheta_dgradn, d2theta_dn2, dn_dtheta_dgradn, dgradn_dtheta_dgradn
    complex(dp), INTENT(OUT) :: gradn_graddeltan

    real(dp) :: y(Nqs), d2P_dq02, dP_dq0, P, const, rho_34, rho_m14
    character(len=70) :: fn

    y = 0.0D0
    y(P_i) = 1.0D0

    !!
    !! P_alpha and derivatives | Num. Recip. Fortran 2nd Ed. p.108
    !!
    d2P_dq02 = a*d2y_dx2(P_i,q_low) + b*d2y_dx2(P_i,q_hi)
    dP_dq0 = (y(q_hi) - y(q_low))/dq - e*d2y_dx2(P_i,q_low) + f*d2y_dx2(P_i,q_hi)
    P = a*y(q_low) + b*y(q_hi) + c*d2y_dx2(P_i,q_low) + d*d2y_dx2(P_i,q_hi)
          
    !!
    !! Thetas
    !!
    const = 1.0D0 / (3.0D0 * b_value**(3.0D0/2.0D0) * pi**(5.0D0/4.0D0) )
    rho_34 = total_rho(i_grid)**(3.0D0/4.0D0)
    rho_m14 = total_rho(i_grid)**(-1.0D0/4.0D0)

    theta     = const*rho_34*P 
 
    dtheta_dn = const*rho_m14*((3.0D0/4.0D0)*P + dP_dq0*dq0_dq(i_grid)*dq_dn_n(i_grid) )
    
    dtheta_dgradn = const*rho_m14*dP_dq0*dq0_dq(i_grid)*dq_dgradn_n_gmod(i_grid)

    d2theta_dn2 = const*rho_m14*( &
                  -(3.0D0/16.0D0)*P &
                  +(1.0D0/2.0D0)*dP_dq0*dq0_dq(i_grid)*dq_dn_n(i_grid) + &
                  d2P_dq02*(dq0_dq(i_grid)**2)*(dq_dn_n(i_grid)**2) + &
                  dP_dq0*d2q0_dq2(i_grid)*(dq_dn_n(i_grid)**2) + &
                  dP_dq0*dq0_dq(i_grid)*dn_dq_dn_n_n(i_grid) )

    dn_dtheta_dgradn = const*rho_m14*( &
                       -(1.0D0/4.0D0)*dP_dq0*dq0_dq(i_grid)*dq_dgradn_n_gmod(i_grid) + &
                       d2P_dq02*(dq0_dq(i_grid)**2)*dq_dn_n(i_grid)*dq_dgradn_n_gmod(i_grid) + &
                       dP_dq0*d2q0_dq2(i_grid)*dq_dn_n(i_grid)*dq_dgradn_n_gmod(i_grid) + &
                       dP_dq0*dq0_dq(i_grid)*dn_dq_dgradn_n_gmod_n(i_grid) )

    dgradn_dtheta_dgradn = const*rho_m14*( &
                           d2P_dq02*(dq0_dq(i_grid)**2)*(dq_dgradn_n_gmod(i_grid)**2) + &
                           dP_dq0*d2q0_dq2(i_grid)*(dq_dgradn_n_gmod(i_grid)**2) + &
                           dP_dq0*dq0_dq(i_grid)*dgradn_dq_dgradn_n_gmod_n_gmod(i_grid))
    !!
    !! Fractions
    !!
    gmod = sqrt(gradient_rho(i_grid,1)**2+gradient_rho(i_grid,2)**2+gradient_rho(i_grid,3)**2)
    gradn_graddeltan = gradient_rho(i_grid,1)*gradient_drho(i_grid,1) + &
                       gradient_rho(i_grid,2)*gradient_drho(i_grid,2) + &
                       gradient_rho(i_grid,3)*gradient_drho(i_grid,3)

  END SUBROUTINE get_thetas_exentended


!! ###############################################################################################################
  !!                                    |                  |
  !!                                    |  GET_Q0_ON_GRID  |
  !!                                    |__________________|

  !! This routine first calculates the q value defined in (DION equations 11 and 12), then
  !! saturates it according to (SOLER equation 7).  
  
  SUBROUTINE fill_q0_extended_on_grid ()
  !!
  !! more specifically it calcultates the following
  !!
  !!     q0(ir) = q0 as defined above
  !!     dq0_dq(ir) =  d q0 /d q
  !!     dq_drho(ir) = total_rho * d q /d rho 
  !!     dq_dgradrho = total_rho / |gradient_rho| * d q / d |gradient_rho|
  !!
    
  USE kernel_table,    ONLY : q_cut, q_min
    
  !                                                                  
  !                                                                        _
  real(dp),   parameter      :: LDA_A  = 0.031091D0, LDA_a1 = 0.2137D0    !
  real(dp),   parameter      :: LDA_b1 = 7.5957D0  , LDA_b2 = 3.5876D0    ! see J.P. Perdew and Yue Wang, Phys. Rev. B 45, 13244 (1992). 
  real(dp),   parameter      :: LDA_b3 = 1.6382D0  , LDA_b4 = 0.49294D0   !_ 
  real(dp)                   :: Z_ab = -0.8491D0                          !! see DION

  integer,    parameter      :: m_cut = 12                                !! How many terms to include in the sum
  real,    parameter         :: C_value = 0.0093                                !! How many terms to include in the sum
  !                                                                       !! of SOLER equation 7
  
  
  real(dp)                   :: kF, r_s, sqrt_r_s, gc                     !! Intermediate variables needed to get q and q0
  real(dp)                   :: LDA_1, LDA_2, exponent, gmod, gmod2       !!
 
  real(dp)                   :: expTemp1, expTemp2
  real(dp)                   :: k, wp2, wg2, w02, w0, alpha
  
  !                                                                       !! Needed by dq0_drho and dq0_dgradrho by the chain rule.
  
  integer                    :: i_grid, index, count=0                    !! Indexing variables
  character(len=70) :: fn

  ! initialize q0-related arrays ... 
  q0(:) = q_cut
  q = 0.0_DP 
  
  dq0_dq(:)           = 0.0_DP ! 
  d2q0_dq2(:)         = 0.0_DP
  dq_dn_n(:)          = 0.0_DP ! total_rho * d q/d rho
  dn_dq_dn_n_n(:)     = 0.0_DP 
  dq_dgradn_n_gmod(:) = 0.0_DP ! total_rho / |gradient_rho| * d q / d |gradient_rho|
  dn_dq_dgradn_n_gmod_n(:) = 0.0_DP ! 
  dgradn_dq_dgradn_n_gmod_n_gmod(:) = 0.0_DP ! 
  
  do i_grid = 1, dfftp%nnr
          
     !! This prevents numerical problems.  If the charge density is negative (an
     !! unphysical situation), we simply treat it as very small.  In that case,
     !! q0 will be very large and will be saturated.  For a saturated q0 the derivative
     !! dq0_dq will be 0 so we set q0 = q_cut and dq0_drho = dq0_dgradrho = 0 and go on
     !! to the next point.
     !! ------------------------------------------------------------------------------------
     
     if (total_rho(i_grid) < epsr)cycle

     gmod2 = gradient_rho(i_grid,1)**2+gradient_rho(i_grid,2)**2+gradient_rho(i_grid,3)**2
     gmod = sqrt(gmod2)


     !! ------------------------------------------------------------------------------------
     
     !! Calculate some intermediate values needed to find q
     !! ------------------------------------------------------------------------------------

     wp2 = 16.0D0*pi*total_rho(i_grid)
     wg2 = 4.0D0*C_value* (gmod/total_rho(i_grid))**4.0D0
     w02 = wg2 + wp2/3.0D0 
     w0 = sqrt( w02 )
     alpha = wg2/w02

     k = b_value*3.0D0* pi* ((total_rho(i_grid)/(9.0D0*pi))**(1.0D0/6.0D0))

     !! ---------------------------------------------------------------
     
     q(i_grid) = w0 / k 
     
     !! ---------------------------------------------------------------
     
     !! Here, we calculate q0 by saturating q according to equation 7 of SOLER.  Also, we find 
     !! the derivative dq0_dq needed for the derivatives dq0_drho and dq0_dgradrh0 discussed below.
     !! ---------------------------------------------------------------------------------------

     exponent = 0.0D0
     dq0_dq(i_grid) = 0.0D0
     expTemp1 = 0.0D0
     expTemp2 = 0.0D0

     do index = 1, m_cut
        
        exponent = exponent + ( (q(i_grid)/q_cut)**index)/index
        dq0_dq(i_grid) = dq0_dq(i_grid) + ( (q(i_grid)/q_cut)**(index-1))
        
        expTemp1 = expTemp1 + ( (q(i_grid)/q_cut)**(index-1))
        expTemp2 = expTemp2 + ( ((index-1)/q_cut)*(q(i_grid)/q_cut)**(index-2))
        
     end do
     

     q0(i_grid) = q_cut*(1.0D0 - exp(-exponent))
     dq0_dq(i_grid) = dq0_dq(i_grid) * exp(-exponent)
     d2q0_dq2(i_grid) = expTemp2*exp(-exponent) - (expTemp1**2)*(1.0D0/q_cut)*exp(-exponent)
   
     
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
  
     dq_dn_n(i_grid)        = q(i_grid)*(1.0D0/3.0D0 -2.5D0 * alpha)
     
     dn_dq_dn_n_n(i_grid)   = q(i_grid) * (1.0D0/9.0D0 + 65.0D0/6.0D0*alpha - 6.25D0*alpha**2.0D0)
 
     dq_dgradn_n_gmod(i_grid) = 8.0D0*q(i_grid)*C_value*(gmod**2.0D0/total_rho(i_grid)**3.0D0)/w02

     dn_dq_dgradn_n_gmod_n(i_grid) = dq_dgradn_n_gmod(i_grid)*(2.5D0*alpha-11.0D0/3.0D0)

     dgradn_dq_dgradn_n_gmod_n_gmod(i_grid) = 16.0D0*q(i_grid)*C_value/total_rho(i_grid)**2.0D0/w02*(1.0D0-alpha) 

  end do

end SUBROUTINE fill_q0_extended_on_grid

!! #####################################################################################################
!!                                        |           |
!!                                        |  delta_u  |
!!                                        |___________|

subroutine get_u_delta_u(u, delta_u, q_point)
  
  USE gvect,           ONLY : nl, nlm, g, gg, ngm, igtongl, gl, ngl, gstart
  USE cell_base,       ONLY : tpiba, omega

  complex(dp), intent(inout)  :: u(dfftp%nnr,Nqs), delta_u(dfftp%nnr,Nqs)
  real(dp), intent(in)   :: q_point(3)
  
  !!
  !! Valirables
  !!
  real(dp), allocatable :: kernel_of_g(:,:), kernel_of_gq(:,:)   
  complex(dp), allocatable :: temp_u(:,:), temp_delta_u(:,:)
  real(dp) :: gmod, gqmod
  integer :: last_g, g_i, q1_i, q2_i, count, i_grid, final_g    !! Index variables
  
  !! -------------------------------------------------------------------------------------------------
  !! Allocate variables
  !!  
  allocate( kernel_of_g(Nqs, Nqs), kernel_of_gq(Nqs, Nqs) )
  allocate( temp_u(dfftp%nnr, Nqs), temp_delta_u(dfftp%nnr, Nqs) )

  temp_u(:,:) = (0.0D0, 0.0D0)
  temp_delta_u(:,:) = (0.0D0, 0.0D0)
  !!
  !! Get argument in reciprocal space
  !!
  call start_clock( 'vdW_ffts')
  do q1_i = 1, Nqs
     CALL fwfft ('Dense', u(:,q1_i), dfftp)
     CALL fwfft ('Dense', delta_u(:,q1_i), dfftp)
  end do
  call stop_clock( 'vdW_ffts')
  
  !!
  !! Integrate in reciprocal space
  !!
  last_g = -1 

  do g_i = 1, ngm
    
     if ( igtongl(g_i) .ne. last_g) then
        gmod = sqrt(gl(igtongl(g_i))) * tpiba
        call interpolate_kernel(gmod, kernel_of_g)
        last_g = igtongl(g_i)

     end if
    
     gqmod = sqrt( (g(1,g_i)+q_point(1))**2 + (g(2,g_i)+q_point(2))**2 + (g(3,g_i)+q_point(3))**2 )*tpiba
     call interpolate_kernel(gqmod, kernel_of_gq)
 
     !! Loop over alpha
     do q2_i = 1, Nqs
        !! Sum over beta
        do q1_i = 1, Nqs
        
           temp_u(nl(g_i), q2_i) = temp_u(nl(g_i), q2_i) + kernel_of_g(q2_i,q1_i)*u(nl(g_i), q1_i)
       
           temp_delta_u(nl(g_i), q2_i) = temp_delta_u(nl(g_i), q2_i) + &
                                        kernel_of_gq(q2_i,q1_i)*delta_u(nl(g_i), q1_i)
        end do
     end do

  end do

  if (gamma_only) then
    temp_u(nlm(:),:) = CONJG(temp_u(nl(:),:))
    temp_delta_u(nlm(:),:) = CONJG(temp_delta_u(nl(:),:))
  endif

  !!
  !! Put everything in real space
  !!
  call start_clock( 'vdW_ffts')
  do q1_i = 1, Nqs
     CALL invfft ('Dense', temp_u(:,q1_i), dfftp)
     CALL invfft ('Dense', temp_delta_u(:,q1_i), dfftp)
  end do
  call stop_clock( 'vdW_ffts')

  u(:,:) = temp_u(:,:)
  delta_u(:,:) = temp_delta_u(:,:)
    
  deallocate(temp_u, temp_delta_u, kernel_of_g, kernel_of_gq)
     
  !! -----------------------------------------------------------------------------------------------
  
end subroutine get_u_delta_u


!! ###############################################################################################################
!!                  |                      |
!!                  |      qgradient       |
!!                  |______________________|


subroutine qgradient (xq, nrxx, a, ngm, g, nl, alat, ga)
  !--------------------------------------------------------------------
  ! Calculates ga = \grad a in R-space (a is also in R-space)
  use control_flags,  ONLY : gamma_only
  USE fft_base,       ONLY: dfftp
  USE fft_interfaces, ONLY: fwfft, invfft
  USE gvect,     ONLY : nlm
  !gamma_only is disregarded for phonon calculations
  USE kinds, only : DP
  USE constants, ONLY: tpi
  implicit none
  integer :: nrxx, ngm, nl (ngm)
  complex(DP) :: a (nrxx), ga (nrxx, 3)
  real(DP) :: g (3, ngm), alat, xq (3)
  integer :: n, ipol
  real(DP) :: tpiba
  complex(DP), allocatable :: aux (:), gaux (:)

  allocate (gaux(  nrxx))
  allocate (aux (  nrxx))

  tpiba = tpi / alat
  ! bring a(r) to G-space, a(G) ...
  aux (:) = a(:)

  CALL fwfft ('Dense', aux, dfftp)
  ! multiply by i(q+G) to get (\grad_ipol a)(q+G) ...
  do ipol = 1, 3
     gaux (:) = (0.d0, 0.d0)
     do n = 1, ngm
        gaux(nl(n)) = CMPLX(0.d0, xq (ipol) + g (ipol, n),kind=DP) * aux (nl(n))
        if (gamma_only) gaux( nlm(n) ) = conjg( gaux( nl(n) ) )
     enddo
     ! bring back to R-space, (\grad_ipol a)(r) ...

     CALL invfft ('Dense', gaux, dfftp)
     ! ...and add the factor 2\pi/a  missing in the definition of q+G
     do n = 1, nrxx
        ga (n, ipol) = gaux (n) * tpiba
     enddo
  enddo
  deallocate (aux)
  deallocate (gaux)
  return

end subroutine qgradient


END MODULE ph_rVV10 

