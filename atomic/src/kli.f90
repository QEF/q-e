!------------------------------------------------------------------------------
! Instituto Balseiro, Condensed Matter Research Group
!------------------------------------------------------------------------------
!
! MODULE: kli
!
!> @author
!> Miguel Carvajal
!
! DESCRIPTION: 
!> Implementation of the KLI exchange potential aproximation.
!
! REVISION HISTORY:
! 18 May 2017 - Initial Version
! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
!------------------------------------------------------------------------------
module kli

use kinds, ONLY : dp
use ld1_parameters, only:  nwfx ! max number of wavefunctions
use ld1inc, only: psi,& ! all electron wavefuntions
                  rho, &   ! all electron density rho(:,1)(up), rho(:,2) (down)
                  num_wave_functions => nwf, &        
                  get_spin => isw , &
                  grid, & ! radial grid
                  orbital_angular_momentum =>  ll, &
                  oc, &
                  nspin,&
                  sl3, &
                  title

use radial_grids, only  : ndmx, hartree
implicit none

private 
public compute_kli_potential ! routine to compute kli potential

real(dp) :: v_x_hf(ndmx,nwfx) !  V^{HF}_i =  \sum \phi_j(r) \int d^3r' \phi_i(r')\phi_j(r')/ | r - r'|
real(dp) :: ux_kli(ndmx,nwfx) !  u_{xi\sigma} = V^{HF}_i / \psi_i
real(dp) :: potential_s(nwfx)! < i| V^{AFA}_\sigma| i> 
real(dp) :: average_ux_kli(nwfx), &  ! < i| u_{xi\sigma}| i> 
            average_kli_potential(nwfx) ! < i| Vx^{KLI}| i> 
integer :: num_up_wavefunctions, num_down_wavefunctions
real(dp) :: mat_m(nwfx, nwfx, 2)  ! M_ij = \int |\psi_i|^2 |\psi_j|^2 / rho
integer :: num_wf(2) ! number of up and down wavefunctions
real(dp) :: oc2(2)
integer :: idx(nwfx,2) ! mapping between the vector of all wawefunctions and the values depending on the spin
real(dp) :: A(nwfx, nwfx), y(nwfx),  ysol(nwfx), AA(nwfx, nwfx) ! linear problem variables \f$ A x = y (with x = V^{KLI}_X \f$
real(dp) :: slater_potential(ndmx,2) ! \f$V^{AFA}_\sigma \sum_i^{N_\sigma} |\phi_{i\sigma}|^2/\rho_sigma * u_{xi\sigma}\f$

contains
    subroutine compute_num_wf(num_wave_functions, num_wf)
    ! determine the number of up and down wavefunctions
        integer , intent(in) :: num_wave_functions
        integer , intent(out) :: num_wf(2) 

        integer :: i, s ! aux variables
        integer :: nu
        num_wf = 0
       

        do i =  1, num_wave_functions
            s = get_spin(i)
            if (oc(i) > 0) then
                num_wf(s) = num_wf(s) + 1
                idx(num_wf(s) , s) = i !index of the jth wf with spin s in the global array
            endif
        enddo   

        if (nspin == 1) then
            if (num_wf(2) /= 0) stop "error"
        endif
    end subroutine compute_num_wf

    function shell_occupancy(i) result(occup)
        integer :: i
        real(dp) :: occup
        
        occup = oc(i)  !* (nspin / 2.0_dp)
    end function shell_occupancy


    subroutine init_module()
        implicit none
        integer :: i
        num_up_wavefunctions = 0
        num_down_wavefunctions = 0
      
       
    end subroutine init_module


    subroutine compute_mat_m(N,grid_size, rho)

        integer,intent(in) :: N
        
        real(dp),intent(in) :: rho(ndmx,2)
        
        real(dp) :: int_0_inf_dr
        
        ! aux variables
        integer :: i, j, k, ii
        integer :: grid_size (2)
        integer ::  s_i, s_j, idx_i, idx_j
        real (dp) :: func(ndmx) , fact1,fact2
        real (dp) :: retval
        integer :: s ! spin interation index
        integer :: nst ! leading behavior for r -> 0 
        real(dp) :: half
        !M_ij = \int |\phi_i|^2 |\phi_j|^2 / rho

        do s = 1,nspin
            do i = 1, num_wf(s)
                idx_i = idx(i,s)
                s_i = get_spin(idx_i) 
                ! s_i should be equal sd

                do j = 1, num_wf(s)
                    ! psi (i) and psi(j) must have the same spin
                    idx_j = idx(j,s)
                    s_j =  get_spin(idx_j)
                    if( s_i == s_j) then

                        func =   psi(:,1,idx_i)**2 * psi(:,1,idx_j)**2
                        do k = 1, grid_size(s_i)
                            ! if (abs(rho(k,s)) > tiny(1.0_dp)) then
                              func(k)  = func(k) /( rho(k,s))
                            ! else if(abs(func(k)) > tiny(1.0_dp)) then
                                ! print *, "Density to small for k", k
                                ! stop 
                            ! end if
                        enddo
                        ! call print_vec(size(func),func)
                        fact1 = (2 * orbital_angular_momentum(idx_i)  + 1)
                        fact2 =  shell_occupancy(idx_j)

                        nst =   2 * orbital_angular_momentum(idx_i) + 2
                        retval = fact1 * fact2 * int_0_inf_dr(func, &
                                            grid, &
                                            grid_size, &
                                            nst) ! this is the asyntotic behavior for r -> 0 of the integrand

                        ! retval goes to the i,j entry of matrix with spin s
                         mat_m(i,j,s) = retval 
#ifdef DEBUG
                         if (retval < tiny(1.0_dp)) then    
                            print *,"Too small"
                            print *,psi(1:10,1,idx_i)
                            print *,psi(1:10,1,idx_j)
                            print *, fact1
                            print *, fact2
                            stop
                         endif
#endif
                    else
                        stop "We got into trouble here"
                    endif ! the elements does not have the same spin
                enddo ! loop over j
            enddo ! lopp over i
            
#ifdef DEBUG      
            ! small check for non nan
            if( mat_m(1,1,s) /= mat_m(1,1,s)) then
                print *, "We got a problem!!!!"
                print *, psi(1:10,1,1)
                stop
            endif

#endif
            
        enddo ! loop over s
        
        

    end subroutine compute_mat_m

    subroutine print_vec(N, v)
        integer:: N,i
        real(dp) :: v(:)
        do i  = 1, N
            print *, "v[",i,"]=", v(i)
        enddo
    end subroutine print_vec


    subroutine  compute_slater_potential(grid_size)



        integer, intent(in) :: grid_size(2)

        integer :: i,j,s
        real(dp) :: work(ndmx)
        real(dp) :: fact
        real(dp) :: int_0_inf_dr
        integer :: nst, spin_i, spin_j
        
        slater_potential = 0.0 ! initialize to zero
        do s = 1, nspin
            do i = 1, num_wf(s)
                fact =  shell_occupancy(idx(i,s))! ok
                slater_potential(1:grid_size(s),s) = slater_potential(1:grid_size(s),s) +  &
                 psi(1:grid_size(s), 1, idx(i,s)) *  v_x_hf(1:grid_size(s), idx(i,s)) * fact 
            enddo
            slater_potential(1:grid_size(s),s) = slater_potential(1:grid_size(s),s)/ (rho(1:grid_size(s),s))
        enddo

!         do j = 1, num_wave_functions
!             spin_j = get_spin(j)
!             nst = 2 * orbital_angular_momentum(j) + 2
!             fact =  2 * orbital_angular_momentum(j) + 1
!             work = psi(:, 1 ,  j ) * psi(:, 1, j) * slater_potential(:,spin_j)
!             ! work = work * ( shell_occupancy(j) * psi(:, 1, j )) / rho(:,1)
!             potential_s(j) = int_0_inf_dr(work, grid, grid_size, nst)  * fact !sshell_occupancy(j)
! #ifdef DEBUG
!             print *, j, potential_s(j)
! #endif

!         enddo
        
    end subroutine compute_slater_potential


    subroutine compute_average_ux_kli(grid_size)

        ! computes  <i| u_xi| i> = \int dr r^2 R_nl(r) V^{HF} (r)

        integer, intent(in) :: grid_size(2)
        integer :: i, fact
        real(dp) :: work(ndmx)
        real(dp) :: int_0_inf_dr
        integer :: nst
        
        work = 0.0_dp
        do i  = 1, num_wave_functions

            work =  psi(:,1,i) *  v_x_hf(:,i) ! shall it be multiplied by the radial function
            
            nst = 2 * orbital_angular_momentum(i) + 2
            average_ux_kli(i) =  int_0_inf_dr(work, grid, grid_size(get_spin(i)), nst)

        enddo
        ! check for nan values
        call   check_nan(num_wave_functions, average_ux_kli)

    end subroutine compute_average_ux_kli


    subroutine solve_linear_problem(num_wave_functions, average_kli_potential)

        integer,intent(in) :: num_wave_functions
        real(dp), intent(out) :: average_kli_potential(nwfx)

        integer :: i, info, N, s,l
        integer :: ipivot(num_wave_functions)
        integer :: lda = nwfx
        
        average_kli_potential = 0
        do s  = 1, nspin
            N = num_wf(s) - 1
            
            if (N > 0) then
                y = 0 
                A = 0
                A = - mat_m(:,:,s)
                ! A = A + I
                do i = 1, N 
                    if(shell_occupancy(idx(i,s)) == 0) then
                        print *, "error"
                        stop
                    endif
                    A(i,i) = A(i,i) + 1.0_dp
                    y(i) = ( potential_s(idx(i,s)) -  average_ux_kli(idx(i,s))) 
                    ! print *,"y", i, y(i)
                enddo

                ! print *,"A", A(1,1:N)
                ! print *,"M", mat_m(1,1,s)
                ! print *, "V^{S}", potential_s(1)
                ! print *, "U_X", average_ux_kli(1)

            ! print *, A(2,1:num_wave_functions)
            ! print *, A(3,1:num_wave_functions)
        
                ! dim(A) = num_wave_functions * num_wave_functions
                ! solve real matrix Ax = b using blas with double precision
                AA = A ! store original
                call DGETRF(N, N, A,lda,ipivot, info)
                if(info /= 0) then
                    print *, "Failed to factorize matrix"
                    stop
                endif
                call DGETRS('N', N, 1, A, lda, ipivot, y, N, info)
                if (info /= 0) then
                    print *, "Failed to solve linear system"
                    stop
                endif
                !    print *, "info", info
                !    print *, y(1:num_wave_functions)
                do l = 1, N

                   if( shell_occupancy(idx(l,s)) == 0) stop 
                   ysol(idx(l,s)) = y(l) !/shell_occupancy(idx(l,s))
                enddo
                ! print *, "y =", y(1:N)
            endif
            
        enddo
        
       
    end subroutine solve_linear_problem

    subroutine compute_average_kli_potential(kli_pot,grid_size, average_kli_potential)
   !---------------------------------------------------------------------------  
   !> @author 
   !> Miguel Carvajal, Instituto Balseiro
   !
   ! DESCRIPTION: 
   !> Computes the average integral of the KLI potential for each KS wavefunction
   !> @brief
   !> The integral that gets solved is $\int dr r^2 R^2_{nl}(r) V^{KLI}(r) $
   !
   ! REVISION HISTORY:
   ! TODO_12_JUL_2017 - Renamed some variables and documented better - Miguel Carvajal
   !
   !> @param[in] kli_pot the radial value of the KLI potential $V^{KLI}(r)$
   !> @param[out] exchange_potential the radial valur of the exchange potential      
   !---------------------------------------------------------------------------   
    real(dp),intent(in) :: kli_pot(ndmx, 2)
    integer ,intent(in) :: grid_size(2)
    real(dp), intent(out) :: average_kli_potential(:)

    ! interval variables
    ! this makes uses of global variables:
    ! `nspin`,  `num_wf`, `orbital_angular_momentum`

    real(dp) :: work(ndmx)
    real(dp) :: int_0_inf_dr
    integer :: s , i, idx_i, nst

    do s = 1, nspin 
        do i = 1, num_wf(s)
            idx_i = idx(i,s)
            nst = 2 * orbital_angular_momentum(idx_i) + 2 ! compute asymtotic power
            work = psi(:,1,idx_i)**2 * kli_pot(:,s) 
            average_kli_potential(idx_i) = int_0_inf_dr(work,grid, grid_size(s), nst) ! integrate for i
        enddo
    enddo

    end subroutine

   subroutine compute_kli_potential(grid_size,  exchange_potential)
   !---------------------------------------------------------------------------  
   !> @author 
   !> Miguel Carvajal, Instituto Balseiro
   !
   ! DESCRIPTION: 
   !> Computes the exchange potential in the KLI aproximation
   !> @brief
   !> calculation of the exchange potential in the KLI ( Kriegger - Li - Ifrate) 
   !> aproximation 
   !> Paper:  Krieger, J. B., Li, Y. & Iafrate, G. J. Phys. Rev. A 46, 5453–5458 (1992).
   !
   ! REVISION HISTORY:
   ! TODO_dd_mmm_yyyy - TODO_describe_appropriate_changes - TODO_name
   !
   !> @param[in] grid_size the size of the radial  grid    
   !> @param[out] exchange_potential the radial valur of the exchange potential      
   !---------------------------------------------------------------------------          
    integer,  intent(in) ::  grid_size  ! number of grid points
    real(dp), intent(out) :: exchange_potential(ndmx, 2) 
    real(dp) :: work(ndmx,2), fact(ndmx), fact1
    real(dp) :: shift
    integer :: i, j, s, last_index(2), idx_i
    logical :: started = .false. !> @var started initialize the function call once

        if (.not. started) then
         average_kli_potential = 0
            started = .true.
        endif

        ! check for valid input pased to the routine
        last_index = 0
        do j = 1, grid_size
            do i =  1, num_wave_functions
                if(psi(j,1,i) /= psi(j,1,i)) then
                    print *, "We got a problem!!!"
                    print *, "Invalid wavefunction passed"
                    print *,psi(j,1,i)
                    stop    
                endif
               
            enddo
            if( rho(j,1 ) > tiny(1.0_dp) .and. last_index(1) == j-1) last_index(1) = j
            if( rho(j,2 ) > tiny(1.0_dp) .and. last_index(2) == j-1) last_index(2) = j 
 
        enddo

        
        ! compute the N_up and N_down wawefunctions
        ! N_up =  num_wf(1)
        ! N_down =  num_wf(2)
        call compute_num_wf(num_wave_functions, num_wf)

        mat_m = 0.0_dp
        ! compute matrix M_sigma 
        ! M(:,:,1)  Matrix up
        ! M(:,:,2)  Matrix down
        ! call compute_mat_m(num_wave_functions, grid_size,rho)    

        do i = 1, num_wave_functions
            call dvex(i,v_x_hf(:, i))    
        enddo


        ! compute the average slater potential
        call compute_slater_potential(last_index)
        call compute_average_ux_kli(last_index)
        ! solve proble A x = b with A =  I - M
        ! call solve_linear_problem(num_wave_functions, average_kli_potential)   
        ! fill the output array for the potential
        exchange_potential  = 0
        do j = 1,10 ! do it ten times 
            do s = 1,nspin
                work = 0

                do i = 1, num_wf(s) - 1
                    idx_i = idx(i,s)
                    fact = average_kli_potential(idx_i) - average_ux_kli(idx_i)  
                    work(:,s) =  work(:,s) + psi(:,1,idx_i)**2 * fact *  shell_occupancy(idx_i)
                enddo   

                if(num_wf(s) > 0) then ! s can be zero as in the hidrogen case
                    work(1:last_index(s),s) =  (work(1:last_index(s),s))/rho(1:last_index(s),s)
                    exchange_potential(1:last_index(s),s) = (slater_potential(1:last_index(s),s)  +  work(1:last_index(s),s))
                    ! all print_vec(last_index(s), exchange_potential(1:last_index(s),s))
                endif
            enddo
            ! compute the average 
            call compute_average_kli_potential(exchange_potential, last_index, average_kli_potential)
        enddo
        ! save exchange potential to file
     
        
    end subroutine compute_kli_potential
    
    subroutine check_nan(N,arr)
        real(dp),intent(in) :: arr(:)
        integer,intent(in) :: N

        ! internal var
        integer :: i
        do i = 1,N
            if(arr(i) /= arr(i)) then 
                print *, "We got a problem" ! NaN value
                print *, "Found NaN value in array, forcing stop"
                stop
            endif
        end do
    end subroutine check_nan
end module kli
