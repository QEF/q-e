!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!this module contains date which defines grids in time and in frequency

  MODULE  times_gw
    USE kinds, only : DP
   
    TYPE times_freqs
       INTEGER :: grid_time!0=Gauss Legendre 1=Gauss Laguerre
       INTEGER :: grid_freq!0=Gauss Legendre 1=Gauss Laguerre
       INTEGER :: n!number of grid points (total of 2n+1 )
       REAL(kind=DP) :: tau!max time
       REAL(kind=DP) :: omega!max frequency
       REAL(kind=DP), POINTER :: times(:)!time grid
       REAL(kind=DP), POINTER :: weights_time(:)!weights on time
       REAL(kind=DP), POINTER :: freqs(:)!frequency grid
       REAL(kind=DP), POINTER :: weights_freq(:)!weights on frequency
       LOGICAL :: l_fft_timefreq!if true uses fft old-style and not grids  
       LOGICAL :: l_fourier_fit_time!if true fits the tails in time
       LOGICAL :: l_fourier_fit_freq!if true fits the tails in freq
       REAL(kind=DP) :: r_tau!ratio for finding outer time point
       REAL(kind=DP) :: r_omega!ratio for finding outer frequency point
       REAL(kind=DP) :: g_tau!ratio for treating bad cases in time
       REAL(kind=DP) :: g_omega!ration for treating bad cases in frequency
       INTEGER  :: grid_fit!grid for self energy ON FREQUENCY: uses the same as for P,W, 1 equally spaced, 2 GL
       REAL(kind=DP) :: omega_fit!max frequency to be considered
       INTEGER :: n_grid_fit!number of grid points on half-axes
       REAL(kind=DP), POINTER :: freqs_fit(:)!frequency grid fot fit
       INTEGER, POINTER :: whois_freq(:)!correspondence for multipoint integration
       REAL (kind=DP), POINTER :: relative_weight(:)!relative weight for multipoint integration
       !options for grid_freq=5
       INTEGER :: second_grid_n=10!sub spacing for second grid
       INTEGER :: second_grid_i=1!max regular step using the second grid

!variable for second frequency grid (for G)
       LOGICAL :: l_g_grid!if true use a dedicated grid on frequency for G
       REAL(kind=DP) :: omega_g! for G: max frequency
       INTEGER :: n_g!for G grid
       REAL(kind=DP), POINTER :: freqs_g(:)!frequency grid
       REAL(kind=DP), POINTER :: weights_freq_g(:)!weights on frequency
       REAL (kind=DP), POINTER :: relative_weight_g(:)!relative weight for multipoint integration
       INTEGER :: grid_freq_g!for G grid
       INTEGER :: second_grid_n_g!for G grid
       INTEGER :: second_grid_i_g!for G grid
       INTEGER, POINTER :: whois_freq_g(:)! for G grid
       REAL(kind=DP), POINTER :: freqs_eff(:)!effective frequency grid
       REAL(kind=DP), POINTER :: freqs_g_eff(:)!effective frequency grid for G
       INTEGER :: grid_levels!for grids of type 4
    END TYPE times_freqs
   

    CONTAINS

      SUBROUTINE free_memory_times_freqs( tf)
        
        implicit none

        TYPE(times_freqs) :: tf

        if(associated(tf%times)) deallocate(tf%times)
        if(associated(tf%weights_time)) deallocate(tf%weights_time)
        if(associated(tf%freqs)) deallocate(tf%freqs)
        if(associated(tf%weights_freq)) deallocate(tf%weights_freq)
        if(associated(tf%freqs_fit)) deallocate(tf%freqs_fit)
        if(associated(tf%whois_freq)) deallocate(tf%whois_freq)
        if(associated(tf%relative_weight)) deallocate(tf%relative_weight)
        if(associated(tf%freqs_g)) deallocate(tf%freqs_g)
        if(associated(tf%relative_weight_g)) deallocate(tf%relative_weight_g) 
        if(associated(tf%weights_freq_g)) deallocate(tf%weights_freq_g)
        if(associated(tf%freqs_eff)) deallocate(tf%freqs_eff)
        if(associated(tf%freqs_g_eff)) deallocate(tf%freqs_g_eff)

        return
      END SUBROUTINE free_memory_times_freqs

      SUBROUTINE setup_timefreq(tf,options)
!sets up and allocates arrays for grids in time and frequency

        USE input_gw,   ONLY : input_options
        USE io_global,  ONLY : stdout
        USE constants,  ONLY : pi
        
        implicit none

        TYPE(input_options) :: options
        TYPE(times_freqs) :: tf


        REAL(kind=DP), ALLOCATABLE :: x(:),w(:)
        INTEGER :: i,j,k,l,ii,nn
        REAL(kind=DP) :: delta

        tf%n = options%n
        tf%grid_time=options%grid_time
        tf%grid_freq=options%grid_freq
        tf%tau=options%tau
        tf%omega=options%omega
        tf%l_fft_timefreq=options%l_fft_timefreq

!fit options available only for Gauss-Legendre grid
        tf%l_fourier_fit_time=options%l_fourier_fit_time
        tf%l_fourier_fit_freq=options%l_fourier_fit_freq
        tf%r_tau=options%r_tau
        tf%r_omega=options%r_omega
        tf%g_tau=options%g_tau
        tf%g_omega=options%g_omega
!options for grid_freq=5
        tf%second_grid_n=options%second_grid_n
        tf%second_grid_i=options%second_grid_i

!options for grid_freq=6
        tf%grid_levels=options%grid_levels

        write(stdout,*) 'DB1',tf%n
       ! allocate(tf%times(-tf%n:tf%n),tf%weights_time(-tf%n:tf%n))
        write(stdout,*) 'DB2'
        if(tf%grid_freq/=5.and.tf%grid_freq/=6) then
           allocate(tf%freqs(-tf%n:tf%n),tf%weights_freq(-tf%n:tf%n),tf%freqs_eff(-tf%n:tf%n))
           nullify(tf%whois_freq)
           !nullify(tf%relative_weight)
           allocate(tf%relative_weight(-nn:nn))
        else if (tf%grid_freq==5) then
           nn=tf%n+tf%second_grid_n*(1+tf%second_grid_i*2)
           allocate(tf%freqs(-nn:nn),tf%weights_freq(-nn:nn))
           allocate(tf%whois_freq(-nn:nn))
           allocate(tf%relative_weight(-nn:nn))
           allocate(tf%freqs_eff(-nn:nn))
        else!grid of type 6
           if(tf%second_grid_i/=0)then
              nn=tf%n-tf%second_grid_i+(tf%second_grid_i*tf%second_grid_n-tf%second_grid_i)*tf%grid_levels+tf%second_grid_i
           else
              nn=tf%n
           endif
           allocate(tf%freqs(-nn:nn),tf%weights_freq(-nn:nn))
           allocate(tf%whois_freq(-nn:nn))
           allocate(tf%relative_weight(-nn:nn))
           allocate(tf%freqs_eff(-nn:nn))
        endif
        write(stdout,*) 'DB3'
        allocate(x(2*tf%n+1),w(2*tf%n+1))

        x(:)=0.d0
        w(:)=0.d0

!frequency grid

        if(tf%grid_freq==0) then!Gauss Legendre
           if(.not.tf%l_fourier_fit_freq) then
              call legzo(tf%n*2+1,x,w)
              tf%freqs(-tf%n:tf%n)=-x(1:2*tf%n+1)*tf%omega
              tf%weights_freq(-tf%n:tf%n)=w(1:2*tf%n+1)*tf%omega
           else
              call legzo(tf%n*2-1,x,w)
              tf%freqs(-tf%n+1:tf%n-1)=-x(1:2*tf%n-1)*tf%omega
              tf%weights_freq(-tf%n+1:tf%n-1)=w(1:2*tf%n-1)*tf%omega
              tf%freqs(-tf%n)=-tf%r_omega*tf%omega
              tf%freqs(tf%n)=tf%r_omega*tf%omega
              tf%weights_freq(-tf%n)=0.d0
              tf%weights_freq(tf%n) =0.d0
           endif
           tf%freqs_eff(-tf%n:tf%n)=tf%freqs(-tf%n:tf%n)
        else if(tf%grid_freq==1) then!Gaus Laguerre
           call lagzo(tf%n,x,w)
           tf%freqs(1:tf%n)=x(1:tf%n)
           do i=1,tf%n
              tf%freqs(-i)=-tf%freqs(i)
           enddo
           tf%freqs(0)=0.d0
           tf%weights_freq(1:tf%n)=w(1:tf%n)*exp(x(1:tf%n))
           do i=1,tf%n
              tf%weights_freq(-i)=tf%weights_freq(i)
           enddo
           tf%weights_freq(0)=0.d0
           tf%freqs_eff(-tf%n:tf%n)=tf%freqs(-tf%n:tf%n)
        else if(tf%grid_freq==2) then
           call legzo(tf%n,x,w)
           tf%freqs(0)=0.d0
           tf%freqs(1:tf%n)=(1.d0-x(1:tf%n))*tf%omega/2.d0
           tf%freqs(-tf%n:-1)=(-1.d0-x(1:tf%n))*tf%omega/2.d0
           tf%weights_freq(0)=0.d0
           tf%weights_freq(1:tf%n)=w(1:tf%n)*tf%omega/2.d0
           tf%weights_freq(-tf%n:-1)=w(1:tf%n)*tf%omega/2.d0
           tf%freqs_eff(-tf%n:tf%n)=tf%freqs(-tf%n:tf%n)
        else if (tf%grid_freq==3) then
           do i=0,tf%n
              tf%freqs(i)=(tf%omega/dble(tf%n))*dble(i)
              tf%freqs(-i)=-tf%freqs(i)
           enddo
           tf%weights_freq(:)=tf%omega/dble(tf%n)
           tf%weights_freq(0)=tf%omega/dble(tf%n)/2.d0
           tf%freqs_eff(-tf%n:tf%n)=tf%freqs(-tf%n:tf%n)
        else if(tf%grid_freq==4) then
           do i=1,tf%n
              tf%freqs(i)=(tf%omega/dble(tf%n))*dble(i)-(0.5d0*tf%omega/dble(tf%n))
              tf%freqs(-i)=-tf%freqs(i)
           enddo
           tf%freqs(0)=0.d0
           tf%weights_freq(:)=(tf%omega/dble(tf%n))
           tf%weights_freq(0)=0.d0
           tf%freqs_eff(-tf%n:tf%n)=tf%freqs(-tf%n:tf%n)
        else if(tf%grid_freq==5) then

           tf%freqs(0)=0.d0
           tf%relative_weight(0)=0.d0
           tf%whois_freq(0)=0

           ii=1
           do i=1,tf%second_grid_n
              tf%freqs(ii)=(tf%omega/dble(2*tf%second_grid_n*tf%n))*dble(i)-0.5d0*tf%omega/dble(2*tf%second_grid_n*tf%n)
              tf%relative_weight(ii)=1.d0/dble(2*tf%second_grid_n)
              tf%whois_freq(ii)=0
              tf%freqs_eff(ii)=0.d0
              ii=ii+1
           enddo
           do j=1,tf%second_grid_i
              do i=1,tf%second_grid_n
                 tf%freqs(ii)=(tf%omega/dble(2*tf%second_grid_n*tf%n))*&
       &dble(i+tf%second_grid_n+2*tf%second_grid_n*(j-1))-0.5d0*tf%omega/dble(2*tf%second_grid_n*tf%n)
                 tf%relative_weight(ii)=1.d0/dble(2*tf%second_grid_n)
                 tf%whois_freq(ii)=j
                 tf%freqs_eff(ii)=tf%omega/dble(tf%n)*dble(j)
                 ii=ii+1
              enddo
              tf%freqs(ii)=tf%omega/dble(tf%n)*dble(j)
              tf%relative_weight(ii)=0.d0
              tf%whois_freq(ii)=j
              tf%freqs_eff(ii)=tf%omega/dble(tf%n)*dble(j)
              ii=ii+1
              do i=1,tf%second_grid_n
                 tf%freqs(ii)=(tf%omega/dble(2*tf%second_grid_n*tf%n))*dble(i+2*tf%second_grid_n*j)-&
                          &0.5d0*tf%omega/dble(2*tf%second_grid_n*tf%n)
                 tf%relative_weight(ii)=1.d0/dble(2*tf%second_grid_n)
                 tf%whois_freq(ii)=j
                 tf%freqs_eff(ii)=tf%omega/dble(tf%n)*dble(j)
                 ii=ii+1
              enddo
           enddo
           do i=tf%second_grid_i+1,tf%n
              tf%freqs(ii)=tf%omega/dble(tf%n)*dble(i)
              tf%relative_weight(ii)=1.d0
              tf%whois_freq(ii)=i
              tf%freqs_eff(ii)=tf%omega/dble(tf%n)*dble(i)
              ii=ii+1
           enddo
           ii=ii-1
           if(ii/=nn) then
              write(stdout,*) 'ERROR ',nn,ii
              stop
           endif
           do i=1,ii
              tf%freqs(-i)=-tf%freqs(i)
              tf%relative_weight(-i)=tf%relative_weight(i)
              tf%whois_freq(-i)=-tf%whois_freq(i)
              tf%freqs_eff(-i)=-tf%freqs_eff(i)
           enddo
           if(.not.options%l_self_time) then
              tf%weights_freq(:)=tf%omega/dble(tf%n)
           else
              tf%weights_freq(0)=0.d0
              ii=1
              do i=1,tf%second_grid_n
                 tf%weights_freq(ii)=tf%omega/dble(tf%n)/dble(2*tf%second_grid_n)
                 ii=ii+1
              enddo
              do j=1,tf%second_grid_i
                 do i=1,tf%second_grid_n
                    tf%weights_freq(ii)=tf%omega/dble(tf%n)/dble(2*tf%second_grid_n)
                    ii=ii+1
                 enddo
                 tf%weights_freq(ii)=0.d0
                 ii=ii+1
                 do i=1,tf%second_grid_n
                    tf%weights_freq(ii)=tf%omega/dble(tf%n)/dble(2*tf%second_grid_n)
                    ii=ii+1
                 enddo
              enddo
              do i=tf%second_grid_i+1,tf%n
                 tf%weights_freq(ii)=tf%omega/dble(tf%n)
                 ii=ii+1
              enddo
              do i=1,nn
                 tf%weights_freq(-i)=tf%weights_freq(i)
                 tf%freqs(-i)=-tf%freqs(i)
              enddo
           endif
           
        else if(tf%grid_freq==6) then

           tf%freqs(0)=0.d0
           tf%weights_freq(0)=0.d0
           tf%relative_weight(0)=0.d0
           tf%whois_freq(0)=0
           ii=1
           do l=1,tf%grid_levels
              if(l==1) then
                 k=1
              else
                 k=tf%second_grid_i+1
              endif
              do j=k,tf%second_grid_n*tf%second_grid_i
                 delta=(tf%omega/dble(tf%n))/(dble(tf%second_grid_n)**(tf%grid_levels-l+1))
                 tf%freqs(ii)=delta*dble(j)-delta/2.d0
                 tf%weights_freq(ii)=delta
                 ii=ii+1
              enddo
           enddo
           delta=(tf%omega/dble(tf%n))
           if(tf%grid_levels==0) then
              j=1
           else
              j=tf%second_grid_i+1
           endif
           do i=j,tf%n
              tf%freqs(ii)=delta*dble(i)-delta/2.d0
              tf%weights_freq(ii)=delta
              ii=ii+1
           enddo
           ii=ii-1
           if(ii/=nn) then
              write(stdout,*) 'ERROR ',nn,ii
              stop
           endif
           do i=1,nn
              tf%weights_freq(-i)=tf%weights_freq(i)
              tf%freqs(-i)=-tf%freqs(i)
           enddo
           tf%freqs_eff(-nn:nn)=tf%freqs(-nn:nn)
           
        else if(tf%grid_freq==7) then
           do i=1,tf%n
              tf%freqs(i)=tf%omega*tan(pi/2.d0/dble(tf%n+1)*dble(i-1)+pi/4.d0/dble(tf%n+1))
              tf%weights_freq(i)=tf%omega*tan(pi/2.d0/dble(tf%n+1)*dble(i))-tf%omega*tan(pi/2.d0/dble(tf%n+1)*dble(i-1))
              tf%freqs(-i)=-tf%freqs(i)
              tf%weights_freq(-i)=tf%weights_freq(i)
           enddo
           tf%freqs(0)=0.d0
           tf%weights_freq(0)=0.d0
           tf%freqs_eff(-tf%n:tf%n)=tf%freqs(-tf%n:tf%n)

        endif
        deallocate(x,w)

!setup frequency grid for fit

        if(.not.(options%l_self_lanczos .and. options%l_lanczos_conv.and. .not.options%l_self_time)) then
           tf%grid_fit=options%grid_fit
           tf%omega_fit=options%omega_fit
           tf%n_grid_fit=options%n_grid_fit
        else
           tf%grid_fit=1
           tf%omega_fit=tf%omega
           tf%n_grid_fit=tf%n
        endif
        if(tf%grid_fit==0) then
           tf%omega_fit=tf%omega
           tf%n_grid_fit=tf%n
        endif
        allocate(tf%freqs_fit(-tf%n_grid_fit:tf%n_grid_fit))
        if(tf%grid_fit==0) then
           tf%freqs_fit(:)=tf%freqs(:)
        else if(tf%grid_fit==1) then
           do i=-tf%n_grid_fit,tf%n_grid_fit
              tf%freqs_fit(i)=(tf%omega_fit/dble(tf%n_grid_fit))*dble(i)
           enddo
        else if(tf%grid_fit==2) then
           allocate(x(2*tf%n_grid_fit+1),w(2*tf%n_grid_fit+1))
           x(:)=0.d0
           w(:)=0.d0
           write(stdout,*) 'CALL LEGZO', tf%n_grid_fit*2+1
           call legzo(tf%n_grid_fit*2+1,x,w)
           write(stdout,*) 'CALLED LEGZO'
           tf%freqs_fit(-tf%n_grid_fit:tf%n_grid_fit)=-x(1:2*tf%n_grid_fit+1)*tf%omega_fit
           deallocate(x,w)
        endif
           

!IN CASE 5 REDEFINE THE TOTAL NUMBER OF FREQUENCIES:
        if(tf%grid_freq==5.or.tf%grid_freq==6) then
           tf%n=nn
           options%n=nn
        endif

!time grid
        allocate(x(2*tf%n+1),w(2*tf%n+1))

        x(:)=0.d0
        w(:)=0.d0


        allocate(tf%times(-tf%n:tf%n),tf%weights_time(-tf%n:tf%n))
        if(tf%grid_time==0) then!Gauss Legendre
           if(.not.tf%l_fourier_fit_time) then
              call legzo(tf%n*2+1,x,w)
              tf%times(-tf%n:tf%n)=-x(1:2*tf%n+1)*tf%tau
              tf%weights_time(-tf%n:tf%n)=w(1:2*tf%n+1)*tf%tau
           else
              call legzo(tf%n*2+1-2,x,w)
              tf%times(-tf%n+1:tf%n-1)=-x(1:2*tf%n-1)*tf%tau
              tf%weights_time(-tf%n+1:tf%n-1)=w(1:2*tf%n-1)*tf%tau
              tf%times(-tf%n)=-tf%r_tau*tf%tau
              tf%times(tf%n)=tf%r_tau*tf%tau
              tf%weights_time(-tf%n)=0.d0
              tf%weights_time(tf%n)=0.d0
           endif
           do i=-tf%n,tf%n
              write(stdout,*) 'TIME:',i, tf%times(i),tf%weights_time(i)
           enddo
        else if(tf%grid_time==1) then!Gaus Laguerre
           call lagzo(tf%n,x,w)
           tf%times(1:tf%n)=x(1:tf%n)
           do i=1,tf%n
              tf%times(-i)=-tf%times(i)
           enddo
           tf%times(0)=0.d0
           tf%weights_time(1:tf%n)=w(1:tf%n)*exp(x(1:tf%n))
           do i=1,tf%n
              tf%weights_time(-i)=tf%weights_time(i)
           enddo
           tf%weights_time(0)=0.d0
        else if(tf%grid_time==2) then
           call legzo(tf%n,x,w)
           tf%times(0)=0.d0
           tf%times(1:tf%n)=(1.d0-x(1:tf%n))*tf%tau/2.d0
           tf%times(-tf%n:-1)=(-1.d0-x(1:tf%n))*tf%tau/2.d0
           tf%weights_time(0)=0.d0
           tf%weights_time(1:tf%n)=w(1:tf%n)*tf%tau/2.d0
           tf%weights_time(-tf%n:-1)=w(1:tf%n)*tf%tau/2.d0

        else if(tf%grid_time==3) then
           do i=0,tf%n
              tf%times(i)=(tf%tau/dble(tf%n))*dble(i)
              tf%times(-i)=-tf%times(i)
           enddo
           tf%weights_time(:)=tf%tau/dble(tf%n)
         else if(tf%grid_time==4) then
           do i=1,tf%n
             tf%times(i)=tf%tau/dble(tf%n)*dble(i)-(0.5d0*tf%tau/dble(tf%n))
              tf%times(-i)=-tf%times(i)
           enddo
           tf%times(0)=0.d0
           tf%weights_time(:)=(tf%tau/dble(tf%n))
           tf%weights_time(0)=0.d0
        endif

        deallocate(x,w)

!options for G grid
        tf%l_g_grid=options%l_g_grid

        if(tf%l_g_grid) then
           tf%n_g=options%n_g
           tf%grid_freq_g=options%grid_freq_g
           tf%second_grid_n_g=options%second_grid_n_g
           tf%second_grid_i_g=options%second_grid_i_g
           tf%omega_g=options%omega_g

 
           if(tf%grid_freq_g/=5) then
              allocate(tf%freqs_g(-tf%n_g:tf%n_g),tf%weights_freq(-tf%n_g:tf%n_g),tf%freqs_g_eff(-tf%n_g:tf%n_g))
              nullify(tf%whois_freq_g)
              nullify(tf%relative_weight_g)
           else
              nn=tf%n_g+tf%second_grid_n_g*(1+tf%second_grid_i_g*2)
              allocate(tf%freqs_g(-nn:nn),tf%weights_freq_g(-nn:nn))
              allocate(tf%whois_freq_g(-nn:nn))
              allocate(tf%relative_weight_g(-nn:nn))
               allocate(tf%freqs_g_eff(-nn:nn))
           endif

           allocate(x(2*tf%n_g+1),w(2*tf%n_g+1))
           x(:)=0.d0
           w(:)=0.d0
           
           if(tf%grid_freq_g==0) then!Gauss Legendre
              call legzo(tf%n_g*2+1,x,w)
              tf%freqs_g(-tf%n_g:tf%n_g)=-x(1:2*tf%n_g+1)*tf%omega_g
              tf%weights_freq_g(-tf%n_g:tf%n_g)=w(1:2*tf%n_g+1)*tf%omega_g
              tf%freqs_g_eff(-tf%n_g:tf%n_g)=tf%freqs_g(-tf%n_g:tf%n_g)
           else if(tf%grid_freq_g==1) then!Gaus Laguerre
              call lagzo(tf%n_g,x,w)
              tf%freqs_g(1:tf%n_g)=x(1:tf%n_g)
              do i=1,tf%n_g
                 tf%freqs_g(-i)=-tf%freqs_g(i)
              enddo
              tf%freqs_g(0)=0.d0
              tf%weights_freq_g(1:tf%n_g)=w(1:tf%n_g)*exp(x(1:tf%n_g))
              do i=1,tf%n_g
                 tf%weights_freq_g(-i)=tf%weights_freq_g(i)
              enddo
              tf%weights_freq_g(0)=0.d0
              tf%freqs_g_eff(-tf%n_g:tf%n_g)=tf%freqs_g(-tf%n_g:tf%n_g)
           else if(tf%grid_freq_g==2) then
              call legzo(tf%n_g,x,w)
              tf%freqs_g(0)=0.d0
              tf%freqs_g(1:tf%n_g)=(1.d0-x(1:tf%n_g))*tf%omega_g/2.d0
              tf%freqs_g(-tf%n_g:-1)=(-1.d0-x(1:tf%n_g))*tf%omega_g/2.d0
              tf%weights_freq_g(0)=0.d0
              tf%weights_freq_g(1:tf%n_g)=w(1:tf%n_g)*tf%omega_g/2.d0
              tf%weights_freq_g(-tf%n_g:1)=w(1:tf%n_g)*tf%omega_g/2.d0
              tf%freqs_g_eff(-tf%n_g:tf%n_g)=tf%freqs_g(-tf%n_g:tf%n_g)
           else if (tf%grid_freq_g==3) then
              do i=0,tf%n_g
                 tf%freqs_g(i)=(tf%omega_g/dble(tf%n_g))*dble(i)
                 tf%freqs_g(-i)=-tf%freqs_g(i)
              enddo
              tf%weights_freq_g(:)=tf%omega_g/dble(tf%n_g)
              tf%freqs_g_eff(-tf%n_g:tf%n_g)=tf%freqs_g(-tf%n_g:tf%n_g)
           else if(tf%grid_freq_g==4) then
              do i=1,tf%n_g
                 tf%freqs_g(i)=(tf%omega_g/dble(tf%n_g))*dble(i)-(0.5d0*tf%omega_g/dble(tf%n_g))
                 tf%freqs_g(-i)=-tf%freqs_g(i)
              enddo
              tf%freqs_g(0)=0.d0
              tf%weights_freq_g(:)=(tf%omega_g/dble(tf%n_g))
              tf%weights_freq_g(0)=0.d0
              tf%freqs_g_eff(-tf%n_g:tf%n_g)=tf%freqs_g(-tf%n_g:tf%n_g)
           else if(tf%grid_freq_g==5) then
              tf%freqs_g(0)=0.d0
              tf%relative_weight_g(0)=0.d0
              tf%whois_freq_g(0)=0
              
              ii=1
              do i=1,tf%second_grid_n_g
                 tf%freqs_g(ii)=(tf%omega_g/dble(2*tf%second_grid_n_g*tf%n_g))*dble(i)-&
                    &0.5d0*tf%omega_g/dble(2*tf%second_grid_n_g*tf%n_g)
                 tf%relative_weight_g(ii)=1.d0/dble(2*tf%second_grid_n_g)
                 tf%whois_freq_g(ii)=0
                 tf%freqs_g_eff(ii)=0.d0
                 ii=ii+1
              enddo
              do j=1,tf%second_grid_i_g
                 do i=1,tf%second_grid_n_g
                    tf%freqs_g(ii)=(tf%omega_g/dble(2*tf%second_grid_n_g*tf%n_g))*&
         &dble(i+tf%second_grid_n_g+2*tf%second_grid_n_g*(j-1))-0.5d0*tf%omega_g/dble(2*tf%second_grid_n_g*tf%n_g)
                    tf%relative_weight_g(ii)=1.d0/dble(2*tf%second_grid_n_g)
                    tf%whois_freq_g(ii)=j
                    tf%freqs_g_eff(ii)=tf%omega_g/dble(tf%n_g)*dble(j)
                    ii=ii+1
                 enddo
                 tf%freqs_g(ii)=tf%omega_g/dble(tf%n_g)*dble(j)
                 tf%relative_weight_g(ii)=0.d0
                 tf%whois_freq_g(ii)=j
                 tf%freqs_g_eff(ii)=tf%omega_g/dble(tf%n_g)*dble(j)
                 ii=ii+1
                 do i=1,tf%second_grid_n_g
                    tf%freqs_g(ii)=(tf%omega_g/dble(2*tf%second_grid_n_g*tf%n_g))*&
     &dble(i+2*tf%second_grid_n_g*j)-0.5d0*tf%omega_g/dble(2*tf%second_grid_n_g*tf%n_g)
                    tf%relative_weight_g(ii)=1.d0/dble(2*tf%second_grid_n_g)
                    tf%whois_freq_g(ii)=j
                    tf%freqs_g_eff(ii)=tf%omega_g/dble(tf%n_g)*dble(j)
                    ii=ii+1
                 enddo
              enddo
              do i=tf%second_grid_i_g+1,tf%n_g
                 tf%freqs_g(ii)=tf%omega_g/dble(tf%n_g)*dble(i)
                 tf%relative_weight_g(ii)=1.d0
                 tf%whois_freq_g(ii)=i
                 tf%freqs_g_eff(ii)=tf%omega_g/dble(tf%n_g)*dble(i)
                 ii=ii+1
              enddo
              ii=ii-1
              if(ii/=nn) then
                 write(stdout,*) 'ERROR ',nn,ii
                 stop
              endif
              do i=1,ii
                 tf%freqs_g(-i)=-tf%freqs_g(i)
                 tf%relative_weight_g(-i)=tf%relative_weight_g(i)
                 tf%whois_freq_g(-i)=-tf%whois_freq_g(i)
                 tf%freqs_g_eff(-i)= tf%freqs_g_eff(i)
              enddo
              if(.not.options%l_self_time) then
                 tf%weights_freq_g(:)=tf%omega_g/dble(tf%n_g)
              else
                 tf%weights_freq_g(0)=0.d0
                 ii=1
                 do i=1,tf%second_grid_n_g
                    tf%weights_freq_g(ii)=tf%omega_g/dble(tf%n_g)/dble(2*tf%second_grid_n_g)
                    ii=ii+1
                 enddo
                 do j=1,tf%second_grid_i_g
                    do i=1,tf%second_grid_n_g
                       tf%weights_freq_g(ii)=tf%omega_g/dble(tf%n_g)/dble(2*tf%second_grid_n_g)
                       ii=ii+1
                    enddo
                    tf%weights_freq_g(ii)=0.d0
                    ii=ii+1
                    do i=1,tf%second_grid_n_g
                       tf%weights_freq_g(ii)=tf%omega_g/dble(tf%n_g)/dble(2*tf%second_grid_n_g)
                       ii=ii+1
                    enddo
                 enddo
                 do i=tf%second_grid_i_g+1,tf%n_g
                    tf%weights_freq_g(ii)=tf%omega_g/dble(tf%n_g)
                    ii=ii+1
                 enddo
                 do i=1,nn
                    tf%weights_freq_g(-i)=tf%weights_freq_g(i)
                    tf%freqs_g(-i)=-tf%freqs_g(i)
                 enddo
              endif
              tf%n_g=nn

           endif

           deallocate(x,w)
           
        else
           allocate(tf%freqs_g(-tf%n:tf%n),tf%weights_freq_g(-tf%n:tf%n),tf%freqs_g_eff(-tf%n:tf%n))
           allocate(tf%whois_freq_g(-tf%n:tf%n))
           allocate(tf%relative_weight_g(-tf%n:tf%n))
           tf%freqs_g(-tf%n:tf%n)= tf%freqs(-tf%n:tf%n)
           tf%freqs_g_eff(-tf%n:tf%n)= tf%freqs_eff(-tf%n:tf%n)
           tf%weights_freq_g(-tf%n:tf%n)=tf%weights_freq(-tf%n:tf%n)
           tf%relative_weight_g(-tf%n:tf%n)=tf%relative_weight(-tf%n:tf%n)
           tf%omega_g=tf%omega
           tf%n_g=tf%n
           tf%grid_freq_g=tf%grid_freq
           tf%second_grid_n_g=tf%second_grid_n
           tf%second_grid_i_g=tf%second_grid_i

        endif
        write(stdout,*) 'N:', tf%n,tf%n_g
        write(stdout,*) 'Omega:', tf%omega,tf%omega_g
        FLUSH(stdout)
        do i=-tf%n,tf%n
           write(stdout,*)'freq:',i, tf%freqs(i),tf%freqs_g(i)
           write(stdout,*)'weight:',i, tf%weights_freq(i),tf%weights_freq_g(i)
           FLUSH(stdout)   
        enddo

        return

      END SUBROUTINE setup_timefreq

    END MODULE times_gw
