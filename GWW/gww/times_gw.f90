!Program GWW  P.Umari
!this module contains date which defines grids in time and in frequency
#ifdef __STD_F95
#define __ALLOCATABLE pointer
#define __ALLOCATED   associated
#else
#define __ALLOCATABLE allocatable
#define __ALLOCATED   allocated
#endif
  MODULE  times_gw
    USE kinds, only : DP

    TYPE times_freqs
       INTEGER :: grid_time!0=Gauss Legendre 1=Gauss Laguerre
       INTEGER :: grid_freq!0=Gauss Legendre 1=Gauss Laguerre
       INTEGER :: n!number of grid points (total of 2n+1 )
       REAL(kind=DP) :: tau!max time
       REAL(kind=DP) :: omega!max frequency
       REAL(kind=DP), __ALLOCATABLE :: times(:)!time grid
       REAL(kind=DP), __ALLOCATABLE :: weights_time(:)!weights on time
       REAL(kind=DP), __ALLOCATABLE :: freqs(:)!frequency grid
       REAL(kind=DP), __ALLOCATABLE :: weights_freq(:)!weights on frequency
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
       REAL(kind=DP), __ALLOCATABLE :: freqs_fit(:)!frequency grid fot fit
    END TYPE times_freqs


    CONTAINS

      SUBROUTINE free_memory_times_freqs( tf)

        implicit none

        TYPE(times_freqs) :: tf

        if(__ALLOCATED (tf%times)) deallocate(tf%times)
        if(__ALLOCATED (tf%weights_time)) deallocate(tf%weights_time)
        if(__ALLOCATED (tf%freqs)) deallocate(tf%freqs)
        if(__ALLOCATED (tf%weights_freq)) deallocate(tf%weights_freq)
        if(__ALLOCATED (tf%freqs_fit)) deallocate(tf%freqs_fit)

        return
      END SUBROUTINE free_memory_times_freqs

      SUBROUTINE setup_timefreq(tf,options)
!sets up and allocates arrays for grids in time and frequency

        USE input_gw,   ONLY : input_options
        USE io_global,  ONLY : stdout

        implicit none

        TYPE(input_options) :: options
        TYPE(times_freqs) :: tf


        REAL(kind=DP), ALLOCATABLE :: x(:),w(:)
        INTEGER :: i

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


        allocate(tf%times(-tf%n:tf%n),tf%weights_time(-tf%n:tf%n))
        allocate(tf%freqs(-tf%n:tf%n),tf%weights_freq(-tf%n:tf%n))

        allocate(x(2*tf%n+1),w(2*tf%n+1))

        x(:)=0.d0
        w(:)=0.d0
!time grid
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

        endif
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
          else if(tf%grid_freq==2) then
             call legzo(tf%n,x,w)
             tf%freqs(0)=0.d0
             tf%freqs(1:tf%n)=(1.d0-x(1:tf%n))*tf%omega/2.d0
             tf%freqs(-tf%n:-1)=(-1.d0-x(1:tf%n))*tf%omega/2.d0
             tf%weights_freq(0)=0.d0
             tf%weights_freq(1:tf%n)=w(1:tf%n)*tf%omega/2.d0
             tf%weights_freq(-tf%n:-1)=w(1:tf%n)*tf%omega/2.d0
        endif


        deallocate(x,w)

!setup frequency grid for fit

        tf%grid_fit=options%grid_fit
        tf%omega_fit=options%omega_fit
        tf%n_grid_fit=options%n_grid_fit
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



        return

      END SUBROUTINE setup_timefreq

    END MODULE times_gw
