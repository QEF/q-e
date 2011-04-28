!
! P.Umari Program GWW
!
  MODULE expansion
!this module conatins descriptions and subroutine for a multipole expansion
!of the self energy
   USE kinds, ONLY : DP


   TYPE self_expansion
!all the parameters for the exapnsion
!the fit is on the POSITIVE imaginary axes
    INTEGER :: max_i !number of states considered
    INTEGER :: i_min!minimum state to be considered
    INTEGER :: i_max!maximum state to be considered
    INTEGER :: n_multipoles!number of multipoles considered
    COMPLEX(kind=DP), DIMENSION(:), POINTER :: a_0!parameters a_0
    COMPLEX(kind=DP), DIMENSION(:,:), POINTER :: a!parameters a (n_multipoles,max_i)
    COMPLEX(kind=DP), DIMENSION(:,:), POINTER :: b!parameters b (n_multipoles,max_i)
  END TYPE self_expansion

  CONTAINS

  SUBROUTINE free_memory_self_expansion(se)
!if allocated deallocates
   implicit none
   TYPE(self_expansion) :: se

   if(associated(se%a_0)) then
     deallocate(se%a_0)
     nullify(se%a_0)
   endif

   if(associated(se%a)) then
     deallocate(se%a)
     nullify(se%a)
   endif

   if(associated(se%b)) then
     deallocate(se%b)
     nullify(se%b)
   endif
   return

  END SUBROUTINE

  SUBROUTINE create_self_energy_fit( tf, se,ss, options)
!this subroutine fit the self energy in the imaginary frequency
!with a multipole complex function
!parallel on states

    USE io_global,  ONLY : stdout
    USE input_gw,   ONLY : input_options
    USE constants,  ONLY : pi
    USE self_energy_storage, ONLY : self_storage
    USE para_gww,   ONLY : is_my_state_range
    USE mp,         ONLY : mp_sum
    USE times_gw,   ONLY : times_freqs

    implicit none

    TYPE(times_freqs), INTENT(in) :: tf!frequency grid
    TYPE(self_expansion) :: se!fit to be created
    TYPE(self_storage) :: ss!data on frequency
    TYPE(input_options) :: options! for number of multipoles

    INTEGER :: ii,jj
    COMPLEX(kind=DP), ALLOCATABLE :: z(:),s(:)
    REAL(kind=DP) :: df,freq, totalperiod, chi, chi0
    INTEGER :: icyc
    COMPLEX(kind=DP) :: a_0_old, a_0_good
    COMPLEX(kind=DP), ALLOCATABLE :: a_old(:), b_old(:), a_good(:), b_good(:)

!sets:
    se%max_i=options%max_i
    se%i_min=options%i_min
    se%i_max=options%i_max
    se%n_multipoles=options%n_multipoles

!allocates:
!    call free_memory_self_expansion(se)
    allocate(se%a_0(se%max_i))
    allocate(se%a(se%n_multipoles,se%max_i))
    allocate(se%b(se%n_multipoles,se%max_i))
    allocate(a_old(se%n_multipoles))
    allocate(b_old(se%n_multipoles))
    allocate(a_good(se%n_multipoles))
    allocate(b_good(se%n_multipoles))


  !  allocate(z(ss%n_grid_fit),s(ss%n_grid_fit))
    allocate(z(options%n_fit),s(options%n_fit))

!allocate and set data arrays
      totalperiod=2.d0*ss%tau+2.d0*ss%tau/real(ss%n)
      df=2.d0*pi/totalperiod

    if(options%offset_fit == 0) then
         do ii=1,options%n_fit-1
            if(tf%l_fft_timefreq) then
               freq=df*real(ii)
            else
               freq=tf%freqs_fit(ii)
            endif
            z(ii+1)=cmplx(0.d0,freq)
         enddo
         z(1)=(0.d0,0.d0)
      else
         do ii=1,options%n_fit
            if(tf%l_fft_timefreq) then
               freq=df*real(ii+options%offset_fit-1)
            else
               freq=tf%freqs_fit(ii+options%offset_fit-1)
            endif
            z(ii)=cmplx(0.d0,freq)
         enddo
      endif

!some checks
    if(ss%ontime .and. tf%grid_fit==0) then
      write(stdout,*) 'Subroutine self_energy_fit: imaginary frequency required'
      stop
    endif

    if(ss%whole_s) then
       write(stdout,*) 'Subroutine self_energy_fit: diagonal terms implemented YET'
       stop
    else
       se%a_0=(0.d0,0.d0)
       se%a(:,:)=(0.d0,0.d0)
       se%b(:,:)=(0.d0,0.d0)


       do ii=se%i_min,se%i_max!loop on states
          chi0=1.d10
!set initial values
          if(is_my_state_range(ii)) then
             if(tf%grid_fit==0) then
                do jj=1,options%n_fit!ss%n
                   s(jj)=ss%diag(ii,jj+ss%n+1)
                enddo
             else
                !do jj=1,ss%n_grid_fit
                do jj=0+options%offset_fit,options%n_fit+options%offset_fit-1!ATTENZIONE
                  s(jj-options%offset_fit+1)=ss%diag_freq_fit(ii,jj+ss%n_grid_fit+1)!ATTENZIONE
                enddo
             endif

             se%a_0(ii)=(0.0,0.0d0)
             do jj=1,options%n_multipoles
                se%a(jj,ii)=cmplx(real(jj)*(0.01d0),0.d0)
                se%b(jj,ii)=cmplx((0.5d0)*real(jj)*(-1.d0)**jj,-0.01d0)
             enddo

             do icyc=1,options%cyc_minpack
                call fit_multipole(options%n_fit,options%n_multipoles,z,s, &
                                   se%a_0(ii),se%a(:,ii),se%b(:,ii),1.d0,  &
                                   options%fit_thres,options%fit_maxiter)
                a_0_old=se%a_0(ii)
                do jj=1,options%n_multipoles
                   a_old(jj)=se%a(jj,ii)
                   b_old(jj)=se%b(jj,ii)
                enddo


                if(options%n_max_minpack /= 0) then
                   write(stdout,*) 'Calling minpack'!ATTENZIONE
                   call fit_multipole_minpack(options%n_fit, &
                        options%n_multipoles,z,s,se%a_0(ii),se%a(:,ii),&
                        se%b(:,ii),options%fit_thres, options%n_max_minpack, chi)
!                   call fit_multipole_verlet2(options%n_fit,options%n_multipoles,z,s, &
!                        se%a_0(ii),se%a(:,ii),se%b(:,ii),options%fit_thres, &
!                        options%n_max_minpack, chi, options%fit_dt, options%fit_frice)

                endif

                if(chi <= chi0) then
                   a_0_good=se%a_0(ii)
                   do jj=1,options%n_multipoles
                      a_good(jj)=se%a(jj,ii)
                      b_good(jj)=se%b(jj,ii)
                   enddo
                   chi0=chi
                endif

                se%a_0(ii)=a_0_old
                do jj=1,options%n_multipoles
                   se%a(jj,ii)=a_old(jj)
                   se%b(jj,ii)=b_old(jj)
                enddo

             enddo
             se%a_0(ii)=a_0_good
             do jj=1,options%n_multipoles
                se%a(jj,ii)=a_good(jj)
                se%b(jj,ii)=b_good(jj)
             enddo

            write(stdout,*) 'FIT state :', ii
            write(stdout,*) 'FIT    a_0:', se%a_0(ii)
            do jj=1,options%n_multipoles
               write(stdout,*) 'FIT    a:',jj,se%a(jj,ii)
               write(stdout,*) 'FIT    b:',jj,se%b(jj,ii)
            enddo
         endif
      enddo
      call mp_sum(se%a_0(:))
      call mp_sum(se%a(:,:))
      call mp_sum(se%b(:,:))
   endif
   deallocate(z,s)
   deallocate(a_old,b_old)
   deallocate(a_good,b_good)
  return
END SUBROUTINE create_self_energy_fit

  SUBROUTINE func_fit(se,z,i,fz)
!this functions returns the value of the fit at z,
!relative to the i-th parameters

   implicit none
   TYPE(self_expansion) :: se!parameters of fits
   COMPLEX(kind=DP) :: z!where
   INTEGER :: i !which set of parameters
   COMPLEX(kind=DP) :: fz
   COMPLEX(kind=DP) :: num, den

   INTEGER :: jj

   fz=se%a_0(i)
   do jj=1,se%n_multipoles
     fz=fz+se%a(jj,i)/(z-se%b(jj,i))
   enddo
!    num=(1.d0,0.d0)
!    do jj=1,se%n_multipoles
!      num=num*(z-se%a(jj,i))
!    enddo
!    den=(1.d0,0.d0)
!    do jj=1,se%n_multipoles
!      den=den*(z-se%b(jj,i))
!    enddo
!    fz=num/den

   return


  END SUBROUTINE

  SUBROUTINE print_fit_onfile(tf, se,ss)
!this subroutines prints the resulta of the fit on file:
!real and imaginary part on imaginary frequency, with results self_energy
!and real and imaginary part on real frequency
!parallel on states
   USE io_files,             ONLY : find_free_unit
   USE self_energy_storage,  ONLY : self_storage
   USE constants,            ONLY : pi
   USE io_global,            ONLY : ionode
   USE para_gww,             ONLY : is_my_state_range
   USE times_gw,             ONLY : times_freqs

   implicit none

   TYPE(times_freqs), INTENT(in) :: tf!for frequency grid
   TYPE(self_expansion) :: se!parameters of fit
   TYPE(self_storage)   :: ss!self energy data


   INTEGER :: ii,jj
   INTEGER :: iun
   CHARACTER(5) :: nfile
   REAL(kind=DP) :: totalperiod,df,freq
   COMPLEX(kind=DP) :: zz
   COMPLEX(kind=DP) :: fz, gz

   do ii=se%i_min,se%i_max!loop on states
      if(is_my_state_range(ii)) then
!set file name

         write(nfile,'(5i1)') &
              & ii/10000,mod(ii,10000)/1000,mod(ii,1000)/100,mod(ii,100)/10,mod(ii,10)

         totalperiod=2.d0*ss%tau+2.d0*ss%tau/real(ss%n)
         df=2.d0*pi/totalperiod

!now real part on imaginary frequency

!openfile
         iun = find_free_unit()
         open( unit=iun, file='re_on_im'// nfile, status='unknown',form='formatted')

         do jj=-ss%n_grid_fit,ss%n_grid_fit
!allocate and set data arrays
            if(tf%l_fft_timefreq) then
               freq=df*real(jj)
            else
               freq=tf%freqs_fit(jj)
            endif
            zz=cmplx(0.d0,freq)
            call value_on_frequency(se,ii,freq,gz)
            !call func_fit(se,zz,ii,fz)
            call value_on_frequency_complex(se,ii,zz,fz)
            if(tf%grid_fit==0) then
               write(iun,'(4f12.6)') freq, real(fz),real(ss%diag(ii,jj+ss%n)),real( gz)
            else
               write(iun,'(4f12.6)') freq, real(fz),real(ss%diag_freq_fit(ii,jj+ss%n_grid_fit+1)),real( gz)
            endif
         enddo

         close(iun)


!now imaginary  part on imaginary frequency

!openfile
         iun = find_free_unit()
         open( unit=iun, file='im_on_im'// nfile, status='unknown',form='formatted')

         do jj=-ss%n_grid_fit,ss%n_grid_fit
!allocate and set data arrays
            if(tf%l_fft_timefreq) then
               freq=df*real(jj)
            else
               freq=tf%freqs_fit(jj)
            endif
            zz=cmplx(0.d0,freq)
            call value_on_frequency(se,ii,freq,gz)
            !call func_fit(se,zz,ii,fz)
            call value_on_frequency_complex(se,ii,zz,fz)
            if(tf%grid_fit==0) then
               write(iun,'(4f12.6)') freq, aimag(fz),aimag(ss%diag(ii,jj+ss%n+1)), aimag(gz)
            else
               write(iun,'(4f12.6)') freq, aimag(fz),aimag(ss%diag_freq_fit(ii,jj+ss%n_grid_fit+1)), aimag(gz)
            endif
         enddo
         close(iun)
      endif
   enddo
   return
  END SUBROUTINE


  SUBROUTINE value_on_frequency(se,is,omega,sigma)
!this subroutine calculates the value of the correlation
!part of the self-energy  on real frequency

   USE io_global, ONLY : stdout

   implicit none

   TYPE(self_expansion),INTENT(in)  :: se!self expansion data
   INTEGER,INTENT(in) :: is!state considered
   REAL(kind=DP), INTENT(in) :: omega!real frequency considered
   COMPLEX(kind=DP), INTENT(out) :: sigma! <\Psi_i|\Sigma_c(w)|\Psi_i>

   INTEGER :: ii
!control is
   if(is>se%max_i) then
     write(stdout,*) 'Routine value_on_frequency is too large'
     stop
   endif

   if(omega >= 0 ) then
     sigma=se%a_0(is)
     do ii=1,se%n_multipoles
       sigma=sigma+se%a(ii,is)/(cmplx(omega,0.d0)-se%b(ii,is))
     enddo
   else
     sigma=conjg(se%a_0(is))
     do ii=1,se%n_multipoles
       sigma=sigma+conjg(se%a(ii,is))/(cmplx(omega,0.d0)-conjg(se%b(ii,is)))
     enddo
   endif

   return

  END SUBROUTINE



  SUBROUTINE derivative_on_frequency(se,is,omega,dsigma)
!this subroutine calculates the value of the correlation
!part of the self-energy  on real frequency

   USE io_global, ONLY : stdout

   implicit none

   TYPE(self_expansion),INTENT(in)  :: se!self expansion data
   INTEGER,INTENT(in) :: is!state considered
   REAL(kind=DP), INTENT(in) :: omega!real frequency considered
   COMPLEX(kind=DP), INTENT(out) :: dsigma! (d<\Psi_i|\Sigma_c(w')|\Psi_i>/dw')_w

   INTEGER :: ii
!control is
   if(is>se%max_i) then
     write(stdout,*) 'Routine value_on_frequency is too large'
     stop
   endif

   if(omega >= 0 ) then
     dsigma=(0.d0,0.d0)
     do ii=1,se%n_multipoles
       dsigma=dsigma-se%a(ii,is)/((cmplx(omega,0.d0)-se%b(ii,is))**2.d0)
     enddo
   else
     dsigma=(0.d0,0.d0)
     do ii=1,se%n_multipoles
       dsigma=dsigma-conjg(se%a(ii,is))/((cmplx(omega,0.d0)-conjg(se%b(ii,is)))**2.d0)
     enddo
   endif

   return

  END SUBROUTINE

    SUBROUTINE value_on_frequency_complex(se,is,omega,sigma)
!this subroutine calculates the value of the correlation
!part of the self-energy  on complex frequency

   USE io_global, ONLY : stdout

   implicit none

   TYPE(self_expansion),INTENT(in)  :: se!self expansion data
   INTEGER,INTENT(in) :: is!state considered
   COMPLEX(kind=DP), INTENT(in) :: omega!real frequency considered
   COMPLEX(kind=DP), INTENT(out) :: sigma! <\Psi_i|\Sigma_c(w)|\Psi_i>

   INTEGER :: ii
!control is
   if(is>se%max_i) then
     write(stdout,*) 'Routine value_on_frequency is too large'
     stop
   endif

   if(real(omega) >= 0 ) then
     sigma=se%a_0(is)
     do ii=1,se%n_multipoles
       sigma=sigma+se%a(ii,is)/(omega-se%b(ii,is))
     enddo
   else
     sigma=conjg(se%a_0(is))
     do ii=1,se%n_multipoles
       sigma=sigma+conjg(se%a(ii,is))/(omega-conjg(se%b(ii,is)))!ATTENZIONE must be checked!!!
     enddo
   endif

   return

  END SUBROUTINE


  END MODULE
