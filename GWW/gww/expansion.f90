!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
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
    INTEGER :: nspin!spin multiplicity
    COMPLEX(kind=DP), DIMENSION(:,:), POINTER :: a_0!parameters a_0
    COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: a!parameters a (n_multipoles,max_i)
    COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: b!parameters b (n_multipoles,max_i)
    LOGICAL :: whole_s!if true consider also off diagonal elements
    INTEGER :: i_min_whole!range for off diagonal elements
    INTEGER :: i_max_whole
    COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: a_0_off!parameters a_0 for off diagonal                     
    COMPLEX(kind=DP), DIMENSION(:,:,:,:), POINTER :: a_off!parameters a (n_multipoles,range,max_i) for  off diagonal
    COMPLEX(kind=DP), DIMENSION(:,:,:,:), POINTER :: b_off!parameters b (n_multipoles,range,max_i)  for off diagonal
  END TYPE self_expansion



  CONTAINS

    SUBROUTINE initialize_self_expansion(se)
      implicit none

      TYPE(self_expansion) :: se
      nullify(se%a_0)
      nullify(se%a)
      nullify(se%b)
      nullify(se%a_0_off)
      nullify(se%a_off)
      nullify(se%b_off)

      return

    END SUBROUTINE initialize_self_expansion
 
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

   if(associated(se%a_0_off)) then
     deallocate(se%a_0_off)
     nullify(se%a_0_off)
   endif

   if(associated(se%a_off)) then
     deallocate(se%a_off)
     nullify(se%a_off)
   endif

   if(associated(se%b_off)) then
     deallocate(se%b_off)
     nullify(se%b_off)
   endif

   return

  END SUBROUTINE

  SUBROUTINE write_self_expansion(se)
!this subroutine writes the multipole expansion on disk
     USE io_global,            ONLY : stdout, ionode
     USE input_gw,             ONLY : input_options
     USE io_files,             ONLY : prefix,tmp_dir

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(self_expansion), INTENT(in) :: se!object to be written
   
    INTEGER :: iun
    if(ionode) then
         iun = find_free_unit()
         open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'self_expansion', status='unknown',form='unformatted')
         write(iun) se%max_i
         write(iun) se%i_min
         write(iun) se%i_max
         write(iun) se%n_multipoles
         write(iun) se%nspin
         write(iun) se%whole_s
         write(iun) se%i_min_whole
         write(iun) se%i_max_whole
         write(iun) se%a_0(1:se%max_i,1:se%nspin)
         write(iun) se%a(1:se%n_multipoles,1:se%max_i,1:se%nspin) 
         write(iun) se%b(1:se%n_multipoles,1:se%max_i,1:se%nspin)  
         if(se%whole_s) then
            write(iun) se%a_0_off(se%i_min_whole:se%i_max_whole,1:se%max_i,1:se%nspin)
            write(iun) se%a_off(1:se%n_multipoles,se%i_min_whole:se%i_max_whole,1:se%max_i,1:se%nspin)
            write(iun) se%b_off(1:se%n_multipoles,se%i_min_whole:se%i_max_whole,1:se%max_i,1:se%nspin)
         endif
         close(iun)
      endif


    return
  END SUBROUTINE write_self_expansion

  SUBROUTINE read_self_expansion(se)
!this subroutine reads the multipole expansion from disk

     USE io_global,            ONLY : stdout, ionode, ionode_id
     USE input_gw,             ONLY : input_options
     USE io_files,             ONLY : prefix,tmp_dir
     USE mp,                   ONLY : mp_bcast
     USE mp_world,             ONLY : world_comm

    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(self_expansion), INTENT(out) :: se!object to be written

    INTEGER :: iun
    if(ionode) then
         iun = find_free_unit()
         open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'self_expansion', status='old',form='unformatted')
         read(iun) se%max_i
         read(iun) se%i_min
         read(iun) se%i_max
         read(iun) se%n_multipoles
         read(iun) se%nspin
         read(iun) se%whole_s
         read(iun) se%i_min_whole
         read(iun) se%i_max_whole
      endif

      call mp_bcast(se%max_i, ionode_id,world_comm)
      call mp_bcast(se%i_min, ionode_id,world_comm)
      call mp_bcast(se%i_max, ionode_id,world_comm)
      call mp_bcast(se%n_multipoles, ionode_id,world_comm)
      call mp_bcast(se%nspin, ionode_id,world_comm)
      call mp_bcast(se%whole_s, ionode_id,world_comm)
      call mp_bcast(se%i_min_whole, ionode_id,world_comm)
      call mp_bcast(se%i_max_whole, ionode_id,world_comm)

      allocate(se%a_0(se%max_i,se%nspin),se%a(se%n_multipoles,se%max_i,se%nspin))
      allocate(se%b(se%n_multipoles,se%max_i,se%nspin))

      if(ionode) then
         read(iun) se%a_0(1:se%max_i,1:se%nspin)
         read(iun) se%a(1:se%n_multipoles,1:se%max_i,1:se%nspin)
         read(iun) se%b(1:se%n_multipoles,1:se%max_i,1:se%nspin)
      endif
      
      call mp_bcast(se%a_0,ionode_id,world_comm)
      call mp_bcast(se%a, ionode_id,world_comm)
      call mp_bcast(se%b, ionode_id,world_comm)
      

      if(se%whole_s) then
         allocate(se%a_0_off(se%i_min_whole:se%i_max_whole,se%max_i,se%nspin))
         allocate( se%a_off(se%n_multipoles,se%i_min_whole:se%i_max_whole,se%max_i,se%nspin))
         allocate(se%b_off(se%n_multipoles,se%i_min_whole:se%i_max_whole,se%max_i,se%nspin))
         if(ionode) then
            read(iun) se%a_0_off(se%i_min_whole:se%i_max_whole,1:se%max_i,1:se%nspin)
            read(iun) se%a_off(1:se%n_multipoles,se%i_min_whole:se%i_max_whole,1:se%max_i,1:se%nspin)
            read(iun) se%b_off(1:se%n_multipoles,se%i_min_whole:se%i_max_whole,1:se%max_i,1:se%nspin)
         endif
         
         call mp_bcast(se%a_0_off,ionode_id,world_comm)
         call mp_bcast(se%a_off, ionode_id,world_comm)
         call mp_bcast(se%b_off, ionode_id,world_comm)

      else
         nullify(se%a_0_off)
         nullify(se%a_off)
         nullify(se%b_off)

      endif
      if(ionode) close(iun)
      return
    END SUBROUTINE read_self_expansion


  SUBROUTINE create_self_energy_fit( tf, se,ss, options,sr,l_real_axis)
!this subroutine fit the self energy in the imaginary frequency
!with a multipole complex function
!parallel on states

    USE io_global,  ONLY : stdout
    USE input_gw,   ONLY : input_options
    USE constants,  ONLY : pi
    USE self_energy_storage, ONLY : self_storage, self_on_real
    USE para_gww,   ONLY : is_my_state_range
    USE mp,         ONLY : mp_sum,mp_barrier
    USE mp_world,   ONLY : world_comm
    USE times_gw,   ONLY : times_freqs

    implicit none

    TYPE(times_freqs), INTENT(in) :: tf!frequency grid
    TYPE(self_expansion) :: se!fit to be created
    TYPE(self_storage) :: ss!data on frequency
    TYPE(input_options) :: options! for number of multipoles
    TYPE(self_on_real), INTENT(in) :: sr!for self energy on real axis
    LOGICAL, INTENT(in) :: l_real_axis

    INTEGER :: ii,jj, kk,is
    COMPLEX(kind=DP), ALLOCATABLE :: z(:),s(:)
    REAL(kind=DP) :: df,freq, totalperiod, chi, chi0
    INTEGER :: icyc
    COMPLEX(kind=DP) :: a_0_old, a_0_good
    COMPLEX(kind=DP), ALLOCATABLE :: a_old(:), b_old(:), a_good(:), b_good(:)
    INTEGER :: n_sample

!sets:
    se%max_i=options%max_i
    se%i_min=options%i_min
    se%i_max=options%i_max
    se%n_multipoles=options%n_multipoles
    se%whole_s=options%whole_s
    se%i_min_whole=options%i_min_whole
    se%i_max_whole=options%i_max_whole
    se%nspin=ss%nspin

!allocates:
!    call free_memory_self_expansion(se)
    allocate(se%a_0(se%max_i,se%nspin))
    allocate(se%a(se%n_multipoles,se%max_i,se%nspin))
    allocate(se%b(se%n_multipoles,se%max_i,se%nspin))
    if(se%whole_s) then
       allocate(se%a_0_off(se%i_min_whole:se%i_max_whole,se%max_i,se%nspin))
       allocate(se%a_off(se%n_multipoles,se%i_min_whole:se%i_max_whole,se%max_i,se%nspin))
       allocate(se%b_off(se%n_multipoles,se%i_min_whole:se%i_max_whole,se%max_i,se%nspin))
    else
       nullify(se%a_0_off)
       nullify(se%a_off)
       nullify(se%b_off)
    endif
    allocate(a_old(se%n_multipoles))
    allocate(b_old(se%n_multipoles))
    allocate(a_good(se%n_multipoles))
    allocate(b_good(se%n_multipoles))

  
  !  allocate(z(ss%n_grid_fit),s(ss%n_grid_fit))
    if(.not.l_real_axis) then
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
    else
       allocate(z(options%n_real_axis),s(options%n_real_axis))
       z(1:options%n_real_axis)=sr%grid(1:options%n_real_axis)
    endif

!some checks
    if(ss%ontime .and. tf%grid_fit==0) then
      write(stdout,*) 'Subroutine self_energy_fit: imaginary frequency required'
      stop
    endif

    if(se%whole_s) then
       se%a_0_off=(0.d0,0.d0)
       se%a_off=(0.d0,0.d0)
       se%b_off=(0.d0,0.d0)
    endif
    se%a_0=(0.d0,0.d0)
    se%a(:,:,:)=(0.d0,0.d0)
    se%b(:,:,:)=(0.d0,0.d0)


 !loop on spin

    do is=1,se%nspin
       do ii=se%i_min,se%i_max!loop on states
          chi0=1.d10
!set initial values
          if(is_my_state_range(ii)) then
             if(.not.l_real_axis) then
                if(tf%grid_fit==0) then
                   do jj=1,options%n_fit!ss%n
                      s(jj)=ss%diag(ii,jj+ss%n+1,is)
                   enddo
                else
             !do jj=1,ss%n_grid_fit
                   do jj=0+options%offset_fit,options%n_fit+options%offset_fit-1!ATTENZIONE
                      s(jj-options%offset_fit+1)=ss%diag_freq_fit(ii,jj+ss%n_grid_fit+1,is)!ATTENZIONE
                   enddo
                endif
             else
                s(1:options%n_real_axis)=sr%diag(1:options%n_real_axis,ii,1)
             endif
             se%a_0(ii,is)=(0.0,0.0d0)
             do jj=1,options%n_multipoles
                se%a(jj,ii,is)=cmplx(real(jj)*(0.01d0),0.d0)
                se%b(jj,ii,is)=cmplx((0.5d0)*real(jj)*(-1.d0)**real(jj),-0.01d0)
             enddo
          
             do icyc=1,options%cyc_minpack
                write(stdout,*) 'Call fit_multipole'                                                               
                FLUSH(stdout)
                if(.not.l_real_axis) then
                   n_sample=options%n_fit
                else
                   n_sample=options%n_real_axis
                endif

                call fit_multipole(n_sample,options%n_multipoles,z,s,se%a_0(ii,is),&
                     &se%a(:,ii,is),se%b(:,ii,is),1.d0,options%fit_thres,options%fit_maxiter)
                write(stdout,*) 'Done'
                FLUSH(stdout)

                a_0_old=se%a_0(ii,is)
                do jj=1,options%n_multipoles
                   a_old(jj)=se%a(jj,ii,is)
                   b_old(jj)=se%b(jj,ii,is)
                enddo
             
                
                if(options%n_max_minpack /= 0) then
                   write(stdout,*) 'Calling minpack'!ATTENZIONE
                   FLUSH(stdout)
                   call fit_multipole_minpack(n_sample,options%n_multipoles,z,s,se%a_0(ii,is),&
                        &se%a(:,ii,is),se%b(:,ii,is),options%fit_thres, options%n_max_minpack, chi)
                               
                   write(stdout,*) 'Done'!ATTENZIONE                                                   
                   FLUSH(stdout)

                endif

                if(chi <= chi0) then
                   a_0_good=se%a_0(ii,is)
                   do jj=1,options%n_multipoles
                      a_good(jj)=se%a(jj,ii,is)
                      b_good(jj)=se%b(jj,ii,is)
                   enddo
                   chi0=chi
                endif
             
                se%a_0(ii,is)=a_0_old
                do jj=1,options%n_multipoles
                   se%a(jj,ii,is)=a_old(jj)
                   se%b(jj,ii,is)=b_old(jj)
                enddo
             
             enddo
             se%a_0(ii,is)=a_0_good
             do jj=1,options%n_multipoles
                se%a(jj,ii,is)=a_good(jj)
                se%b(jj,ii,is)=b_good(jj)
             enddo
          
             write(stdout,*) 'FIT state :', ii,is
             write(stdout,*) 'FIT    a_0:', se%a_0(ii,is)
             do jj=1,options%n_multipoles
                write(stdout,*) 'FIT    a:',jj,se%a(jj,ii,is)
                write(stdout,*) 'FIT    b:',jj,se%b(jj,ii,is)
             enddo
             FLUSH(stdout)
          endif
       enddo
       call mp_sum(se%a_0(:,is),world_comm)
       call mp_sum(se%a(:,:,is),world_comm)
       call mp_sum(se%b(:,:,is),world_comm)
   

!!!!!!!!!!!now  off diagonal part
       if(se%whole_s) then

          do kk=se%i_min_whole,se%i_max_whole
             do ii=se%i_min,se%i_max!lo                                           
                chi0=1.d10
             !set initial values                                    
                if(is_my_state_range(ii)) then
                   if(tf%grid_fit==0) then
                      do jj=1,options%n_fit!ss%n
                         s(jj)=ss%whole(kk,ii,jj+ss%n+1,is)
                      enddo
                   else
                      do jj=0+options%offset_fit,options%n_fit+options%offset_fit-1!ATTENZIONE 
                         s(jj-options%offset_fit+1)=ss%whole_freq_fit(kk,ii,jj+ss%n_grid_fit+1,is)!ATTENZIONE 
                      enddo
                   endif
                   
                   se%a_0_off(kk,ii,is)=(0.0,0.0d0)
                   do jj=1,options%n_multipoles
                      se%a_off(jj,kk,ii,is)=cmplx(real(jj)*(0.01d0),0.d0)
                      se%b_off(jj,kk,ii,is)=cmplx((0.5d0)*real(jj)*(-1.d0)**real(jj),-0.01d0)
                   enddo

                   do icyc=1,options%cyc_minpack
                      call fit_multipole(options%n_fit,options%n_multipoles,z,s,&
                           &se%a_0_off(kk,ii,is),se%a_off(:,kk,ii,is),se%b_off(:,kk,ii,is),&
                           &1.d0,options%fit_thres,options%fit_maxiter)
                      a_0_old=se%a_0_off(kk,ii,is)
                      do jj=1,options%n_multipoles
                         a_old(jj)=se%a_off(jj,kk,ii,is)
                         b_old(jj)=se%b_off(jj,kk,ii,is)
                      enddo


                      if(options%n_max_minpack /= 0) then
                         write(stdout,*) 'Calling minpack'!ATTENZIONE
                         FLUSH(stdout)
                         call fit_multipole_minpack(options%n_fit,options%n_multipoles,z,s,&
    &se%a_0_off(kk,ii,is),se%a_off(:,kk,ii,is),se%b_off(:,kk,ii,is),options%fit_thres, options%n_max_minpack, chi)

                      endif

                      if(chi <= chi0) then
                         a_0_good=se%a_0_off(kk,ii,is)
                         do jj=1,options%n_multipoles
                            a_good(jj)=se%a_off(jj,kk,ii,is)
                            b_good(jj)=se%b_off(jj,kk,ii,is)
                         enddo
                         chi0=chi
                      endif
                      
                      se%a_0_off(kk,ii,is)=a_0_old
                      do jj=1,options%n_multipoles
                         se%a_off(jj,kk,ii,is)=a_old(jj)
                         se%b_off(jj,kk,ii,is)=b_old(jj)
                      enddo

                   enddo
                   se%a_0_off(kk,ii,is)=a_0_good
                   do jj=1,options%n_multipoles
                      se%a_off(jj,kk,ii,is)=a_good(jj)
                      se%b_off(jj,kk,ii,is)=b_good(jj)
                   enddo
                   
                   write(stdout,*) 'FIT off diagonal :', ii, kk,is
                   write(stdout,*) 'FIT    a_0:', se%a_0_off(kk,ii,is)
                   do jj=1,options%n_multipoles
                      write(stdout,*) 'FIT    a:',jj,se%a_off(jj,kk,ii,is)
                      write(stdout,*) 'FIT    b:',jj,se%b_off(jj,kk,ii,is)
                   enddo
                   FLUSH(stdout)
                endif
             enddo


          enddo
          call mp_sum(se%a_0_off(:,:,is),world_comm)
          call mp_sum(se%a_off(:,:,:,is),world_comm)
          call mp_sum(se%b_off(:,:,:,is),world_comm)
          


       endif

    enddo!nspin
    deallocate(z,s)
    deallocate(a_old,b_old)
    deallocate(a_good,b_good)

    call mp_barrier( world_comm )
    write(stdout,*) 'Out of create_self_energy_fit'
    FLUSH(stdout)

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

   fz=se%a_0(i,1)
   do jj=1,se%n_multipoles
     fz=fz+se%a(jj,i,1)/(z-se%b(jj,i,1))
   enddo

   return


  END SUBROUTINE 

  SUBROUTINE print_fit_onfile(tf, se,ss)
!this subroutines prints the resulta of the fit on file:
!real and imaginary part on imaginary frequency, with results self_energy
!and real and imaginary part on real frequency
!parallel on states
    USE self_energy_storage,  ONLY : self_storage
   USE constants,            ONLY : pi
   USE io_global,            ONLY : ionode
   USE para_gww,             ONLY : is_my_state_range
   USE times_gw,             ONLY : times_freqs
      USE io_files,  ONLY : prefix, tmp_dir

   implicit none
   INTEGER, EXTERNAL :: find_free_unit
   TYPE(times_freqs), INTENT(in) :: tf!for frequency grid
   TYPE(self_expansion) :: se!parameters of fit
   TYPE(self_storage)   :: ss!self energy data
   

   INTEGER :: ii,jj,kk,is
   INTEGER :: iun
   CHARACTER(5) :: nfile,mfile
   REAL(kind=DP) :: totalperiod,df,freq
   COMPLEX(kind=DP) :: zz
   COMPLEX(kind=DP) :: fz, gz


   do is=1,se%nspin
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
            if(is==1) then
               open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'re_on_im'// nfile, status='unknown',form='formatted')
            else
               open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'re_on_im2'// nfile, status='unknown',form='formatted')
            endif

            do jj=-ss%n_grid_fit,ss%n_grid_fit
!allocate and set data arrays
               if(tf%l_fft_timefreq) then
                  freq=df*real(jj)
               else
                  freq=tf%freqs_fit(jj)
               endif
               zz=cmplx(0.d0,freq)
               call value_on_frequency(se,ii,freq,gz,is)
            !call func_fit(se,zz,ii,fz)
               call value_on_frequency_complex(se,ii,zz,fz,is)            
               if(tf%grid_fit==0) then
                  write(iun,'(4f14.8)') freq, real(fz),real(ss%diag(ii,jj+ss%n,is)),real( gz)     
               else
                  write(iun,'(4f14.8)') freq, real(fz),real(ss%diag_freq_fit(ii,jj+ss%n_grid_fit+1,is)),real( gz)
               endif
            enddo
         
            close(iun)


!now imaginary  part on imaginary frequency

!openfile
            iun = find_free_unit()
            if(is==1) then
               open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'im_on_im'// nfile, status='unknown',form='formatted')
            else
               open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'im_on_im2'// nfile, status='unknown',form='formatted')
            endif

            do jj=-ss%n_grid_fit,ss%n_grid_fit
!allocate and set data arrays
               if(tf%l_fft_timefreq) then
                  freq=df*real(jj)
               else
                  freq=tf%freqs_fit(jj)
               endif
               zz=cmplx(0.d0,freq)
               call value_on_frequency(se,ii,freq,gz,is)
               !call func_fit(se,zz,ii,fz)
               call value_on_frequency_complex(se,ii,zz,fz,is)
               if(tf%grid_fit==0) then
                  write(iun,'(4f14.8)') freq, aimag(fz),aimag(ss%diag(ii,jj+ss%n+1,is)), aimag(gz)
               else
                  write(iun,'(4f14.8)') freq, aimag(fz),aimag(ss%diag_freq_fit(ii,jj+ss%n_grid_fit+1,is)), aimag(gz)
               endif
            enddo
            close(iun)
!now off diagonal terms
            if(se%whole_s) then

!set file name
            
               do kk=se%i_min_whole,se%i_max_whole

                  write(nfile,'(5i1)') &
                       & kk/10000,mod(kk,10000)/1000,mod(kk,1000)/100,mod(kk,100)/10,mod(kk,10)
                  write(mfile,'(5i1)') &
                       & ii/10000,mod(ii,10000)/1000,mod(ii,1000)/100,mod(ii,100)/10,mod(ii,10)

                  totalperiod=2.d0*ss%tau+2.d0*ss%tau/real(ss%n)
                  df=2.d0*pi/totalperiod

!now real part on imaginary frequency
!openfile                            
                  iun = find_free_unit()
                  if(is==1) then
                     open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'off_re_on_im'// nfile // '_' // mfile, &
                          &status='unknown',form='formatted')
                  else
                     open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'off_re_on_im2'// nfile // '_' // mfile, &
                          &status='unknown',form='formatted')

                  endif
                  do jj=-ss%n_grid_fit,ss%n_grid_fit
!allocate and set data arrays
                     if(tf%l_fft_timefreq) then
                        freq=df*real(jj)
                     else
                        freq=tf%freqs_fit(jj)
                     endif
                     zz=cmplx(0.d0,freq)
                     call value_on_frequency_off(se,kk,ii,freq,gz,is)
            !call func_fit(se,zz,ii,fz)
                     call value_on_frequency_complex_off(se,kk,ii,zz,fz,is)
                     if(tf%grid_fit==0) then
                        write(iun,'(4f12.6)') freq, real(fz),real(ss%whole(kk,ii,jj+ss%n,is)),real( gz)
                     else
                        write(iun,'(4f12.6)') freq, real(fz),real(ss%whole_freq_fit(kk,ii,jj+ss%n_grid_fit+1,is)),real( gz)
                     endif
                  enddo

                  close(iun)


!now imaginary  part on imaginary frequency
!openfile                                  
                  iun = find_free_unit()
                  if(is==1) then
                     open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'off_im_on_im'// nfile // '_' // mfile, &
                          &status='unknown',form='formatted')
                  else
                     open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'off_im_on_im2'// nfile // '_' // mfile, &
                          &status='unknown',form='formatted')
                  endif
                  do jj=-ss%n_grid_fit,ss%n_grid_fit
!allocate and set data arrays
                     if(tf%l_fft_timefreq) then
                        freq=df*real(jj)
                     else
                        freq=tf%freqs_fit(jj)
                     endif
                     zz=cmplx(0.d0,freq)
                     call value_on_frequency_off(se,kk,ii,freq,gz,is)
!call func_fit(se,zz,ii,fz)
                     call value_on_frequency_complex_off(se,kk,ii,zz,fz,is)
                     if(tf%grid_fit==0) then
                        write(iun,'(4f12.6)') freq, aimag(fz),aimag(ss%whole(kk,ii,jj+ss%n+1,is)), aimag(gz)
                     else
                        write(iun,'(4f12.6)') freq, aimag(fz),aimag(ss%whole_freq_fit(kk,ii,jj+ss%n_grid_fit+1,is)), aimag(gz)
                     endif
                  enddo
                  close(iun)

               enddo
            endif
            
         endif
      enddo
   enddo!loop on spin
   return
 END SUBROUTINE print_fit_onfile


  SUBROUTINE value_on_frequency(se,is,omega,sigma,ispin)
!this subroutine calculates the value of the correlation
!part of the self-energy  on real frequency

   USE io_global, ONLY : stdout

   implicit none

   TYPE(self_expansion),INTENT(in)  :: se!self expansion data
   INTEGER,INTENT(in) :: is!state considered
   REAL(kind=DP), INTENT(in) :: omega!real frequency considered
   COMPLEX(kind=DP), INTENT(out) :: sigma! <\Psi_i|\Sigma_c(w)|\Psi_i>
   INTEGER, INTENT(in) :: ispin!spin channel

   INTEGER :: ii 
!control is
   if(is>se%max_i) then
     write(stdout,*) 'Routine value_on_frequency is too large'
     stop
   endif

   if(omega >= 0 ) then
     sigma=se%a_0(is,ispin)
     do ii=1,se%n_multipoles
       sigma=sigma+se%a(ii,is,ispin)/(cmplx(omega,0.d0)-se%b(ii,is,ispin))
     enddo
   else
     sigma=conjg(se%a_0(is,ispin))
     do ii=1,se%n_multipoles
       sigma=sigma+conjg(se%a(ii,is,ispin))/(cmplx(omega,0.d0)-conjg(se%b(ii,is,ispin)))
     enddo
   endif  

   return

  END SUBROUTINE
 

  SUBROUTINE value_on_frequency_off(se,is,js,omega,sigma,ispin)
!this subroutine calculates the value of the correlation                                                                                                   
!part of the self-energy  on real frequency                                                                                                                

   USE io_global, ONLY : stdout

   implicit none

   TYPE(self_expansion),INTENT(in)  :: se!self expansion data                                                                                              
   INTEGER,INTENT(in) :: is,js!state considered                                                                        
   REAL(kind=DP), INTENT(in) :: omega!real frequency considered
   COMPLEX(kind=DP), INTENT(out) :: sigma! <\Psi_i|\Sigma_c(w)|\Psi_i>                                                                                     
   INTEGER, INTENT(in) :: ispin!spin channel

   INTEGER :: ii 
  
   if(omega >= 0 ) then
     sigma=se%a_0_off(is,js,ispin)
     do ii=1,se%n_multipoles
       sigma=sigma+se%a_off(ii,is,js,ispin)/(cmplx(omega,0.d0)-se%b_off(ii,is,js,ispin))
     enddo
   else
     sigma=conjg(se%a_0_off(is,js,ispin))
     do ii=1,se%n_multipoles
       sigma=sigma+conjg(se%a_off(ii,is,js,ispin))/(cmplx(omega,0.d0)-conjg(se%b_off(ii,is,js,ispin)))
     enddo
   endif

   return

 END SUBROUTINE value_on_frequency_off

  SUBROUTINE derivative_on_frequency(se,is,omega,dsigma,ispin)
!this subroutine calculates the value of the correlation
!part of the self-energy  on real frequency

   USE io_global, ONLY : stdout

   implicit none

   TYPE(self_expansion),INTENT(in)  :: se!self expansion data
   INTEGER,INTENT(in) :: is!state considered
   REAL(kind=DP), INTENT(in) :: omega!real frequency considered
   COMPLEX(kind=DP), INTENT(out) :: dsigma! (d<\Psi_i|\Sigma_c(w')|\Psi_i>/dw')_w
   INTEGER, INTENT(in) :: ispin!spin channel

   INTEGER :: ii
!control is
   if(is>se%max_i) then
     write(stdout,*) 'Routine value_on_frequency is too large'
     stop
   endif

   if(omega >= 0 ) then
     dsigma=(0.d0,0.d0)
     do ii=1,se%n_multipoles
       dsigma=dsigma-se%a(ii,is,ispin)/((cmplx(omega,0.d0)-se%b(ii,is,ispin))**2.d0)
     enddo
   else
     dsigma=(0.d0,0.d0)
     do ii=1,se%n_multipoles
       dsigma=dsigma-conjg(se%a(ii,is,ispin))/((cmplx(omega,0.d0)-conjg(se%b(ii,is,ispin)))**2.d0)
     enddo
   endif

   return

  END SUBROUTINE

    SUBROUTINE value_on_frequency_complex(se,is,omega,sigma,ispin)
!this subroutine calculates the value of the correlation
!part of the self-energy  on complex frequency

   USE io_global, ONLY : stdout

   implicit none

   TYPE(self_expansion),INTENT(in)  :: se!self expansion data
   INTEGER,INTENT(in) :: is!state considered
   COMPLEX(kind=DP), INTENT(in) :: omega!real frequency considered
   COMPLEX(kind=DP), INTENT(out) :: sigma! <\Psi_i|\Sigma_c(w)|\Psi_i>
   INTEGER, INTENT(in) :: ispin!spin channel

   INTEGER :: ii
!control is
   if(is>se%max_i) then
     write(stdout,*) 'Routine value_on_frequency is too large'
     stop
   endif

   if(real(omega) >= 0 ) then
     sigma=se%a_0(is,ispin)
     do ii=1,se%n_multipoles
       sigma=sigma+se%a(ii,is,ispin)/(omega-se%b(ii,is,ispin))
     enddo
   else
     sigma=conjg(se%a_0(is,ispin))
     do ii=1,se%n_multipoles
       sigma=sigma+conjg(se%a(ii,is,ispin))/(omega-conjg(se%b(ii,is,ispin)))!ATTENZIONE must be checked!!!
     enddo
   endif

   return

  END SUBROUTINE

 SUBROUTINE value_on_frequency_complex_off(se,is,js,omega,sigma,ispin)
!this subroutine calculates the value of the correlation
!part of the self-energy  on complex frequency

   USE io_global, ONLY : stdout

   implicit none

   TYPE(self_expansion),INTENT(in)  :: se!self expansion data
                                                                                                                        
   INTEGER,INTENT(in) :: is,js!state considered 
                                                                                                                                        
   COMPLEX(kind=DP), INTENT(in) :: omega!real frequency considered 
   COMPLEX(kind=DP), INTENT(out) :: sigma! <\Psi_i|\Sigma_c(w)|\Psi_i>
   INTEGER, INTENT(in) :: ispin!spin channel

   INTEGER :: ii

   if(real(omega) >= 0 ) then
     sigma=se%a_0_off(is,js,ispin)
     do ii=1,se%n_multipoles
       sigma=sigma+se%a_off(ii,is,js,ispin)/(omega-se%b_off(ii,is,js,ispin))
     enddo
   else
     sigma=conjg(se%a_0_off(is,js,ispin))
     do ii=1,se%n_multipoles
       sigma=sigma+conjg(se%a_off(ii,is,js,ispin))/(omega-conjg(se%b_off(ii,is,js,ispin)))
     enddo
   endif

   return

 END SUBROUTINE value_on_frequency_complex_off


  END MODULE
