!
! P.Umari Program GWW
! Modified by G. Stenuit
!
MODULE self_energy_storage
!this modules contains the structure and subroutines
!to store the expectation values of the self-energy
!and to perform ffts and fits
!in parallel version the calculations on times are parallelized

  USE kinds, ONLY  : DP

  TYPE self_storage
!descriptor of <Psi_i|\Sigma|\Psi_j>
     LOGICAL :: ontime!if .true. data is on imaginary time , otherwise imaginary frequency
     LOGICAL :: whole_s!if .true. also the off-diagonal elements are considered
     INTEGER :: n!number of sample on positive and on negative times (total of 2*n+1 samples)
     INTEGER :: n_grid_fit!number of sample on positive and on negative frequencies for fit (total of 2*n+1 samples)
     INTEGER :: max_i!number of states considered
     INTEGER :: i_min!minimum state to be calculated
     INTEGER :: i_max!maximum state to be calculated
     REAL(kind=DP) :: tau!max time (on imaginary axes)
     COMPLEX(kind=DP), DIMENSION(:,:), POINTER :: diag !values <Psi_i|\Sigma|\Psi_i>,time_j
     COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: whole !values  <Psi_i|\Sigma|\Psi_j>,time_k
     COMPLEX(kind=DP), DIMENSION(:,:), POINTER :: diag_freq_fit !values <Psi_i|\Sigma|\Psi_i>,on frequency for fit
     COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: whole_freq_fit !values  <Psi_i|\Sigma|\Psi_j>, on frequency for fit
     REAL(kind=DP), POINTER, DIMENSION(:) :: ene_remainder!for storing remainders
  END  TYPE self_storage

CONTAINS

  SUBROUTINE write_self_storage_ondisk(ss, options)
!this subroutine writes the green function on disk
!the file name is taken from the label

    USE io_files,           ONLY : find_free_unit
    USE io_global,          ONLY : stdout, ionode
    USE input_gw,           ONLY : input_options
    USE mp,                 ONLY : mp_barrier
    implicit none

    TYPE(self_storage) :: ss!the self_energy descriptor to be written on file
    TYPE(input_options) :: options!for debug flag

    INTEGER :: iw, jw, kw, iun

    if(ionode) then
       iun = find_free_unit()
       if(.not. options%debug) then
          open( unit=iun, file='storage', status='unknown',form='unformatted')
       else
          open( unit=iun, file='storage', status='unknown',form='formatted')
       endif

       if(.not.options%debug) then
          write(iun) ss%ontime
          write(iun) ss%whole_s
          write(iun) ss%n
          write(iun) ss%max_i
          write(iun) ss%i_min
          write(iun) ss%i_max
          write(iun) ss%tau
          write(iun) ss%n_grid_fit
          if(.not.ss%whole_s) then
             do iw=1,2*ss%n+1
                write(iun) ss%diag(1:ss%max_i,iw)
             end do
          else
             do iw=1,2*ss%n+1
                write(iun) ss%whole(1:ss%max_i,1:ss%max_i,iw)
             end do
          endif
          write(iun) ss%ene_remainder(1:ss%max_i)
          if(.not.ss%whole_s) then
             do iw=1,2*ss%n_grid_fit+1
                write(iun) ss%diag_freq_fit(1:ss%max_i,iw)
             end do
          else
             do iw=1,2*ss%n_grid_fit+1
                write(iun) ss%whole_freq_fit(1:ss%max_i,1:ss%max_i,iw)
             end do
          endif

       else
          write(iun,*) ss%ontime
          write(iun,*) ss%whole_s
          write(iun,*) ss%n
          write(iun,*) ss%max_i
          write(iun,*) ss%i_min
          write(iun,*) ss%i_max
          write(iun,*) ss%tau
          write(iun,*) ss%n_grid_fit
          if(.not.ss%whole_s) then
             do iw=1,2*ss%n+1
                do jw=1,ss%max_i
                   write(iun,*) ss%diag(jw,iw)
                enddo
             end do
          else
             do iw=1,2*ss%n+1
                do jw=1,ss%max_i
                   do kw=1,ss%max_i
                      write(iun,*) ss%whole(jw,kw,iw)
                   enddo
                enddo
             end do
          endif
          do iw=1,ss%max_i
             write(iun,*) ss%ene_remainder(iw)
          enddo
          if(.not.ss%whole_s) then
             do iw=1,2*ss%n_grid_fit+1
                do jw=1,ss%max_i
                   write(iun,*) ss%diag_freq_fit(jw,iw)
                enddo
             end do
          else
             do iw=1,2*ss%n_grid_fit+1
                do jw=1,ss%max_i
                   do kw=1,ss%max_i
                      write(iun,*) ss%whole_freq_fit(jw,kw,iw)
                   enddo
                enddo
             end do
          endif
       endif
       close(iun)
    endif
    call mp_barrier
  END SUBROUTINE write_self_storage_ondisk

  SUBROUTINE read_self_storage_ondisk(ss, options)
!this subroutine writes the green function on disk
!the file name is taken from the label

    USE io_files,           ONLY : find_free_unit
    USE io_global,          ONLY : stdout, ionode, ionode_id
    USE input_gw,           ONLY : input_options
    USE mp,                 ONLY : mp_barrier, mp_bcast
    implicit none

    TYPE(self_storage) :: ss!the self_energy descriptor to be read from file
    TYPE(input_options) :: options!for debug flag

    INTEGER :: iw, jw, kw, iun

    if(ionode) then
       iun = find_free_unit()
       if(.not. options%debug) then
          open( unit=iun, file='storage', status='old',form='unformatted')
       else
          open( unit=iun, file='storage', status='old',form='formatted')
       endif
    endif

  !  call free_memory_self_storage(ss)

    if(ionode) then
       if(.not.options%debug) then
          read(iun) ss%ontime
          read(iun) ss%whole_s
          read(iun) ss%n
          read(iun) ss%max_i
          read(iun) ss%i_min
          read(iun) ss%i_max
          read(iun) ss%tau
          read(iun) ss%n_grid_fit
       else
          read(iun,*) ss%ontime
          read(iun,*) ss%whole_s
          read(iun,*) ss%n
          read(iun,*) ss%max_i
          read(iun,*) ss%i_min
          read(iun,*) ss%i_max
          read(iun,*) ss%tau
          read(iun,*) ss%n_grid_fit
       endif
    endif

    call mp_bcast(ss%ontime, ionode_id)
    call mp_bcast(ss%whole_s, ionode_id)
    call mp_bcast(ss%n, ionode_id)
    call mp_bcast(ss%max_i, ionode_id)
    call mp_bcast(ss%i_min, ionode_id)
    call mp_bcast(ss%i_max, ionode_id)
    call mp_bcast(ss%tau, ionode_id)
    call mp_bcast(ss%n_grid_fit, ionode_id)


!check for consistency
    if(ss%max_i/=options%max_i) then
      write(stdout,*) 'Routine  read_self_storage_ondisk max_i wrong'
      stop
    endif

!allocates
    if(ss%whole_s) then
      allocate(ss%whole(ss%max_i,ss%max_i,2*ss%n+1))
      nullify(ss%diag)
    else
      allocate(ss%diag(ss%max_i,2*ss%n+1))
      nullify(ss%whole)
    endif



    allocate(ss%ene_remainder(ss%max_i))

    if(ss%whole_s) then
      allocate(ss%whole_freq_fit(ss%max_i,ss%max_i,2*ss%n_grid_fit+1))
      nullify(ss%diag_freq_fit)
    else
      allocate(ss%diag_freq_fit(ss%max_i,2*ss%n_grid_fit+1))
      nullify(ss%whole_freq_fit)
    endif



    if(ionode) then
       if(.not.options%debug) then
          if(.not.ss%whole_s) then
             do iw=1,2*ss%n+1
                read(iun) ss%diag(1:ss%max_i,iw)
             end do
          else
             do iw=1,2*ss%n+1
                read(iun) ss%whole(1:ss%max_i,1:ss%max_i,iw)
             end do
          endif
          read(iun) ss%ene_remainder(1:ss%max_i)
           if(.not.ss%whole_s) then
             do iw=1,2*ss%n_grid_fit+1
                read(iun) ss%diag_freq_fit(1:ss%max_i,iw)
             end do
          else
             do iw=1,2*ss%n_grid_fit+1
                read(iun) ss%whole_freq_fit(1:ss%max_i,1:ss%max_i,iw)
             end do
          endif
       else
          if(.not.ss%whole_s) then
             do iw=1,2*ss%n+1
                do jw=1,ss%max_i
                   read(iun,*) ss%diag(jw,iw)
                enddo
             end do
          else
             do iw=1,2*ss%n+1
                do jw=1,ss%max_i
                   do kw=1,ss%max_i
                      read(iun) ss%whole(jw,kw,iw)
                   enddo
                enddo
             end do
          endif
          do iw=1,ss%max_i
             read(iun,*) ss%ene_remainder(iw)
          enddo
          if(.not.ss%whole_s) then
             do iw=1,2*ss%n_grid_fit+1
                do jw=1,ss%max_i
                   read(iun,*) ss%diag_freq_fit(jw,iw)
                enddo
             end do
          else
             do iw=1,2*ss%n_grid_fit+1
                do jw=1,ss%max_i
                   do kw=1,ss%max_i
                      read(iun) ss%whole_freq_fit(jw,kw,iw)
                   enddo
                enddo
             end do
          endif

       endif
       close(iun)
    endif
    if(.not.ss%whole_s) then
       call mp_bcast(ss%diag(:,:), ionode_id)
    else
       call mp_bcast(ss%whole(:,:,:), ionode_id)
    endif
    call mp_bcast(ss%ene_remainder(:), ionode_id)
    if(.not.ss%whole_s) then
       call mp_bcast(ss%diag_freq_fit(:,:), ionode_id)
    else
       call mp_bcast(ss%whole_freq_fit(:,:,:), ionode_id)
    endif

    return
   END SUBROUTINE read_self_storage_ondisk


  SUBROUTINE free_memory_self_storage(ss)
!deallocate if allocated
    implicit none
    TYPE(self_storage) :: ss
    if(associated(ss%diag)) deallocate(ss%diag)
    nullify(ss%diag)
    if(associated(ss%whole)) deallocate(ss%whole)
    nullify(ss%whole)
    if(associated(ss%ene_remainder)) deallocate(ss%ene_remainder)
    nullify(ss%ene_remainder)
    if(associated(ss%diag_freq_fit)) deallocate(ss%diag_freq_fit)
    nullify(ss%diag_freq_fit)
    if(associated(ss%whole_freq_fit)) deallocate(ss%whole_freq_fit)
    nullify(ss%whole_freq_fit)

  END SUBROUTINE

  SUBROUTINE set_remainder(ss, qp)
!this subroutine simply copy the self-energy remainders
!from ss to qp, in order to allow restarting
    USE energies_gww, ONLY : quasi_particles

    implicit none

    TYPE(self_storage) :: ss
    TYPE(quasi_particles) :: qp

    if(.not.associated(qp%ene_remainder))  allocate(qp%ene_remainder(ss%max_i))
    qp%ene_remainder(:)=ss%ene_remainder(:)

    return
  END SUBROUTINE set_remainder


  SUBROUTINE create_self_ontime(tf, ss,options,qp)
!this subroutine creates the structure self_storege
!on imaginary time
   USE  constants,         ONLY : eps8
   USE io_global,          ONLY : stdout, ionode, ionode_id
   USE input_gw,           ONLY : input_options
   USE basic_structures,   ONLY : q_mat, wannier_u, wp_psi,v_pot,free_memory
   USE green_function,     ONLY : green,read_green,free_memory_green, initialize_green
   USE polarization,       ONLY : polaw,free_memory_polaw,read_polaw, initialize_polaw, &
                                   &invert_v_pot,distribute_v_pot, collect_v_pot
   USE compact_product
   USE mp,                 ONLY : mp_sum, mp_barrier, mp_bcast
   USE para_gww,           ONLY : is_my_time, is_my_pola
   USE energies_gww,       ONLY : quasi_particles
   USE times_gw,           ONLY : times_freqs
   USE w_divergence
   USE mp_global,          ONLY : nproc,mpime
   USE io_files,           ONLY : find_free_unit, prefix

   implicit none


   TYPE(times_freqs), INTENT(in) :: tf!for times grid
   TYPE(input_options), INTENT(in) :: options! for imaginary time range and number of samples
   TYPE(self_storage) :: ss!
   TYPE(quasi_particles), INTENT(in) :: qp!for the HF energies if required

   TYPE(green)     :: gg,gm!green function
   TYPE(q_mat)     :: qm!overlap of orthonormalized wannier products with wannier products
   TYPE(polaw)     :: ww!dressed interaction
   TYPE(wannier_u) :: uu!transformation matrix ks to wannier
   TYPE(contraction) :: cr!to speed up calculation
   TYPE(contraction_index) :: cri! index of contraction
   TYPE(contraction_state) :: crs!state contraction data
   TYPE(wp_psi) :: wp!for remainder calculations
   TYPE(gv_time) :: gt!for the treatment of the G=0,G=0 divergence of W
   TYPE(v_pot) :: vp,vpi,vpid



   REAL(kind=DP) :: time,dt
   INTEGER       :: iw,ii,jj,it
   REAL(kind=DP) :: offset
   COMPLEX(kind=DP) :: sene
   INTEGER :: l_blk, nbegin,nend, iun
   REAL(kind=DP), ALLOCATABLE :: wtemp(:,:)
   LOGICAL :: exst

   nullify(vp%vmat)
   nullify(vpi%vmat)
   nullify(vpid%vmat)

   if(options%l_self_from_pola .or. options%l_self_beta) then
      if(options%w_divergence == 2) then
         call read_data_pw_v(vp,options%prefix,options%debug,0,.true.)
      else
         call read_data_pw_v(vp,options%prefix,options%debug,0,.false.)
      endif
      call invert_v_pot(vp,vpi)
      call free_memory(vp)
      call distribute_v_pot(vpi,vpid)
      call free_memory(vpi)
   endif

!set self_energy descriptor

!  call free_memory_self_storage(ss)
  ss%ontime=.true.
  ss%max_i=options%max_i
  ss%i_min=options%i_min
  ss%i_max=options%i_max
  ss%n=options%n
  ss%tau=options%tau
  ss%whole_s=options%whole_s

  if(tf%grid_fit/=0) then
     ss%n_grid_fit=tf%n_grid_fit
  else
     ss%n_grid_fit=tf%n
  endif

  if(ss%whole_s) then
   allocate(ss%whole(ss%max_i,ss%max_i,2*ss%n+1))
   ss%whole(:,:,:)=(0.d0,0.d0)
   allocate(ss%whole_freq_fit(ss%max_i,ss%max_i,2*ss%n_grid_fit+1))
   ss%whole_freq_fit(:,:,:)=(0.d0,0.d0)
   nullify(ss%diag)
   nullify(ss%diag_freq_fit)
  else
   allocate(ss%diag(ss%max_i,2*ss%n+1))
   ss%diag(:,:)=(0.d0,0.d0)
   nullify(ss%whole)
   allocate(ss%diag_freq_fit(ss%max_i,2*ss%n_grid_fit+1))
   ss%diag_freq_fit(:,:)=(0.d0,0.d0)
   nullify(ss%whole_freq_fit)
  endif

!set up self-energy remainders

  allocate(ss%ene_remainder(ss%max_i))
  if(options%remainder == 3 .or. options%remainder == 4) then
     ss%ene_remainder(:)=qp%ene_remainder(:)
  else
     ss%ene_remainder(:)=0.d0
  endif


  if(.not.options%lvcprim_file .and. .not.options%l_self_beta) then

!read U matrix
     call read_data_pw_u(uu,options%prefix)
!read overlap matrix Q
     call read_data_pw_q(qm,options%prefix, options%l_self_from_pola)
     dt = ss%tau/real(ss%n)
     if(options%use_contractions) then
        if(.not.options%l_contraction_single_state) then
           write(stdout,*) 'call do_contraction'!ATTENZIONE
           call  do_contraction(qm,uu,cr, options%max_i)
           write(stdout,*) 'done do_contraction'!ATTENZIONE
           call  write_contraction(cr,options)
           write(stdout,*) 'done do_contraction'!ATTENZIONE
        else
!contraction index and states already available on disk
           call read_contraction_index(cri, options)
        endif
     endif
     !loop
     call initialize_green(gg)
     call initialize_polaw(ww)

     l_blk= (2*ss%n+1)/nproc
     if(l_blk*nproc < (2*ss%n+1)) l_blk = l_blk+1
     nbegin=mpime*l_blk+1 -(ss%n+1)
     nend=nbegin+l_blk-1
     if(nend > ss%n) nend = ss%n

!     do iw=-ss%n,ss%n
!        if(is_my_time(iw)) then

     do iw=nbegin,nbegin+l_blk-1
        if(iw <= ss%n)  then

           write(stdout,*) 'Time :',iw!ATTENZIONE
           call flush_unit(stdout)
           time=dt*real(iw)
!read dressed interaction
!we take care of the symmetry t ==> -t
           call read_polaw(abs(iw),ww,options%debug)

!some controls
           if(.not. ww%ontime) then
              write(stdout,*) 'Routine create_self_ontime: imaginary time required'
              stop
           endif
           if(tf%l_fft_timefreq) then
              if(abs(time-ww%time) >= eps8) then
                 write(stdout,*) 'Routine create_self_ontime: imaginary time does not correspond'
                 stop
              endif
           endif

           if(options%l_self_from_pola) then
!if required  obtains the dressed polarization
              call collect_v_pot(vpi,vpid)
              allocate(wtemp(ww%numpw,ww%numpw))
              call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,vpi%vmat,ww%numpw,ww%pw,ww%numpw,&
                    &0.d0, wtemp,ww%numpw)
              call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,wtemp,ww%numpw,vpi%vmat,ww%numpw,&
                   &0.d0,ww%pw,ww%numpw)
              call free_memory(vpi)
              deallocate(wtemp)
           endif
           call read_green(iw,gg,options%debug,.false.)



!some controls
           if(.not. gg%ontime) then
              write(*,*) 'Routine create_self_ontime: imaginary time required'
              stop
           endif
           if(tf%l_fft_timefreq) then
              if(abs(time-gg%time) >= eps8) then
                 write(*,*) 'Routine create_self_ontime: imaginary time does not correspond'
                 stop
              endif
           endif
       !calculate elements
           if(ss%whole_s) then
              do ii=ss%i_min,ss%i_max
                 do jj=ss%i_min,ss%i_max
                    if(.not.options%use_contractions) then
                       call self_energy(ii,jj,ss%whole(ii,jj,iw+ss%n+1),time,qm,uu,gg,ww)
                    else
                       call self_energy_contraction(ii,jj,ss%whole(ii,jj,iw+ss%n+1),time,cr,gg,ww)
                    endif
                 enddo
              enddo
           else
              do ii=ss%i_min,ss%i_max
                 write(stdout,*) 'State:', ii
                 call flush_unit(stdout)
                 if(.not.options%use_contractions) then
                    call self_energy(ii,ii,ss%diag(ii,iw+ss%n+1),time,qm,uu,gg,ww)
                 else
                    if(.not.options%l_contraction_single_state) then
                       call self_energy_contraction(ii,ii,ss%diag(ii,iw+ss%n+1),time,cr,gg,ww)
                    else
                       crs%state=ii
                       call read_contraction_state(cri,crs,options)
                       call self_energy_contraction_state(ii,ii,ss%diag(ii,iw+ss%n+1),time,cri,crs,gg,ww)
                       call free_memory_contraction_state(crs)
                    endif
                 endif
              enddo
           endif

!at zero time 1/2 positive G and 1/2 negative

           if(iw==0) then

              do ii=1,ss%max_i
                 ss%diag(ii,iw+ss%n+1)=0.5d0*ss%diag(ii,iw+ss%n+1)
              enddo
              call read_green(iw,gg,options%debug,.true.)
              do ii=ss%i_min,ss%i_max
                  write(stdout,*) 'State:', ii
                 call flush_unit(stdout)
                 if(.not.options%use_contractions) then
                    call self_energy(ii,ii,sene,time,qm,uu,gg,ww)
                 else
                    if(.not.options%l_contraction_single_state) then
                       call self_energy_contraction(ii,ii,sene,time,cr,gg,ww)
                    else
                       crs%state=ii
                       call read_contraction_state(cri,crs,options)
                       call self_energy_contraction_state(ii,ii,sene,time,cri,crs,gg,ww)
                       call free_memory_contraction_state(crs)
                    endif
                 endif
                 ss%diag(ii,iw+ss%n+1)=ss%diag(ii,iw+ss%n+1)+0.5d0*sene
              enddo
           endif
        else
           if(options%l_self_from_pola) then
              call collect_v_pot(vpi,vpid)
              call free_memory(vpi)
           endif


        endif!on is_my_time
     enddo

     call free_memory(vpid)
!if required add remainder time to negative ones
     if(options%remainder==1 .or. options%remainder==2) then

        write(stdout,*) 'enter remainder'

        call read_data_pw_wp_psi(wp,options%prefix)



        if(.not.options%l_hf_energies) then
           if(uu%nums > uu%nums_occ) then
              offset=-(uu%ene(uu%nums_occ+1)+uu%ene(uu%nums_occ))/2.d0
           else
              offset=-uu%ene(uu%nums_occ)
           endif
        else
           if(uu%nums > uu%nums_occ) then
              offset=-(qp%ene_hf(uu%nums_occ+1)+qp%ene_hf(uu%nums_occ))/2.d0
           else
              offset=-qp%ene_hf(uu%nums_occ)
           endif
        endif

        call read_green(0,gg,options%debug,.false.)
        call read_green(0,gm,options%debug,.true.)

        do iw=-ss%n,0
           if(is_my_pola(-iw)) then
              write(stdout,*) 'Remainder time:', iw
              if(ss%whole_s) then
                 write(stdout,*) 'Routine create_self_ontime: remainder and whole matrix not implemented YET'
                 stop
              else
                 if(tf%l_fft_timefreq) then
                    time=dt*real(iw)
                 else
                    time=tf%times(iw)
                 endif
                 call read_polaw(iw,ww,options%debug)

                 do ii=ss%i_min,ss%i_max
                    if(options%l_contraction_single_state) then
                       crs%state=ii
                       call read_contraction_state(cri,crs,options)
                    endif
                    if(.not.options%use_contractions) then
                       call self_energy(ii,ii,sene,time,qm,uu,gg,ww)
                    else
                       if(.not.options%l_contraction_single_state) then
                          call self_energy_contraction(ii,ii,sene,time,cr,gg,ww)
                       else
                          call self_energy_contraction_state(ii,ii,sene,time,cri,crs,gg,ww)
                       endif
                    endif
                    if(options%remainder==1) then
                       if(.not.options%l_hf_energies) then
                          sene=sene*exp((uu%ene(uu%nums)+offset)*time)
                       else
                          sene=sene*exp((qp%ene_hf(qp%max_i)+offset)*time)
                       endif
                    endif
!sene changes sign because we are on the negative axes!!

                    if(iw==0) sene=sene*0.5d0

                    ss%diag(ii,iw+ss%n+1)=ss%diag(ii,iw+ss%n+1)+sene

                    write(stdout,*) 'SENE 0', iw, sene

                    if(.not.options%use_contractions) then
                       call self_energy(ii,ii,sene,time,qm,uu,gm,ww)
                    else
                      if(.not.options%l_contraction_single_state) then
                         call self_energy_contraction(ii,ii,sene,time,cr,gg,ww)
                      else
                         call self_energy_contraction_state(ii,ii,sene,time,cri,crs,gg,ww)
                      endif
                   endif
                   if(options%remainder==1) then
                      if(.not.options%l_hf_energies) then
                         sene=sene*exp((uu%ene(uu%nums)+offset)*time)
                      else
                         sene=sene*exp((qp%ene_hf(qp%max_i)+offset)*time)
                      endif
                   endif
                   if(iw==0) sene=sene*0.5d0


                   ss%diag(ii,iw+ss%n+1)=ss%diag(ii,iw+ss%n+1)-sene
                   write(stdout,*) 'SENE 1', iw, sene

                   call self_energy_remainder(ii,sene,time,wp,ww)
                   if(options%remainder==1) then
                      if(.not.options%l_hf_energies) then
                         sene=sene*exp((uu%ene(uu%nums)+offset)*time)
                      else
                         sene=sene*exp((qp%ene_hf(qp%max_i)+offset)*time)
                      endif
                   endif
                   if(iw==0) sene=sene*0.5d0

                   ss%diag(ii,iw+ss%n+1)=ss%diag(ii,iw+ss%n+1)+sene

                   write(stdout,*) 'SENE 2', iw, sene
                   if(options%l_contraction_single_state) &
                        & call free_memory_contraction_state(crs)
                enddo
             endif
          endif
       enddo
       call free_memory(wp)
    endif


    if(ss%whole_s) then
       call mp_sum(ss%whole(:,:,:))
    else
       call mp_sum(ss%diag(:,:))
    end if

    call free_memory(uu)
    call free_memory(qm)
    call free_memory_polaw(ww)
    call free_memory_green(gg)
    if(.not.options%l_contraction_single_state) &
         & call free_memory_contraction(cr)

 else
   !FROM VCPRIM FILE
   !!! call selfenergy_ontime_file(ss,tf,options)
   !!! NEW VERSION
   if (ionode) then
      write(stdout,*) 'check if the file prefix.self_energie_negative_time exist'
      call flush_unit(stdout)
      inquire(file=trim(options%prefix)//'.self_energie_negative_time', exist = exst)
      write(stdout,*) 'exist=', exst
      call flush_unit(stdout)
   endif

   call mp_bcast(exst, ionode_id)

   if (.not. exst) then
       write(stdout,*) 'Do not exist => start first with the negative times and then do the positive times'
       call flush_unit(stdout)
       ! nt= negative_time, pt=positive_time
       ! name shortened to prevent trouble on some compilers - PG
       call selfenergy_on_split_nt_file(ss, tf, options)
       call selfenergy_on_split_pt_file(ss, tf, options)
    else
       write(stdout,*) 'Exist => Read the file from the negative times and then do the positive times'
       call flush_unit(stdout)
       if (ionode) then
          iun = find_free_unit()
          open( unit=iun, file=trim(options%prefix)//'.self_energie_negative_time', status='unknown',form='unformatted')
          write(stdout,*) 'Double loop from ss%i_min=', ss%i_min, ' to ss%i_max=', ss%i_max
          write(stdout,*) 'from it=-ss%n=', -ss%n, ' to 0'
          call flush_unit(stdout)
          do ii=ss%i_min,ss%i_max
             do it=-ss%n,0
                 read(iun) ss%diag(ii,it+ss%n+1)
             enddo
          enddo
          close (iun)
       endif
       write(stdout,*) 'Before mp_bcast and ionode_id=', ionode_id
       call flush_unit(stdout)
       do ii=ss%i_min,ss%i_max
          write(stdout,*) 'ss%diag(ii,1)=',ss%diag(ii,1)
          write(stdout,*) 'ss%diag(ii,ss%n+1)=',ss%diag(ii,ss%n+1)
       enddo
       call flush_unit(stdout)

       call mp_bcast(ss%diag(:,:), ionode_id)

       write(stdout,*) 'After mp_bcast'
       call flush_unit(stdout)
       do ii=ss%i_min,ss%i_max
          write(stdout,*) 'ss%diag(ii,1)=',ss%diag(ii,1)
          write(stdout,*) 'ss%diag(ii,ss%n+1)=',ss%diag(ii,ss%n+1)
       enddo
       call flush_unit(stdout)


       write(stdout,*) 'File read, do positive times'
       call flush_unit(stdout)
       call selfenergy_on_split_pt_file(ss, tf, options)
       write(stdout,*) 'Do positive times'
       call flush_unit(stdout)
    endif
 endif

!if required add coulomb-like term for the treatment of the (G=0,G=0) divergence of W
 if(options%w_divergence == 2 ) then
    call initialize_gv_time(gt)
    call read_gv_time(gt)
    !consistency check
     if(options%max_i /= gt%max_i) then
        write(stdout,*) 'max_i not correct'
        stop
     endif
     call setup_gv_time(gt)

     do iw=1,2*gt%n+1
        ss%diag(:,iw)=ss%diag(:,iw)+gt%ex(:,iw)
     enddo
     call free_memory_gv_time(gt)
  endif

  return

END SUBROUTINE create_self_ontime


  SUBROUTINE write_storage(tf,ss)
!this subroutine write on standard output
!the values of write_storage
    USE io_global,      ONLY : stdout, ionode
    USE constants,      ONLY : pi
    USE mp,             ONLY : mp_barrier
    USE times_gw,       ONLY : times_freqs

    implicit none

    TYPE(times_freqs), INTENT(in) :: tf!for time grid
    TYPE(self_storage), INTENT(in) :: ss
    INTEGER :: iw,ii,jj
    REAL(kind=DP) :: time,dt,totalfrequency,totalperiod,omega

    if(ionode) then
       if(ss%ontime) then
          write(stdout,*) '--------Sigma on imaginary time----------'
          dt=ss%tau/real(ss%n)
          do iw=-ss%n,ss%n
             if(tf%l_fft_timefreq) then
                time=dt*real(iw)
             else
                time=tf%times(iw)
             endif
             if(ss%whole_s) then
                do ii=1,ss%max_i
                   do jj=1,ss%max_i
                      write(stdout,*) time,ii,jj,ss%whole(ii,jj,iw+ss%n+1)
                   enddo
                enddo
             else
                do ii=1,ss%max_i
                   write(stdout,*) iw, time,ii, ss%diag(ii,iw+ss%n+1)
                enddo
             endif
          enddo
       else
          write(stdout,*) '--------Sigma on imaginary frequency----------'
          totalperiod=2.d0*ss%tau+2.d0*ss%tau/real(ss%n)
          totalfrequency=(2.d0*pi/totalperiod)

          do iw=-ss%n,ss%n
             if(tf%l_fft_timefreq) then
                omega=totalfrequency*real(iw)
             else
                omega=tf%freqs(iw)
             endif
             if(ss%whole_s) then
                do ii=1,ss%max_i
                   do jj=1,ss%max_i
                      write(stdout,*) omega,ii,jj,ss%whole(ii,jj,iw+ss%n+1)
                   enddo
                enddo
             else
                do ii=1,ss%max_i
                   write(stdout,*) omega,ii, ss%diag(ii,iw+ss%n+1)
                enddo
             endif
          enddo
       endif
    endif
    return

  END SUBROUTINE write_storage

  SUBROUTINE fft_storage_grid(tf,ss)
!this subroutine performs a FFT on the storage data
!inverse or direct determined by ontime
!uses grid

   USE io_global, ONLY : stdout
   USE constants, ONLY : pi
   USE times_gw,  ONLY : times_freqs


   implicit none

   TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids
   TYPE(self_storage), INTENT(inout) :: ss!input data

   COMPLEX(kind=DP), DIMENSION(:), ALLOCATABLE :: ss_old, tmpc
   COMPLEX(kind=DP),  DIMENSION(:,:), ALLOCATABLE :: factors
   INTEGER :: ii,jj, is,js,kk


   INTEGER, PARAMETER :: nmesh=30

   REAL(kind=DP) :: b_p,b_m,r_p,r_m
   COMPLEX(kind=DP) :: a_p,a_m, cor_1,cor_2
   REAL(kind=DP), ALLOCATABLE :: x(:),w(:)
   COMPLEX(kind=DP), ALLOCATABLE ::  fij(:,:), fp(:),fm(:)


   allocate(ss_old(2*tf%n+1), tmpc(2*tf%n+1))
   allocate(factors(-tf%n:tf%n, -tf%n:tf%n))

!setup factors for every time posistion
   do ii=-tf%n, tf%n
      if(ss%ontime) then!time to frequency transform
         do jj=-tf%n,tf%n
             factors(jj,ii)=tf%weights_time(jj)*exp((0.d0,-1.d0)*tf%freqs(ii)*tf%times(jj))
         enddo
         factors(:,ii)=factors(:,ii)*(0.d0,-1.d0)
      else!frequency to time transform
         do jj=-tf%n,tf%n
            factors(jj,ii)=tf%weights_freq(jj)*exp((0.d0,1.d0)*tf%times(ii)*tf%freqs(jj))
         enddo
         factors(:,ii)=factors(:,ii)*(0.d0,1.d0)/(2.d0*pi)
      endif
   enddo
   if(ss%whole_s) then!full matrix
      do is=1,ss%max_i
         do js=1,ss%max_i
!copy array to be transformed
            ss_old(:)=ss%whole(is,js,:)
!transform
            do ii=-tf%n,tf%n
               do kk=-tf%n,tf%n
                  tmpc(kk+tf%n+1)=ss_old(kk+tf%n+1)*factors(kk,ii)
               enddo
               ss%whole(is,js,ii+tf%n+1)=sum(tmpc(1:2*tf%n+1))
            enddo
         enddo
      enddo
   else
      if(tf%l_fourier_fit_time .and. ss%ontime) then
         allocate(fij(-tf%n:tf%n,nmesh))
         allocate(fp(nmesh),fm(nmesh))
         allocate(x(nmesh),w(nmesh))
         x(:)=0.d0
         w(:)=0.d0
         call legzo(nmesh,x,w)

         !x(:)=x(:)*tf%tau/2.d0
         !x(:)=x(:)+tf%tau/2.d0
         !w(:)=w(:)*tf%tau/2.d0

         x(:)=x(:)*(tf%times(tf%n)-tf%tau)/2.d0
         x(:)=x(:)+(tf%times(tf%n)-tf%tau)/2.d0+tf%tau
         w(:)=w(:)*(tf%times(tf%n)-tf%tau)/2.d0

         do ii=-tf%n,tf%n
            do jj=1,nmesh
               fij(ii,jj)=exp((0.d0,-1.d0)*tf%freqs(ii)*x(jj))
            enddo
         enddo
      endif
      do is=1,ss%max_i
 !copy array to be transformed
         ss_old(:)=ss%diag(is,:)
!transform
         do ii=-tf%n,tf%n
            do kk=-tf%n,tf%n
               tmpc(kk+tf%n+1)=ss_old(kk+tf%n+1)*factors(kk,ii)
            enddo
            ss%diag(is,ii+tf%n+1)=sum(tmpc(1:2*tf%n+1))
         enddo

         if(tf%l_fourier_fit_time .and. ss%ontime) then
            r_p=dble(ss_old(2*tf%n)/ss_old(2*tf%n+1))
            write(stdout,*) 'RP',ss_old(2*tf%n),ss_old(2*tf%n-5)
            if(r_p <= 1.d0) r_p = tf%g_tau
            b_p=log(r_p)/(tf%times(tf%n)-tf%times(tf%n-1))
            a_p=ss_old(2*tf%n)/(exp(-b_p*tf%times(tf%n-1)))
            if(r_p == tf%g_tau) a_p=0.d0

            if(abs(ss_old(2)) > 1.d-10 .and. abs(ss_old(1)) > 1.d-10) then
               r_m=dble(ss_old(2)/ss_old(1))
               if(r_m <= 1.d0) r_m = tf%g_tau
               b_m=log(r_m)/(tf%times(-tf%n+1)-tf%times(-tf%n))
               a_m=ss_old(2)/(exp(b_m*tf%times(-tf%n+1)))
               if(r_m == tf%g_tau) a_m=0.d0
            else
               r_m=0.d0
               a_m=(0.d0,0.d0)
               b_m=0.d0
            endif


            do jj=1,nmesh
               fp(jj)=a_p*exp(-b_p*x(jj))*w(jj)
            enddo
            if(r_m /=0.d0) then
               do jj=1,nmesh
                  fm(jj)=a_m*exp(-b_m*x(jj))*w(jj)
               enddo
            endif

            do ii=-tf%n,tf%n
             !  cor_1=(0.d0,-1.d0)*(a_p/(b_p+(0.d0,1.d0)*tf%freqs(ii)))
             !  if(r_m /= 0.d0) then
             !     cor_1=cor_1+(0.d0,-1.d0)*(a_m/(b_m-(0.d0,1.d0)*tf%freqs(ii)))
             !  endif
               cor_2=0.d0
               do jj=1,nmesh
                  cor_2=cor_2-fij(ii,jj)*fp(jj)
                  if(r_m /=0.d0) then
                     cor_2=cor_2-conjg(fij(ii,jj))*fm(jj)
                  endif
               enddo
               cor_2=cor_2*(0.d0,-1.d0)
               ss%diag(is,ii+tf%n+1)=ss%diag(is,ii+tf%n+1)!-cor_2!+cor_1+cor_2
               write(stdout,*) 'COR2' , cor_2
            enddo
         endif
      enddo
      if(tf%l_fourier_fit_time .and.ss%ontime) deallocate(fij,fp,fm,x,w)
   endif
   if(ss%ontime) then
      ss%ontime=.false.
   else
      ss%ontime=.true.
   endif


   deallocate(ss_old,tmpc)
   deallocate(factors)

   return
 END SUBROUTINE fft_storage_grid

  SUBROUTINE fft_storage(ss)
!this subroutine performs a FFT on the storage data
!inverse or direct determined by ontime
   USE io_global, ONLY : stdout
   USE constants, ONLY : pi
   USE fft_scalar, ONLY : cft_1z

   implicit none

   TYPE(self_storage) :: ss!input data



   REAL(kind=DP) :: totalperiod,omega,time,totalfrequency
   INTEGER :: iw,ii,ipos
   COMPLEX(kind=DP), ALLOCATABLE :: inz(:),outz(:)
   COMPLEX(kind=DP) :: fact
   INTEGER*8 :: plan

   totalperiod=2.d0*ss%tau+2.d0*ss%tau/real(ss%n)
   totalfrequency=(2.d0*pi/totalperiod)*real(2*ss%n+2)
   allocate(inz(2*ss%n+2),outz(2*ss%n+2))

   if(.not.ss%whole_s) then
     if(ss%ontime) then!time to frequency transformation
       ss%ontime=.false.
!loop on states
       do ii=1,ss%max_i

         inz(:)=(0.d0,0.d0)
         do iw=-ss%n,ss%n
           ipos=iw+ss%n+2
           inz(ipos)=ss%diag(ii,iw+ss%n+1)
         enddo
         inz(1)=inz(2)

         call cft_1z(inz,1,2*ss%n+2,2*ss%n+2, -1,outz)
         outz(:)=outz(:)*dble(2*ss%n+2)

         do iw=0,2*ss%n+2-1
           if(iw <= (2*ss%n+1)) then
             omega=(2.d0*pi/totalperiod)*real(iw)
           else
             omega=(2.d0*pi/totalperiod)*real(iw-2*ss%n-2)
           endif
           fact=exp((0.d0,-1.d0)*omega*totalperiod/2.d0)*(0.d0,-1.d0)*(ss%tau/real(ss%n))
           outz(iw+1)=outz(iw+1)*fact
         enddo

         do iw=0,2*ss%n+1
           if(iw/=(ss%n+1)) then
             if(iw < (ss%n+1)) then
               ss%diag(ii,ss%n+iw+1)=outz(iw+1)
             else
               ss%diag(ii,iw-ss%n-2+1)=outz(iw+1)
             endif
           endif
         enddo
         write(*,*) 'ELIMINATO:', outz(ss%n+1)
       enddo

   else !frequency to time transform
       ss%ontime=.true.
!loop on states
       do ii=1,ss%max_i

         inz(:)=(0.d0,0.d0)
         do iw=-ss%n,ss%n
           ipos=iw+ss%n+2
           inz(ipos)=ss%diag(ii,iw+ss%n+1)
         enddo



         call cft_1z(inz,1,2*ss%n+2,2*ss%n+2, 1,outz)

         do iw=0,2*ss%n+2-1
           if(iw <= (2*ss%n+1)) then
             time=(ss%tau/real(ss%n))*real(iw)
           else
             time=(ss%tau/real(ss%n))*real(iw-2*ss%n-2)
           endif
           fact=exp((0.d0,+1.d0)*time*totalfrequency/2.d0)*(0.d0,+1.d0)/totalperiod
           outz(iw+1)=outz(iw+1)*fact
         enddo


         do iw=0,2*ss%n+1
           if(iw/=(ss%n+1)) then
             if(iw < (ss%n+1)) then
               ss%diag(ii,ss%n+iw+1)=outz(iw+1)
             else
               ss%diag(ii,iw-ss%n-2+1)=outz(iw+1)
             endif
           endif
         enddo

       enddo


   endif

   else
     write(stdout,*) 'ENTIRE SIGMA NOT IMPLEMENTED YET'
   endif
  deallocate(inz,outz)
  return

  END SUBROUTINE

  SUBROUTINE test_fft(tf)
!just a fft test

  USE times_gw, ONLY : times_freqs

  implicit none

  TYPE(times_freqs), INTENT(in) :: tf

  TYPE(self_storage) :: ss
  INTEGER :: n,iw
  REAL(kind=DP) :: tau, lambda

  n=100
  tau=25.
  lambda=2.


  ss%ontime=.true.
  ss%whole_s=.false.
  ss%n=tf%n
  ss%tau=tf%tau
  ss%max_i=1
  n=tf%n
  allocate(ss%diag(1,2*n+1))
  nullify(ss%whole)

  do iw=-n,n
    ss%diag(1,iw+n+1)=exp(-(real(iw)*tau/real(n)/lambda)**2.)
  enddo

  call write_storage(tf,ss)

  call fft_storage(ss)

  call write_storage(tf,ss)

  call fft_storage(ss)

  call write_storage(tf,ss)

  call free_memory_self_storage(ss)
  return
  END SUBROUTINE

  SUBROUTINE addconduction_self_ontime(ss, options)
!this subroutine adds to the self_energy of conduction states
!on negative imaginary times, the part due to terms \Psi_c'\Psic\w_P

    USE io_global,     ONLY : stdout, ionode
    USE input_gw,      ONLY : input_options
    USE basic_structures,  ONLY : v_pot,wannier_u_prim, v_pot_prim,free_memory, ortho_polaw
    USE green_function,    ONLY : green, read_green, free_memory_green, initialize_green
    USE polarization,      ONLY : polaw, free_memory_polaw, read_polaw, invert_v_pot, invert_ortho_polaw,&
         & orthonormalize_inverse, initialize_polaw, orthonormalize_vpot, distribute_ortho_polaw, collect_ortho_polaw,&
         &  distribute_v_pot, collect_v_pot
    USE mp,                ONLY : mp_sum
    USE para_gww,          ONLY : is_my_pola
    USE mp_global,            ONLY : nproc,mpime

    implicit none

    TYPE(input_options) :: options
    TYPE(self_storage)  :: ss

    TYPE(v_pot) :: vp,vpi,vpid
    TYPE(ortho_polaw) :: op,opi, opd, opid
    TYPE(polaw) :: ww!dressed interaction
    TYPE(wannier_u_prim) :: wup
    TYPE(v_pot_prim) :: vpp
    TYPE(green) :: gg

    INTEGER iw,jw,kw,it,ii
    REAL(kind=DP), ALLOCATABLE :: wtemp(:,:)
    REAL(kind=DP), ALLOCATABLE :: cp(:,:,:) !arrys for contraction c',c, numpw
    REAL(kind=DP), ALLOCATABLE :: qg(:,:)
    COMPLEX(kind=DP), ALLOCATABLE :: sene(:,:)
    REAL(kind=DP), ALLOCATABLE :: gf_t(:,:)
    REAL(kind=DP), ALLOCATABLE :: pwcp_t(:,:)
    REAL(kind=DP), EXTERNAL :: ddot

    INTEGER :: l_blk, nbegin,nend
    INTEGER  :: i_first

    nullify(vp%vmat)
    nullify(vpi%vmat)
    nullify(op%on_mat)
    nullify(opi%on_mat)
    nullify(wup%umat)
    nullify(vpp%ij)
    nullify(vpp%vmat)
    nullify(opd%on_mat)
    nullify(opid%on_mat)
    nullify(vpid%vmat)

    call initialize_green(gg)
    call initialize_polaw(ww)
    write(stdout,*) 'addconduction_self_ontime OLD 1'!ATTENZIONE
    call flush_unit(stdout)


!read coulombian potential and calculate inverse

    if(ss%whole_s) then
       write(stdout,*) 'Whole s not implemented YET'
       stop
    endif




    call read_data_pw_u_prim(wup,options%prefix)

    i_first=max(ss%i_min,wup%nums_occ+1)

    write(stdout,*) 'addconduction_self_ontime1_2'!ATTENZIONE
    call flush_unit(stdout)

    if(options%w_divergence==2) then
       call read_data_pw_v_pot_prim(vpp, options%prefix,.true.)
    else
       call read_data_pw_v_pot_prim(vpp, options%prefix,.false.)
    endif


    write(stdout,*) 'addconduction_self_ontime1_3'!ATTENZIONE
     call flush_unit(stdout)

   allocate(sene(-ss%n:0,ss%i_max-wup%nums_occ))
   sene(:,:)=(0.d0,0.d0)


!set up contraction array \sum_j U^{C'}_ij Vjkl

  !  allocate(cp(vpp%numpw, wup%nums-wup%nums_occ,options%max_i-wup%nums_occ))
   allocate(cp(vpp%numpw, wup%nums-wup%nums_occ,i_first:ss%i_max))
   cp(:,:,:)=0.d0
    do iw=1,vpp%numpw_prim
       do ii=i_first,ss%i_max
          do kw=1,vpp%numpw
             cp(kw,vpp%ij(2,iw)-wup%nums_occ,ii)=cp(kw,vpp%ij(2,iw)-wup%nums_occ,ii)+&
                  &dble(wup%umat(ii-wup%nums_occ,vpp%ij(1,iw)))*vpp%vmat(iw,kw)
          enddo
       enddo
    enddo

    call free_memory(vpp)
    call free_memory(wup)!in this way only the data is deallocated

     write(stdout,*) 'addconduction_self_ontime1_4'!ATTENZIONE
     call flush_unit(stdout)

    if(options%w_divergence == 2) then
       call read_data_pw_v(vp,options%prefix,options%debug,0,.true.)
    else
       call read_data_pw_v(vp,options%prefix,options%debug,0,.false.)
    endif
    if(options%lnonorthogonal) then
       call read_data_pw_ortho_polaw(op,options%prefix)
       call orthonormalize_vpot(op,vp)
    endif
    call invert_v_pot(vp,vpi)
    call free_memory(vp)
    write(stdout,*) 'addconduction_self_ontime1_45'
    call distribute_v_pot(vpi,vpid)
    call free_memory(vpi)
    if(options%lnonorthogonal) then
       call invert_ortho_polaw(op,opi)
       write(stdout,*) 'addconduction_self_ontime1_5 op',op%numpw!ATTENZIONE
       call distribute_ortho_polaw(op,opd)
       call free_memory(op)
       write(stdout,*) 'addconduction_self_ontime1_6 opd',opd%numpw!ATTENZIONE

       call distribute_ortho_polaw(opi,opid)
       call free_memory(opi)
    endif


    l_blk= (ss%n+1)/nproc
     if(l_blk*nproc < (ss%n+1)) l_blk = l_blk+1
     nbegin=mpime*l_blk+1 -(ss%n+1)
     nend=nbegin+l_blk-1
     if(nend > 0) nend = 0




    write(stdout,*) 'addconduction_self_ontime5',nbegin,l_blk!ATTENZIONE
    call flush_unit(stdout)
!loop on negative imaginary times
    do it=nbegin,nbegin+l_blk-1
       if(it <= 0)  then
          write(stdout,*) 'addconduction_self_ontime time', it!ATTENZIONE
          call flush_unit(stdout)
!we take care of the symmetru t ==> -t

          call read_polaw(abs(it),ww,options%debug)
          write(stdout,*) 'addconduction_self_ontime6 ww', ww%numpw!ATTENZIONE
           if(options%lnonorthogonal) then
              call collect_ortho_polaw(opi,opid)
              write(stdout,*) 'dimensions', opi%numpw, opid%numpw
              call orthonormalize_inverse(opi,ww)
              call free_memory(opi)
           endif
           write(stdout,*) 'addconduction_self_ontime7'!ATTENZIONE
          allocate(wtemp(ww%numpw,ww%numpw))

          call collect_v_pot(vpi,vpid)

          call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,vpi%vmat,ww%numpw,ww%pw,ww%numpw,&
               &0.d0, wtemp,ww%numpw)



          call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,wtemp,ww%numpw,vpi%vmat,ww%numpw,&
              &0.d0,ww%pw,ww%numpw)

          call free_memory(vpi)

          deallocate(wtemp)

          if(options%lnonorthogonal) then
             call collect_ortho_polaw(op,opd)
             call orthonormalize_inverse(op,ww)
             call free_memory(op)
          endif
           write(stdout,*) 'addconduction_self_ontime8'!ATTENZIONE
           call flush_unit(stdout)
          call read_green(it,gg,options%debug,.true.)

          allocate(gf_t(wup%nums-wup%nums_occ,wup%nums-wup%nums_occ))
          do iw=1,(wup%nums-wup%nums_occ)
             do jw=1,(wup%nums-wup%nums_occ)
                gf_t(jw,iw)=gg%gf_p(jw+wup%nums_occ, iw+wup%nums_occ)
             enddo
          enddo

          do ii=i_first,ss%i_max
             write(stdout,*) 'II' , ii
             call flush_unit(stdout)
             allocate(qg(ww%numpw,wup%nums-wup%nums_occ))
             call dgemm('N','N',ww%numpw,wup%nums-wup%nums_occ,wup%nums-wup%nums_occ,1.d0,&
                & cp(:,:,ii),ww%numpw,gf_t,wup%nums-wup%nums_occ,0.d0,qg,ww%numpw)


             allocate(pwcp_t(ww%numpw,wup%nums-wup%nums_occ))
             call dgemm('N','N',ww%numpw,wup%nums-wup%nums_occ,ww%numpw,1.d0,&
              &ww%pw,ww%numpw,cp(:,:,ii),ww%numpw,0.d0,pwcp_t,ww%numpw)

             do iw=1,(wup%nums-wup%nums_occ)
                sene(it,ii-wup%nums_occ)=sene(it,ii-wup%nums_occ)+&
                     &ddot(ww%numpw,qg(:,iw),1,pwcp_t(:,iw),1)*gg%factor*ww%factor
             enddo

             deallocate(pwcp_t)
             deallocate(qg)

             sene(it,ii-wup%nums_occ)=sene(it,ii-wup%nums_occ)*(0.d0,1.d0)
             if(it==0) sene(it,ii-wup%nums_occ)=sene(it,ii-wup%nums_occ)*0.5d0
             write(stdout,*) 'Conduction contribution', it,ii, sene(it,ii-wup%nums_occ)
          enddo
          deallocate(gf_t)
       else
          if(options%lnonorthogonal) then
             call collect_ortho_polaw(opi,opid)
             call free_memory(opi)
          endif
          call collect_v_pot(vpi,vpid)
          call free_memory(vpi)
          if(options%lnonorthogonal) then
             call collect_ortho_polaw(op,opd)
             call free_memory(op)
          endif
       endif
    enddo
    call mp_sum(sene(-ss%n:0,:))
    do ii=1,ss%i_max-wup%nums_occ
       do it=-ss%n,0
          ss%diag(ii+wup%nums_occ,it+ss%n+1)=ss%diag(ii+wup%nums_occ, it+ss%n+1)+sene(it,ii)
       enddo
    enddo
!!!!!!!!!!!
    call free_memory(vpid)
    if(options%lnonorthogonal) then
       call free_memory(opd)
       call free_memory(opi)
       call free_memory(opid)
    endif
    call free_memory_polaw(ww)
    call free_memory_green(gg)
    deallocate(cp)
    deallocate(sene)
    return
  END SUBROUTINE addconduction_self_ontime

  SUBROUTINE fft_storage_grid_fit(tf,ss)
!this subroutine performs a FFT from time to frequency on the storage data
!from W,P grid to fit grid

   USE io_global, ONLY : stdout
   USE constants, ONLY : pi
   USE times_gw,  ONLY : times_freqs


   implicit none

   TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids
   TYPE(self_storage), INTENT(inout) :: ss!input data

   COMPLEX(kind=DP), DIMENSION(:), ALLOCATABLE :: tmpc
   COMPLEX(kind=DP),  DIMENSION(:,:), ALLOCATABLE :: factors
   INTEGER :: ii,jj, is,js,kk


   INTEGER, PARAMETER :: nmesh=30

   REAL(kind=DP) :: b_p,b_m,r_p,r_m
   COMPLEX(kind=DP) :: a_p,a_m, cor_1,cor_2
   REAL(kind=DP), ALLOCATABLE :: x(:),w(:)
   COMPLEX(kind=DP), ALLOCATABLE ::  fij(:,:), fp(:),fm(:)





   allocate(tmpc(2*tf%n+1))
   allocate(factors(-tf%n:tf%n, -tf%n_grid_fit:tf%n_grid_fit))




!setup factors for every time position
   do ii=-tf%n_grid_fit, tf%n_grid_fit
      if(ss%ontime) then!time to frequency transform
         do jj=-tf%n,tf%n
             factors(jj,ii)=tf%weights_time(jj)*exp((0.d0,-1.d0)*tf%freqs_fit(ii)*tf%times(jj))
         enddo
         factors(:,ii)=factors(:,ii)*(0.d0,-1.d0)
      else
         write(stdout,*) 'ONLY TIME TO FREQUENCY FIT GRID IMPLEMENTED'
         stop
      endif
   enddo
   if(ss%whole_s) then!full matrix
      do is=1,ss%max_i
         do js=1,ss%max_i
!copy array to be transformed
            !transform
            do ii=-tf%n_grid_fit,tf%n_grid_fit
               do kk=-tf%n,tf%n
                  tmpc(kk+tf%n+1)=ss%whole(is,js,kk+tf%n+1)*factors(kk,ii)
               enddo
               ss%whole_freq_fit(is,js,ii+tf%n_grid_fit+1)=sum(tmpc(1:2*tf%n+1))
            enddo
         enddo
      enddo
   else
      if(tf%l_fourier_fit_time .and. ss%ontime) then
         allocate(fij(-tf%n_grid_fit:tf%n_grid_fit,nmesh))
         allocate(fp(nmesh),fm(nmesh))
         allocate(x(nmesh),w(nmesh))
         x(:)=0.d0
         w(:)=0.d0
         call legzo(nmesh,x,w)

         !x(:)=x(:)*tf%tau/2.d0
         !x(:)=x(:)+tf%tau/2.d0
         !w(:)=w(:)*tf%tau/2.d0

         x(:)=x(:)*(tf%times(tf%n)-tf%tau)/2.d0
         x(:)=x(:)+(tf%times(tf%n)-tf%tau)/2.d0+tf%tau
         w(:)=w(:)*(tf%times(tf%n)-tf%tau)/2.d0

         do ii=-tf%n_grid_fit,tf%n_grid_fit
            do jj=1,nmesh
               fij(ii,jj)=exp((0.d0,-1.d0)*tf%freqs_fit(ii)*x(jj))
            enddo
         enddo
      endif
      do is=1,ss%max_i

!transform
         do ii=-tf%n_grid_fit,tf%n_grid_fit
            do kk=-tf%n,tf%n
               tmpc(kk+tf%n+1)=ss%diag(is,kk+tf%n+1)*factors(kk,ii)
            enddo
            ss%diag_freq_fit(is,ii+tf%n_grid_fit+1)=sum(tmpc(1:2*tf%n+1))
         enddo


         if(tf%l_fourier_fit_time .and. ss%ontime) then
            r_p=dble(ss%diag(is,2*tf%n)/ss%diag(is,2*tf%n+1))
            if(r_p <= 1.d0) r_p = tf%g_tau
            b_p=log(r_p)/(tf%times(tf%n)-tf%times(tf%n-1))
            a_p=ss%diag(is,2*tf%n)/(exp(-b_p*tf%times(tf%n-1)))
            if(r_p == tf%g_tau) a_p=0.d0

            if(abs(ss%diag(is,2)) > 1.d-10 .and. abs(ss%diag(is,1)) > 1.d-10) then
               r_m=dble(ss%diag(is,2)/ss%diag(is,1))
               if(r_m <= 1.d0) r_m = tf%g_tau
               b_m=log(r_m)/(tf%times(-tf%n+1)-tf%times(-tf%n))
               a_m=ss%diag(is,2)/(exp(b_m*tf%times(-tf%n+1)))
               if(r_m == tf%g_tau) a_m=0.d0
            else
               r_m=0.d0
               a_m=(0.d0,0.d0)
               b_m=0.d0
            endif


            do jj=1,nmesh
               fp(jj)=a_p*exp(-b_p*x(jj))*w(jj)
            enddo
            if(r_m /=0.d0) then
               do jj=1,nmesh
                  fm(jj)=a_m*exp(-b_m*x(jj))*w(jj)
               enddo
            endif

            do ii=-tf%n_grid_fit,tf%n_grid_fit
             !  cor_1=(0.d0,-1.d0)*(a_p/(b_p+(0.d0,1.d0)*tf%freqs_fit(ii)))
             !  if(r_m /= 0.d0) then
             !     cor_1=cor_1+(0.d0,-1.d0)*(a_m/(b_m-(0.d0,1.d0)*tf%freqs_fit(ii)))
             !  endif
               cor_2=0.d0
               do jj=1,nmesh
                  cor_2=cor_2-fij(ii,jj)*fp(jj)
                  if(r_m /=0.d0) then
                     cor_2=cor_2-conjg(fij(ii,jj))*fm(jj)
                  endif
               enddo
               cor_2=cor_2*(0.d0,-1.d0)
               ss%diag_freq_fit(is,ii+tf%n_grid_fit+1)=ss%diag_freq_fit(is,ii+tf%n_grid_fit+1)!-cor_2!+cor_1+cor_2
            enddo
         endif
      enddo
      if(tf%l_fourier_fit_time .and.ss%ontime) deallocate(fij,fp,fm,x,w)
   endif

   deallocate(tmpc)
   deallocate(factors)

   return
 END SUBROUTINE fft_storage_grid_fit

  SUBROUTINE addconduct_self_ontime_file(ss, tf ,options)
!this subroutine adds to the self_energy of conduction states
!on negative imaginary times, the part due to terms \Psi_c'\Psic\w_P
!using terms from file

    USE io_global,     ONLY : stdout, ionode
    USE input_gw,      ONLY : input_options
    USE basic_structures,  ONLY : v_pot,wannier_u,free_memory, ortho_polaw,initialize_memory,cprim_prod
    USE green_function,    ONLY : green, read_green, free_memory_green, initialize_green
    USE polarization,      ONLY : polaw, free_memory_polaw, read_polaw, invert_v_pot, invert_ortho_polaw,&
         & orthonormalize_inverse, initialize_polaw, orthonormalize_vpot, distribute_ortho_polaw, collect_ortho_polaw,&
         &  distribute_v_pot, collect_v_pot
    USE mp,                ONLY : mp_sum
    USE para_gww,          ONLY : is_my_pola
    USE mp_global,            ONLY : nproc,mpime
    USE times_gw,  ONLY : times_freqs


   implicit none

   TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids
   TYPE(input_options) :: options
   TYPE(self_storage)  :: ss

   TYPE(v_pot) :: vp,vpi,vpid
   TYPE(ortho_polaw) :: op,opi, opd, opid
   TYPE(polaw) :: ww!dressed interaction
   TYPE(wannier_u) :: wu!structure to be read and initialized
   TYPE(cprim_prod) :: cpp!the producs c' c' v  wp

   INTEGER iw,jw,kw,it,ii,jj
   INTEGER :: l_blk, nbegin,nend
   REAL(kind=DP) :: offset
   COMPLEX(kind=DP), ALLOCATABLE :: sene(:,:)
   REAL(kind=DP), ALLOCATABLE:: wtemp(:,:), vtemp(:)
   REAL(kind=DP), EXTERNAL :: ddot
   LOGICAL :: ok_read

   nullify(vp%vmat)
   nullify(vpi%vmat)
   nullify(op%on_mat)
   nullify(opi%on_mat)
   nullify(opd%on_mat)
   nullify(opid%on_mat)
   nullify(vpid%vmat)
   nullify(wu%umat)

   call initialize_memory(cpp)

   call read_data_pw_u(wu, options%prefix)
   deallocate(wu%umat)

   if(.not.options%l_hf_energies) then
      if(wu%nums > wu%nums_occ) then
         offset=-(wu%ene(wu%nums_occ+1)+wu%ene(wu%nums_occ))/2.d0
      else
         offset=-wu%ene(wu%nums_occ)
      endif
   else
      write(stdout,*) 'HF energies to be implemented YET'
      stop
      !if(wu%nums > wu%nums_occ) then
      !   offset=-(ene_hf(wu%nums_occ+1)+ene_hf(wu%nums_occ))/2.d0
      !else
      !   offset=-ene_hf(wu%nums_occ)
      !endif
   endif



   call initialize_polaw(ww)
   write(stdout,*) 'addconduction_self_ontime_file1'!ATTENZIONE

!read coulombian potential and calculate inverse

   if(ss%whole_s) then
      write(stdout,*) 'Whole s not implemented YET'
      stop
   endif


   if(options%w_divergence == 2) then
      call read_data_pw_v(vp,options%prefix,options%debug,0,.true.)
   else
      call read_data_pw_v(vp,options%prefix,options%debug,0,.false.)
   endif
   call read_data_pw_ortho_polaw(op,options%prefix)
   call orthonormalize_vpot(op,vp)
   call invert_v_pot(vp,vpi)
   call free_memory(vp)
   write(stdout,*) 'addconduction_self_ontime1_45'
   call distribute_v_pot(vpi,vpid)
   call free_memory(vpi)
   call invert_ortho_polaw(op,opi)

   write(stdout,*) 'addconduction_self_ontime1_5 op',op%numpw!ATTENZIONE
   call distribute_ortho_polaw(op,opd)
   call free_memory(op)
   write(stdout,*) 'addconduction_self_ontime1_6 opd',opd%numpw!ATTENZIONE
   call distribute_ortho_polaw(opi,opid)
   call free_memory(opi)

   l_blk= (ss%n+1)/nproc
   if(l_blk*nproc < (ss%n+1)) l_blk = l_blk+1
   nbegin=mpime*l_blk+1 -(ss%n+1)
   nend=nbegin+l_blk-1
   if(nend > 0) nend = 0


   write(stdout,*) 'addconduction_self_ontime5',nbegin,l_blk!ATTENZIONE
   call flush_unit(stdout)

   allocate(sene(-ss%n:0,options%i_min:options%i_max))
   sene(:,:)=(0.d0,0.d0)



!loop on negative imaginary times
   do it=nbegin,nbegin+l_blk-1
      if(it <= 0)  then
         write(stdout,*) 'addconduction_self_ontime time', it!ATTENZIONE
!we take care of the symmetru t ==> -t

         call read_polaw(abs(it),ww,options%debug)
         write(stdout,*) 'addconduction_self_ontime6'!ATTENZIONE
         call flush_unit(stdout)
         call collect_ortho_polaw(opi,opid)
         write(stdout,*) 'addconduction_self_ontime6.1'!ATTENZIONE
         call orthonormalize_inverse(opi,ww)
         write(stdout,*) 'addconduction_self_ontime6.2'!ATTENZIONE
         call free_memory(opi)
         write(stdout,*) 'addconduction_self_ontime7'!ATTENZIONE
         call flush_unit(stdout)
         allocate(wtemp(ww%numpw,ww%numpw))

         call collect_v_pot(vpi,vpid)

         call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,vpi%vmat,ww%numpw,ww%pw,ww%numpw,&
                 &0.d0, wtemp,ww%numpw)



         call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,wtemp,ww%numpw,vpi%vmat,ww%numpw,&
              &0.d0,ww%pw,ww%numpw)

         call free_memory(vpi)

         deallocate(wtemp)

         call collect_ortho_polaw(op,opd)
         call orthonormalize_inverse(op,ww)
         call free_memory(op)

!!now ww contains \tilde{ww}
         write(stdout,*) 'addconduction_self_ontime8'!ATTENZIONE
         call flush_unit(stdout)

!read in cprim_prod
!first multiplication
!second multiplication
!copy on sene
            !loop on c' states
         do ii=max(options%i_min,wu%nums_occ+1),options%i_max
            cpp%cprim=ii
            call read_data_pw_cprim_prod(cpp, options%prefix,.false.,ok_read,.false.,.false.)

!loop on c
            allocate(vtemp(cpp%numpw))
            do jj=1,cpp%nums_cond

!multiply W_ijS_jc =T_ic
               call dgemv('N',ww%numpw,ww%numpw,1.d0,ww%pw,ww%numpw,cpp%cpmat(:,jj),1,0.d0,vtemp,1)
!multiply S_icTi_c
               sene(it,ii)=sene(it,ii)+ddot(cpp%numpw,vtemp,1,cpp%cpmat(:,jj),1)*&
                    & exp((wu%ene(jj+wu%nums_occ)+offset)*tf%times(it))*ww%factor*(0.d0,-1.d0)

            enddo
            sene(it,ii)=sene(it,ii)*(0.d0,1.d0)
            if(it==0) sene(it,ii)=sene(it,ii)*0.5d0
            write(stdout,*) 'Conduction contribution', it,ii, sene(it,ii)
            call flush_unit(stdout)
            deallocate(vtemp)
         enddo
      else
         call collect_ortho_polaw(opi,opid)
         call free_memory(opi)
         call collect_v_pot(vpi,vpid)
         call free_memory(vpi)
         call collect_ortho_polaw(op,opd)
         call free_memory(op)
         do ii=max(options%i_min,wu%nums_occ+1),options%i_max
            cpp%cprim=ii
            call read_data_pw_cprim_prod(cpp, options%prefix,.false.,ok_read,.false.,.false.)
         enddo
      endif
   enddo
   call mp_sum(sene(-ss%n:0,:))
   do ii=max(options%i_min,wu%nums_occ+1),options%i_max
       do it=-ss%n,0
          ss%diag(ii,it+ss%n+1)=ss%diag(ii, it+ss%n+1)+sene(it,ii)
       enddo
    enddo

 !copy sene results on ss with opportune factors
!!!!!!!!!!!
   call free_memory(vpid)
   call free_memory(opd)
   call free_memory(opi)
   call free_memory(opid)
   call free_memory_polaw(ww)
   call free_memory( cpp)
   deallocate(sene)
   return
 END SUBROUTINE addconduct_self_ontime_file



  SUBROUTINE selfenergy_ontime_file(ss, tf ,options)
!this subroutine calculates the self_energy of selected states
!using terms from file or from strategy BETA

    USE io_global,     ONLY : stdout, ionode
    USE input_gw,      ONLY : input_options
    USE basic_structures,  ONLY : v_pot,wannier_u,free_memory, ortho_polaw,initialize_memory,cprim_prod,q_mat,&
                                        & wannier_u_prim,v_pot_prim
    USE green_function,    ONLY : green, read_green, free_memory_green, initialize_green
    USE polarization,      ONLY : polaw, free_memory_polaw, read_polaw, invert_v_pot, invert_ortho_polaw,&
         & orthonormalize_inverse, initialize_polaw, orthonormalize_vpot, distribute_ortho_polaw, collect_ortho_polaw,&
         &  distribute_v_pot, collect_v_pot
    USE mp,                ONLY : mp_sum, mp_barrier
    USE para_gww,          ONLY : is_my_pola
    USE mp_global,            ONLY : nproc,mpime
    USE times_gw,  ONLY : times_freqs


   implicit none

   TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids
   TYPE(input_options) :: options
   TYPE(self_storage)  :: ss

   TYPE(v_pot) :: vp,vpi,vpid
   TYPE(ortho_polaw) :: op,opi, opd, opid
   TYPE(polaw) :: ww!dressed interaction
   TYPE(wannier_u) :: wu!structure to be read and initialized
   TYPE(cprim_prod) :: cpp,cppd!the producs c' c' v  wp
   TYPE(q_mat) :: qm, qmd!for strategy beta
   TYPE(wannier_u_prim) :: wup!for strategy beta
   TYPE(v_pot_prim) :: vpp,vppd!for strategy beta


   INTEGER iw,jw,kw,it,ii,jj
   INTEGER :: l_blk, nbegin,nend
   REAL(kind=DP) :: offset
   COMPLEX(kind=DP), ALLOCATABLE :: sene(:,:)
   REAL(kind=DP), ALLOCATABLE:: wtemp(:,:), vtemp(:,:)
   REAL(kind=DP), EXTERNAL :: ddot
   LOGICAL :: ok_read

   nullify(vp%vmat)
   nullify(vpi%vmat)
   nullify(op%on_mat)
   nullify(opi%on_mat)
   nullify(opd%on_mat)
   nullify(opid%on_mat)
   nullify(vpid%vmat)
   nullify(wu%umat)
   nullify(wup%umat)
   nullify(vpp%vmat)

   call initialize_memory(cpp)
   call initialize_memory(cppd)


   if(options%l_self_beta) ok_read=.true.
   call read_data_pw_u(wu, options%prefix)
   if(.not.options%l_self_beta) deallocate(wu%umat)

   if(.not.options%l_hf_energies) then
      if(wu%nums > wu%nums_occ) then
         offset=-(wu%ene(wu%nums_occ+1)+wu%ene(wu%nums_occ))/2.d0
      else
         offset=-wu%ene(wu%nums_occ)
      endif
   else
      write(stdout,*) 'HF energies to be implemented YET'
      stop
      !if(wu%nums > wu%nums_occ) then
      !   offset=-(ene_hf(wu%nums_occ+1)+ene_hf(wu%nums_occ))/2.d0
      !else
      !   offset=-ene_hf(wu%nums_occ)
      !endif
   endif



   call initialize_polaw(ww)
   write(stdout,*) 'addconduction_self_ontime1'!ATTENZIONE
   call flush_unit(stdout)

!read coulombian potential and calculate inverse

   if(ss%whole_s) then
      write(stdout,*) 'Whole s not implemented YET'
      stop
   endif


   if(options%w_divergence == 2) then
      call read_data_pw_v(vp,options%prefix,options%debug,0,.true.)
   else
      call read_data_pw_v(vp,options%prefix,options%debug,0,.false.)
   endif
   if(options%lnonorthogonal) then
      call read_data_pw_ortho_polaw(op,options%prefix)
      call orthonormalize_vpot(op,vp)
   endif
   call invert_v_pot(vp,vpi)
   call free_memory(vp)
   write(stdout,*) 'addconduction_self_ontime1_45'
   call flush_unit(stdout)
   call distribute_v_pot(vpi,vpid)
   call free_memory(vpi)
   if(options%lnonorthogonal) then
      call invert_ortho_polaw(op,opi)
   endif

   write(stdout,*) 'addconduction_self_ontime1_5 op',op%numpw!ATTENZIONE
   call flush_unit(stdout)
   if(options%lnonorthogonal) then
      call distribute_ortho_polaw(op,opd)
      call free_memory(op)

      write(stdout,*) 'addconduction_self_ontime1_6 opd',opd%numpw!ATTENZIONE
      call flush_unit(stdout)
      call distribute_ortho_polaw(opi,opid)
      call free_memory(opi)
   endif


   l_blk= (2*ss%n+1)/nproc
   if(l_blk*nproc < (2*ss%n+1)) l_blk = l_blk+1
   nbegin=mpime*l_blk+1 -(ss%n+1)
   nend=nbegin+l_blk-1
   if(nend > ss%n) nend = ss%n


   write(stdout,*) 'addconduction_self_ontime5',nbegin,l_blk!ATTENZIONE
   call flush_unit(stdout)
   allocate(sene(-ss%n:ss%n,options%i_min:options%i_max))
   sene(:,:)=(0.d0,0.d0)

!if required read and distribute q_mat
   if(options%l_self_beta) then
      call read_data_pw_q(qm,options%prefix,.true.)
      call distribute_qmat(qm,qmd)
      call free_memory(qm)
      if(options%i_max > wu%nums_occ) then
         if(options%w_divergence == 2) then
            call read_data_pw_v_pot_prim(vpp,options%prefix, .true.)
         else
            call read_data_pw_v_pot_prim(vpp,options%prefix, .false.)
         endif
         call distribute_v_pot_prim(vpp,vppd)
         call free_memory(vpp)
         call read_data_pw_u_prim(wup,options%prefix)
      endif
   endif



!loop on negative imaginary times
   do it=nbegin,nbegin+l_blk-1
      if(it <= ss%n)  then
         write(stdout,*) 'addconduction_self_ontime time', it!ATTENZIONE
         call flush_unit(stdout)
!we take care of the symmetru t ==> -t
         call read_polaw(abs(it),ww,options%debug)
         write(stdout,*) 'addconduction_self_ontime6'!ATTENZIONE
         call flush_unit(stdout)
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(opi,opid)
            write(stdout,*) 'addconduction_self_ontime6.1'!ATTENZIONE
            call orthonormalize_inverse(opi,ww)
            write(stdout,*) 'addconduction_self_ontime6.2'!ATTENZIONE
            call free_memory(opi)
         endif
         write(stdout,*) 'addconduction_self_ontime7'!ATTENZIONE
         call flush_unit(stdout)
         allocate(wtemp(ww%numpw,ww%numpw))
         call collect_v_pot(vpi,vpid)
         call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,vpi%vmat,ww%numpw,ww%pw,ww%numpw,&
                 &0.d0, wtemp,ww%numpw)
         call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,wtemp,ww%numpw,vpi%vmat,ww%numpw,&
              &0.d0,ww%pw,ww%numpw)
         call free_memory(vpi)
         deallocate(wtemp)
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(op,opd)
            call orthonormalize_inverse(op,ww)
            call free_memory(op)
         endif
!!now ww contains \tilde{ww}
         write(stdout,*) 'addconduction_self_ontime8'!ATTENZIONE
         call flush_unit(stdout)
!read in cprim_prod
!first multiplication
!second multiplication
!copy on sene
            !loop on c' states
         do ii=options%i_min,options%i_max
            cpp%cprim=ii
            call mp_barrier
            if(.not.options%l_self_beta) then
               call read_data_pw_cprim_prod(cpp, options%prefix,.true.,ok_read,.false.,.false.)
            else
!read qmat
               call  create_vcprim(cppd, ii ,wu, qmd)
               if(ii>wu%nums_occ) then
!if required adds the conduction term
                  call add_vcprim_conduction(cppd, wu, wup, vppd)
               end if
               call collect_cprim_prod(cpp,cppd)
               call free_memory(cppd)
            endif
            if(ok_read) then
!loop on c
               allocate(vtemp(cpp%numpw,max(cpp%nums_occ,cpp%nums-options%nc_minus)))
               !!allocate(vtemp(cpp%numpw,max(cpp%nums_occ,cpp%nums-cpp%nums_occ-options%nc_minus)))
               !!!!!!!!!! negative time
               if(it <= 0) then
                  call dgemm('N','N',ww%numpw,cpp%nums-cpp%nums_occ-options%nc_minus,ww%numpw,1.d0,ww%pw,ww%numpw,&
                       & cpp%cpmat(:,cpp%nums_occ+1:cpp%nums-options%nc_minus),cpp%lda,0.d0,vtemp,ww%numpw)
                  do jj=cpp%nums_occ+1,cpp%nums-options%nc_minus
!multiply W_ijS_jc =T_ic
!                     call dgemv('N',ww%numpw,ww%numpw,1.d0,ww%pw,ww%numpw,cpp%cpmat(:,jj),1,0.d0,vtemp,1)

!multiply S_icTi_c
                     sene(it,ii)=sene(it,ii)+ddot(cpp%numpw,vtemp(:,jj-cpp%nums_occ),1,cpp%cpmat(:,jj),1)*&
                          & exp((wu%ene(jj)+offset)*tf%times(it))*ww%factor*(0.d0,-1.d0)

                  enddo
                  sene(it,ii)=sene(it,ii)*(0.d0,1.d0)
                  write(stdout,*) 'Conduction contribution', it,ii, sene(it,ii)
                  call flush_unit(stdout)
               endif
               !!!!!!!!!! end negative time
               !!!!!!!!!! positive time
               if(it >= 0) then
                   call dgemm('N','N',ww%numpw,cpp%nums_occ,ww%numpw,1.d0,ww%pw,ww%numpw,&
                       &cpp%cpmat(:,1:cpp%nums_occ),cpp%lda,0.d0,vtemp,ww%numpw)

                  do jj=1,cpp%nums_occ

!multiply W_ijS_jc =T_ic
!                     call dgemv('N',ww%numpw,ww%numpw,1.d0,ww%pw,ww%numpw,cpp%cpmat(:,jj),1,0.d0,vtemp,1)
!multiply S_icTi_c
                     sene(it,ii)=sene(it,ii)+ddot(cpp%numpw,vtemp(:,jj),1,cpp%cpmat(:,jj),1)*&
                          & exp((wu%ene(jj)+offset)*tf%times(it))*ww%factor*(0.d0,+1.d0)

                  enddo
                  sene(it,ii)=sene(it,ii)*(0.d0,1.d0)
                  write(stdout,*) 'Conduction contribution', it,ii, sene(it,ii)
                  call flush_unit(stdout)
               endif
               !!!!!!! end positive time
               if(it==0) sene(it,ii)=sene(it,ii)*0.5d0
               deallocate(vtemp)
            endif
         enddo
      else
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(opi,opid)
            call free_memory(opi)
         endif
         call collect_v_pot(vpi,vpid)
         call free_memory(vpi)
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(op,opd)
            call free_memory(op)
         endif
         do ii=options%i_min,options%i_max
            cpp%cprim=ii
            call mp_barrier
            if(.not.options%l_self_beta) then
               call read_data_pw_cprim_prod(cpp, options%prefix,.true.,ok_read,.false.,.false.)
            else
!read qmat
               call  create_vcprim(cppd, ii ,wu, qmd)
               if(ii>wu%nums_occ) then
!if required adds the conduction term
                  call add_vcprim_conduction(cppd, wu, wup, vppd)
               end if
               call collect_cprim_prod(cpp,cppd)
               call free_memory(cppd)
            endif
         enddo
      endif
   enddo
   call mp_sum(sene(-ss%n:ss%n,:))
   do ii=options%i_min,options%i_max
       do it=-ss%n,ss%n
          ss%diag(ii,it+ss%n+1)=ss%diag(ii, it+ss%n+1)+sene(it,ii)
       enddo
    enddo

 !copy sene results on ss with opportune factors
!!!!!!!!!!!
    if(options%l_self_beta) call free_memory(qmd)
    if(options%l_self_beta .and. options%i_max > wu%nums_occ ) then
       call free_memory(vppd)
       call free_memory(wup)
    endif
   call free_memory(vpid)
   if(options%lnonorthogonal) then
      call free_memory(opd)
      call free_memory(opi)
      call free_memory(opid)
   endif
   call free_memory_polaw(ww)
   call free_memory( cpp)
   call free_memory(cppd)
   call free_memory( wu)
   deallocate(sene)
  return
 END SUBROUTINE selfenergy_ontime_file

!!!!!! split in two the calculation of the self: 1: negative times, and 2: positive times
 SUBROUTINE selfenergy_on_split_nt_file(ss, tf, options)
!this subroutine calculates the self_energy of selected states
!using terms from file or from strategy BETA
!!! FOR THE NEGATIVE TIME !!!

    USE io_global,         ONLY : stdout, ionode
    USE input_gw,          ONLY : input_options
    USE basic_structures,  ONLY : v_pot,wannier_u,free_memory, ortho_polaw,&
                                  & initialize_memory, cprim_prod,q_mat,&
                                  & wannier_u_prim,v_pot_prim
    USE green_function,    ONLY : green, read_green, free_memory_green, initialize_green
    USE polarization,      ONLY : polaw, free_memory_polaw, read_polaw, &
                                  & invert_v_pot, invert_ortho_polaw, &
                                  & orthonormalize_inverse,&
                                  & initialize_polaw, orthonormalize_vpot,&
                                  & distribute_ortho_polaw, collect_ortho_polaw,&
                                  & distribute_v_pot, collect_v_pot
    USE mp,                ONLY : mp_sum, mp_barrier
    USE para_gww,          ONLY : is_my_pola
    USE mp_global,         ONLY : nproc,mpime
    USE times_gw,          ONLY : times_freqs
    USE io_files,          ONLY : find_free_unit, prefix

   implicit none

   TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids
   TYPE(input_options) :: options
   TYPE(self_storage)  :: ss

   TYPE(v_pot) :: vp,vpi,vpid
   TYPE(ortho_polaw) :: op,opi, opd, opid
   TYPE(polaw) :: ww!dressed interaction
   TYPE(wannier_u) :: wu!structure to be read and initialized
   TYPE(cprim_prod) :: cpp,cppd!the producs c' c' v  wp
   TYPE(q_mat) :: qm, qmd!for strategy beta
   TYPE(wannier_u_prim) :: wup!for strategy beta
   TYPE(v_pot_prim) :: vpp,vppd!for strategy beta


   INTEGER iw,jw,kw,it,ii,jj, aa, bb, iun
   INTEGER :: l_blk, nbegin,nend
   REAL(kind=DP) :: offset
   COMPLEX(kind=DP), ALLOCATABLE :: sene(:,:)
   REAL(kind=DP), ALLOCATABLE:: wtemp(:,:), vtemp(:,:)
   REAL(kind=DP), EXTERNAL :: ddot
   LOGICAL :: ok_read

   nullify(vp%vmat)
   nullify(vpi%vmat)
   nullify(op%on_mat)
   nullify(opi%on_mat)
   nullify(opd%on_mat)
   nullify(opid%on_mat)
   nullify(vpid%vmat)
   nullify(wu%umat)
   nullify(wup%umat)
   nullify(vpp%vmat)

   call initialize_memory(cpp)
   call initialize_memory(cppd)

   if(options%l_self_beta) ok_read=.true.
   call read_data_pw_u(wu, options%prefix)
   if(.not.options%l_self_beta) deallocate(wu%umat)

   if(.not.options%l_hf_energies) then
      if(wu%nums > wu%nums_occ) then
         offset=-(wu%ene(wu%nums_occ+1)+wu%ene(wu%nums_occ))/2.d0
      else
         offset=-wu%ene(wu%nums_occ)
      endif
   else
      write(stdout,*) 'HF energies to be implemented YET'
      stop
   endif

   call initialize_polaw(ww)
   write(stdout,*) 'selfenergy_on_split_time_file_1'!ATTENZIONE
   call flush_unit(stdout)
!read coulombian potential and calculate inverse
   if(ss%whole_s) then
      write(stdout,*) 'Whole s not implemented YET'
      stop
   endif

   if(options%w_divergence == 2) then
      call read_data_pw_v(vp,options%prefix,options%debug,0,.true.)
   else
      call read_data_pw_v(vp,options%prefix,options%debug,0,.false.)
   endif
   if(options%lnonorthogonal) then
      call read_data_pw_ortho_polaw(op,options%prefix)
      call orthonormalize_vpot(op,vp)
   endif
   call invert_v_pot(vp,vpi)
   call free_memory(vp)
   write(stdout,*) 'selfenergy_on_split_time_file_2'
   call flush_unit(stdout)
   call distribute_v_pot(vpi,vpid)
   call free_memory(vpi)
   if(options%lnonorthogonal) then
      call invert_ortho_polaw(op,opi)
   endif

   write(stdout,*) 'selfenergy_on_split_time_file_3 with op%numpw=',op%numpw!ATTENZIONE
   call flush_unit(stdout)
   if(options%lnonorthogonal) then
      call distribute_ortho_polaw(op,opd)
      call free_memory(op)

      write(stdout,*) 'selfenergy_on_split_time_file_4 with opd%numpw=',opd%numpw!ATTENZIONE
      call flush_unit(stdout)
      call distribute_ortho_polaw(opi,opid)
      call free_memory(opi)
   endif

   !!! here, we first cover the negative time on spread them
   !! onto the nproc processors
   l_blk= (ss%n+1)/nproc
   if(l_blk*nproc < (ss%n+1)) l_blk = l_blk+1
   nbegin=mpime*l_blk+1 -(ss%n+1)
   nend=nbegin+l_blk-1
   if(nend > 0) nend = 0
   !!!!!!!!

   !!!!!!!! CHECK
   write(stdout,*) 'selfenergy_on_split_time_file_5 with nbegin=',nbegin, ' and l_blk=', l_blk
   call flush_unit(stdout)
   !!!!!!!!

   !!!!!!!! allocate sene !!!!!!
   allocate(sene(-ss%n:0,options%i_min:options%i_max))
   sene(:,:)=(0.d0,0.d0)

   !if required read and distribute q_mat
   if(options%l_self_beta) then
      write(stdout,*) 'selfenergy_on_split_time_file_6 with options%l_self_beta=', options%l_self_beta
      call flush_unit(stdout)
      call read_data_pw_q(qm,options%prefix,.true.)
      call distribute_qmat(qm,qmd)
      call free_memory(qm)
      if(options%i_max > wu%nums_occ) then
         if(options%w_divergence == 2) then
            call read_data_pw_v_pot_prim(vpp,options%prefix, .true.)
         else
            call read_data_pw_v_pot_prim(vpp,options%prefix, .false.)
         endif
         call distribute_v_pot_prim(vpp,vppd)
         call free_memory(vpp)
         call read_data_pw_u_prim(wup,options%prefix)
      endif
   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !loop on negative imaginary times
   do it=nbegin,nbegin+l_blk-1
      if(it <= 0)  then
         write(stdout,*) 'Negative imaginary time with it=', it
         call flush_unit(stdout)
         !we take care of the symmetru t ==> -t
         call read_polaw(abs(it),ww,options%debug)
         write(stdout,*) 'after read_polaw'
         call flush_unit(stdout)
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(opi,opid)
            write(stdout,*) 'after collect_ortho_polaw'
            call orthonormalize_inverse(opi,ww)
            write(stdout,*) 'after orthonormalize_inverse'
            call free_memory(opi)
         endif
         allocate(wtemp(ww%numpw,ww%numpw))
         call collect_v_pot(vpi,vpid)
         call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,vpi%vmat,ww%numpw,ww%pw,ww%numpw,&
                 &0.d0, wtemp,ww%numpw)
         call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,wtemp,ww%numpw,vpi%vmat,ww%numpw,&
              &0.d0,ww%pw,ww%numpw)
         call free_memory(vpi)
         deallocate(wtemp)
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(op,opd)
            call orthonormalize_inverse(op,ww)
            call free_memory(op)
         endif
         write(stdout,*) 'CHECK 1'
         call flush_unit(stdout)
         !!now ww contains \tilde{ww}
         !read in cprim_prod
         !first multiplication
         !second multiplication
         !copy on sene
         !!!!!!!!! !loop on c' states !!!!!!!!!!!
         do ii=options%i_min,options%i_max
            cpp%cprim=ii
            call mp_barrier
            if(.not.options%l_self_beta) then
               call read_data_pw_cprim_prod(cpp, options%prefix,.true.,ok_read,.false.,.false.)
            else
            !read qmat
               call  create_vcprim(cppd, ii ,wu, qmd)
               if(ii>wu%nums_occ) then
                  !if required adds the conduction term
                  call add_vcprim_conduction(cppd, wu, wup, vppd)
               end if
               call collect_cprim_prod(cpp,cppd)
               call free_memory(cppd)
            endif !! if(.not.options%l_self_beta)
            !!!call mp_barrier
            write(stdout,*) 'CHECK 2 and ok_read=', ok_read
            call flush_unit(stdout)
            if(ok_read) then
               !loop on c
               write(stdout,*) 'CHECK 2.1'
               call flush_unit(stdout)
               !allocate(vtemp(cpp%numpw,max(cpp%nums_occ,cpp%nums-cpp%nums_occ-options%nc_minus)))
               !!allocate(vtemp(cpp%numpw,cpp%nums_occ+1:max(cpp%nums_occ,cpp%nums-options%nc_minus)))
               allocate(vtemp(cpp%numpw,max(cpp%nums_occ,cpp%nums-options%nc_minus)))
               write(stdout,*) 'CHECK 2.1 cpp%numpw=', cpp%numpw, ' ww%numpw=', ww%numpw
               write(stdout,*) 'CHECK 2.1 cpp%nums=', cpp%nums
               write(stdout,*) 'CHECK 2.2: after allocate vtemp'
               call flush_unit(stdout)
               call dgemm('N','N',ww%numpw,cpp%nums-cpp%nums_occ-options%nc_minus,ww%numpw,1.d0,ww%pw,ww%numpw,&
                   & cpp%cpmat(:,cpp%nums_occ+1:cpp%nums-options%nc_minus),cpp%lda,0.d0,vtemp,ww%numpw)
               !!call dgemm('N','N',ww%numpw,cpp%nums-cpp%nums_occ-options%nc_minus,ww%numpw,1.d0,ww%pw,ww%numpw,&
               !!     & cpp%cpmat(1:ww%numpw,cpp%nums_occ+1:cpp%nums-options%nc_minus),cpp%lda,0.d0,vtemp,ww%numpw)
               write(stdout,*) 'CHECK 3; after dgemm'
               call flush_unit(stdout)
               do jj=cpp%nums_occ+1,cpp%nums-options%nc_minus
                      !multiply W_ijS_jc =T_ic
                      !multiply S_icTi_c
                      sene(it,ii)=sene(it,ii)+ddot(cpp%numpw,vtemp(:,jj-cpp%nums_occ),&
                 & 1,cpp%cpmat(:,jj),1)*exp((wu%ene(jj)+offset)*tf%times(it))*&
                 & ww%factor*(0.d0,-1.d0)

               enddo
               sene(it,ii)=sene(it,ii)*(0.d0,1.d0)
               write(stdout,*) 'sene(', it,',', ii, ') = ', sene(it,ii)
               call flush_unit(stdout)
               if(it == 0) then !! i add also the contribution from the VB
                  !!call dgemm('N','N',ww%numpw,cpp%nums_occ,ww%numpw,1.d0,ww%pw,ww%numpw,&
                  !!     &cpp%cpmat(:,1:cpp%nums_occ),cpp%lda,0.d0,vtemp,ww%numpw)
                  call dgemm('N','N',ww%numpw,cpp%nums_occ,ww%numpw,1.d0,ww%pw,ww%numpw,&
                       &cpp%cpmat(1:ww%numpw,1:cpp%nums_occ),cpp%lda,0.d0,vtemp,ww%numpw)
                  write(stdout,*) 'CHECK 3; after dgemm at it=0'
                  call flush_unit(stdout)
                  do jj=1,cpp%nums_occ
                     sene(it,ii)=sene(it,ii)+ddot(cpp%numpw,vtemp(:,jj),1,cpp%cpmat(:,jj),1)*&
                          & exp((wu%ene(jj)+offset)*tf%times(it))*ww%factor*(0.d0,+1.d0)

                  enddo
                  sene(it,ii)=sene(it,ii)*(0.d0,1.d0)
                  sene(it,ii)=sene(it,ii)*0.5d0
                  write(stdout,*) 'sene(', it,',', ii, ') = ', sene(it,ii)
                  call flush_unit(stdout)
               endif !! (it == 0)

               write(stdout,*) 'Conduction contribution', it,ii, sene(it,ii)
               call flush_unit(stdout)

               deallocate(vtemp)
            endif !! ok_read
         enddo !! loop on c' : ii=options%i_min,options%i_max
      else !! if(it <= 0)
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(opi,opid)
            call free_memory(opi)
         endif
         call collect_v_pot(vpi,vpid)
         call free_memory(vpi)
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(op,opd)
            call free_memory(op)
         endif
         do ii=options%i_min,options%i_max
            cpp%cprim=ii
            call mp_barrier
            if(.not.options%l_self_beta) then
               call read_data_pw_cprim_prod(cpp, options%prefix,.true.,ok_read,.false.,.false.)
            else
               !read qmat
               call  create_vcprim(cppd, ii ,wu, qmd)
               if(ii>wu%nums_occ) then
                  !if required adds the conduction term
                  call add_vcprim_conduction(cppd, wu, wup, vppd)
               end if
               call collect_cprim_prod(cpp,cppd)
               call free_memory(cppd)
            endif
            !!!call mp_barrier
         enddo
      endif !! if(it <= 0)
   enddo !!loop on negative imaginary times: do it=nbegin,nbegin+l_blk-1
   call mp_sum(sene(-ss%n:0,:))

   do ii=options%i_min,options%i_max
       do it=-ss%n,0
          ss%diag(ii,it+ss%n+1)=ss%diag(ii, it+ss%n+1)+sene(it,ii)
       enddo
   enddo

   if(ionode) then
       iun = find_free_unit()
       if(.not. options%debug) then
          open( unit=iun, file=trim(options%prefix)//'.self_energie_negative_time', status='unknown',form='unformatted')
       else
          open( unit=iun, file=trim(options%prefix)//'.self_energie_negative_time', status='unknown',form='formatted')
       endif
       do ii=options%i_min,options%i_max
          do it=-ss%n,0
             write(iun) ss%diag(ii,it+ss%n+1)
          enddo
       enddo
       close(iun)
   endif

   write(stdout,*) 'File written'
   call flush_unit(stdout)

 !copy sene results on ss with opportune factors
!!!!!!!!!!!
   if(options%l_self_beta) call free_memory(qmd)
   if(options%l_self_beta .and. options%i_max > wu%nums_occ ) then
       call free_memory(vppd)
       call free_memory(wup)
   endif
   call free_memory(vpid)
   if(options%lnonorthogonal) then
      call free_memory(opd)
      call free_memory(opi)
      call free_memory(opid)
   endif
   call free_memory_polaw(ww)
   call free_memory( cpp)
   call free_memory(cppd)
   call free_memory( wu)
   deallocate(sene)
  return

 END SUBROUTINE selfenergy_on_split_nt_file

 SUBROUTINE selfenergy_on_split_pt_file(ss, tf, options)
!this subroutine calculates the self_energy of selected states
!using terms from file or from strategy BETA
!!! FOR THE POSITIVE TIME !!!

    USE io_global,         ONLY : stdout, ionode
    USE input_gw,          ONLY : input_options
    USE basic_structures,  ONLY : v_pot,wannier_u,free_memory, ortho_polaw,&
                                  & initialize_memory, cprim_prod,q_mat,&
                                  & wannier_u_prim,v_pot_prim
    USE green_function,    ONLY : green, read_green, free_memory_green, initialize_green
    USE polarization,      ONLY : polaw, free_memory_polaw, read_polaw, invert_v_pot,&
                                  & invert_ortho_polaw, orthonormalize_inverse,&
                                  & initialize_polaw, orthonormalize_vpot,&
                                  & distribute_ortho_polaw, collect_ortho_polaw,&
                                  & distribute_v_pot, collect_v_pot
    USE mp,                ONLY : mp_sum, mp_barrier
    USE para_gww,          ONLY : is_my_pola
    USE mp_global,         ONLY : nproc,mpime
    USE times_gw,          ONLY : times_freqs
    USE io_files,          ONLY : find_free_unit, prefix

   implicit none

   TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids
   TYPE(input_options) :: options
   TYPE(self_storage)  :: ss

   TYPE(v_pot) :: vp,vpi,vpid
   TYPE(ortho_polaw) :: op,opi, opd, opid
   TYPE(polaw) :: ww!dressed interaction
   TYPE(wannier_u) :: wu!structure to be read and initialized
   TYPE(cprim_prod) :: cpp,cppd!the producs c' c' v  wp
   TYPE(q_mat) :: qm, qmd!for strategy beta
   TYPE(wannier_u_prim) :: wup!for strategy beta
   TYPE(v_pot_prim) :: vpp,vppd!for strategy beta


   INTEGER iw,jw,kw,it,ii,jj
   INTEGER :: l_blk, nbegin,nend
   REAL(kind=DP) :: offset
   COMPLEX(kind=DP), ALLOCATABLE :: sene(:,:)
   REAL(kind=DP), ALLOCATABLE:: wtemp(:,:), vtemp(:,:)
   REAL(kind=DP), EXTERNAL :: ddot
   LOGICAL :: ok_read
   INTEGER :: iunet

   nullify(vp%vmat)
   nullify(vpi%vmat)
   nullify(op%on_mat)
   nullify(opi%on_mat)
   nullify(opd%on_mat)
   nullify(opid%on_mat)
   nullify(vpid%vmat)
   nullify(wu%umat)
   nullify(wup%umat)
   nullify(vpp%vmat)

   call initialize_memory(cpp)
   call initialize_memory(cppd)


   !!iunet=find_free_unit()
   !!open( unit=iunet, file=trim(options%prefix)//'.self_energie_negatve_time', status='old',form='unformatted')
   !!do ii=options%i_min,options%i_max
   !!       do it=-ss%n,0
   !!          read (iunet) ss%diag(ii,it+ss%n+1)
   !!       enddo
   !!enddo
   !!close(iunet)

   if(options%l_self_beta) ok_read=.true.
   call read_data_pw_u(wu, options%prefix)
   if(.not.options%l_self_beta) deallocate(wu%umat)

   if(.not.options%l_hf_energies) then
      if(wu%nums > wu%nums_occ) then
         offset=-(wu%ene(wu%nums_occ+1)+wu%ene(wu%nums_occ))/2.d0
      else
         offset=-wu%ene(wu%nums_occ)
      endif
   else
      write(stdout,*) 'HF energies to be implemented YET'
      stop
   endif

   call initialize_polaw(ww)
   write(stdout,*) 'selfenergy_on_split_time_file_1'!ATTENZIONE
   call flush_unit(stdout)
!read coulombian potential and calculate inverse
   if(ss%whole_s) then
      write(stdout,*) 'Whole s not implemented YET'
      stop
   endif

   if(options%w_divergence == 2) then
      call read_data_pw_v(vp,options%prefix,options%debug,0,.true.)
   else
      call read_data_pw_v(vp,options%prefix,options%debug,0,.false.)
   endif
   if(options%lnonorthogonal) then
      call read_data_pw_ortho_polaw(op,options%prefix)
      call orthonormalize_vpot(op,vp)
   endif
   call invert_v_pot(vp,vpi)
   call free_memory(vp)
   write(stdout,*) 'selfenergy_on_split_time_file_2'
   call flush_unit(stdout)
   call distribute_v_pot(vpi,vpid)
   call free_memory(vpi)
   if(options%lnonorthogonal) then
      call invert_ortho_polaw(op,opi)
   endif

   write(stdout,*) 'selfenergy_on_split_time_file_3 with op%numpw=',op%numpw!ATTENZIONE
   call flush_unit(stdout)
   if(options%lnonorthogonal) then
      call distribute_ortho_polaw(op,opd)
      call free_memory(op)

      write(stdout,*) 'selfenergy_on_split_time_file_4 with opd%numpw=',opd%numpw!ATTENZIONE
      call flush_unit(stdout)
      call distribute_ortho_polaw(opi,opid)
      call free_memory(opi)
   endif

   !!! here, we cover now the positive time and spread them
   !! onto the nproc processors
   l_blk= (ss%n+1)/nproc
   if(l_blk*nproc < (ss%n+1)) l_blk = l_blk+1
   nbegin=mpime*l_blk+1
   nend=nbegin+l_blk-1
   if(nend > ss%n) nend = ss%n
   !!!!!!!!

   !!!!!!!! CHECK
   write(stdout,*) 'selfenergy_on_split_time_file_5 with nbegin=',nbegin, ' and l_blk=', l_blk
   call flush_unit(stdout)
   !!!!!!!!

   !!!!!!!! allocate sene !!!!!!
   allocate(sene(1:ss%n,options%i_min:options%i_max))
   sene(:,:)=(0.d0,0.d0)

   !if required read and distribute q_mat
   if(options%l_self_beta) then
      write(stdout,*) 'selfenergy_on_split_time_file_6 with options%l_self_beta=', options%l_self_beta
      call flush_unit(stdout)
      call read_data_pw_q(qm,options%prefix,.true.)
      call distribute_qmat(qm,qmd)
      call free_memory(qm)
      if(options%i_max > wu%nums_occ) then
         if(options%w_divergence == 2) then
            call read_data_pw_v_pot_prim(vpp,options%prefix, .true.)
         else
            call read_data_pw_v_pot_prim(vpp,options%prefix, .false.)
         endif
         call distribute_v_pot_prim(vpp,vppd)
         call free_memory(vpp)
         call read_data_pw_u_prim(wup,options%prefix)
      endif
   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !loop on negative imaginary times
   do it=nbegin,nbegin+l_blk-1
      if(it <= ss%n)  then
         write(stdout,*) 'Positive imaginary time with it=', it
         call flush_unit(stdout)
         !we take care of the symmetru t ==> -t
         call read_polaw(abs(it),ww,options%debug)
         write(stdout,*) 'after read_polaw'
         call flush_unit(stdout)
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(opi,opid)
            write(stdout,*) 'after collect_ortho_polaw'
            call orthonormalize_inverse(opi,ww)
            write(stdout,*) 'after orthonormalize_inverse'
            call free_memory(opi)
         endif
         write(stdout,*) 'before allocate wtemp'
         call flush_unit(stdout)
         allocate(wtemp(ww%numpw,ww%numpw))
         write(stdout,*) 'before collect_v_pot'
         call flush_unit(stdout)
         call collect_v_pot(vpi,vpid)
         write(stdout,*) 'before dgemm1'
         call flush_unit(stdout)
         call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,vpi%vmat,ww%numpw,ww%pw,ww%numpw,&
                 &0.d0, wtemp,ww%numpw)
         write(stdout,*) 'before dgemm2'
         call flush_unit(stdout)
         call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,wtemp,ww%numpw,vpi%vmat,ww%numpw,&
              &0.d0,ww%pw,ww%numpw)
         write(stdout,*) 'before free_memory'
         call flush_unit(stdout)
         call free_memory(vpi)
         write(stdout,*) 'deallocate'
         call flush_unit(stdout)
         deallocate(wtemp)
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(op,opd)
            call orthonormalize_inverse(op,ww)
            call free_memory(op)
         endif
         write(stdout,*) 'CHECK 1'
         call flush_unit(stdout)
         !!now ww contains \tilde{ww}
         !read in cprim_prod
         !first multiplication
         !second multiplication
         !copy on sene
         !!!!!!!!! !loop on c' states !!!!!!!!!!!
         do ii=options%i_min,options%i_max
            cpp%cprim=ii
            call mp_barrier
            write(stdout,*) 'START ii=', ii
            call flush_unit(stdout)
            if(.not.options%l_self_beta) then
               call read_data_pw_cprim_prod(cpp, options%prefix,.true.,ok_read,.false.,.false.)
            else
            !read qmat
               call  create_vcprim(cppd, ii ,wu, qmd)
               if(ii>wu%nums_occ) then
                  !if required adds the conduction term
                  call add_vcprim_conduction(cppd, wu, wup, vppd)
               end if
               call collect_cprim_prod(cpp,cppd)
               call free_memory(cppd)
            endif !! if(.not.options%l_self_beta)
            write(stdout,*) 'ok_read=', ok_read
            call flush_unit(stdout)
            !!call mp_barrier
            if(ok_read) then
       !loop on c
               allocate(vtemp(cpp%numpw,max(cpp%nums_occ,cpp%nums-cpp%nums_occ)))
               call dgemm('N','N',ww%numpw,cpp%nums_occ,ww%numpw,1.d0,ww%pw,ww%numpw,&
                    &cpp%cpmat(:,1:cpp%nums_occ),cpp%lda,0.d0,vtemp,ww%numpw)

               do jj=1,cpp%nums_occ
                      !multiply W_ijS_jc =T_ic
                      !multiply S_icTi_c
                      sene(it,ii)=sene(it,ii)+ddot(cpp%numpw,vtemp(:,jj),1,cpp%cpmat(:,jj),1)*&
                                 & exp((wu%ene(jj)+offset)*tf%times(it))*ww%factor*(0.d0,+1.d0)

               enddo
               sene(it,ii)=sene(it,ii)*(0.d0,1.d0)
               write(stdout,*) 'sene(', it,',', ii, ') = ', sene(it,ii)
               call flush_unit(stdout)
               deallocate(vtemp)
            endif !! ok_read
         enddo !! loop on c' : ii=options%i_min,options%i_max
      else !! if(it <= ss%n)
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(opi,opid)
            call free_memory(opi)
         endif
         call collect_v_pot(vpi,vpid)
         call free_memory(vpi)
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(op,opd)
            call free_memory(op)
         endif
         do ii=options%i_min,options%i_max
            cpp%cprim=ii
            call mp_barrier
            if(.not.options%l_self_beta) then
               call read_data_pw_cprim_prod(cpp, options%prefix,.true.,ok_read,.false.,.false.)
            else
               !read qmat
               call  create_vcprim(cppd, ii ,wu, qmd)
               if(ii>wu%nums_occ) then
                  !if required adds the conduction term
                  call add_vcprim_conduction(cppd, wu, wup, vppd)
               end if
               call collect_cprim_prod(cpp,cppd)
               call free_memory(cppd)
            endif
            !!call mp_barrier
         enddo
      endif !! if(it <= 0)
   enddo !!loop on negative imaginary times: do it=nbegin,nbegin+l_blk-1
   call mp_sum(sene(1:ss%n,:))

   do ii=options%i_min,options%i_max
       do it=1,ss%n
          ss%diag(ii,it+ss%n+1)=ss%diag(ii, it+ss%n+1)+sene(it,ii)
       enddo
   enddo
 !copy sene results on ss with opportune factors
!!!!!!!!!!!
   if(options%l_self_beta) call free_memory(qmd)
   if(options%l_self_beta .and. options%i_max > wu%nums_occ ) then
       call free_memory(vppd)
       call free_memory(wup)
   endif
   call free_memory(vpid)
   if(options%lnonorthogonal) then
      call free_memory(opd)
      call free_memory(opi)
      call free_memory(opid)
   endif
   call free_memory_polaw(ww)
   call free_memory( cpp)
   call free_memory(cppd)
   call free_memory( wu)

   if (ionode) then
      iunet=find_free_unit()
      open( unit=iunet, file=trim(options%prefix)//'.self_energie_negative_time', status='old',form='unformatted')
      close(iunet, status='delete')
   endif

   deallocate(sene)
  return

 END SUBROUTINE selfenergy_on_split_pt_file

!!!!!!!!!!!!
  SUBROUTINE selfenergy_ontime_upper(ss, tf ,options)
!this subroutine calculates the self_energy of selected states


    USE io_global,     ONLY : stdout, ionode
    USE input_gw,      ONLY : input_options
    USE basic_structures,  ONLY : v_pot,wannier_u,free_memory, ortho_polaw,initialize_memory,cprim_prod,q_mat,&
                                        & wannier_u_prim,v_pot_prim, upper_states
    USE green_function,    ONLY : green, read_green, free_memory_green, initialize_green
    USE polarization,      ONLY : polaw, free_memory_polaw, read_polaw, invert_v_pot, invert_ortho_polaw,&
         & orthonormalize_inverse, initialize_polaw, orthonormalize_vpot, distribute_ortho_polaw, collect_ortho_polaw,&
         &  distribute_v_pot, collect_v_pot
    USE mp,                ONLY : mp_sum
    USE para_gww,          ONLY : is_my_pola
    USE mp_global,            ONLY : nproc,mpime
    USE times_gw,  ONLY : times_freqs


   implicit none

   TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids
   TYPE(input_options) :: options
   TYPE(self_storage)  :: ss

   TYPE(v_pot) :: vp,vpi,vpid
   TYPE(ortho_polaw) :: op,opi, opd, opid
   TYPE(polaw) :: ww!dressed interaction
   TYPE(wannier_u) :: wu!structure to be read and initialized
   TYPE(cprim_prod) :: cpp,cppd!the producs c' c' v  wp
   TYPE(upper_states) :: us


   INTEGER iw,jw,kw,it,ii,jj
   INTEGER :: l_blk, nbegin,nend
   REAL(kind=DP) :: offset
   COMPLEX(kind=DP), ALLOCATABLE :: sene(:,:)
   REAL(kind=DP), ALLOCATABLE:: wtemp(:,:), vtemp(:,:)
   REAL(kind=DP), EXTERNAL :: ddot
   LOGICAL :: ok_read

   nullify(vp%vmat)
   nullify(vpi%vmat)
   nullify(op%on_mat)
   nullify(opi%on_mat)
   nullify(opd%on_mat)
   nullify(opid%on_mat)
   nullify(vpid%vmat)
   nullify(wu%umat)



   call initialize_memory(cpp)
   call initialize_memory(cppd)


   call read_data_pw_u(wu, options%prefix)
   deallocate(wu%umat)
   call initialize_memory(us)
   call read_data_pw_upper_states(us,options%prefix)


   if(.not.options%l_hf_energies) then
      if(wu%nums > wu%nums_occ) then
         offset=-(wu%ene(wu%nums_occ+1)+wu%ene(wu%nums_occ))/2.d0
      else
         offset=-wu%ene(wu%nums_occ)
      endif
   else
      write(stdout,*) 'HF energies to be implemented YET'
      stop
      !if(wu%nums > wu%nums_occ) then
      !   offset=-(ene_hf(wu%nums_occ+1)+ene_hf(wu%nums_occ))/2.d0
      !else
      !   offset=-ene_hf(wu%nums_occ)
      !endif
   endif



   call initialize_polaw(ww)
   write(stdout,*) 'addconduction_self_upper1'!ATTENZIONE
   call flush_unit(stdout)

!read coulombian potential and calculate inverse

   if(ss%whole_s) then
      write(stdout,*) 'Whole s not implemented YET'
      stop
   endif


   if(options%w_divergence == 2) then
      call read_data_pw_v(vp,options%prefix,options%debug,0,.true.)
   else
      call read_data_pw_v(vp,options%prefix,options%debug,0,.false.)
   endif
   if(options%lnonorthogonal) then
      call read_data_pw_ortho_polaw(op,options%prefix)
      call orthonormalize_vpot(op,vp)
   endif
   call invert_v_pot(vp,vpi)
   call free_memory(vp)
   write(stdout,*) 'addconduction_self_upper1_45'
   call flush_unit(stdout)
   call distribute_v_pot(vpi,vpid)
   call free_memory(vpi)
   if(options%lnonorthogonal) then
      call invert_ortho_polaw(op,opi)
   endif

   write(stdout,*) 'addconduction_self_upper1_5 op',op%numpw!ATTENZIONE
   call flush_unit(stdout)
   if(options%lnonorthogonal) then
      call distribute_ortho_polaw(op,opd)
      call free_memory(op)

      write(stdout,*) 'addconduction_self_upper_6 opd',opd%numpw!ATTENZIONE
      call flush_unit(stdout)
      call distribute_ortho_polaw(opi,opid)
      call free_memory(opi)
   endif

   l_blk= (ss%n+1)/nproc
   if(l_blk*nproc < (ss%n+1)) l_blk = l_blk+1
   nbegin=mpime*l_blk + 1 - (ss%n+1)
   nend=nbegin+l_blk-1
   if(nend > 0) nend = 0


   write(stdout,*) 'addconduction_self_upper5',nbegin,l_blk!ATTENZIONE
   call flush_unit(stdout)
   allocate(sene(-ss%n:ss%n,options%i_min:options%i_max))
   sene(:,:)=(0.d0,0.d0)




!loop on negative imaginary times
   do it=nbegin,nbegin+l_blk-1
      if(it <= ss%n)  then
         write(stdout,*) 'addconduction_self_ontime time', it!ATTENZIONE
         call flush_unit(stdout)
!we take care of the symmetru t ==> -t

         call read_polaw(abs(it),ww,options%debug)
         write(stdout,*) 'addconduction_self_upper6'!ATTENZIONE
         call flush_unit(stdout)
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(opi,opid)
            write(stdout,*) 'addconduction_self_ontime6.1'!ATTENZIONE
            call orthonormalize_inverse(opi,ww)
            write(stdout,*) 'addconduction_self_ontime6.2'!ATTENZIONE
            call free_memory(opi)
         endif
         write(stdout,*) 'addconduction_self_upper7'!ATTENZIONE
         call flush_unit(stdout)
         allocate(wtemp(ww%numpw,ww%numpw))

         call collect_v_pot(vpi,vpid)

         call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,vpi%vmat,ww%numpw,ww%pw,ww%numpw,&
                 &0.d0, wtemp,ww%numpw)



         call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,wtemp,ww%numpw,vpi%vmat,ww%numpw,&
              &0.d0,ww%pw,ww%numpw)

         call free_memory(vpi)

         deallocate(wtemp)

         if(options%lnonorthogonal) then
            call collect_ortho_polaw(op,opd)
            call orthonormalize_inverse(op,ww)
            call free_memory(op)
         endif

!!now ww contains \tilde{ww}
         write(stdout,*) 'addconduction_self_ontime8'!ATTENZIONE
         call flush_unit(stdout)

!read in cprim_prod
!first multiplication
!second multiplication
!copy on sene
            !loop on c' states
         do ii=options%i_min,options%i_max
            cpp%cprim=ii
            call read_data_pw_cprim_prod(cpp, options%prefix,.true.,ok_read,.false.,.true.)

            if(ok_read) then

!loop on c
               allocate(vtemp(cpp%numpw,us%nums_reduced))
               if(it <= 0) then
                  call dgemm('N','N',ww%numpw,us%nums_reduced,ww%numpw,1.d0,ww%pw,ww%numpw,&
                       & cpp%cpmat,cpp%lda,0.d0,vtemp,ww%numpw)

                  do jj=1,us%nums_reduced

!multiply W_ijS_jc =T_ic
!                     call dgemv('N',ww%numpw,ww%numpw,1.d0,ww%pw,ww%numpw,cpp%cpmat(:,jj),1,0.d0,vtemp,1)

!multiply S_icTi_c
                     sene(it,ii)=sene(it,ii)+ddot(cpp%numpw,vtemp(:,jj),1,cpp%cpmat(:,jj),1)*&
                          & exp((us%ene(jj)+offset)*tf%times(it))*ww%factor*(0.d0,-1.d0)

                  enddo
                  sene(it,ii)=sene(it,ii)*(0.d0,1.d0)
                  write(stdout,*) 'Conduction contribution', it,ii, sene(it,ii)
                  call flush_unit(stdout)
               endif
               if(it >= 0) then
               endif
               if(it==0) sene(it,ii)=sene(it,ii)*0.5d0
               deallocate(vtemp)
            endif
         enddo
      else
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(opi,opid)
            call free_memory(opi)
         endif
         call collect_v_pot(vpi,vpid)
         call free_memory(vpi)
         if(options%lnonorthogonal) then
            call collect_ortho_polaw(op,opd)
            call free_memory(op)
         endif
         do ii=options%i_min,options%i_max
            cpp%cprim=ii
            call read_data_pw_cprim_prod(cpp, options%prefix,.true.,ok_read,.false.,.true.)
         enddo
      endif
   enddo
   call mp_sum(sene(-ss%n:ss%n,:))
   do ii=options%i_min,options%i_max
       do it=-ss%n,ss%n
          ss%diag(ii,it+ss%n+1)=ss%diag(ii, it+ss%n+1)+sene(it,ii)
       enddo
    enddo

 !copy sene results on ss with opportune factors
!!!!!!!!!!!


   call free_memory(vpid)
   if(options%lnonorthogonal) then
      call free_memory(opd)
      call free_memory(opi)
      call free_memory(opid)
   endif
   call free_memory_polaw(ww)
   call free_memory( cpp)
   call free_memory(cppd)
   call free_memory( wu)
   call free_memory(us)
   deallocate(sene)
  return
 END SUBROUTINE selfenergy_ontime_upper


 END MODULE
