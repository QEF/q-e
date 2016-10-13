!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!


!this is the main part of the GWW program

   PROGRAM gww

   USE  input_gw,            ONLY : input_options, read_input_gww
   USE  io_global,           ONLY : stdout, ionode
   USE  self_energy_storage
   USE  expansion
   USE  energies_gww
   USE  start_end
   USE  mp_world,           ONLY : mpime, world_comm
   USE  para_gww
   USE  times_gw
   USE  w_divergence
   USE  mp,                  ONLY : mp_barrier
   USE  contour
      USE io_files,  ONLY : prefix, tmp_dir

   implicit none

   TYPE(input_options)   :: options
   TYPE(self_storage)    :: ss
   TYPE(self_expansion)  :: se
   TYPE(self_on_real)    :: sr
   TYPE(quasi_particles) :: qp
   TYPE(times_freqs)     :: tf
   TYPE(gv_time)         :: gt
   TYPE(w_expectation)   :: we
   TYPE(w_poles)         :: wp

   INTEGER :: ispin
   CHARACTER(5) :: name_proc
   INTEGER :: ie
   REAL(kind=DP) :: energy

#if defined(__OPENMP)
   INTEGER :: omp_get_num_threads, omp_get_max_threads
   EXTERNAL omp_set_num_threads, omp_get_num_threads, omp_get_max_threads
#endif

   tmp_dir=' '

!setup MPI environment

   call startup

   !CALL remove_stack_limit ( )

#if defined(__OPENMP)
     ntids=omp_get_max_threads()
     ! call omp_set_num_threads(1)
#endif



#if defined(__OPENMP)
   write(stdout,*)  'ntids = ', ntids
#endif



!initialize arrays
   call initialize_quasi_particle(qp)


!  read in input structure

   call read_input_gww(options)
#if defined(__MPI)
   if(options%l_verbose) then
      write(name_proc,'(5i1)') &
           & (mpime+1)/10000,mod(mpime+1,10000)/1000,mod(mpime+1,1000)/100,mod(mpime+1,100)/10,mod(mpime+1,10)
      OPEN( UNIT = stdout, FILE = trim(tmp_dir)//trim(prefix)//'-out_'//name_proc, STATUS = 'UNKNOWN' )
   endif
#endif

   FLUSH(stdout)
   if(options%grid_freq/=5.and.options%grid_freq/=6) then
      call setup_para_gww(options%n, options%max_i, options%i_min, options%i_max)
   else
      call setup_para_gww(options%n+(1+2*options%second_grid_i)*options%second_grid_n, options%max_i, options%i_min, options%i_max)
   endif
   FLUSH(stdout)
! setup time/frequency grid if required
   call setup_timefreq(tf,options)

!Step 0
!calculates the exchange energies



   if(options%starting_point <=1) then
      call  go_exchange_main( options, qp)
      call write_quasi_particles(qp, options,.false.)
   else
      call read_quasi_particles(qp,options,.false.)
   endif


!Step 1
!create the Green function G_0 in imaginary time and save on file
!it also calculates here the exchage energies


   if(options%starting_point <= 1 .and. options%ending_point >= 1) then
      if(.not.options%lpola_file .and. .not. options%lvcprim_file) then
         call go_green(tf,options, qp)
      endif
   endif



!Step 2
!create the polarization in imaginary time and save on  file


!loop on spin

   do ispin=1,options%nspin

      if(options%starting_point <= 2 .and. options%ending_point >=2 ) then
         if(options%l_t_wannier) then
            call calculate_compact_pola_lanczos(options,ispin)
         endif
         
      endif

      if(options%starting_point <= 3 .and. options%ending_point >= 3 ) then
         write(stdout,*) "*******************************"
         write(stdout,*) "     RESTART FROM POINT 3"
         write(stdout,*) "*******************************"
!Step 3
!FFT of polarization to imaginary frequency and save on file

     
         call do_polarization_lanczos(tf,options,ispin)
     
      endif
   enddo

  if(options%starting_point<=4  .and. options%ending_point >= 4) then


!Step 3.1
!calculate dresses interaction W, and save on file



     write(stdout,*) 'Call go_dressed_w'


     call go_dressed_w(options)


  
     call write_quasi_particles(qp,options,.false.)
  


  endif


!Step 3.2
!FFT of W to imaginary time and save on file
  if(options%starting_point<=5  .and. options%ending_point >= 5) then


     
     call read_quasi_particles(qp,options,.false.)
     
     if(.not.  options%l_self_lanczos) then
        write(stdout,*) 'Call FFT'
        call go_fft_para2(tf, options)
!if required do fft of gt structure
        if(options%w_divergence==2) then
           write(stdout,*) 'Go fft gt'
           call initialize_gv_time(gt)
           write(stdout,*) 'Go fft gt 1'
           call read_gv_time(gt)
           write(stdout,*) 'Go fft gt 1.5'
           call fft_gv_time(gt,tf)
           write(stdout,*) 'Go fft gt2'
           call write_gv_time(gt)
           write(stdout,*) 'Go fft gt3'
           call free_memory_gv_time(gt)
        endif
     else
        call do_reducible_pola(tf ,options)

     endif
  endif


   if(options%starting_point <= 6  .and. options%ending_point >= 6) then
!Step 4
      write(stdout,*) '*******************************'
      write(stdout,*) '     RESTART FROM POINT 6'
      write(stdout,*) '*******************************'


      if(options%n_real_axis>0) then
         call initialize_w_expectation(we)
         call create_w_expectation(we, tf, options)
         call write_w_expectation(we)
         call free_memory_w_expectation(we)
      endif


      if(.not.  options%l_self_lanczos) then
!calculate the expectation value of Sigma in imaginary time and save on file
         call create_self_ontime(tf, ss,options,qp)

         if(options%lconduction.and. .not.options%lvcprim_file .and. .not.options%l_self_beta) then
            if(.not.options%lcprim_file) then
               call addconduction_self_ontime(ss, options)
            else
               call addconduction_self_ontime_file(ss, tf, options)
            endif
         endif
         if(options%l_self_upper) then
            call selfenergy_ontime_upper(ss, tf ,options)
         endif
         
         if(options%debug) call write_storage(tf,ss)
         
         if(options%l_fft_timefreq) then
            call fft_storage(ss)
         else
            if(tf%grid_fit==0) then
               call fft_storage_grid(tf,ss)
            else
               call fft_storage_grid_fit(tf, ss)
            endif
         endif
         if(options%debug) call write_storage(tf,ss)
         call write_self_storage_ondisk(ss, options)
      else
!lanczos calculation of self-energy
         if(options%n_real_axis==0) then
            if(.not.options%l_self_time) then
               call do_self_lanczos(ss, tf ,options)
            else
               if(.not.options%l_full) then
                  call do_self_lanczos_time(ss, tf ,options,.false.,0.d0)
               else
                  call do_self_lanczos_full(ss, tf ,options,.false.,0.d0)
               endif
               call fft_storage_grid_fit(tf, ss)
            endif
            call write_self_storage_ondisk(ss, options)
         else
            call do_self_on_real(options,tf,ss,sr)
            call write_self_on_real(sr,0)
            call free_memory_self_on_real(sr)
         endif
         call write_self_storage_ondisk(ss, options)
      endif


   endif





   if(options%starting_point <= 7 .and. options%ending_point >= 7) then
!Step 7
! fit self_energy with a multipole expansion
      call read_self_storage_ondisk(ss, options)
      call create_self_energy_fit( tf, se, ss, options,sr,.false.)
      call mp_barrier( world_comm )
      call print_fit_onfile(tf, se,ss)
      call mp_barrier( world_comm )
      call free_memory_self_storage(ss)
      call mp_barrier( world_comm )
      call create_quasi_particles(options,qp,se)
      call mp_barrier( world_comm )
      call write_self_expansion(se)
      call free_memory_self_expansion(se)
      call mp_barrier( world_comm )
      call printout_quasi_particles(qp)
   endif


   if(options%starting_point <= 8 .and. options%ending_point >= 8) then
!if the whole self_energy matrix has been calculate do use it for obtaining QPEs and QPAs
      if(options%whole_s) then
         call initialize_self_expansion(se)
         call read_self_expansion(se)

         call create_quasi_particles_off(options,qp,se)

         call printout_quasi_particles_off(qp)
         call free_memory_self_expansion(se)
      endif
   endif

    if(options%starting_point <= 9 .and. options%ending_point >= 9) then
!here does analytic continuation for contour integration
       call initialize_w_expectation(we)
       call initialize_w_poles(wp)
       call read_w_expectation(we)
       call create_w_poles(we,wp,options)
       call write_w_poles(wp)
       call free_memory_w_expectation(we)
       call free_memory_w_poles(wp)
    endif

    if(options%starting_point <= 10 .and. options%ending_point >= 10) then
!adds poles
       call initialize_w_poles(wp)
       call initialize_self_on_real(sr)
       call read_w_poles(wp)
       call read_self_on_real(sr,0)
       !call self_on_real_print(sr)
!NOT_TO_BE_INCLUDED_START
       call do_contour(sr,wp,options)
!NOT_TO_BE_INCLUDED_END
       call write_self_on_real(sr,1)
       call self_on_real_print(sr)
       call free_memory_w_poles(wp)
       call free_memory_self_on_real(sr)
    endif
    if(options%starting_point <= 11 .and. options%ending_point >= 11) then
       call initialize_self_on_real(sr)
       call read_self_on_real(sr,1)
       call create_quasi_particle_on_real(options,qp,sr)
       call printout_quasi_particles(qp)
       call free_memory_self_on_real(sr)
    endif
!stops MPI
    call free_memory_times_freqs(tf)
    call free_memory_para_gww

    call stop_run

   stop
   END PROGRAM gww

