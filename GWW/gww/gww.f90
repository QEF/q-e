!
! P.Umari  program GWW
! Modified by G. Stenuit
!
! This is the main part of the GWW program

   PROGRAM gww

   USE  input_gw,            ONLY : input_options, read_input_gww
   USE  io_global,           ONLY : stdout, ionode
   USE  self_energy_storage
   USE  expansion
   USE  energies_gww
   USE  start_end
   USE  mp_global,           ONLY : mpime
   USE  para_gww
   USE  times_gw
   USE  w_divergence


   implicit none

   TYPE(input_options)   :: options
   TYPE(self_storage)    :: ss
   TYPE(self_expansion)  :: se
   TYPE(quasi_particles) :: qp
   TYPE(times_freqs)     :: tf
   TYPE(gv_time)         :: gt
   real*8 :: tini, t1, t2, t3, t4, t5, t6, t7, tfinal

!setup MPI environment

   call startup



   write(stdout,*)  'PROGRAM GWW: Version 0.3'
   call cpu_time(tini)


!initialize arrays
   call initialize_quasi_particle(qp)


!  read in input structure

   call read_input_gww(options)
   call flush_unit(stdout)
   call setup_para_gww(options%n, options%max_i, options%i_min, options%i_max)
   call flush_unit(stdout)
! setup time/frequency grid if required
   write(stdout,*) 'ATT1'
   call flush_unit(stdout)
   call setup_timefreq(tf,options)
   write(stdout,*) 'ATT2'
   call flush_unit(stdout)

   call flush_unit(stdout)
!  call test_fft(tf)

!Step 0
!calculates the exchange energies

   if(options%starting_point <=1) then
      call cpu_time(t1)
      write(stdout,*) 'starting_point <=1'
      call flush_unit(stdout)
      call  go_exchange_main( options, qp)
      call write_quasi_particles(qp, options,.false.)
   else
      if(options%starting_point >= 6 .and. (options%remainder==3 .or. options%remainder==4)) then
         call read_quasi_particles(qp,options,.true.)
      else
         call read_quasi_particles(qp,options,.false.)
      endif
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
   if(options%starting_point <= 2 .and. options%ending_point >=2 ) then
      call cpu_time(t2)
      write(stdout,*) 'starting_point <=2'
      call flush_unit(stdout)
     call go_polarization(tf,options)
   endif


   if(options%starting_point <= 3 .and. options%ending_point >= 3 ) then
      call cpu_time(t3)
      write(stdout,*) 'starting_point <=3'
      call flush_unit(stdout)
!Step 3
!FFT of polarization to imaginary frequency and save on file
      write(stdout,*) 'Call FFT'
      call go_fft_para2(tf, options)  ! TEST new fft
  endif

  if(options%starting_point<=4  .and. options%ending_point >= 4) then
      call cpu_time(t4)
      write(stdout,*) 'starting_point <=4'
      call flush_unit(stdout)


!Step 3.1
!calculate dresses interaction W, and save on file
     write(stdout,*) 'Call go_dressed_w'
     call go_dressed_w(options)


     if(options%remainder==4) call create_dressed_polarization( options)

     if(options%remainder==3 .or. options%remainder==4) then
        call remainder(options, qp)
        call write_quasi_particles(qp,options,.true.)
     else
        call write_quasi_particles(qp,options,.false.)
     endif
  endif


!Step 3.2
!FFT of W to imaginary time and save on file
  if(options%starting_point<=5  .and. options%ending_point >= 5) then
      call cpu_time(t5)
      write(stdout,*) 'starting_point <=5'
      call flush_unit(stdout)
     if(options%remainder==3.or.options%remainder==4) then
        call read_quasi_particles(qp,options,.true.)
     else
        call read_quasi_particles(qp,options,.false.)
     endif
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
   endif

 WRITE(STDOUT,*) 'ENE_HF2' , QP%ENE_HF(2)!ATTENZIONE


   if(options%starting_point <= 6  .and. options%ending_point >= 6) then
      call cpu_time(t6)
      write(stdout,*) 'starting_point <=6'
      call flush_unit(stdout)
!Step 4
!calculate the expectation value of Sigma in imaginary time and save on file
     call create_self_ontime(tf, ss,options,qp)

     if(options%lconduction.and. .not.options%lvcprim_file .and. .not.options%l_self_beta) then
        if(.not.options%lcprim_file) then
           call addconduction_self_ontime(ss, options)
        else
           ! name shortened to fit into max allowed length - PG
           call addconduct_self_ontime_file(ss, tf, options)
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
   endif






    WRITE(STDOUT,*) 'ENE_HF' , QP%ENE_HF(2)!ATTENZIONE
   if(options%starting_point <= 7 .and. options%ending_point >= 7) then
      call cpu_time(t7)
      write(stdout,*) 'starting_point <=7'
      call flush_unit(stdout)
!Step 7
! fit self_energy with a multipole expansion
      call read_self_storage_ondisk(ss, options)
      call set_remainder(ss, qp)
      call create_self_energy_fit( tf, se, ss, options)
      call print_fit_onfile(tf, se,ss)
      call free_memory_self_storage(ss)
      call create_quasi_particles(options,qp,se)
      call free_memory_self_expansion(se)

      call printout_quasi_particles(qp)
   endif



    call cpu_time(tfinal)
    write(stdout,*) '---------------------------------------------------------'
    write(stdout,*) 'TIME REPORT : '
    write(stdout,*) '---------------------------------------------------------'
    write(stdout,*) 'time between starting_point=0 and starting_point=1: TIME=', t1-tini
    write(stdout,*) 'time between starting_point=1 and starting_point=2: TIME=', t2-t1
    write(stdout,*) 'time between starting_point=2 and starting_point=3: TIME=', t3-t2
    write(stdout,*) 'time between starting_point=3 and starting_point=4: TIME=', t4-t3
    write(stdout,*) 'time between starting_point=4 and starting_point=5: TIME=', t5-t4
    write(stdout,*) 'time between starting_point=5 and starting_point=6: TIME=', t6-t5
    write(stdout,*) 'time between starting_point=6 and starting_point=7: TIME=', t7-t6
    write(stdout,*) '---------------------------------------------------------'
    write(stdout,*) 'GWW :: FULL  TIME  =  ', tfinal-tini
    write(stdout,*) '---------------------------------------------------------'
    call flush_unit(stdout)

!stops MPI
    call free_memory_times_freqs(tf)
    call free_memory_para_gww
    call stop_run

   stop
   END PROGRAM gww

