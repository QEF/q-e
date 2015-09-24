!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

 SUBROUTINE go_fft(tf, options) 
!this subroutine perform FFT on polarization written on disk 

   USE kinds,              ONLY : DP
   USE io_global,          ONLY : stdout, ionode
   USE input_gw,           ONLY : input_options
   USE polarization,       ONLY : polaw,free_memory_polaw,read_polaw
   USE fft_gw
   USE times_gw,           ONLY : times_freqs

   implicit none

   TYPE(times_freqs) :: tf
   TYPE(input_options) :: options! for imaginary time range,number of samples and number of rows

   TYPE(fft_data) :: fftd
   TYPE(polaw) :: pw
   INTEGER :: iw
   INTEGER :: numpw
   INTEGER :: firstr,lastr

!read in polarization for checks

   if(ionode) then !ATTENZIONE provvisorio

      call read_polaw(options%n,pw,options%debug,options%l_verbose)
      numpw=pw%numpw

      if(numpw < options%num_rows) then
         write(stdout,*) 'Routine go_fft: num_rows too big'
         stop
      endif
      call free_memory_polaw(pw)

!first loop, construct fftd and save on disk

      do iw=1,ceiling(real(numpw)/real(options%num_rows))
!set limits
         firstr=(iw-1)*options%num_rows+1
         lastr=iw*options%num_rows
         if(lastr > numpw) lastr=numpw

!create fft descriptor
 
         call create_fft_data(tf,firstr,lastr,options%tau,options%n,iw,fftd,options%debug)

!write on file

         call write_fft_data(fftd,options%debug)
      enddo


!sencond loop, read fftd, perform fftd, and update polaw's

      do iw=1,ceiling(real(numpw)/real(options%num_rows))
!set limits
         firstr=(iw-1)*options%num_rows+1
         lastr=iw*options%num_rows
         if(lastr > numpw) lastr=numpw

! read fftd descriptor

         call read_fft_data(iw,fftd,options%debug)

!do fft

         call transform_fft_data(fftd)

!update polaw's

         call save_fft_data(tf,fftd,options%debug)
         
         call write_fft_data(fftd,options%debug)

      enddo

   endif
   call free_memory_fft_data(fftd)
   return
 END SUBROUTINE go_fft
 

 SUBROUTINE go_fft_para(tf, options)
!this subroutine perform FFT on polarization written on disk
 
   USE kinds,              ONLY : DP
   USE io_global,          ONLY : stdout, ionode
   USE input_gw,           ONLY : input_options
   USE polarization,       ONLY : polaw,free_memory_polaw,read_polaw_range, initialize_polaw,&
                                    &fake_polarization_io
   USE fft_gw
   USE mp_world,           ONLY : world_comm, mpime, nproc
   USE mp,                 ONLY : mp_barrier
   USE times_gw,           ONLY : times_freqs


   implicit none
 
   TYPE(times_freqs), INTENT(in) :: tf!for time and frequency grids
   TYPE(input_options), INTENT(in) :: options! for imaginary time range,number of samples and number of rows
 
   TYPE(fft_data) :: fftd
   TYPE(polaw) :: pw
   INTEGER :: iw
   INTEGER :: numpw
   INTEGER :: firstr,lastr
 
   INTEGER :: number_fft
   LOGICAL, ALLOCATABLE :: is_my_fft(:)
   INTEGER  number_my_fft, done_fft, fft_first
   INTEGER :: iqq,ndelta, it, ip
   LOGICAL :: ontime


!read in polarization for checks
   
      call initialize_polaw(pw)

      call read_polaw_range(options%n,pw,options%debug,1,1,.true.)

      numpw  = pw%numpw
      ontime = pw%ontime
 
      if(numpw < options%num_rows) then
         write(stdout,*) 'Routine go_fft: num_rows too big'
         stop
      endif
      write(stdout,*) 'Prima'!ATTENZIONE
      call free_memory_polaw(pw)
 
      number_fft = ceiling(real(numpw)/real(options%num_rows))
      number_my_fft=0
      ndelta = number_fft / nproc
      if( nproc*ndelta < number_fft ) ndelta = ndelta + 1
      allocate( is_my_fft( 0 : nproc*ndelta ) )
      is_my_fft(:)=.false.
      iqq = 1
      fft_first = -1
      do ip = 0, nproc - 1
         do it = 1, ndelta
           ! write(stdout,*) 'fft -> ', ip, it, iqq, number_fft
           if( iqq <= number_fft .and. ( mpime == ip ) ) then
              if(fft_first  == -1) fft_first = iqq
              is_my_fft( iqq ) = .true.
              number_my_fft = number_my_fft + 1
           endif
           iqq = iqq + 1
         enddo
       enddo

      call mp_barrier( world_comm )
!sencond loop, read fftd, perform fftd, and update polaw's
 
      done_fft=0
      do it = 1, ndelta
       iw = it + fft_first - 1
       if( is_my_fft( iw ) .and. fft_first /= -1 ) then
           done_fft=done_fft+1
!set limits
           firstr=(iw-1)*options%num_rows+1
           lastr=iw*options%num_rows

           ! write(stdout,*) 'Continue FFT: ',iw , firstr, lastr, numpw

           if(lastr > numpw) lastr=numpw
 
! read fftd descriptor

           write(stdout,*) 'Create iw = ', iw ! ATTENZIONE 
           call create_fft_data(tf,firstr,lastr,options%tau,options%n,iw,fftd,options%debug)

           ! avoid reading time/freq parameter from matrix,
           ! use the one read at the beginning.

           fftd%ontime = ontime  
 
!do fft
           write(stdout,*) 'Transform'!ATTENZIONE
           if(options%l_fft_timefreq) then
              call transform_fft_data(fftd)
           else
              call transform_fft_data_grid(tf,fftd)
           endif
 
!update polaw's
           write(stdout,*) 'Save'!ATTENZIONE
           call save_fft_data(tf,fftd,options%debug)
          
        else
           call fake_polarization_io(options%n)
        endif
      enddo
 
   call free_memory_fft_data(fftd)
   deallocate(is_my_fft)
   return
 END SUBROUTINE go_fft_para



 SUBROUTINE go_fft_para2(tf, options)
!this subroutine perform FFT on polarization written on disk
 
   USE kinds,              ONLY : DP
   USE io_global,          ONLY : stdout, ionode
   USE input_gw,           ONLY : input_options
   USE polarization,       ONLY : polaw,free_memory_polaw,read_polaw_range, initialize_polaw,&
                                    &fake_polarization_io
   USE fft_gw
   USE mp_world,           ONLY : world_comm, mpime, nproc, root
   USE mp,                 ONLY : mp_barrier, mp_bcast
   USE times_gw,           ONLY : times_freqs


   implicit none
 
   TYPE(times_freqs), INTENT(in) :: tf!for time and frequency grids
   TYPE(input_options), INTENT(in) :: options! for imaginary time range,number of samples and number of rows
 
   TYPE(fft_data) :: fftd
   TYPE(polaw) :: pw
   INTEGER :: iw
   INTEGER :: numpw
   INTEGER :: firstr,lastr
 
   INTEGER :: number_fft
   INTEGER  number_my_fft, done_fft, fft_first
   INTEGER :: iqq,ndelta, it, ip, iblk, ipown
   LOGICAL :: ontime
   COMPLEX(kind=DP) :: factor

!read in polarization for checks

   
      call initialize_polaw( pw )

       write(stdout,*) 'Routine go_fft_para2'!ATTENZIONE
       FLUSH(stdout)

      ipown = 0
      if( mpime == ipown ) then
          write(stdout,*) 'Prima read_polaw_range'!ATTENZIONE
         call read_polaw_range( 0, pw, options%debug, 1, 1, .true. )
         write(stdout,*) 'Dopo read_polaw_range'!ATTENZIONE
      end if
      call mp_bcast( pw%numpw,  ipown, world_comm )
      call mp_bcast( pw%ontime, ipown, world_comm )
      call mp_bcast( pw%factor, ipown, world_comm )

      numpw  = pw%numpw
      ontime = pw%ontime
      factor = pw%factor
 
      if( numpw < options%num_rows ) then
         write(stdout,*) 'Routine go_fft: num_rows too big'
         stop
      endif
      write(stdout,*) 'Prima',factor!ATTENZIONE
      call free_memory_polaw(pw)

      !
      ! options%num_rows, number of rows to be read by each proc.
      ! from its own files
      ! number_fft, number of collective read cycles
      ! 
 
      number_fft = ceiling( real(numpw)/ real(options%num_rows) )

      call mp_barrier( world_comm )

      do iblk = 1, number_fft
         !
         firstr = (iblk-1) * options%num_rows + 1
         lastr  = (iblk  ) * options%num_rows
         if( lastr > numpw ) lastr = numpw
         !
         fftd%ontime = ontime  
         fftd%numpw  = numpw  
         fftd%factor = factor

         CALL create_fft_data2( tf, firstr, lastr, options%tau, options%n, fftd, options%debug )

         write(stdout,*) 'Transform'!ATTENZIONE
         FLUSH(stdout)

         if( options%l_fft_timefreq ) then
            call transform_fft_data( fftd )
         else
            call transform_fft_data_grid( tf, fftd )
         endif

         write(stdout,*) 'Save'!ATTENZIONE
         FLUSH(stdout)
         call save_fft_data2( tf, fftd, options%debug )
         
      end do
      


   call free_memory_fft_data(fftd)
   return
 END SUBROUTINE go_fft_para2


