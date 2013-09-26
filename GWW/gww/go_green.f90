!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!


 SUBROUTINE go_green(tf, options, qp)
!this subroutine at every imaginary time, calculate the green function
!and save it on file

   USE kinds,              ONLY : DP
   USE input_gw,           ONLY : input_options
   USE basic_structures,   ONLY : wannier_u, free_memory
   USE green_function,     ONLY : green,create_green_part,write_green,free_memory_green,initialize_green
   USE para_gww,           ONLY : is_my_time, is_my_last
   USE mp,                 ONLY : mp_barrier
   USE mp_world,           ONLY : world_comm
   USE io_global,          ONLY : stdout
   USE energies_gww,           ONLY : quasi_particles
   USE times_gw,        ONLY : times_freqs
   
   implicit none

   TYPE(times_freqs), INTENT(in)     :: tf!time grid
   TYPE(input_options), INTENT(in)   :: options! for imaginary time range and number of samples
   TYPE(quasi_particles), INTENT(in) :: qp!for the HF energies if required

   TYPE(wannier_u)  :: wu
   TYPE(green) :: gr
   INTEGER :: iw
   REAL(kind=DP) :: time, dt


   call initialize_green(gr)

!read in U tranformation matrix and KS eneregies

   call read_data_pw_u(wu,options%prefix)

!loop on samples

   dt=options%tau/real(options%n)

   do iw=-options%n,options%n
      if(is_my_time(iw)) then
         write(stdout,*) 'Green: ', iw, time
         if(options%l_fft_timefreq) then
            time=dt*real(iw)
         else
            time=tf%times(iw)
         endif
         call create_green_part(gr,wu,time,options%debug,.false.,options%l_hf_energies, qp%ene_hf(:,1))
         gr%label=iw
         write(stdout,*) 'Green created: ', iw, time
         call write_green(gr,options%debug)
      endif
   enddo
  
!now insert the zero time negative one
   
   if(is_my_last) then
      write(stdout,*) 'green 0'
      call create_green_part(gr,wu,0.d0,options%debug,.true.,options%l_hf_energies, qp%ene_hf(:,1))
      gr%label=0
      call write_green(gr,options%debug)
      write(stdout,*) 'green 0 created'
   endif

   call mp_barrier( world_comm )

   call free_memory_green(gr)
   call free_memory(wu)
   return

 END SUBROUTINE
 
 
