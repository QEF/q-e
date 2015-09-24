!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!


 SUBROUTINE go_polarization(tf, options, qp)
!this subroutines reads in the green functions and calculates
!the polarization at every imaginary time
!only positive imaginary times are required!!
   USE kinds,              ONLY : DP
   USE input_gw,           ONLY : input_options
   USE basic_structures,   ONLY : q_mat, wannier_u,free_memory
   USE green_function,     ONLY : green,read_green,free_memory_green, initialize_green
   USE polarization,       ONLY : polaw,free_memory_polaw,write_polaw,create_polarization,&
                              &create_polarization_contraction,create_polarization_contraction_state,&
                              &create_polarization_file,create_polarization_beta, create_polarization_upper
   USE io_global,          ONLY : stdout, ionode
   USE compact_product,    ONLY : contraction_pola,do_contraction_pola,free_memory_contraction_pola, &
                                           &do_contraction_pola_state
   USE para_gww,           ONLY : is_my_pola
   USE mp,                 ONLY : mp_barrier
   USE mp_world,           ONLY : world_comm
   USE energies_gww,           ONLY : quasi_particles
   USE times_gw,           ONLY : times_freqs


   implicit none

   TYPE(times_freqs) ,  INTENT(in) :: tf!for time grid
   TYPE(input_options), INTENT(in) :: options! for imaginary time range and number of samples
   TYPE(quasi_particles),INTENT(in) :: qp!for HF energies

   TYPE(green) :: gp,gm
   TYPE(q_mat) :: qm
   TYPE(polaw) :: pp
   TYPE(wannier_u) :: uu
   TYPE(contraction_pola) :: cp
   REAL(kind=DP) :: time,dt
   INTEGER :: iw

!read in overlap matrix

   write(stdout,*) 'GO POLARIZATION' !ATTENZIONE
   FLUSH(stdout)

   if(.not. options%lpola_file) then

      call read_data_pw_q(qm,options%prefix,.false.)
      call initialize_green(gp)
      call initialize_green(gm)

!loop on time samples
      dt=options%tau/real(options%n)
      
      write(stdout,*) 'GO POLARIZATION1' !ATTENZIONE
      
      if(options%use_contractions .and. .not.options%l_pola_beta) then
         call read_data_pw_u(uu,options%prefix)
         write(stdout,*) 'Calculates contraction of polarization'
         if(.not.options%l_contraction_single_state) then
            call do_contraction_pola(qm,uu,cp)
         else
            call do_contraction_pola_state(qm, uu, options)
         endif
      endif
      
      if(options%l_pola_beta) then
         call read_data_pw_u(uu,options%prefix)
      endif

      do iw=0,options%n
         if(is_my_pola(iw)) then
            write(stdout,*) 'Calculate Polarization:', iw
            if(options%l_fft_timefreq) then
               time=dt*real(iw)
            else
               time=tf%times(iw)
            endif
            if(.not.options%l_pola_beta) then
               if(.not.options%use_contractions) then
                  write(*,*) 'CREATE POLARIZATION AND NO CONTRACTIONS'
                  stop
                  call read_green(iw,gp,options%debug,.false.)
                  call read_green(-iw,gm,options%debug,.true.)
                  call create_polarization(time,pp,gp,gm,qm,options%debug)
               else
                  if(.not.options%l_contraction_single_state) then
                     call create_polarization_contraction(time,pp,cp,uu,options%l_hf_energies,qp%ene_hf(:,1))
                  else
                     call create_polarization_contraction_state(time,pp,uu,options%l_hf_energies,qp%ene_hf(:,1),options)
                  endif
               endif
            else
               call  create_polarization_beta(time, pp, uu, qm)
            endif
            pp%label=iw
            call write_polaw(pp,options%debug)
!we take advantage of the P(t)=P(-t) symmetry
!         if(iw /= 0) then
!            pp%time=-time
!            pp%label=-iw
!            call write_polaw(pp,options%debug)
!         endif
         endif
      enddo

      call mp_barrier( world_comm )
      
      call free_memory(qm)
      call free_memory_green(gp)
      call free_memory_green(gm)
      call free_memory_polaw(pp)
      call free_memory_contraction_pola(cp)
   else
      call read_data_pw_u(uu,options%prefix)
      call create_polarization_file(uu, tf, options%prefix)
   endif
   if(options%l_pola_upper) then
      call create_polarization_upper(uu, tf, options%prefix)
   endif
   call free_memory(uu)
   return
 END SUBROUTINE
