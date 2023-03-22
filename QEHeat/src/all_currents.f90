!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
! by Riccardo Bertossa (SISSA)
! written during the 2020 and 2021 pandemic years
!
! implements Marcolongo, Umari, and Baroni, Nat. Phys. 12, 80 (2016)
! details of the original implementation are described in
!
! Marcolongo, Bertossa, Tisi and Baroni, Computer Physics Communications, 269, 108090 (2021)
! https://doi.org/10.1016/j.cpc.2021.108090
! https://arxiv.org/abs/2104.06383
!-----------------------------------------------------------------------

program all_currents
   ! code that orchestrate all the calculation of the energy current

   ! module to compute the kohn_sham part
   use kohn_sham_mod, only: init_kohn_sham, current_kohn_sham
   ! module to compute the hartree and exchange-correlation part
   use hartree_xc_mod, only: current_hartree_xc
   ! code to compute the pseudopotential part
   use zero_mod, only: allocate_zero, deallocate_zero, init_zero, current_zero
   ! code to compute the coulomb part
   use ionic_mod, only: init_ionic, ionic_init_type, current_ionic
   ! utilities to hold multiple wavefunctions and atomic position/velocities data
   use scf_result_mod, only: multiple_scf_result_allocate, &
                             scf_result_set_from_global_variables, multiple_scf_result_deallocate, &
                             multiple_scf_result, scf_result_set_global_variables

   USE environment, ONLY: environment_start, environment_end
   use io_global, ONLY: ionode
   use wavefunctions, only: evc
   use wvfct, only: nbnd, npwx, npw
   use kinds, only: dp

   ! utilities to read the cp.x produced trajectory
   use cpv_traj, only: cpv_trajectory, &
                       cpv_trajectory_initialize, cpv_trajectory_deallocate

   use dynamics_module, only: vel
   use ions_base, ONLY: tau, tau_format, nat
   USE control_flags, ONLY: ethr, use_gpu
   USE extrapolation, ONLY: update_pot

!from ../PW/src/pwscf.f90
   USE mp_global, ONLY: mp_startup, ionode_id
   USE mp_world, ONLY: world_comm
   use mp, ONLY: mp_bcast, mp_barrier
   USE mp_pools, ONLY: intra_pool_comm
   USE mp_bands, ONLY: intra_bgrp_comm, inter_bgrp_comm
   USE read_input, ONLY: read_input_file
   USE command_line_options, ONLY: input_file_, command_line, ndiag_, nimage_
   USE check_stop, ONLY: check_stop_init

!from ../Modules/read_input.f90
   USE read_namelists_module, ONLY: read_namelists

   use input_parameters, only : outdir
   USE read_cards_module, ONLY: read_cards
   use averages, only: online_average, online_average_init

   use wvfct, only: nbnd, npw, npwx
   use wavefunctions, only: psic, evc
   use gvect, only: g, ngm, gstart, gg, igtongl, gl, ngl
   USE cell_base, ONLY: tpiba, omega, tpiba2, alat, at, bg
   use ions_base, only: tau, nsp, zv, nat, ityp, amass
   use uspp, ONLY: vkb, nkb, deeq
   use uspp_param, ONLY: upf, nh, nbetam
   use uspp_data, only: dq, nqxq
   use klist, only: xk, igk_k
   use wvfct, ONLY: g2kin, et
   use fft_base, only: dffts
   use atom, only: rgrid
   ! testing only!!!
   use test_h_psi, only: init_test

   implicit none

   type J_all
      !holds all the result calculated by this program
      real(dp) :: i_current(3), i_current_a(3), i_current_b(3), i_current_c(3), i_current_d(3), i_current_e(3)
      real(dp) ::z_current(3)
      real(kind=DP) ::J_kohn(3), J_kohn_a(3), J_kohn_b(3), J_hartree(3), J_xc(3), J_electron(3)
      real(kind=DP), allocatable :: v_cm(:, :)
   end type
   type(J_all) :: j

   integer :: exit_status, ios, irepeat, n_max
   logical :: print_stat
   type(cpv_trajectory) :: traj
   type(ionic_init_type) :: ionic_data
   real(kind=dp), allocatable :: tau_save(:, :), &
                                 tabr(:, :, :, :), H_g(:, :, :, :) ! current zero
   type(multiple_scf_result) :: scf_all ! kohn_sham & hartree
   complex(kind=DP), allocatable :: dvpsi_save(:, :, :) ! to save the solution of the system between iterations (kohn_sham
   real(kind=dp) :: vel_factor

   real(kind=dp) :: eta
   logical :: save_dvpsi = .true. ! if true dvpsi_save is allocated and used
   CHARACTER(len=256) :: file_output, trajdir = ''
   type(online_average) :: ave_cur
   real(kind=DP) ::delta_t, ethr_small_step, ethr_big_step
   integer :: first_step, last_step, step_mul, step_rem, n_workers, worker_id, &
              n_repeat_every_step, n_digit, &
              first_s_chunk, last_s_chunk, steps_per_chunk
   logical :: restart ! if true try to read last calculated step from output and set first_step
   logical :: subtract_cm_vel ! if true do velocity renormalization
   logical :: re_init_wfc_1 = .false., re_init_wfc_2 = .false. ! initialize again evc before scf step number 1 or 2
   logical :: re_init_wfc_3 = .false. ! initialize again evc before scf step number 1 or 2
   logical :: three_point_derivative = .true. ! compute hartree derivative with 3 points
   logical :: add_i_current_b = .false.  ! if true adds i_current_b to the final result, in the output. 
   ! note: i_current_b is proportional to the ionic velocities. In principle is not needed to calculate the thermal
   ! conductivity since it does not influence the final result. It is implemented only for a cubic cell.

   character(len=256) :: vel_input_units = 'PW', worker_id_char, format_string
   logical :: ec_test, hpsi_test ! activates tests for debugging purposes
   logical :: continue_not_converged ! don't stop the calculation if a step does not converge
   LOGICAL,EXTERNAL :: check_gpu_support
   !from ../PW/src/pwscf.f90
   include 'laxlib.fh'

!from ../PW/src/pwscf.f90
   use_gpu = check_gpu_support()
   if(use_gpu) Call errore('QEHeat', 'QEHeat with GPU NYI.', 1)
   CALL mp_startup( images_only=.true. )
   CALL environment_start('QEHeat')
   call start_clock('all_currents')
   IF (ionode) THEN
      write (*,*) ''
      write (*,*) '============================================================'
      write (*,*) '============================================================'
      write (*,*) ' This code implements Marcolongo, A., Umari, P. and Baroni, S'
      write (*,*) '  Nature Phys 12, 80-84 (2016). https://doi.org/10.1038/nphys3509'
      write (*,*) ''
      write (*,*) ' The details of the implementation are described in'
      write (*,*) '  Marcolongo, Bertossa, Tisi and Baroni'
      write (*,*) '  Computer Physics Communications, 269, 108090 (2021)'
      write (*,*) '  https://doi.org/10.1016/j.cpc.2021.108090'
      write (*,*) '  https://arxiv.org/abs/2104.06383'
      write (*,*) '============================================================'
      write (*,*) '============================================================'
      write (*,*) ''
      CALL input_from_file()
      ! all_currents input
      call read_all_currents_namelists(5, &
                                       delta_t, &
                                       file_output, trajdir, vel_input_units, &
                                       eta, n_max, first_step, last_step, &
                                       ethr_small_step, ethr_big_step, &
                                       restart, subtract_cm_vel, step_mul, &
                                       step_rem, ec_test, add_i_current_b, &
                                       save_dvpsi, re_init_wfc_1, re_init_wfc_2, &
                                       re_init_wfc_3, three_point_derivative, &
                                       n_repeat_every_step, hpsi_test, &
                                       n_workers, worker_id, &
                                       continue_not_converged)
     !if the problem is parallelized simply by running many times the code over the same trajectory
     !with a different starting and ending timestep you can use the n_worker and the worker_id variables
     !
     ! if n_workers > 0, append worker_id to file_output
     ! set first/last step accordingly
     ! append worker_id to outdir
     if (n_workers>0 ) then
        if (worker_id >= n_workers .or. worker_id<0) then
           call errore ('all_currents', 'worker_id must be one of 0, 1, ..., n_workers-1')
        end if
        n_digit = floor(log10(real(n_workers+1)))
        write (format_string, '(A2,I1,A1)') "(I",n_digit, ")"

        write (worker_id_char, format_string) worker_id
        file_output=trim(file_output) // '.'//trim(worker_id_char)
        ! calculate first step / last step for the chunk
        steps_per_chunk = (last_step-first_step + 1)/n_workers
        if (steps_per_chunk < 0) call errore('all_currents', 'last_step must be greater than first_step',1)
        if (steps_per_chunk == 0 ) then
           steps_per_chunk = 1
           write(*,*) 'WARNING: n_workers is too high: some chunks will have no work' 
        end if
        first_s_chunk = first_step + steps_per_chunk*worker_id
        if (worker_id < n_workers - 1 ) then
           last_s_chunk = first_step + steps_per_chunk*(worker_id+1)
        else
           last_s_chunk = last_step
        end if

        write (*,*) 'This worker has steps from ', first_s_chunk, ' to ', last_s_chunk
        first_step=first_s_chunk
        last_step=last_s_chunk - 1

     end if

   endif
   call bcast_all_current_namelist( &
      delta_t, &
      file_output, trajdir, vel_input_units, &
      eta, n_max, first_step, last_step, &
      ethr_small_step, ethr_big_step, &
      restart, subtract_cm_vel, step_mul, &
      step_rem, ec_test, add_i_current_b, &
      save_dvpsi, re_init_wfc_1, re_init_wfc_2, &
      re_init_wfc_3, three_point_derivative, &
      n_repeat_every_step, hpsi_test, &
      n_workers, worker_id, &
      continue_not_converged)
   ! PW input
   call read_namelists('PW', 5)
   if (n_workers>0 ) then
      CALL mp_bcast(worker_id_char, ionode_id, world_comm)
      outdir = trim(outdir) //trim(worker_id_char)
   endif
   
   call read_cards('PW', 5)

   call check_input()

   call mp_barrier(intra_pool_comm)
   if (vel_input_units == 'CP') then ! atomic units of cp are different
      vel_factor = 2.0_dp
      if (ionode) &
         write (*, *) 'Using CP units for velocities'
   else if (vel_input_units == 'PW') then
      if (ionode) &
         write (*, *) 'Using PW units for velocities'
      vel_factor = 1.0_dp
   else
      call errore('all_currents', 'error: unknown vel_input_units', 1)
   endif

   ! try to see if there was in interrupted calculation
   call set_first_step_restart(restart, file_output, first_step)
   call iosys()    ! ../PW/src/input.f90    save in internal variables
   call check_stop_init() ! ../PW/src/input.f90

   !eventually read new positions and velocities from trajectory
   if (ionode) then
      !initialize trajectory reading
      call cpv_trajectory_initialize(traj, trajdir, nat, 1.0_dp, vel_factor, 1.0_dp, ios=ios, circular=.true.)
      if (ios == 0) then
         write (*, *) 'After first step from input file, I will read from the CPV trajectory ', trim(trajdir)
      else
         write (*, *) 'Cannot open trajectory file', trim(trajdir), '. I can calculate only a single step from input file'
      endif
   endif
   if (first_step > 0) then
      if (ionode) &
         write (*, *) 'SKIPPING STEP FROM INPUT FILE'
      if (.not. read_next_step(traj, first_step, last_step, step_mul, step_rem)) then
         if (ionode) &
            write (*, *) 'NOTHING TO DO IN TRAJECTORY FILE'
         goto 100 ! skip everything and exit
      endif
   end if

   call setup()    ! ../PW/src/setup.f90    setup the calculation
   call init_run() ! ../PW/src/init_run.f90 allocate stuff
   ! now scf is ready to start, but I first initialize energy current stuff

   ! pseudopotential (zero) current initialization
   call allocate_zero() ! only once per all trajectory
   ! current zero (pseudopotential part) quantities that do not depend on the positions/velocities but only on the cell and atomic types
   allocate (H_g(ngm, 3, 3, nsp))
   allocate (tabr(nqxq, nbetam, nsp, 3))
   call init_zero(tabr, H_g, nsp, zv, tpiba2, tpiba, omega, at, alat, &
                  ngm, gg, gstart, g, igtongl, gl, ngl, dq, &
                  upf, rgrid, nqxq, intra_bgrp_comm, nat, ityp) ! only once per all trajectory
   ! coulomb (ionic) current initialization
   call init_ionic(ionic_data, eta, n_max, ngm, gstart, at, alat, omega, gg, g, tpiba2)
   ! kohn-sham initialization
   call init_kohn_sham()
   if (save_dvpsi) then ! to use the previous result as initial guess of the solution of the system solved in project.f90
      if (.not. allocated(dvpsi_save)) then
         allocate (dvpsi_save(npwx, nbnd, 3))
         dvpsi_save = (0.d0, 0.d0)
      end if
   end if

   call setup_nbnd_occ() ! only once per all trajectory

   if (ionode .and. first_step == 0) then !set velocities factor also in the input file step
      vel = vel*vel_factor
      call convert_tau(tau_format, nat, vel)
   end if
   CALL mp_bcast(vel, ionode_id, world_comm)

   !initialize scf results object. Later it will saves evc, tau and vel arrays for t-dt, t and t+dt
   call multiple_scf_result_allocate(scf_all, .true.)
   if (n_repeat_every_step > 1) then
      allocate (tau_save(3, nat))
   end if

   !loop over input trajectory steps starts here
   do
      if (ionode) then
         if (subtract_cm_vel) then
            !calculates center of mass velocity for each atomic species and subtract it
            call cm_vel(j%v_cm, vel)
         else
            !calculates center of mass velocity only, without modifying velocities
            call cm_vel(j%v_cm)
         end if
      endif
      if (n_repeat_every_step > 1) then
         tau_save = tau
         !initialize average object. Average is computed in write_results sub
         call online_average_init(ave_cur, .true.)
      end if
      !the same step can be calculated more than one time, to estimate stability or variance of the result
      do irepeat = 1, n_repeat_every_step
         if (irepeat > 1) then
            if (ionode) &
               write (*, *) 'REPETITION ', irepeat - 1
            tau = tau_save
         end if
         if (n_repeat_every_step > 1 .and. irepeat .eq. n_repeat_every_step) then
            print_stat = .true.
         else
            print_stat = .false.
         end if

         call prepare_next_step(-1, delta_t, ethr_small_step, three_point_derivative) !-1 goes back by dt, so we are in t-dt. Inside, after setting tau, we call hinit1 and update_pot
         if (re_init_wfc_1) then
            call init_wfc(1)
            call sum_band()
         end if
         ethr = ethr_big_step
         call run_electrons(exit_status, continue_not_converged)
         if (exit_status == 2 .and. continue_not_converged) then 
              continue
         else if (exit_status /= 0) then
             goto 100 !shutdown everything and exit
         end if
         !save evc, tau and vel for t-dt
         call scf_result_set_from_global_variables(scf_all%t_minus)
         if (three_point_derivative) then
            call prepare_next_step(1, delta_t, ethr_small_step, three_point_derivative) !1 advance by dt (so we are in the original positions)
            if (re_init_wfc_2) then ! eventually, to set a random initial evc to do statistical tests
               call init_wfc(1)
               call sum_band()
            end if
            ethr = ethr_small_step
            call run_electrons(exit_status, continue_not_converged)
            !evc_due = evc
            if (exit_status == 2 .and. continue_not_converged) then 
               continue
            else if (exit_status /= 0) then
               goto 100 !shutdown everything and exit
            end if
            !save evc, tau and vel for t
            call scf_result_set_from_global_variables(scf_all%t_zero)
            if (hpsi_test) &
               call init_test(evc) ! TESTING ONLY
            call current_zero(j%z_current, tabr, H_g, &
                              nbnd, npwx, npw, dffts, nsp, zv, nat, ityp, amass, tau, &
                              vel, tpiba, tpiba2, at, alat, omega, psic, evc, ngm, gg, g, gstart, &
                              nkb, vkb, deeq, upf, nh, xk, igk_k, bg, ec_test) ! routine zero should be called in t
            call current_ionic(ionic_data, &
                          j%i_current, j%i_current_a, j%i_current_b, j%i_current_c, j%i_current_d, j%i_current_e, add_i_current_b, &
                               nat, tau, vel, zv, ityp, alat, at, bg, tpiba, gstart, g, gg, npw, amass)
         else
            call scf_result_set_from_global_variables(scf_all%t_zero) !if we don't have 3pt derivative, zero and minus are equal
         end if

         call prepare_next_step(1, delta_t, ethr_small_step, three_point_derivative) !1 advances by dt, so we are in t+dt
         !if we don't do 3pt we are in t now

         if (re_init_wfc_3) then
            call init_wfc(1)
            call sum_band()
         end if
         ethr = ethr_small_step
         call run_electrons(exit_status, continue_not_converged)
         if (exit_status == 2 .and. continue_not_converged) then 
              continue
         else if (exit_status /= 0) then 
              goto 100 !shutdown everything and exit
         end if
         !save evc, tau and vel for t+dt
         call scf_result_set_from_global_variables(scf_all%t_plus)

         if (three_point_derivative) then
            ! restore wfct and potentials for t=0 (needed only if last point was t+dt)
            call scf_result_set_global_variables(scf_all%t_zero)
         else
            ! we are in t in this case, so we call here routines that do not compute numerical derivatives
            if (hpsi_test) &
               call init_test(evc) ! TESTING ONLY
            call current_zero(j%z_current, tabr, H_g, &
                              nbnd, npwx, npw, dffts, nsp, zv, nat, ityp, amass, tau, &
                              vel, tpiba, tpiba2, at, alat, omega, psic, evc, ngm, gg, g, gstart, &
                              nkb, vkb, deeq, upf, nh, xk, igk_k, bg, ec_test)
            call current_ionic(ionic_data, &
                          j%i_current, j%i_current_a, j%i_current_b, j%i_current_c, j%i_current_d, j%i_current_e, add_i_current_b, &
                               nat, tau, vel, zv, ityp, alat, at, bg, tpiba, gstart, g, gg, npw, amass)
         end if
         !calculate second part of energy current

         call current_hartree_xc(three_point_derivative, delta_t, scf_all, &
                                 j%j_hartree, j%j_xc, nbnd, npw, npwx, dffts, psic, g, ngm, gstart, &
                                 tpiba, omega, tpiba2)
         if (.not. three_point_derivative) &
            scf_all%t_zero%evc = scf_all%t_plus%evc
         call current_kohn_sham(j%J_kohn, j%J_kohn_a, j%J_kohn_b, j%J_electron, delta_t, scf_all, &
                                dvpsi_save, save_dvpsi, &
                                nbnd, npw, npwx, dffts, evc, g, ngm, gstart, &
                                tpiba2, at, vkb, nkb, xk, igk_k, g2kin, et, hpsi_test, &
                                omega, gg, intra_bgrp_comm, nat, ityp)
         call write_results(traj, print_stat, j, ave_cur)
      end do
      !read new velocities and positions and continue, or exit the loop
      if (.not. read_next_step(traj, first_step, last_step, step_mul, step_rem)) exit
   end do

   call stop_clock('all_currents')
   call print_clock('all_currents')
   call start_clock('PWSCF')
   ! shutdown stuff
100 call laxlib_end()
   call cpv_trajectory_deallocate(traj)
   call deallocate_zero()
   call multiple_scf_result_deallocate(scf_all)
   if (allocated(tau_save)) deallocate (tau_save)
   if (allocated(dvpsi_save)) deallocate (dvpsi_save)
   call stop_run(exit_status)
   call do_stop(exit_status)
   stop

contains

   subroutine write_results(traj, write_ave, j, ave_cur)
      use kinds, only: dp
      use ions_base, ONLY: nsp, na
      use cell_base, only: alat
      use io_global, ONLY: ionode
      use cpv_traj, only: cpv_trajectory, cpv_trajectory_get_last_step
      use traj_object, only: timestep
      use averages
      implicit none
      type(cpv_trajectory), intent(in)  :: traj
      logical, intent(in) :: write_ave
      type(j_all), intent(in) :: j
      type(online_average), intent(inout) :: ave_cur
      type(timestep) :: ts
      integer :: iun, icoord, step, itype
      integer, external :: find_free_unit
      real(dp) :: time, J_tot(3)
      logical :: file_exists

      INQUIRE (FILE=trim(file_output)//'.dat', EXIST=file_exists)
      if (traj%traj%nsteps > 0) then
         call cpv_trajectory_get_last_step(traj, ts)
         step = ts%nstep
         time = ts%tps
      else
         step = 0
         time = 0.d0
      endif
      if (ionode) then
         iun = find_free_unit()
         open (iun, file=trim(file_output), position='append')
         write (iun, *) 'step: ', step, ' dt: ', delta_t
         write (iun, '(A,10E20.12)') 'h&K-XC', j%J_xc(:)
         write (iun, '(A,10E20.12)') 'h&K-H', j%J_hartree(:)
         write (iun, '(A,3E20.12)') 'h&K-K', j%J_kohn(1:3)
         write (iun, '(A,3E20.12)') 'h&K-K_a', j%J_kohn_a(1:3)
         write (iun, '(A,3E20.12)') 'h&K-K_b', j%J_kohn_b(1:3)
         write (iun, '(A,3E20.12)') 'h&K-ELE', j%J_electron(1:3)
         write (iun, '(A,3E20.12)') 'ionic:', j%i_current(:)
         write (iun, '(A,3E20.12)') 'ionic_a:', j%i_current_a(:)
         write (iun, '(A,3E20.12)') 'ionic_b:', j%i_current_b(:)
         write (iun, '(A,3E20.12)') 'ionic_c:', j%i_current_c(:)
         write (iun, '(A,3E20.12)') 'ionic_d:', j%i_current_d(:)
         write (iun, '(A,3E20.12)') 'ionic_e:', j%i_current_e(:)
         write (iun, '(A,3E20.12)') 'zero:', j%z_current(:)
         J_tot = j%J_xc + j%J_hartree + j%J_kohn + j%i_current + j%z_current
         write (iun, '(A,3E20.12)') 'total: ', J_tot
         write (*, '(A,3E20.12)') 'total energy current: ', J_tot
         close (iun)
         !WARNING: if you modify the following lines
         !remember to modify the set_first_step_restart() subroutine,
         !so we can read the file that here we are writing in the correct way
         open (iun, file=trim(file_output)//'.dat', position='append')
         if (.not. file_exists) then
            write (iun, '(A)') '#units:'
            write (iun, '(A)') '# Ry = 2.1799e-18J = 13.606eV'
            write (iun, '(A)') '# a0 = 5.2918e-11m,'
            write (iun, '(A)') '# tau = 4.8378e-17s'
            write (iun, '(A)') '# TIME: same as input trajectory'
            write (iun, '(A)') '# J: Ry*a0/tau'
            write (iun, '(A)') '# J_el: number * a0/tau'
            write (iun, '(A)') '# J_cm: number * a0/tau'
            write (iun, '(A)', advance='no') '# STEP TIME J[1] J[2] J[3] J_el[1] J_el[2] J_el[3] '
            do itype = 1, nsp
               do icoord = 1, 3
                  write (iun, '(A,I1,A,I1,A)', advance='no') 'J_cm', itype, '[', icoord, '] '
               end do
            end do
            write (iun, '(A)') ''
         end if
         write (iun, '(1I7,1E14.6,3E20.12,3E20.12)', advance='no') step, time, &
            J_tot, j%J_electron(1:3)
         do itype = 1, nsp
            write (iun, '(3E20.12)', advance='no') real(na(itype), dp)*alat*j%v_cm(:, itype)
            write (*, '(A,1I3,A,3E20.12)') 'center of mass velocity of type ', itype, ': ', alat*j%v_cm(:, itype)
         end do
         write (iun, '(A)') ''
         close (iun)
         call online_average_do(ave_cur, J_tot)
         if (write_ave) then
            INQUIRE (FILE=trim(file_output)//'.stat', EXIST=file_exists)
            open (iun, file=trim(file_output)//'.stat', position='append')
            if (.not. file_exists) then
               write (iun, '(A)') '#units:'
               write (iun, '(A)') '# Ry = 2.1799e-18J = 13.606eV'
               write (iun, '(A)') '# a0 = 5.2918e-11m,'
               write (iun, '(A)') '# tau = 4.8378e-17s'
               write (iun, '(A)') '# TIME: same as input trajectory'
               write (iun, '(A)') '# J: Ry*a0/tau (average of the results)'
               write (iun, '(A)') '# sigma_J: Ry*a0/tau (standard deviation of the results)'
               write (iun, '(A)') '# STEP TIME J[1] J[2] J[3] sigma_J[1] sigma_J[2] sigma_J[3] '
            end if
            write (iun, '(1I7,1E14.6)', advance='no') step, time
            call online_average_print(ave_cur, iun)
            write (iun, '(A)') ''
            close (iun)
         end if

      end if

   end subroutine

   subroutine read_all_currents_namelists(iunit, &
                                          delta_t, &
                                          file_output, trajdir, vel_input_units, &
                                          eta, n_max, first_step, last_step, &
                                          ethr_small_step, ethr_big_step, &
                                          restart, subtract_cm_vel, step_mul, &
                                          step_rem, ec_test, add_i_current_b, &
                                          save_dvpsi, re_init_wfc_1, re_init_wfc_2, &
                                          re_init_wfc_3, three_point_derivative, &
                                          n_repeat_every_step, hpsi_test, &
                                          n_workers, worker_id, &
                                          continue_not_converged)
      use io_global, ONLY: stdout, ionode, ionode_id
      implicit none
      integer, intent(in) :: iunit
      real(kind=dp), intent(inout) :: eta
      logical, intent(inout) :: save_dvpsi
      CHARACTER(len=256), intent(inout) :: file_output, trajdir
      real(kind=DP), intent(inout) ::delta_t, ethr_small_step, ethr_big_step
      integer, intent(inout) :: first_step, last_step, step_mul, step_rem, n_repeat_every_step, &
                                n_workers, worker_id
      logical, intent(inout) :: restart
      logical, intent(inout) :: subtract_cm_vel
      logical, intent(inout) :: re_init_wfc_1, re_init_wfc_2
      logical, intent(inout) :: re_init_wfc_3
      logical, intent(inout) :: three_point_derivative
      logical, intent(inout) :: add_i_current_b
      character(len=256), intent(inout) :: vel_input_units
      logical, intent(inout) :: ec_test, hpsi_test ! activates tests for debugging purposes
      integer, intent(out) :: n_max
      logical, intent(inout) :: continue_not_converged
      integer :: ios

      NAMELIST /energy_current/ delta_t, &
         file_output, trajdir, vel_input_units, &
         eta, n_max, first_step, last_step, &
         ethr_small_step, ethr_big_step, &
         restart, subtract_cm_vel, step_mul, &
         step_rem, ec_test, add_i_current_b, &
         save_dvpsi, re_init_wfc_1, re_init_wfc_2, &
         re_init_wfc_3, three_point_derivative, &
         n_repeat_every_step, hpsi_test, &
         n_workers, worker_id, continue_not_converged
      !
      !   set default values for variables in namelist
      !
      delta_t = 1.d0
      n_max = 5 ! number of periodic cells in each direction used to sum stuff in zero current
      eta = 1.d0 ! ewald sum convergence parameter
      file_output = "current_hz"
      trajdir = ''
      ethr_small_step = 1.d-7
      ethr_big_step = 1.d-3
      first_step = 0
      last_step = 0
      step_mul = 1
      step_rem = 0
      restart = .false.
      subtract_cm_vel = .false.
      ec_test = .false.
      add_i_current_b = .false.
      save_dvpsi = .true.
      re_init_wfc_1 = .false.
      re_init_wfc_2 = .false.
      re_init_wfc_3 = .false.
      three_point_derivative = .true.
      vel_input_units = 'PW'
      n_repeat_every_step = 1
      hpsi_test = .false.
      n_workers = 0
      worker_id = 0
      continue_not_converged = .false.
      READ (iunit, energy_current, IOSTAT=ios)
      IF (ios /= 0) CALL errore('main', 'reading energy_current namelist', ABS(ios))

   end subroutine

   subroutine bcast_all_current_namelist( &
      delta_t, &
      file_output, trajdir, vel_input_units, &
      eta, n_max, first_step, last_step, &
      ethr_small_step, ethr_big_step, &
      restart, subtract_cm_vel, step_mul, &
      step_rem, ec_test, add_i_current_b, &
      save_dvpsi, re_init_wfc_1, re_init_wfc_2, &
      re_init_wfc_3, three_point_derivative, &
      n_repeat_every_step, hpsi_test, &
      n_workers, worker_id, &
      continue_not_converged)
      use io_global, ONLY: stdout, ionode, ionode_id
      use mp_world, ONLY: mpime, world_comm
      use mp, ONLY: mp_bcast
      implicit none
      real(kind=dp), intent(inout) :: eta
      logical, intent(inout) :: save_dvpsi
      CHARACTER(len=256), intent(inout) :: file_output, trajdir
      real(kind=DP), intent(inout) ::delta_t, ethr_small_step, ethr_big_step
      integer, intent(inout) :: first_step, last_step, step_mul, step_rem, n_repeat_every_step, &
                                n_workers, worker_id
      logical, intent(inout) :: restart
      logical, intent(inout) :: subtract_cm_vel
      logical, intent(inout) :: re_init_wfc_1, re_init_wfc_2
      logical, intent(inout) :: re_init_wfc_3
      logical, intent(inout) :: three_point_derivative
      logical, intent(inout) :: add_i_current_b
      character(len=256), intent(inout) :: vel_input_units
      logical, intent(inout) :: ec_test, hpsi_test
      integer, intent(out) :: n_max
      logical, intent(inout) :: continue_not_converged
      CALL mp_bcast(trajdir, ionode_id, world_comm)
      CALL mp_bcast(delta_t, ionode_id, world_comm)
      CALL mp_bcast(eta, ionode_id, world_comm)
      CALL mp_bcast(restart, ionode_id, world_comm)
      CALL mp_bcast(subtract_cm_vel, ionode_id, world_comm)
      CALL mp_bcast(first_step, ionode_id, world_comm)
      CALL mp_bcast(last_step, ionode_id, world_comm)
      CALL mp_bcast(ethr_small_step, ionode_id, world_comm)
      CALL mp_bcast(ethr_big_step, ionode_id, world_comm)
      CALL mp_bcast(n_max, ionode_id, world_comm)
      CALL mp_bcast(save_dvpsi, ionode_id, world_comm)
      CALL mp_bcast(file_output, ionode_id, world_comm)
      CALL mp_bcast(step_mul, ionode_id, world_comm)
      CALL mp_bcast(step_rem, ionode_id, world_comm)
      CALL mp_bcast(ec_test, ionode_id, world_comm)
      CALL mp_bcast(add_i_current_b, ionode_id, world_comm)
      CALL mp_bcast(re_init_wfc_1, ionode_id, world_comm)
      CALL mp_bcast(re_init_wfc_2, ionode_id, world_comm)
      CALL mp_bcast(re_init_wfc_3, ionode_id, world_comm)
      CALL mp_bcast(three_point_derivative, ionode_id, world_comm)
      CALL mp_bcast(n_repeat_every_step, ionode_id, world_comm)
      CALL mp_bcast(hpsi_test, ionode_id, world_comm)
      CALL mp_bcast(n_workers, ionode_id, world_comm)
      CALL mp_bcast(worker_id, ionode_id, world_comm)
      CALL mp_bcast(vel_input_units, ionode_id, world_comm)
      CALL mp_bcast(continue_not_converged, ionode_id, world_comm)
   end subroutine

   subroutine set_first_step_restart(restart, file_output, first_step)
      use io_global, ONLY: ionode, ionode_id
      use mp_world, ONLY: mpime, world_comm
      use mp, ONLY: mp_bcast
      use ions_base, ONLY: nsp
      implicit none
      logical, intent(in) :: restart
      CHARACTER(len=256), intent(in) :: file_output
      integer, intent(inout) :: first_step
      integer :: iun, step, ios, i, step1
      integer, external :: find_free_unit
      real(dp) :: time, J(3), time1, J1(3), Jdummy(3)

      if (.not. restart) return
      if (ionode) then
         iun = find_free_unit()
         ios = 0
         step1 = -1
         open (iun, file=trim(file_output)//'.dat')
         write (*, *) 'reading file ', trim(file_output)//'.dat'
         do while (ios == 0)
            read (iun, *, iostat=ios) step, time, J, (Jdummy, i=1, nsp + 1)
            if (ios == 0) then
               write (*, *) 'found: ', step, time, J
               step1 = step
               time1 = time
               J1 = J
            end if
         enddo
         if (step1 /= -1) then
            first_step = step1 + 1
            write (*, *) 'RESTARTING AFTER STEP ', step1, time1, J1
         else
            write (*, *) 'NOTHING TO RESTART'
         endif
         close (iun)
      endif
      !broadcast first_step
      CALL mp_bcast(first_step, ionode_id, world_comm)

   end subroutine

   subroutine check_input()
      use input_parameters, only: rd_pos, tapos, rd_vel, tavel, atomic_positions, ion_velocities
      use ions_base, ONLY: tau, tau_format, nat
      implicit none
      if (first_step == 0) then
         if (.not. tavel) &
            call errore('read_vel', 'error: must provide velocities in input', 1)
         if (ion_velocities /= 'from_input') &
            call errore('read_vel', 'error: atomic_velocities must be "from_input"', 1)
      else
         if (tavel) &
            write (*, *) 'WARNING: VELOCITIES FROM INPUT FILE WILL BE IGNORED'
      end if
   end subroutine

   subroutine run_electrons(exit_status, continue_not_converged)
      USE control_flags, ONLY: conv_elec, gamma_only, ethr, lscf, treinit_gvecs
      USE check_stop, ONLY: check_stop_init, check_stop_now
      USE qexsd_module, ONLY: qexsd_set_status
      implicit none
      INTEGER, INTENT(OUT) :: exit_status
      logical, intent(in) :: continue_not_converged
      exit_status = 0
      call start_clock('PWSCF')
      IF (.NOT. lscf) THEN
         CALL non_scf()
      ELSE
         CALL electrons()
      END IF
      call stop_clock('PWSCF')
      IF (.NOT. conv_elec) exit_status = 2
      IF (check_stop_now() .OR. & 
             ( (.NOT. conv_elec ) .and. ( .not. continue_not_converged)) ) THEN
         IF (check_stop_now()) exit_status = 255
         CALL qexsd_set_status(exit_status)
         CALL punch('config')
         RETURN
      ENDIF
   end subroutine

   subroutine cm_vel(v_cm, vel_cm)
      !calculates center of mass velocities for each atomic type and eventually subtract it from the velocity of each atom
      use kinds, only: dp
      use ions_base, ONLY: nat, nsp, ityp, na
      use dynamics_module, only: vel
      implicit none
      real(dp), allocatable, intent(inout) :: v_cm(:, :)
      real(dp), intent(out), optional :: vel_cm(:, :)
      integer :: iatom, itype
      integer, allocatable :: counter(:)
      real(dp) :: delta(3), mean(3)
      if (.not. allocated(v_cm)) &
         allocate (v_cm(3, nsp))
      allocate (counter(nsp))
      counter = 0
      v_cm = 0.0_dp

      do iatom = 1, nat
         itype = ityp(iatom)
         counter(itype) = counter(itype) + 1
         delta = (vel(:, iatom) - v_cm(:, itype))/real(counter(itype), dp)
         v_cm(:, itype) = v_cm(:, itype) + delta
      end do
      if (na(1) .eq. 0) &
         na(1:nsp) = counter
      if (present(vel_cm)) then
         do iatom = 1, nat
            itype = ityp(iatom)
            vel_cm(:, iatom) = vel_cm(:, iatom) - v_cm(:, itype)
         end do
      end if

      deallocate (counter)
   end subroutine

   subroutine prepare_next_step(ipm, delta_t, ethr_small_step, three_point_derivative)
      ! advance/go back by dt according to ipm = -1 1 and if dt/2 has to be used (three point derivative)
      USE extrapolation, ONLY: update_pot
      USE control_flags, ONLY: ethr
      use ions_base, ONLY: tau, tau_format, nat
      use cell_base, only: alat
      use dynamics_module, only: vel
      use io_global, ONLY: ionode, ionode_id
      USE mp_world, ONLY: world_comm
      use mp, ONLY: mp_bcast, mp_barrier
      use wavefunctions, only: evc
      implicit none
      real(kind=DP), intent(in) ::delta_t, ethr_small_step
      logical, intent(in) :: three_point_derivative
      integer, intent(in) :: ipm

      if (ipm /= 0 .and. ipm /= -1 .and. ipm /= 1) &
         call errore('prepare_next_step', 'error: unknown timestep type')
      !broadcast
      CALL mp_bcast(tau, ionode_id, world_comm)
      CALL mp_bcast(vel, ionode_id, world_comm)
      !set new positions
      if (three_point_derivative) then
         tau = tau + delta_t*vel*real(ipm, dp)/2.0_dp
      else
         tau = tau + delta_t*vel*real(ipm, dp)
      end if
      call mp_barrier(world_comm)
      call update_pot()
      call hinit1()

   end subroutine

   function read_next_step(t, first_step, last_step, step_mul, step_rem) result(res)
      ! read next trajectory step. if no next step exists, return false
      USE extrapolation, ONLY: update_pot
      use cpv_traj, only: cpv_trajectory, cpv_trajectory_initialize, cpv_trajectory_deallocate, &
                          cpv_trajectory_read_step, cpv_trajectory_get_step
      use traj_object, only: timestep ! type for timestep data
      use kinds, only: dp
      use ions_base, ONLY: tau, tau_format, nat
      use cell_base, only: alat
      use dynamics_module, only: vel
      use io_global, ONLY: ionode, ionode_id
      USE mp_world, ONLY: world_comm
      use mp, ONLY: mp_bcast, mp_barrier
      implicit none
      integer, intent(in) :: first_step, last_step, step_mul, step_rem
      type(cpv_trajectory), intent(inout) :: t
      type(timestep) :: ts
      logical :: res
      integer, save :: nstep = 0
      integer, save :: step_idx = 0
      logical, save :: first = .true.
      if (ionode) then
         do while (cpv_trajectory_read_step(t))
            step_idx = step_idx + 1
            call cpv_trajectory_get_step(t, step_idx, ts)
            !if we are done exit the loop
            nstep = ts%nstep
            if (ts%nstep > last_step .and. last_step > 0) then
               write (*, *) 'STEP ', ts%nstep, ts%tps, ' > ', last_step
               exit
            end if
            !if needed go on in the reading of trajectory
            if (ts%nstep < first_step .or. mod(ts%nstep, step_mul) /= step_rem) then
               write (*, *) 'SKIPPED STEP ', ts%nstep, ts%tps
               cycle
            end if
            write (*, *) 'STEP ', ts%nstep, ts%tps
            vel = ts%vel
            tau = ts%tau
            CALL convert_tau(tau_format, nat, tau)
            call convert_tau(tau_format, nat, vel)
            res = .true.
            goto 200 ! exit the loop skipping 'finish': we will do an other calculation
         enddo
         !else
         write (*, *) 'Finished reading trajectory ', trim(t%fname), ' at step ', nstep
         res = .false.
200      continue
      endif
      CALL mp_bcast(res, ionode_id, world_comm)
      if (res) then
         CALL mp_bcast(vel, ionode_id, world_comm)
         CALL mp_bcast(tau, ionode_id, world_comm)
         if (first_step > 0 .and. first) then
            first = .false.
            return
         end if
         first = .false.

      end if

   end function

end program all_currents
