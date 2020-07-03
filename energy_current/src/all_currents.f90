program all_currents
   use hartree_mod, only: evc_uno, evc_due, trajdir, first_step, dvpsi_save
   USE environment, ONLY: environment_start, environment_end
   use io_global, ONLY: ionode
   use wavefunctions, only: evc
   use kinds, only: dp
   !trajectory reading stuff
   use ions_base, only: nat
   use cpv_traj, only: cpv_trajectory, &
                       cpv_trajectory_initialize, cpv_trajectory_deallocate

!from ../PW/src/pwscf.f90
   USE mp_global, ONLY: mp_startup
   USE mp_world, ONLY: world_comm
   use mp, ONLY: mp_bcast, mp_barrier
   USE mp_pools, ONLY: intra_pool_comm
   USE mp_bands, ONLY: intra_bgrp_comm, inter_bgrp_comm
   USE read_input, ONLY: read_input_file
   USE command_line_options, ONLY: input_file_, command_line, ndiag_, nimage_
   USE check_stop, ONLY: check_stop_init
!from ../Modules/read_input.f90
   USE read_namelists_module, ONLY: read_namelists
   USE read_cards_module, ONLY: read_cards

   implicit none
   integer :: exit_status, ios
   type(cpv_trajectory) :: traj

!from ../PW/src/pwscf.f90
   include 'laxlib.fh'

!from ../PW/src/pwscf.f90
   CALL mp_startup()
   CALL laxlib_start(ndiag_, world_comm, intra_bgrp_comm, &
                     do_distr_diag_inside_bgrp_=.TRUE.)
   CALL set_mpi_comm_4_solvers(intra_pool_comm, intra_bgrp_comm, &
                               inter_bgrp_comm)
   CALL environment_start('PWSCF')

   IF (ionode) THEN
      CALL input_from_file()
      ! all_currents input
      call read_all_currents_namelists(5)
   endif
   ! PW input
   call read_namelists('PW', 5)
   call read_cards('PW', 5)

   call check_input()

   call mp_barrier(intra_pool_comm)
   call bcast_all_current_namelist()
   call set_first_step_restart()
   call iosys()    ! ../PW/src/input.f90    save in internal variables
   call check_stop_init() ! ../PW/src/input.f90

   !eventually read new positions and velocities from trajectory
   if (ionode) then
      !initialize trajectory reading
      call cpv_trajectory_initialize(traj, trajdir, nat, 1.0_dp, 1.0_dp, 1.0_dp, ios=ios, circular=.true.)
      if (ios == 0) then
         write (*, *) 'After first step from input file, I will read from the CPV trajectory ', trim(trajdir)
      else
         write (*, *) 'Cannot open trajectory file', trim(trajdir), '. I can calculate only a single step from input file'
      endif
   endif
   if (first_step > 0) then
      if (ionode) &
         write (*, *) 'SKIPPING STEP FROM INPUT FILE'
      if (.not. read_next_step(traj)) then
         if (ionode) &
            write (*, *) 'NOTHING TO DO IN TRAJECTORY FILE'
         goto 100 ! skip everything and exit
      endif
   end if

   call setup()    ! ../PW/src/setup.f90    setup the calculation
   call init_run() ! ../PW/src/init_run.f90 allocate stuff
   ! now scf is ready to start, but I first initialize energy current stuff
   call allocate_zero() ! only once per all trajectory
   call init_zero() ! only once per all trajectory
   call setup_nbnd_occ() ! only once per all trajectory

   do
      call run_pwscf(exit_status)
      if (exit_status /= 0) exit

      call prepare_next_step() ! this stores value of evc and setup tau and ion_vel

      call run_pwscf(exit_status)
      if (exit_status /= 0) goto 100 !shutdown everything and exit
      if (allocated(evc_uno)) then
         evc_uno = evc
      else
         allocate (evc_uno, source=evc)
      end if

      !calculate energy current
      call routine_hartree()
      call routine_zero()
      call write_results(traj)
      !read new velocities and positions and continue, or exit the loop
      if (.not. read_next_step(traj)) exit
   end do

   ! shutdown stuff
100 call laxlib_end()
   call cpv_trajectory_deallocate(traj)
   call deallocate_zero()
   if (allocated(evc_uno)) deallocate (evc_uno)
   if (allocated(evc_due)) deallocate (evc_due)
   if (allocated(dvpsi_save)) deallocate(dvpsi_save)
   call stop_run(exit_status)
   call do_stop(exit_status)
   stop

contains

   subroutine write_results(traj)
      use kinds, only: dp
      use ions_base, ONLY: nsp
      use hartree_mod
      use zero_mod
      use io_global, ONLY: ionode
      use cpv_traj, only: cpv_trajectory, cpv_trajectory_get_last_step
      use traj_object, only: timestep
      implicit none
      type(cpv_trajectory), intent(in)  :: traj
      type(timestep) :: ts
      integer :: iun, step, itype
      integer, external :: find_free_unit
      real(dp) :: time

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
         write (iun, *) 'Passo: ', step
         write (iun, '(A,10E20.12)') 'h&K-XC', J_xc(:)
         write (iun, '(A,10E20.12)') 'h&K-H', J_hartree(:)
         write (iun, '(A,1F15.7,9E20.12)') 'h&K-K', delta_t, J_kohn(1:3), J_kohn_a(1:3), J_kohn_b(1:3)
         write (iun, '(A,3E20.12)') 'h&K-ELE', J_electron(1:3)
         write (iun, '(A,3E20.12)') 'ionic:', i_current(:)
         write (iun, '(A,3E20.12)') 'ionic_a:', i_current_a(:)
         write (iun, '(A,3E20.12)') 'ionic_b:', i_current_b(:)
         write (iun, '(A,3E20.12)') 'ionic_c:', i_current_c(:)
         write (iun, '(A,3E20.12)') 'ionic_d:', i_current_d(:)
         write (iun, '(A,3E20.12)') 'ionic_e:', i_current_e(:)
         write (iun, '(A,3E20.12)') 'zero:', z_current(:)
         write (iun, '(A,3E20.12)') 'total: ', J_xc + J_hartree + J_kohn + i_current + z_current
         write (*, '(A,3E20.12)') 'total energy current: ', J_xc + J_hartree + J_kohn + i_current + z_current
         close (iun)
         !WARNING: if you modify the following lines
         !remember to modify the set_first_step_restart() subroutine,
         !so we can read the file that here we are writing in the correct way
         open (iun, file=trim(file_output)//'.dat', position='append')
         write (iun, '(1I7,1E14.6,3E20.12,3E20.12)', advance='no') step, time, &
            J_xc + J_hartree + J_kohn + i_current + z_current, J_electron(1:3)
         do itype = 1, nsp
            write (iun, '(3E20.12)', advance='no') v_cm(:, itype)
            write (*, '(A,1I3,A,3E20.12)') 'center of mass velocity of type ', itype, ': ', v_cm(:, itype)
         end do
         write (iun, '(A)') ''
         close (iun)
      end if

   end subroutine

   subroutine read_all_currents_namelists(iunit)
      use zero_mod
      use hartree_mod
      use io_global, ONLY: stdout, ionode, ionode_id
      implicit none
      integer, intent(in) :: iunit
      integer :: ios
      CHARACTER(LEN=256), EXTERNAL :: trimcheck

      NAMELIST /energy_current/ delta_t, &
         file_output, trajdir, vel_input_units, &
         eta, n_max, first_step, last_step, &
         ethr_small_step, ethr_big_step, &
         restart, subtract_cm_vel, step_mul, &
         step_rem, ec_test, add_i_current_b
      !
      !   set default values for variables in namelist
      !
      delta_t = 1.d0
      n_max = 5 ! number of periodic cells in each direction used to sum stuff in zero current
      eta = 1.d0 ! ewald sum convergence parameter
!      init_linear = "nothing" ! 'scratch' or 'restart'. If 'scratch', saves a restart file in project routine. If 'restart', it starts from the saved restart file, and then save again it.
      file_output = "current_hz"
      ethr_small_step = 1.d-7
      ethr_big_step = 1.d-3
      first_step = 0
      last_step = 0
      restart = .false.
      subtract_cm_vel = .false.
      step_mul = 1
      step_rem = 0
      ec_test = .false.
      add_i_current_b = .false.
      READ (iunit, energy_current, IOSTAT=ios)
      IF (ios /= 0) CALL errore('main', 'reading energy_current namelist', ABS(ios))

   end subroutine

   subroutine bcast_all_current_namelist()
      use zero_mod
      use hartree_mod
      use io_global, ONLY: stdout, ionode, ionode_id
      use mp_world, ONLY: mpime, world_comm
      use mp, ONLY: mp_bcast
      implicit none
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
 !     CALL mp_bcast(init_linear, ionode_id, world_comm)
      CALL mp_bcast(file_output, ionode_id, world_comm)
      CALL mp_bcast(step_mul, ionode_id, world_comm)
      CALL mp_bcast(step_rem, ionode_id, world_comm)
      CALL mp_bcast(ec_test, ionode_id, world_comm)
      CALL mp_bcast(add_i_current_b, ionode_id, world_comm)
   end subroutine

   subroutine set_first_step_restart()
      use hartree_mod, only: restart, file_output, first_step
      use io_global, ONLY: ionode, ionode_id
      use mp_world, ONLY: mpime, world_comm
      use mp, ONLY: mp_bcast
      use ions_base, ONLY: nsp
      implicit none
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

   subroutine run_pwscf(exit_status)
      USE control_flags, ONLY: conv_elec, gamma_only, ethr, lscf, treinit_gvecs
      USE check_stop, ONLY: check_stop_init, check_stop_now
      USE qexsd_module, ONLY: qexsd_set_status
      implicit none
      INTEGER, INTENT(OUT) :: exit_status
      exit_status = 0
      IF (.NOT. lscf) THEN
         CALL non_scf()
      ELSE
         CALL electrons()
      END IF
      IF (check_stop_now() .OR. .NOT. conv_elec) THEN
         IF (check_stop_now()) exit_status = 255
         IF (.NOT. conv_elec) exit_status = 2
         CALL qexsd_set_status(exit_status)
         CALL punch('config')
         RETURN
      ENDIF
   end subroutine

   subroutine cm_vel(vel_cm)
      use kinds, only: dp
      use ions_base, ONLY: nat, nsp, ityp
      use dynamics_module, only: vel
      use hartree_mod, only: v_cm
      implicit none
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
      if (present(vel_cm)) then
         do iatom = 1, nat
            itype = ityp(iatom)
            vel_cm(:, iatom) = vel_cm(:, iatom) - v_cm(:, itype)
         end do
      end if

      deallocate (counter)
   end subroutine

   subroutine prepare_next_step()
      USE extrapolation, ONLY: update_pot
      USE control_flags, ONLY: ethr
      use ions_base, ONLY: tau, tau_format, nat
      use cell_base, only: alat
      use dynamics_module, only: vel
      use io_global, ONLY: ionode, ionode_id
      USE mp_world, ONLY: world_comm
      use mp, ONLY: mp_bcast, mp_barrier
      use hartree_mod, only: evc_due, delta_t, ethr_small_step, &
                             subtract_cm_vel
      use zero_mod, only: vel_input_units, ion_vel
      use wavefunctions, only: evc
      implicit none
      !save old evc
      if (allocated(evc_due)) then
         evc_due = evc
      else
         allocate (evc_due, source=evc)
      end if
      !set new positions
      if (ionode) then
         if (vel_input_units == 'CP') then ! atomic units of cp are different
            vel = 2.d0*vel
         else if (vel_input_units == 'PW') then
            !do nothing
         else
            call errore('read_vel', 'error: unknown vel_input_units', 1)
         endif
         if (subtract_cm_vel) then
            !calculate center of mass velocity for each atomic species and subtract it
            call cm_vel(vel)
         else
            call cm_vel()
         end if
      endif
      !broadcast
      CALL mp_bcast(tau, ionode_id, world_comm)
      CALL mp_bcast(vel, ionode_id, world_comm)
      if (.not. allocated(ion_vel)) then
         allocate (ion_vel, source=vel)
      else
         ion_vel = vel
      endif
      call convert_tau(tau_format, nat, vel)
      tau = tau + delta_t*vel
      call mp_barrier(world_comm)
      call update_pot()
      call hinit1()
      ethr = ethr_small_step
   end subroutine

   function read_next_step(t) result(res)
      USE extrapolation, ONLY: update_pot
      USE control_flags, ONLY: ethr
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
      use zero_mod, only: vel_input_units, ion_vel
      use hartree_mod, only: first_step, last_step, ethr_big_step, &
                             step_mul, step_rem
      implicit none
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
         call update_pot()
         call hinit1()
         ethr = ethr_big_step

      end if

   end function

end program all_currents
