MODULE oscdft_functions
#if defined (__OSCDFT)
   USE kinds,                    ONLY : DP
   USE io_global,                ONLY : ionode, ionode_id, stdout
   USE mp,                       ONLY : mp_bcast
   USE mp_images,                ONLY : intra_image_comm
   USE oscdft_enums
   USE oscdft_occupations,       ONLY : oscdft_new_ns,&
                                        oscdft_get_occupation_numbers
   USE oscdft_context,           ONLY : oscdft_context_type, oscdft_ns_type
   USE oscdft_input,             ONLY : oscdft_input_type
   USE oscdft_indices,           ONLY : oscdft_indices_type,&
                                        oscdft_orbital_indices_type,&
                                        oscdft_constr_indices_type
   USE oscdft_wavefunction,      ONLY : oscdft_wavefunction_type,&
                                        check_bec_type_unallocated
   USE oscdft_wavefunction_subs, ONLY : oscdft_get_buffer

   PRIVATE
   PUBLIC oscdft_run_pwscf,&
          oscdft_electrons,&
          oscdft_electrons_restart,&
          oscdft_scf_energy,&
          oscdft_print_energies,&
          oscdft_close_files,&
          oscdft_h_diag,&
          oscdft_h_psi,&
          oscdft_mix_rho,&
          oscdft_sum_band,&
          oscdft_print_ns,&
          oscdft_init_constraint,&
          oscdft_constrain_ns,&
          oscdft_v_constraint,&
          oscdft_ns_set
   CONTAINS
      SUBROUTINE oscdft_get_deriv(ctx)
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         INTEGER                                          :: iconstr, ioscdft
         TYPE(oscdft_input_type),   POINTER               :: inp
         TYPE(oscdft_indices_type), POINTER               :: idx
         TYPE(oscdft_ns_type),      POINTER               :: nst
         LOGICAL :: have_clamped_gradient

         inp => ctx%inp
         idx => ctx%idx
         nst => ctx%nst

         DO iconstr=1,idx%nconstr
            ioscdft = idx%iconstr2ioscdft(iconstr)
            nst%gradient(iconstr) = nst%occup_numbers(iconstr) - inp%target_occup(ioscdft)
         END DO

         WRITE(stdout, 100)
         WRITE(stdout, 200) nst%gradient

         have_clamped_gradient = .false.
         DO iconstr=1,idx%nconstr
            ioscdft = idx%iconstr2ioscdft(iconstr)
            IF (inp%constraint_applied(ioscdft) == CONSTR_GE2) THEN
               have_clamped_gradient = .true.
               nst%gradient(iconstr) = MIN(nst%gradient(iconstr), 0.D0)
            ELSE IF (inp%constraint_applied(ioscdft) == CONSTR_LE2) THEN
               have_clamped_gradient = .true.
               nst%gradient(iconstr) = MAX(nst%gradient(iconstr), 0.D0)
            ELSE IF (inp%constraint_applied(ioscdft) == CONSTR_GE3) THEN
               have_clamped_gradient = .true.
               nst%gradient(iconstr) = MIN(nst%gradient(iconstr),  0.9D0 * ctx%conv_thr)
            ELSE IF (inp%constraint_applied(ioscdft) == CONSTR_LE3) THEN
               have_clamped_gradient = .true.
               nst%gradient(iconstr) = MAX(nst%gradient(iconstr), -0.9D0 * ctx%conv_thr)
            END IF
         END DO

         IF (have_clamped_gradient) THEN
            WRITE(stdout, 101)
            WRITE(stdout, 200) nst%gradient
         END IF

         100 FORMAT("OSCDFT: gradient")
         101 FORMAT("OSCDFT: clamped gradient")
         200 FORMAT("OSCDFT:   ", *(ES14.7, :, " "))
      END SUBROUTINE oscdft_get_deriv

      SUBROUTINE oscdft_optimize_multipliers(ctx)
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),   POINTER               :: inp
         TYPE(oscdft_indices_type), POINTER               :: idx
         TYPE(oscdft_ns_type),      POINTER               :: nst
         REAL(DP) :: delta_mul(ctx%idx%nconstr),&
                     mask(ctx%idx%nconstr)

         inp => ctx%inp
         idx => ctx%idx
         nst => ctx%nst

         mask = MERGE(1.0, 0.0, inp%start_index(idx%iconstr2ioscdft) <= ctx%global_start_index)
         WRITE(stdout, 300) "global_start_index mask"
         WRITE(stdout, 400) mask

         IF (idx%nconstr > 0) THEN
            SELECT CASE (inp%optimization_method)
               CASE (OPT_GRADIENT_DESCENT)
                  WRITE(stdout, 102) "gamma", inp%min_gamma_n
                  delta_mul = inp%min_gamma_n * nst%gradient
                  ctx%multipliers = ctx%multipliers + delta_mul * mask
               CASE (OPT_GRADIENT_DESCENT2)
                  WRITE(stdout, 102) "min_gamma_n", inp%min_gamma_n
                  WRITE(stdout, 300) "gamma_val"
                  WRITE(stdout, 400) inp%gamma_val(idx%iconstr2ioscdft)
                  WRITE(stdout, 300) "gamma"
                  WRITE(stdout, 400) inp%min_gamma_n * inp%gamma_val(idx%iconstr2ioscdft)
                  delta_mul = inp%min_gamma_n *&
                              inp%gamma_val(idx%iconstr2ioscdft) *&
                              nst%gradient
                  ctx%multipliers = ctx%multipliers + delta_mul * mask
               CASE DEFAULT
                  CALL errore("oscdft_optimize_multipliers",&
                              "invalid optimization method",&
                              ABS(inp%optimization_method))
            END SELECT
         END IF
         102 FORMAT("OSCDFT: ", A, ": ", ES14.7)
         300 FORMAT("OSCDFT: ", A)
         400 FORMAT("OSCDFT:   ", *(ES14.7, :, " "))
      END SUBROUTINE oscdft_optimize_multipliers

      SUBROUTINE oscdft_clamp_multipliers(ctx)
         IMPLICIT NONE
         TYPE(oscdft_context_type), INTENT(INOUT) :: ctx

         IF (ctx%inp%has_max_multiplier) ctx%multipliers = MIN(ctx%multipliers, ctx%inp%max_multiplier)
         IF (ctx%inp%has_min_multiplier) ctx%multipliers = MAX(ctx%multipliers, ctx%inp%min_multiplier)
      END SUBROUTINE oscdft_clamp_multipliers

      FUNCTION oscdft_array_to_scalar(ctx, array) result(scalar)
         IMPLICIT NONE
         TYPE(oscdft_context_type), INTENT(IN) :: ctx
         REAL(DP),                  INTENT(IN) :: array(:)
         REAL(DP)                              :: scalar
         INTEGER                               :: length, idx
         REAL(DP)                              :: mask(ctx%idx%nconstr)

         length = ctx%idx%nconstr
         mask = MERGE(1.0, 0.0, ctx%inp%start_index(ctx%idx%iconstr2ioscdft) <= ctx%global_start_index)

         SELECT CASE (ctx%inp%array_convergence_func)
            CASE (CONV_FUNC_MAXVAL)
               scalar = MAXVAL(ABS(array * mask))
            CASE (CONV_FUNC_NORM)
               scalar = 0.D0
               DO idx=1,length
                  scalar = scalar + array(idx) * array(idx) * mask(idx)
               ENDDO
               scalar = SQRT(scalar)
            CASE (CONV_FUNC_NORM_AVERAGE)
               scalar = 0.D0
               DO idx=1,length
                  scalar = scalar + array(idx) * array(idx) * mask(idx)
               ENDDO
               scalar = SQRT(scalar / length)
            CASE DEFAULT
               CALL errore("oscdft_array_to_scalar", "invalid array_convergence_func", 1)
         END SELECT
      END FUNCTION oscdft_array_to_scalar

      FUNCTION oscdft_multipliers_convergence_test(ctx, threshold, threshold_label) RESULT(conv)
         IMPLICIT NONE
         TYPE(oscdft_context_type), INTENT(IN), TARGET    :: ctx
         REAL(DP), INTENT(IN)                             :: threshold
         CHARACTER(LEN=*), INTENT(IN)                     :: threshold_label
         LOGICAL                                          :: conv
         REAL(DP)                                         :: test
         REAL(DP)                                         :: temp(ctx%idx%nconstr)
         TYPE(oscdft_input_type),   POINTER               :: inp
         TYPE(oscdft_indices_type), POINTER               :: idx
         TYPE(oscdft_ns_type),      POINTER               :: nst
         INTEGER                                          :: iconstr, ioscdft

         inp => ctx%inp
         idx => ctx%idx
         nst => ctx%nst

         test = 0.D0
         SELECT CASE (ctx%inp%convergence_type)
            CASE (CONV_MULTIPLIERS)
               test = oscdft_array_to_scalar(ctx, ctx%multipliers - ctx%old_multipliers)
               conv = (test < threshold)
            CASE (CONV_GRADIENT)
               temp = nst%gradient
               DO iconstr=1,idx%nconstr
                  ioscdft = idx%iconstr2ioscdft(iconstr)
                  IF (inp%constraint_applied(ioscdft) == CONSTR_GE) THEN
                     temp(iconstr) = MIN(temp(iconstr), 0.0D0)
                  ELSE IF (inp%constraint_applied(ioscdft) == CONSTR_LE) THEN
                     temp(iconstr) = MAX(temp(iconstr), 0.0D0)
                  END IF
               END DO
               test = oscdft_array_to_scalar(ctx, temp)
               conv = (test < threshold)
            CASE (CONV_ENERGY)
               test = ABS(SUM(ctx%nst%gradient * ctx%multipliers))
               conv = (test < threshold).AND..NOT.ALL(ctx%multipliers.EQ.0.D0)
            CASE (CONV_ALWAYS_FALSE)
               conv = .false.
            CASE (CONV_ALWAYS_TRUE)
               conv = .true.
            CASE DEFAULT
               CALL errore("oscdft_multipliers_convergence_test",&
                           "invalid convergence",&
                           ABS(ctx%inp%convergence_type))
         END SELECT

         WRITE(stdout, 100) test, TRIM(threshold_label), threshold
         IF (conv) THEN
            WRITE(stdout, 101) TRIM(threshold_label)
         ELSE
            WRITE(stdout, 102) TRIM(threshold_label)
         END IF
         100 FORMAT("OSCDFT: convergence test of", ES14.7, " vs. ", A, " threshold of", ES14.7)
         101 FORMAT("OSCDFT: ", A, " test: PASSED")
         102 FORMAT("OSCDFT: ", A, " test: FAILED")
      END FUNCTION oscdft_multipliers_convergence_test

      SUBROUTINE oscdft_run_pwscf(ctx)
         USE control_flags,    ONLY : conv_elec, restart
         USE oscdft_wfcO,      ONLY : oscdft_init_wfcO
         USE wavefunctions,    ONLY : evc
         USE klist,            ONLY : nks
         USE io_files,         ONLY : iunwfc, nwordwfc
         USE buffers,          ONLY : get_buffer, save_buffer
         USE input_parameters, ONLY : startingwfc
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),   POINTER               :: inp
         TYPE(oscdft_indices_type), POINTER               :: idx
         TYPE(oscdft_ns_type),      POINTER               :: nst
         LOGICAL  :: multipliers_converge
         INTEGER  :: ioscdft, iconstr, ik
         REAL(DP) :: change_in_multipliers(ctx%idx%nconstr)

         inp => ctx%inp
         idx => ctx%idx
         nst => ctx%nst

         ctx%inp%hpsi_sum_band = .false.
         CALL oscdft_init_wfcO(ctx)

         SELECT CASE (inp%iteration_type)
            CASE (ITER_MULTIPLIERS_RHO)
               IF (inp%get_ground_state_first) THEN
                  CALL electrons()
                  CALL oscdft_new_ns(ctx, "oscdft_run_pwscf")
                  CALL oscdft_get_occupation_numbers(ctx, .false.)
               END IF
               inp%get_ground_state_first = .false.

               CALL electrons()
            CASE (ITER_RHO_MULTIPLIERS)
               ctx%warming_up = .true.

               IF (ctx%inp%get_ground_state_first) THEN
                  CALL electrons()
                  CALL oscdft_new_ns(ctx)
                  CALL oscdft_get_occupation_numbers(ctx, .false.)
               ELSE IF (startingwfc == 'file') THEN
                  CALL oscdft_new_ns(ctx)
                  CALL oscdft_get_occupation_numbers(ctx, .false.)
               END IF

               multipliers_converge = .false.
               ctx%warming_up = .false.
               CALL mp_bcast(ctx%warming_up, ionode_id, intra_image_comm)
               IF (.NOT.restart.OR.ctx%preserve_iter_when_restart_done) THEN
                  ctx%multiplier_iter = 0
               END IF
               ctx%preserve_iter_when_restart_done = .TRUE.

               nst%ns(:,:,:) = 0.D0
               nst%gradient(:) = 0.D0

               DO WHILE (.NOT.multipliers_converge)
                  ctx%multiplier_iter = ctx%multiplier_iter + 1
                  WRITE(stdout, 500) ctx%multiplier_iter

                  CALL electrons()

                  IF (conv_elec) THEN
                     conv_elec = .false.
                     CALL oscdft_new_ns(ctx)
                     CALL oscdft_get_occupation_numbers(ctx, .false.)
                     ctx%old_multipliers = ctx%multipliers
                     CALL oscdft_get_deriv(ctx)

                     IF (idx%nconstr > 0) THEN
                        CALL oscdft_optimize_multipliers(ctx)

                        IF (ctx%has_debug) THEN
                           DO iconstr=1,idx%nconstr
                              ioscdft = idx%iconstr2ioscdft(iconstr)
                              SELECT CASE (inp%constraint_applied(ioscdft))
                                 CASE (CONSTR_D1)
                                    ctx%multipliers(iconstr) = ctx%old_multipliers(iconstr)+&
                                                               inp%debug_step(ioscdft)
                              END SELECT
                           END DO
                        END IF

                        CALL oscdft_clamp_multipliers(ctx)

                        WRITE(stdout, 300) "current multipliers"
                        WRITE(stdout, 400) ctx%old_multipliers
                        WRITE(stdout, 300) "updated multipliers"
                        WRITE(stdout, 400) ctx%multipliers
                        multipliers_converge = oscdft_multipliers_convergence_test(ctx, ctx%conv_thr, "convergence")

                        IF (inp%rhoiter) THEN
                           IF (ctx%iter > 2) THEN
                              multipliers_converge = .false.
                           END IF
                        END IF
                        IF (inp%total_energy_threshold /= 0.D0) THEN
                           IF (ABS(ctx%total_energy - ctx%old_total_energy) >= inp%total_energy_threshold) THEN
                              multipliers_converge = .false.
                           END IF
                        END IF

                        IF (multipliers_converge) THEN
                           ctx%global_start_index = ctx%global_start_index + 1
                           WRITE(stdout, 101) "global_start_index", ctx%global_start_index
                           IF (ctx%global_start_index > MAXVAL(inp%start_index)) THEN
                              multipliers_converge = .true.
                              WRITE(stdout, 600)
                           ELSE
                              multipliers_converge = .false.
                           END IF
                        END IF
                     ELSE
                        multipliers_converge = .true.
                     END IF
                     IF (multipliers_converge) conv_elec = .true.
                  ELSE
                     multipliers_converge = .true.
                  END IF

                  IF (ctx%multiplier_iter > inp%maxiter) multipliers_converge = .true.
                  IF (ctx%multiplier_iter < inp%miniter) multipliers_converge = .false.
               END DO
               IF (multipliers_converge) THEN
                  ctx%multipliers = ctx%old_multipliers
               END IF
         END SELECT
         ctx%conv_thr = MIN(inp%max_conv_thr,&
                            MAX(inp%min_conv_thr,&
                                ctx%conv_thr * inp%conv_thr_multiplier))

         101 FORMAT("OSCDFT: ", A, ": ", I5)
         300 FORMAT("OSCDFT: ", A)
         400 FORMAT("OSCDFT:   ", *(ES14.7, :, " "))
         500 FORMAT("OSCDFT: oscdft loop #", I5)
         600 FORMAT("OSCDFT: All multipliers have converged")
      END SUBROUTINE oscdft_run_pwscf

      SUBROUTINE oscdft_electrons(ctx, iter, et, nbnd, nkstot, nks)
         USE mp,               ONLY : mp_barrier, mp_max
         USE control_flags,    ONLY : isolve, restart
         USE mp_pools,         ONLY : inter_pool_comm, intra_pool_comm
         USE check_stop,       ONLY : check_stop_now
         USE parallel_include
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),   POINTER               :: inp
         TYPE(oscdft_indices_type), POINTER               :: idx
         TYPE(oscdft_ns_type),      POINTER               :: nst
         INTEGER,                   INTENT(IN)            :: iter, nbnd, nkstot, nks
         REAL(DP),                  INTENT(IN)            :: et(nbnd,nkstot)
         REAL(DP),                  EXTERNAL              :: get_clock
         LOGICAL                                          :: multipliers_converge
         INTEGER                                          :: iconstr, ioscdft

         inp => ctx%inp
         idx => ctx%idx
         nst => ctx%nst

         SELECT CASE (inp%iteration_type)
            CASE (ITER_RHO_MULTIPLIERS)
               CALL c_bands(iter)
               CALL mp_barrier(intra_image_comm)
               ctx%iter = iter
            CASE (ITER_MULTIPLIERS_RHO)
               ! DISABLE THE START INDEX HERE
               ctx%global_start_index = MAXVAL(inp%start_index) + 1
               ctx%iter = iter
               IF (inp%get_ground_state_first) THEN
                  CALL c_bands(iter)
                  CALL mp_barrier(intra_image_comm)
               ELSE
                  IF (ctx%warm_up_iter >= inp%warm_up_niter) ctx%warming_up = .false.
                  IF (ctx%warming_up) THEN
                     ctx%warm_up_iter = ctx%warm_up_iter + 1
                     WRITE(stdout, 500) ctx%warm_up_iter
                     CALL c_bands(iter)
                     CALL mp_barrier(intra_image_comm)
                  ELSE IF (idx%nconstr == 0) THEN
                     CALL c_bands(iter)
                     CALL mp_barrier(intra_image_comm)
                     CALL oscdft_new_ns(ctx, "oscdft_electrons")
                     CALL oscdft_get_occupation_numbers(ctx, .false.)
                  ELSE
                     WRITE(stdout, 600) iter
                     IF (.NOT.restart.OR.ctx%preserve_iter_when_restart_done) THEN
                        ctx%multiplier_iter = 0
                     END IF
                     ctx%preserve_iter_when_restart_done = .TRUE.
                     multipliers_converge = .false.
                     DO WHILE (.NOT.multipliers_converge)
                        ctx%multiplier_iter = ctx%multiplier_iter + 1
                        WRITE(stdout, 601) ctx%multiplier_iter
                        ! this is just for debugging
                        IF ((inp%test_exit_rho_iter >= 0)    .AND.&
                            (inp%test_exit_oscdft_iter >= 0) .AND.&
                            (inp%test_exit_rho_iter <= iter).AND.&
                            (inp%test_exit_oscdft_iter <= ctx%multiplier_iter)) THEN
                            CALL oscdft_debug_write_exit_file()
                        END IF
                        WRITE(stdout, 900) get_clock('PWSCF')
                        CALL c_bands(iter)
                        IF (check_stop_now()) THEN
                           RETURN
                        ENDIF

                        CALL oscdft_new_ns(ctx, "oscdft_electrons")
                        CALL oscdft_get_occupation_numbers(ctx, .false.)
                        CALL oscdft_get_deriv(ctx)

                        ctx%old_multipliers = ctx%multipliers
                        CALL oscdft_optimize_multipliers(ctx)

                        IF (ctx%has_debug) THEN
                           DO iconstr=1,idx%nconstr
                              ioscdft = idx%iconstr2ioscdft(iconstr)
                              IF (inp%constraint_applied(ioscdft) == CONSTR_D1) THEN
                                 ctx%multipliers(iconstr) = ctx%old_multipliers(iconstr) +&
                                                            inp%debug_step(ioscdft)
                              END IF
                           END DO
                        END IF

                        CALL oscdft_clamp_multipliers(ctx)

                        WRITE(stdout, 300) "current multipliers"
                        WRITE(stdout, 400) ctx%old_multipliers(1:idx%nconstr)
                        WRITE(stdout, 300) "updated multipliers"
                        WRITE(stdout, 400) ctx%multipliers(1:idx%nconstr) 
                        multipliers_converge = oscdft_multipliers_convergence_test(ctx, ctx%conv_thr, "inner convergence")

                        IF (multipliers_converge) THEN
                           ! see comment CONV_THR_REDUCE
                           ctx%conv_thr_not_minimum = ctx%conv_thr > ctx%inp%min_conv_thr
                           IF (ctx%conv_thr_not_minimum) THEN
                              WRITE(stdout, 300) "multipliers converge, decreasing multiplier conv_thr"
                           ELSE IF (ctx%inp%final_conv_thr > 0.D0) THEN
                              WRITE(stdout, 300) "multipliers converge, checking for final_conv_thr convergence"
                           ELSE
                              WRITE(stdout, 300) "multipliers converge, looping"
                           END IF
                        END IF
                     END DO
                     ! CONV_THR_REDUCE
                     ! mix rho will be called after this line,
                     ! so the ctx%conv_thr > ctx%inp%min_conv_thr in mix_rho
                     ! should be calculated before this 
                     ctx%conv_thr = MIN(inp%max_conv_thr,&
                                        MAX(inp%min_conv_thr,&
                                            ctx%conv_thr * inp%conv_thr_multiplier))
                  END IF
               END IF
         END SELECT
         300 FORMAT("OSCDFT: ", A)
         400 FORMAT("OSCDFT:   ", *(ES14.7, :, " "))
         500 FORMAT("OSCDFT: warm up loop #", I5, "; constraints are NOT applied")
         600 FORMAT("OSCDFT: main loop #", I5, "; constraints are applied")
         601 FORMAT("OSCDFT: oscdft loop #", I5)
         900 FORMAT(/'     total cpu time spent up to now is ',F10.1,' secs' )
      END SUBROUTINE oscdft_electrons

      SUBROUTINE oscdft_electrons_restart(ctx)
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),           POINTER       :: inp

         inp => ctx%inp

         IF (.NOT.(inp%oscdft_type==1)) RETURN

         ! alternatively, the weights could be saved during exit and loaded here
         CALL weights()
      END SUBROUTINE oscdft_electrons_restart

      ! this is just for debugging
      SUBROUTINE oscdft_debug_write_exit_file
         USE io_files, ONLY : exit_file
         IMPLICIT NONE

         INTEGER :: iun
         INTEGER, EXTERNAL :: find_free_unit

         IF (ionode) THEN
            iun = find_free_unit()
            OPEN(unit=iun, FILE=exit_file, FORM="FORMATTED", STATUS="REPLACE")
            CLOSE(iun)
         END IF
      END SUBROUTINE oscdft_debug_write_exit_file

      SUBROUTINE oscdft_scf_energy(ctx, plugin_etot)
         IMPLICIT NONE
         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),   POINTER       :: inp
         REAL(DP),                  INTENT(INOUT) :: plugin_etot

         inp => ctx%inp

         IF (.NOT.(inp%oscdft_type==1)) RETURN

         IF (ctx%initialized.AND.(.NOT.ctx%warming_up)) THEN
            CALL oscdft_calc_energy_offset(ctx)
            plugin_etot = plugin_etot + ctx%energy
         END IF
      END SUBROUTINE oscdft_scf_energy

      SUBROUTINE oscdft_calc_energy_offset(ctx)
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT) :: ctx
         INTEGER                                  :: iconstr, ioscdft, oidx, osum

         IF (ctx%initialized.AND.(.NOT.ctx%warming_up)) THEN
            IF (ctx%recalculate_ns) THEN
               CALL oscdft_new_ns(ctx, "oscdft_calc_energy_offset")
               CALL oscdft_get_occupation_numbers(ctx, .true.)
            END IF
            ctx%energy = 0.D0
            DO iconstr=1,ctx%idx%nconstr
               ioscdft = ctx%idx%iconstr2ioscdft(iconstr)
               oidx = ctx%inp%occup_index(ioscdft)
               IF (oidx == OCCUP_SUM) THEN
                  osum = ctx%inp%occup_index_sum(2,ioscdft)
                  ctx%energy = ctx%energy - (ctx%multipliers(iconstr) * ctx%nst%occup_eigvals(osum,iconstr))
               ELSE IF (ctx%inp%constraint_applied(ioscdft) == CONSTR_GE  .OR. &
                        ctx%inp%constraint_applied(ioscdft) == CONSTR_LE  .OR. &
                        ctx%inp%constraint_applied(ioscdft) == CONSTR_GE2 .OR. &
                        ctx%inp%constraint_applied(ioscdft) == CONSTR_LE2) THEN
                  ctx%energy = ctx%energy - (ctx%multipliers(iconstr) * ctx%nst%occup_numbers(iconstr))
               ELSE
                  ctx%energy = ctx%energy - (ctx%multipliers(iconstr) * ctx%inp%target_occup(ioscdft))
               END IF
            END DO
         END IF
      END SUBROUTINE oscdft_calc_energy_offset

      SUBROUTINE oscdft_save_total_energy(ctx, etot)
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT) :: ctx
         REAL(DP),                  INTENT(IN)    :: etot

         ctx%old_total_energy = ctx%total_energy
         ctx%total_energy = etot
      END SUBROUTINE oscdft_save_total_energy

      SUBROUTINE oscdft_print_energies(ctx)
         USE ener,             ONLY : etot
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),   POINTER       :: inp

         inp => ctx%inp

         IF (.NOT.(inp%oscdft_type==1)) RETURN

         CALL oscdft_save_total_energy(ctx, etot)

         IF (.NOT.ctx%initialized) RETURN
         IF (ctx%warming_up) RETURN
         IF (ALL(ctx%inp%constraint_applied == CONSTR_FALSE)) RETURN

         WRITE(stdout, 100) ctx%energy
         100 FORMAT( '     OSCDFT energy             =',F17.8,' Ry' )
      END SUBROUTINE oscdft_print_energies

      SUBROUTINE oscdft_save(ctx)
         USE io_files, ONLY : restart_dir
         USE klist,    ONLY : nelup, neldw
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT) :: ctx
         INTEGER                                  :: iun, curr_dim, row, ioscdft
         INTEGER, EXTERNAL                        :: find_free_unit

         IF (ionode) THEN
            iun = find_free_unit()
            OPEN(UNIT=iun, FILE=TRIM(restart_dir())//'/oscdft_save',&
               FORM="FORMATTED", STATUS="REPLACE")
            WRITE(iun, 101) nelup, neldw
            WRITE(iun, 101) ctx%total_energy
            IF (ctx%idx%nconstr > 0) THEN
               WRITE(iun, 101) ctx%old_multipliers
            ELSE
               WRITE(iun, 100) 0
            END IF
            WRITE(iun, 101) ctx%conv_thr
            WRITE(iun, 100) ctx%warm_up_iter, ctx%global_start_index, ctx%multiplier_iter
            WRITE(iun, 102) ctx%inp%get_ground_state_first

            DO ioscdft=1,ctx%inp%noscdft
               curr_dim = ctx%idx%ns_dim(ioscdft)
               DO row=1,curr_dim
                  WRITE(iun, 101) ctx%nst%ns(row,1:curr_dim,ioscdft)
               END DO
            END DO

            CLOSE(iun)
         END IF
         100 FORMAT(*(I5))
         101 FORMAT(*(ES14.7))
         102 FORMAT(L)
      END SUBROUTINE oscdft_save

      SUBROUTINE oscdft_close_files(ctx)
         USE io_files,                     ONLY : restart_dir
         USE oscdft_wavefunction_subs,     ONLY : oscdft_destroy_wavefunctions
         USE control_flags,                ONLY : use_gpu
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),   POINTER       :: inp
         INTEGER                                  :: file_size, iun
         INTEGER, EXTERNAL                        :: find_free_unit
         CHARACTER(LEN=:), ALLOCATABLE            :: str

         inp => ctx%inp

         IF (.NOT.(inp%oscdft_type==1)) RETURN

         CALL oscdft_destroy_wavefunctions(ctx%wfcO, "DELETE")
         CALL oscdft_destroy_wavefunctions(ctx%wfcS, "DELETE")
         CALL oscdft_destroy_wavefunctions(ctx%forces%wfcO, "DELETE")
         IF (ionode) THEN
            INQUIRE(file="oscdft.in", SIZE=file_size)
            iun = find_free_unit()
            IF (file_size > 0) THEN
               OPEN(UNIT=iun, FILE="oscdft.in", status="old",&
                    FORM="UNFORMATTED", ACTION="READ", ACCESS="STREAM")
               ALLOCATE(CHARACTER(LEN=file_size) :: str)
               READ(UNIT=iun) str
               CLOSE(iun)
               OPEN(UNIT=iun, FILE=TRIM(restart_dir())//'/oscdft.in',&
                    FORM="FORMATTED", STATUS="REPLACE")
               WRITE(iun, 100, advance="no") str
               CLOSE(iun)
               DEALLOCATE(str)
            END IF
         END IF
         CALL oscdft_save(ctx)
         ctx%warming_up = .true.
         100 FORMAT(A)
      END SUBROUTINE oscdft_close_files

      SUBROUTINE oscdft_h_diag(ctx, h_diag)
         USE lsda_mod,         ONLY : isk, current_spin
         USE klist,            ONLY : ngk
         USE klist,            ONLY : nks
         USE buffers,          ONLY : get_buffer
         USE wvfct,            ONLY : current_k, npwx
         USE noncollin_module, ONLY : npol

         IMPLICIT NONE

         REAL(dp), INTENT(INOUT) :: h_diag(npwx,npol)
         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),           POINTER       :: inp
         TYPE(oscdft_indices_type),         POINTER       :: idx
         TYPE(oscdft_ns_type),              POINTER       :: nst
         TYPE(oscdft_wavefunction_type),    POINTER       :: wfcO

         REAL(DP)                                         :: alpha
         INTEGER                                          :: curr_dim, h, k, i, ik,&
                                                             iconstr, ioscdft, oidx,&
                                                             isum, osum, h_off, k_off,&
                                                             npw
         COMPLEX(DP), POINTER :: wfcO_wfc(:,:)

         inp => ctx%inp
         idx => ctx%idx
         nst => ctx%nst
         wfcO => ctx%wfcO

         IF (.NOT.(inp%oscdft_type==1)) RETURN
         IF (.NOT.ctx%initialized) RETURN
         IF (ctx%warming_up) RETURN
         IF (idx%nconstr == 0) RETURN
         IF (.NOT.ANY(inp%spin_index(idx%iconstr2ioscdft(:)) == current_spin)) RETURN

         CALL start_clock("oscdft_hdiag")
         npw = ngk(current_k)

         IF (.NOT.ctx%wfc_allocated) THEN
            CALL errore("oscdft_h_diag", "wfc not allocated", 1)
         END IF

         ik = current_k
         IF (isk(ik) == current_spin) THEN
            IF (nks > 1) THEN
               CALL oscdft_get_buffer(wfcO, ik)
            END IF
            wfcO_wfc => wfcO%wfc
            !$acc data copyin(wfcO_wfc)

            DO iconstr=1,idx%nconstr
               ioscdft = idx%iconstr2ioscdft(iconstr)
               oidx = inp%occup_index(ioscdft)
               IF (inp%start_index(ioscdft) > ctx%global_start_index) CYCLE
               IF (inp%spin_index(ioscdft) == current_spin) THEN
                  curr_dim = idx%ns_dim(ioscdft)
                  DO h=1,curr_dim
                     h_off = wfcO%get_offset(idx%constr, h, iconstr, -1)
                     DO k=1,curr_dim
                        k_off = wfcO%get_offset(idx%constr, k, iconstr, -1)
                        IF (oidx == OCCUP_TRACE) THEN
                           alpha = MERGE(ctx%multipliers(iconstr), 0.D0, h == k)
                        ELSE IF (oidx == OCCUP_SUM) THEN
                           osum = inp%occup_index_sum(2,ioscdft)
                           alpha = ctx%multipliers(iconstr)*&
                                   ctx%nst%occup_eigvects(h,osum,iconstr)*&
                                   ctx%nst%occup_eigvects(k,osum,iconstr)
                        ELSE
                           alpha = ctx%multipliers(iconstr)*&
                                   ctx%nst%occup_eigvects(h,oidx,iconstr)*&
                                   ctx%nst%occup_eigvects(k,oidx,iconstr)
                        END IF
                        !$acc data present(h_diag)
                        !$acc data present(wfcO_wfc)
                        !$acc parallel loop
                        DO i=1,npw
                           h_diag(i,1) = h_diag(i,1)+&
                                         alpha*DBLE(wfcO_wfc(i,k_off)*CONJG(wfcO_wfc(i,h_off)))
                        END DO
                        !$acc end data
                        !$acc end data
                     END DO
                  END DO
               END IF
            END DO
            !$acc end data
         END IF
         CALL stop_clock("oscdft_hdiag")
      END SUBROUTINE oscdft_h_diag

      SUBROUTINE oscdft_h_psi(ctx, lda, n, m, psi, hpsi)
         USE noncollin_module,    ONLY : npol
         USE becmod,              ONLY : bec_type, calbec, allocate_bec_type, deallocate_bec_type
         USE klist,               ONLY : nks, nelup, neldw
         USE buffers,             ONLY : get_buffer
         USE lsda_mod,            ONLY : isk, current_spin
         USE io_files,            ONLY : nwordwfc
         USE gvect,               ONLY : gstart
         USE control_flags,       ONLY : gamma_only
         USE wvfct,               ONLY : btype, current_k, nbnd, wg
         USE mp_bands,            ONLY : intra_bgrp_comm
         USE mp,                  ONLY : mp_barrier
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),        POINTER          :: inp
         TYPE(oscdft_indices_type),      POINTER          :: idx
         TYPE(oscdft_wavefunction_type), POINTER          :: wfcO
         INTEGER,                   INTENT(IN)            :: lda, n, m
         COMPLEX(DP),               INTENT(IN)            :: psi(lda*npol,m)
         COMPLEX(DP),               INTENT(INOUT)         :: hpsi(lda*npol,m)
         TYPE(bec_type)                                   :: proj
         LOGICAL                                          :: done_calbec
         INTEGER                                          :: iconstr, ioscdft,&
                                                             ik, h, k, curr_dim,&
                                                             ibnd, i, oidx, isum, osum, jsum,&
                                                             h_off, k_off
         REAL(DP)                                         :: occ_const
         COMPLEX(DP)                                      :: alpha! , deriv_f(lda*npol,m)

         inp => ctx%inp
         idx => ctx%idx
         wfcO => ctx%wfcO

         IF (.NOT.(inp%oscdft_type==1)) RETURN
         IF (.NOT.ctx%initialized) RETURN
         IF (ctx%warming_up) RETURN
         IF (idx%nconstr == 0) RETURN

         CALL start_clock("oscdft_hpsi")
         ! sum_hk u_h^I u_k^I |phi_k> <phi_h|psi>
         CALL check_bec_type_unallocated(proj)
         ik = current_k
         CALL allocate_bec_type(m, wfcO%n, proj)
         IF (nks > 1) CALL oscdft_get_buffer(wfcO, ik)
         ! gets <psi|phi_h>, then call ZGERC where this term is turned to its complex conjugate (<phi_h|psi>)
         CALL calbec(n, psi, wfcO%wfc, proj)
         ! deriv_f = (0.D0, 0.D0)

         DO iconstr=1,idx%nconstr
            ioscdft = idx%iconstr2ioscdft(iconstr)
            oidx = inp%occup_index(ioscdft)
            IF (inp%start_index(ioscdft) > ctx%global_start_index) CYCLE
            IF ((inp%constraint_applied(ioscdft) /= CONSTR_FALSE).AND.&
                (inp%spin_index(ioscdft) == current_spin)) THEN
               curr_dim = idx%ns_dim(ioscdft)
               DO k=1,curr_dim
                  k_off = wfcO%get_offset(idx%constr, k, iconstr, -1)
                  DO h=1,curr_dim
                     h_off = wfcO%get_offset(idx%constr, h, iconstr, -1)
                     IF (oidx == OCCUP_TRACE) THEN
                        occ_const = MERGE(1.D0, 0.D0, h == k)
                     ELSE IF (oidx == OCCUP_SUM) THEN
                        osum = inp%occup_index_sum(2,ioscdft)
                        occ_const = ctx%nst%occup_eigvects(h,osum,iconstr)*&
                                    ctx%nst%occup_eigvects(k,osum,iconstr)
                     ELSE
                        occ_const = ctx%nst%occup_eigvects(h,oidx,iconstr)*&
                                    ctx%nst%occup_eigvects(k,oidx,iconstr)
                     ENDIF
                     alpha = CMPLX(occ_const * ctx%multipliers(iconstr), 0.D0, kind=DP)
                     IF (gamma_only) THEN
                        CALL DGER(2*n,m,DBLE(alpha),&
                                  wfcO%wfc(1,h_off),1,&
                                  proj%r(1:m,k_off),1,&
                                  hpsi,2*lda)
                     ELSE
                        CALL ZGERC(n,m,alpha,&
                                   wfcO%wfc(1,h_off),1,&
                                   proj%k(1:m,k_off),1,&
                                   hpsi,lda)
                     END IF
                  END DO
               END DO
            END IF
         END DO
         CALL deallocate_bec_type(proj)
         CALL mp_barrier(intra_bgrp_comm)
         CALL stop_clock("oscdft_hpsi")
      END SUBROUTINE oscdft_h_psi

      SUBROUTINE oscdft_mix_rho(ctx, conv)
         IMPLICIT NONE

         TYPE(oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),   POINTER       :: inp
         LOGICAL,                   INTENT(INOUT) :: conv
         LOGICAL                                  :: test

         inp => ctx%inp

         IF (.NOT.(inp%oscdft_type==1)) RETURN

         IF (ctx%inp%iteration_type /= ITER_MULTIPLIERS_RHO) RETURN
         IF (ctx%inp%get_ground_state_first) RETURN
         IF (ctx%warming_up) THEN
            conv = .false.
            RETURN
         END IF
         IF (ctx%iter < ctx%inp%miniter) THEN
            WRITE(stdout, 300) "iter < miniter"
            conv = .false.
            RETURN
         END IF
         IF (ctx%idx%nconstr > 0) THEN
            IF (ctx%conv_thr_not_minimum) THEN
               conv = .false.
               RETURN
            END IF
            IF (ctx%inp%final_conv_thr > 0.0) THEN
               test = oscdft_multipliers_convergence_test(ctx, ctx%inp%final_conv_thr, "outer convergence")
               IF (.NOT.test) THEN
                  conv = .false.
                  RETURN
               END IF
            END IF
         END IF
         300 FORMAT("OSCDFT: ", A)
      END SUBROUTINE oscdft_mix_rho

      SUBROUTINE oscdft_sum_band(ctx)
         IMPLICIT NONE

         TYPE (oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE (oscdft_input_type),           POINTER       :: inp

         inp => ctx%inp

         IF (.NOT.(inp%oscdft_type==1)) RETURN

         IF (ctx%initialized.AND.(.NOT.ctx%warming_up).AND.ctx%inp%sum_band) THEN
            CALL oscdft_new_ns(ctx, "oscdft_sum_band")
            CALL oscdft_get_occupation_numbers(ctx, .true.)
         END IF
      END SUBROUTINE oscdft_sum_band

      SUBROUTINE oscdft_print_ns(ctx)
         IMPLICIT NONE

         TYPE (oscdft_context_type), INTENT(INOUT), TARGET :: ctx
         TYPE(oscdft_input_type),    POINTER       :: inp

         inp => ctx%inp

         IF (.NOT.(inp%oscdft_type==1)) RETURN

         IF (ctx%initialized.AND.(.NOT.ctx%warming_up)) THEN
            CALL oscdft_new_ns(ctx, "oscdft_print_ns")
            CALL oscdft_get_occupation_numbers(ctx, .false.)
         END IF
      END SUBROUTINE oscdft_print_ns

      SUBROUTINE oscdft_init_constraint (ctx, ldmx, starting_ns, lda_plus_u, lda_plus_u_kind, is_hubbard)        
         !
         !! This routine allocates and initializes some arrays
         !! needed for constrained calculations
         !
         USE kinds,               ONLY : DP
         USE ions_base,           ONLY : nat, ityp
         USE parameters,          ONLY : ntypx
         USE upf_params,          ONLY : lqmax
         USE lsda_mod,            ONLY : nspin
         USE noncollin_module,    ONLY : noncolin
         USE oscdft_context,      ONLY : oscdft_context_type

         IMPLICIT NONE
         CLASS(oscdft_context_type), INTENT(INOUT) :: ctx
         INTEGER, INTENT(IN) :: ldmx
         REAL(DP), INTENT(INOUT) :: starting_ns(lqmax,2,ntypx)
         LOGICAL, INTENT(IN) :: lda_plus_u
         LOGICAL, INTENT(IN) :: is_hubbard(ntypx)
         INTEGER, INTENT(IN) :: lda_plus_u_kind
         INTEGER :: na, nt, is, m1, m2, ldim
         !
         IF (.NOT.(ctx%inp%oscdft_type==2)) RETURN
         !
         ctx%is_constraint = .FALSE.
         !
         ALLOCATE(ctx%constraining(nat))
         ctx%constraining(:) = .FALSE.
         !
         ctx%conv = .TRUE.
         !
         ! Sanity checks
         !
         IF (.NOT.lda_plus_u) CALL errore( 'oscdft_init_constraint', &
             'Constrained calculation with oscdft_type=2 can be done only for DFT+Hubbard', 1 ) 
         !
         IF (lda_plus_u_kind/=0 .AND. lda_plus_u_kind/=2) CALL errore( 'oscdft_init_constraint', &
             'Constrained calculation with oscdft_type=2 can be done only for lda_plus_u_kind=0 or 2', 1 )
         !
         IF (noncolin) CALL errore( 'oscdft_init_constraint', &
             'Constrained calculation with oscdft_type=2 and noncolin=true is not supported', 1 ) 
         !
         ! Check whether the user specified starting_ns_eigenvalue in the input
         IF (ANY(starting_ns /= -1.0_dp)) THEN
            WRITE( stdout, '(/5X,"WARNING!!! starting_ns_eigenvalue is not allowed when OSCDFT is enabled.")')
            WRITE( stdout, '( 5X,"Instead, specify the target eigenvalues of the occupation matrix")')
            WRITE( stdout, '( 5X,"in the OSCDFT input, by specifying constraint_diag = .true.")')
            WRITE( stdout, '( 5X,"See the OSCDFT documentation for more details.")')
            CALL errore('dftu_iosys',' starting_ns_eigenvalue and OSCDFT are not compatible', 1)
         ENDIF
         !
         IF (ANY(ctx%inp%occupation /= -2.0_DP)) THEN
            !
            ALLOCATE (ctx%constraint(ldmx, ldmx, nspin, nat))
            ctx%constraint(:,:,:,:) = 0.0_DP
            !
            ctx%conv = .FALSE.
            !
            DO na = 1, nat
               !
               nt = ityp(na)
               !
               ctx%constraining(na) = ANY(ctx%inp%occupation(:,:,:,na) /= -2.0_DP)
               !
               IF (ctx%constraining(na) .AND. .NOT.is_hubbard(nt)) CALL errore( 'oscdft_init_constraint', &
                  'Constraints can be applied only to Hubbard atoms', 1 )
               !
               IF (ctx%inp%constraint_diag) THEN
                  ! In this case, the array ctx%inp%occupation has been used as
                  ! a workspace to store the eigenvalues of the occupation matrix. 
                  ! Now we copy this data to starting_ns that will be used to 
                  ! build the occupation matrix in the routine ns_adj, which in turn
                  ! will write the result in the array ctx%inp%occupation for the 
                  ! constrained calculation.
                  DO is = 1, nspin
                     DO m1 = 1, ldmx
                        starting_ns(m1, is, nt) = ctx%inp%occupation(m1, m1, is, na)
                     ENDDO
                  ENDDO
               ENDIF
               !  
            ENDDO 
            !
            IF (ANY(ctx%constraining(:))) ctx%is_constraint = .TRUE. 
            ! 
         ENDIF
         !
         RETURN
         !
      END SUBROUTINE oscdft_init_constraint

      SUBROUTINE oscdft_constrain_ns (ctx, lmax, l, nscur, conv, dr2, tr2, iter)
      !
      !! This routine updates the constraint term based on the occupation matrix ns
      !! It computes the mean absolute error (MAE) and the root mean square error (RMSE),
      !! based on the current and target occupation matrices. When MAE is smaller than 
      !! the threshold constraint_conv_thr, the convergence on the constraint is reached
      !! and it is released, and the SCF calculation continues without the constraint.
      !! If the system "likes" the energy minimum with the target occupations, then 
      !! the SCF calculation will be stabilized there after the releasing of the constraint;
      !! otherwise, the system will drift away from this target minim.
      !
      USE kinds,           ONLY : DP
      USE io_global,       ONLY : stdout
      USE parameters,      ONLY : ntypx
      USE ions_base,       ONLY : nat, ityp
      USE lsda_mod,        ONLY : nspin
      USE oscdft_context,  ONLY : oscdft_context_type

      IMPLICIT NONE
      !
      CLASS(oscdft_context_type), INTENT(INOUT) :: ctx
      INTEGER, INTENT(IN) :: lmax
      INTEGER, INTENT(IN) :: l(ntypx)
      LOGICAL, INTENT(INOUT) :: conv
      REAL(DP), INTENT(IN) ::dr2, tr2
      INTEGER, INTENT(IN) :: iter
      REAL(DP), INTENT(IN) :: nscur(2*lmax+1,2*lmax+1,nspin,nat)
      !
      INTEGER, PARAMETER :: ldmx = 7
      INTEGER :: m1, m2, ldim, na, is, nt, counter
      REAL(DP) :: ns_diff, totdist, totdist2, &
                  mae, & ! Mean absolute error (MAE)
                         ! MAE = \sum_{I,s,m1,m2=1}^{N} |n^Is_m1m2 - n^Is_m1m2_target| / N
                  rmse   ! Root mean square error (RMSE)
                         ! RMSE = SQRT{ \sum_{I,s,m1,m2=1}^{N} (n^Is_m1m2 - n^Is_m1m2_target)^2 / N }
      !
      IF (.NOT.(ctx%inp%oscdft_type==2)) RETURN
      !
      totdist = 0.0_DP
      totdist2 = 0.0_DP
      mae = 0.0_DP
      rmse = 0.0_DP
      counter = 0
      !
      IF (.NOT. conv) THEN
         !
         WRITE(stdout,'(/5x,19("="), " APPLYING CONSTRAINTS ",19("="))')   
         WRITE(stdout, '(4x," Constraint strength: ", f7.3,/)') ctx%inp%constraint_strength
         !
         DO na = 1, nat
            !
            IF (ctx%constraining(na)) THEN
               !
               nt = ityp(na)
               ldim = 2 * l(nt) + 1
               !
               DO is = 1, nspin
                  !
                  DO m1 = 1, ldim
                     DO m2 = 1, ldim
                        ! 
                        ns_diff = nscur(m1, m2, is, na) - ctx%inp%occupation(m1, m2, is, na)
                        !
                        ctx%constraint(m1, m2, is, na) = ctx%constraint(m1, m2, is, na) + &
                                                         ctx%inp%constraint_mixing_beta * ns_diff
                        !
                        totdist  = totdist  + ABS(ns_diff)
                        totdist2 = totdist2 + (ns_diff)**2
                        !
                     ENDDO
                  ENDDO
                  !
                  WRITE( stdout,'(5x,"Constraint matrix for atom ",1x,i3,1x," and spin ",1x,i2)') na, is
                  DO m1 = 1, ldim
                    WRITE( stdout,'(5x,7f7.3)') ( DBLE(ctx%constraint(m1,m2,is,na)), m2=1, ldim )
                  ENDDO
                  !
               ENDDO
               !
               counter = counter + nspin*ldim*ldim
               !
            ENDIF
            !
         ENDDO
         !
         ! Compute the MAE and RMSE
         IF (counter>0) THEN
            mae = totdist/counter
            rmse = SQRT(totdist2/counter)
         ENDIF
         WRITE(stdout, '(/5x,"MAE = ", f10.6, 2x," RMSE = ", f10.6)') mae, rmse
         !
         ! Check whether the convergence threshold on MAE is reached 
         IF (mae < ctx%inp%constraint_conv_thr .AND. iter > 1) THEN
             conv = .true.
             WRITE(stdout,'(/5x,"CONSTRAINT CONVERGENCE THRESHOLD REACHED")')
             WRITE(stdout,'(/5x,"SWITCHING OFF THE CONSTRAINT...")')
         ELSEIF (iter > ctx%inp%constraint_maxstep) THEN
             WRITE(stdout,'("constraint_maxstep reached")')
             CALL errore( 'oscdft_constrain_ns', &
                  'CONSTRAINT CONVERGENCE THRESHOLD HAS NOT BEEN REACHED', 1 ) 
         ENDIF
         !
      ELSEIF (dr2 < tr2) THEN
          conv = .true.
          WRITE(stdout,'("SCF convergence reached")')
      ENDIF
      !
      RETURN
      !
   END SUBROUTINE oscdft_constrain_ns

   SUBROUTINE oscdft_v_constraint (ctx, lmax, l, ns, vtot, etot)
      !
      !! Computes the contribution to the potential from the Lagrange multipliers used
      !! to constrain the occupation matrix to the target (DFT+U case).
      !
      USE kinds,           ONLY : DP
      USE parameters,      ONLY : ntypx
      USE ions_base,       ONLY : nat, ityp
      USE lsda_mod,        ONLY : nspin
      USE io_global,       ONLY : stdout
      USE control_flags,   ONLY : iverbosity
      USE oscdft_context,  ONLY : oscdft_context_type  
 
      IMPLICIT NONE

      CLASS(oscdft_context_type), INTENT(IN) :: ctx
      INTEGER, INTENT(IN) :: lmax
      INTEGER, INTENT(IN) :: l(ntypx)
      REAL(DP), INTENT(IN)    :: ns(2 * lmax + 1, 2 * lmax + 1, nspin, nat)
      REAL(DP), INTENT(INOUT) :: vtot(2 * lmax + 1, 2 * lmax + 1, nspin, nat)
      REAL(DP), INTENT(INOUT) :: etot
      REAL(DP) :: ec
      INTEGER :: is, na, nt, ldim, m1, m2
      !
      IF (.NOT.(ctx%inp%oscdft_type==2) .OR. .NOT.ctx%is_constraint .OR. ctx%conv) RETURN
      !
      ec = 0
      !
      DO na = 1, nat
         IF (ctx%constraining(na)) THEN
            nt = ityp(na)
            ldim = 2 * l(nt) + 1
            DO is = 1, nspin
               DO m1 = 1, ldim
                  DO m2 = 1, ldim
                     vtot(m1,m2,is,na) = vtot(m1, m2, is, na) + &
                                         ctx%inp%constraint_strength * &
                                         ctx%constraint(m1,m2,is,na)
                     ec = ec + ctx%inp%constraint_strength * &
                               ctx%constraint(m1,m2,is,na) * &
                               (ns(m1,m2,is,na) - ctx%inp%occupation(m1,m2,is,na))
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      !
      IF (nspin == 1) THEN
           ec = 2.d0 * ec
      ENDIF
      etot = etot + ec
      !
      IF (iverbosity > 0) THEN
         WRITE(stdout, '(/5x,"CONSTRAINT ENERGY = ", f9.4, 1x," (Ry)")') ec
      ENDIF
      ! 
      RETURN       
      !
   END SUBROUTINE oscdft_v_constraint

   SUBROUTINE oscdft_ns_set (ctx, lmax, l, ns, lflag)
      !
      !! This routine adjusts (modifies) the ns based on constraints
      !! If lflag=1, then copy ctx%inp%occupation to ns
      !! If lflag=2, then copy ns to ctx%inp%occupation
      !
      USE kinds,           ONLY : DP
      USE parameters,      ONLY : ntypx
      USE io_global,       ONLY : stdout
      USE ions_base,       ONLY : nat, ityp
      USE lsda_mod,        ONLY : nspin
      USE oscdft_context,  ONLY : oscdft_context_type
      !
      IMPLICIT NONE
      CLASS(oscdft_context_type), INTENT(INOUT) :: ctx
      INTEGER, INTENT(IN) :: lmax
      INTEGER, INTENT(IN) :: l(ntypx)
      REAL(DP), INTENT(INOUT) :: ns(2 * lmax + 1, 2 * lmax + 1, nspin, nat)
      INTEGER, INTENT(IN) :: lflag
      !
      INTEGER :: na, nt, ldim, is, m1, m2
      LOGICAL :: found
      !
      IF (.NOT.(ctx%inp%oscdft_type==2) .OR. .NOT.ctx%is_constraint) RETURN
      !
      found = .true.
      !
      WRITE(stdout, '(/5x,"Setting starting occupation matrices according to input constrained values")')
      !
      DO na = 1, nat
         IF (ctx%constraining(na)) THEN
            nt = ityp(na)
            ldim = 2*l(nt)+1
            DO is = 1, nspin
               DO m1 = 1, ldim
                  DO m2 = 1, ldim
                     IF (lflag==1) THEN
                        IF (ctx%inp%occupation(m1,m2,is,na)==-2.d0) THEN
                           WRITE(stdout, '(/5x,"Warning!!! Missing element: ",4(1x,i4))') na, is, m1, m2 
                           found = .false.
                        ELSE
                           ns(m1,m2,is,na) = ctx%inp%occupation(m1,m2,is,na)
                        ENDIF
                     ELSEIF (lflag==2) THEN
                        ctx%inp%occupation(m1,m2,is,na) = ns(m1,m2,is,na)
                     ENDIF
                  ENDDO
               ENDDO       
            ENDDO
         ENDIF
      ENDDO
      !
      IF (.NOT.found) CALL errore( 'oscdft_ns_set', 'Missing target occupation matrix element', 1 )
      !
      RETURN
      !
   END SUBROUTINE oscdft_ns_set

#endif
END MODULE oscdft_functions
