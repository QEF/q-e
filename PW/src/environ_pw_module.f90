!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE environ_pw_module
    !------------------------------------------------------------------------------------
#if defined (__ENVIRON)
    !
    USE environ_base_module, ONLY: environ
    !
    USE global_version, ONLY: version_number
    !
    USE mp, ONLY: mp_bcast
    USE mp_images, ONLY: intra_image_comm
    USE io_global, ONLY: stdout, ionode_id
    USE kinds, ONLY: DP
    !
    USE control_flags, ONLY: lscf, lbfgs, conv_ions, istep, nstep, lforce => tprnfor
    !
    USE fft_base, ONLY: dfftp
    !
    USE scf, ONLY: scf_type
    USE klist, ONLY: nelec, tot_charge
    USE lsda_mod, ONLY: nspin
    !
    USE ener, ONLY: ef
    USE constants, ONLY: rytoev
    USE ions_base, ONLY: nat, ityp, zv
    USE extrapolation, ONLY: update_pot
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: update_environ_potential, calc_environ_potential
    !
    PUBLIC :: is_ms_gcs, init_ms_gcs, run_ms_gcs
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_potential(vltot)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: vltot(dfftp%nnr)
        !
        CHARACTER(LEN=80) :: sub_name = 'update_environ_potential'
        !
        !--------------------------------------------------------------------------------
        !
        CALL environ%main%update_potential(dfftp%nnr, vltot)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_potential
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_environ_potential(rhoin, converged, dr2, vltot)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(scf_type), INTENT(IN) :: rhoin
        LOGICAL, INTENT(IN) :: converged
        REAL(DP), INTENT(IN) :: dr2
        !
        REAL(DP), INTENT(INOUT) :: vltot(dfftp%nnr)
        !
        LOGICAL :: update_venviron = .FALSE.
        INTEGER :: local_verbose = 0
        !
        REAL(DP) :: rhoaux(dfftp%nnr)
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_environ_potential'
        !
        !--------------------------------------------------------------------------------
        ! Reduce output at each scf iteration
        !
        IF (.NOT. lscf .OR. converged) local_verbose = 1
        !
        !--------------------------------------------------------------------------------
        ! update electrons-related quantities in environ
        !
        rhoaux = rhoin%of_r(:, 1)
        !
        IF (version_number == '6.3') THEN
            IF (nspin == 2) rhoaux = rhoaux + rhoin%of_r(:, 2)
        END IF
        !
        CALL environ%main%update_electrons(dfftp%nnr, rhoaux, nelec)
        !
        !--------------------------------------------------------------------------------
        ! Environ contribution to the local potential
        !
        IF (dr2 > 0.0_DP) THEN
            update_venviron = .NOT. converged .AND. dr2 < environ%setup%get_threshold()
        ELSE
            update_venviron = environ%setup%is_restart() .OR. environ%setup%is_tddfpt()
            !
            !----------------------------------------------------------------------------
            ! For subsequent steps of optimization or dynamics, compute Environ
            ! contribution during initialization
            !
            CALL environ%setup%set_restart(.TRUE.)
            !
        END IF
        !
        IF (update_venviron) WRITE (stdout, 1000)
        !
        CALL environ%calc%potential(update_venviron, local_verbose)
        !
        vltot = environ%main%get_vzero(dfftp%nnr) + environ%main%get_dvtot(dfftp%nnr)
        !
        IF (.NOT. lscf .OR. converged) CALL environ%main%print_potential_shift()
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 5X, "add environment contribution to local potential")
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_environ_potential
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    LOGICAL FUNCTION is_ms_gcs()
        !--------------------------------------------------------------------------------
        !
        is_ms_gcs = environ%setup%lmsgcs
        !
        !--------------------------------------------------------------------------------
    END FUNCTION is_ms_gcs
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_ms_gcs()
        !--------------------------------------------------------------------------------
        !
        CALL start_clock("semiconductor")
        !
        lforce = .TRUE.
        lbfgs = .FALSE.
        nstep = 100
        tot_charge = 0.0
        !
        WRITE (stdout, 1001)
        WRITE (stdout, 1002) tot_charge
        !
        CALL stop_clock("semiconductor")
        !
        !--------------------------------------------------------------------------------
        !
1001    FORMAT(5X, //"*******************************************"//, &
               "  Please cite                              "//, &
               "  Q. Campbell, D. Fisher and I. Dabo, Phys. Rev. Mat. 3, 015404 (2019)."//, &
               "  doi: 10.1103/PhysRevMaterials.3.015404   "//, &
               "  In any publications resulting from this work.")
        !
1002    FORMAT(5X, //"*******************************************"//, &
               "     Running initial calculation for flatband."// &
               "     Using charge of: ", F14.8, // &
               "*******************************************")
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_ms_gcs
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE run_ms_gcs()
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP) :: cur_chg
        REAL(DP) :: prev_chg, prev_chg2
        REAL(DP) :: cur_dchg
        REAL(DP) :: cur_fermi
        REAL(DP) :: prev_dchg
        REAL(DP) :: gamma_mult
        REAL(DP) :: prev_step_size
        REAL(DP) :: ss_chg, charge
        INTEGER :: chg_step, na
        REAL(DP) :: surf_area
        REAL(DP) :: chg_per_area
        REAL(DP) :: ss_chg_per_area
        REAL(DP) :: ss_potential, total_potential
        REAL(DP) :: dft_chg_max, dft_chg_min
        REAL(DP) :: change_vec
        REAL(DP) :: v_cut, bulk_potential
        REAL(DP) :: ionic_charge
        LOGICAL :: converge
        !
        CHARACTER(LEN=80) :: sub_name = 'run_ms_gcs'
        !
        !--------------------------------------------------------------------------------
        !
        gamma_mult = 0.15
        converge = .TRUE.
        !
        !--------------------------------------------------------------------------------
        ! Calculating ionic charge ! TODO why is this before the start_clock?
        !
        ionic_charge = 0._DP
        !
        DO na = 1, nat
            ionic_charge = ionic_charge + zv(ityp(na))
        END DO
        !
        !--------------------------------------------------------------------------------
        !
        CALL start_clock("semiconductor")
        !
        !--------------------------------------------------------------------------------
        !
        ASSOCIATE (sc => environ%main%semiconductor%base)
            !
            !----------------------------------------------------------------------------
            ! Initializing the constraints of possible DFT charges
            !
            chg_step = istep
            !
            IF (chg_step == 1) THEN
                !
                IF (sc%electrode_charge > 0.0) THEN
                    dft_chg_max = sc%electrode_charge
                    dft_chg_min = 0.0
                ELSE
                    dft_chg_min = sc%electrode_charge
                    dft_chg_max = 0.0
                END IF
                !
            END IF
            !
            !----------------------------------------------------------------------------
            ! After first scf step, extract fermi level and set charge
            !
            IF (chg_step == 0) THEN
                sc%flatband_fermi = ef
                tot_charge = 0.7 * sc%electrode_charge
                sc%slab_charge = tot_charge
                conv_ions = .FALSE.
                nelec = ionic_charge - tot_charge
                !
                WRITE (stdout, 1003) &
                    sc%flatband_fermi * rytoev, tot_charge
                !
                !------------------------------------------------------------------------
                ! Trigger next SCF cycle
                !
                istep = istep + 1
                !
                CALL update_pot()
                !
                CALL hinit1()
                !
            ELSE
                !
                !------------------------------------------------------------------------
                ! Calculate steepest descent for changing charge at each scf step
                !
                cur_fermi = ef
                cur_dchg = sc%bulk_sc_fermi - cur_fermi
                bulk_potential = (sc%bulk_sc_fermi - sc%flatband_fermi) * rytoev
                ss_chg = tot_charge
                !
                !------------------------------------------------------------------------
                ! Making sure constraints are updated - stop if...
                ! 1) electrode chg positive and dft is negative
                ! 2) electrode chg negative and dft is positive
                !
                IF (sc%electrode_charge > 0) THEN
                    !
                    IF (ss_chg < 0.0) THEN
                        dft_chg_min = tot_charge
                        converge = .FALSE.
                    ELSE
                        prev_chg2 = tot_charge
                    END IF
                    !
                ELSE
                    !
                    IF (ss_chg > 0.0) THEN
                        dft_chg_max = tot_charge
                        converge = .FALSE.
                    ELSE
                        prev_chg2 = tot_charge
                    END IF
                    !
                END IF
                !
                CALL mp_bcast(dft_chg_min, ionode_id, intra_image_comm)
                !
                CALL mp_bcast(dft_chg_max, ionode_id, intra_image_comm)
                !
                !------------------------------------------------------------------------
                ! Updating the steepest descent parameter, gamma_mult
                !
                IF (chg_step > 1) THEN
                    gamma_mult = (cur_chg - prev_chg) / (cur_dchg - prev_dchg)
                END IF
                !
                change_vec = -gamma_mult * cur_dchg
                prev_chg = tot_charge
                !
                !------------------------------------------------------------------------
                ! Updating tot_charge within constraints
                !
                IF ((tot_charge + change_vec) > dft_chg_max) THEN
                    !
                    IF (tot_charge >= dft_chg_max) THEN
                        tot_charge = prev_chg2 + 0.7 * (dft_chg_max - prev_chg2)
                    ELSE
                        tot_charge = tot_charge + 0.7 * (dft_chg_max - tot_charge)
                    END IF
                    !
                ELSE IF ((tot_charge + change_vec) < dft_chg_min) THEN
                    !
                    IF (tot_charge <= dft_chg_min) THEN
                        tot_charge = prev_chg2 - 0.7 * (prev_chg2 - dft_chg_min)
                    ELSE
                        tot_charge = tot_charge - 0.7 * (tot_charge - dft_chg_min)
                    END IF
                    !
                ELSE
                    tot_charge = tot_charge + change_vec
                END IF
                !
                CALL mp_bcast(tot_charge, ionode_id, intra_image_comm)
                !
                !------------------------------------------------------------------------
                ! Updating variables based on new_tot_charge
                !
                cur_chg = tot_charge
                prev_step_size = ABS(cur_chg - prev_chg)
                prev_dchg = cur_dchg
                !
                CALL mp_bcast(converge, ionode_id, intra_image_comm)
                !
                CALL mp_bcast(prev_step_size, ionode_id, intra_image_comm)
                !
                !------------------------------------------------------------------------
                ! Check charge convergence
                !
                IF (((prev_step_size > sc%charge_threshold) .OR. (.NOT. converge)) .AND. &
                    (chg_step < nstep - 1)) THEN
                    !
                    !--------------------------------------------------------------------
                    ! Not converged
                    !
                    conv_ions = .FALSE.
                    !
                    WRITE (STDOUT, 1004) &
                        chg_step, cur_fermi * rytoev, &
                        ss_chg, &
                        prev_step_size, &
                        cur_dchg, &
                        tot_charge
                    !
                    sc%slab_charge = tot_charge
                    nelec = ionic_charge - tot_charge
                    !
                    CALL mp_bcast(nelec, ionode_id, intra_image_comm)
                    !
                    !--------------------------------------------------------------------
                    ! Trigger next SCF cycle
                    !
                    istep = istep + 1
                    !
                    CALL update_pot()
                    !
                    CALL hinit1()
                    !
                ELSE
                    !
                    !--------------------------------------------------------------------
                    ! Converged
                    !
                    IF (chg_step == nstep - 1) THEN
                        !
                        !----------------------------------------------------------------
                        ! Case where about to exceed max number of steps
                        !
                        WRITE (stdout, *) &
                            NEW_LINE("a")//"   Exceeded Max number steps!"// &
                            NEW_LINE("a")//"   Results probably out of accurate range"// &
                            NEW_LINE("a")//"   Smaller chg_thr recommended."// &
                            NEW_LINE("a")//"   Writing current step to q-v.dat."
                        !
                    END IF
                    !
                    !--------------------------------------------------------------------
                    ! Writing output for successful charge convergance
                    !
                    WRITE (stdout, 1005) &
                        chg_step, prev_step_size, ss_chg, cur_dchg, bulk_potential
                    !
                    OPEN (21, file="q-v.dat", status="unknown")
                    !
                    WRITE (21, *) &
                        & "Potential (V-V_fb)  Surface State Potential (V-V_cut)", &
                        & "  Electrode Charge (e)", &
                        & "  Surface States Charge (e)    ", &
                        & "Electrode Charge per surface area (e/cm^2)     ", &
                        & "Surface State Charge per surface area (e/cm^2)"
                    !
                    !--------------------------------------------------------------------
                    ! TODO documentation
                    !
                    surf_area = sc%surf_area_per_sq_cm
                    chg_per_area = sc%electrode_charge / surf_area
                    ss_chg_per_area = ss_chg / surf_area
                    ss_potential = -bulk_potential
                    !
                    CALL mp_bcast(ss_potential, ionode_id, intra_image_comm)
                    !
                    WRITE (21, 1006) &
                        total_potential, ss_potential, &
                        sc%electrode_charge, ss_chg, &
                        chg_per_area, ss_chg_per_area
                    !
                    CLOSE (21)
                    !
                END IF
                !
            END IF
            !
        END ASSOCIATE
        !
        !--------------------------------------------------------------------------------
        !
        CALL stop_clock("semiconductor")
        !
        !--------------------------------------------------------------------------------
        !
1003    FORMAT(5X, //"***************************************************", // &
               "     Flatband potential calculated as ", F14.8, // &
               "     Now using initial charge of:  ", F14.8, // &
               "***************************************************")
!
1004    FORMAT(5X, //"***************************************************", // &
               "     Finished Charge convergence step : ", I3, // &
               "     DFT Fermi level calculated as ", F14.8, // &
               "     Charge trapped in surface states: ", F14.8, " e", // &
               "     Charge Accuracy < ", F14.8, // &
               "     Difference between bulk and DFT fermi: ", F14.8, // &
               "     Now using DFT charge of:  ", F14.8, // &
               "***************************************************")
1005    FORMAT(5X, //"***************************************************", // &
               "     Finished charge convergence step : ", I3, // &
               "     Convergence of charge with accuracy < ", F14.8, " e", // &
               "     Charge trapped in surface states: ", F14.8, // &
               "     Difference between bulk and DFT fermi: ", F14.8, // &
               "     Final Potential: ", F14.8, " V", // &
               "     Output written to q-v.dat       ", // &
               "***************************************************")
1006    FORMAT(1X, 4F14.8, 2ES12.5)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE run_ms_gcs
    !------------------------------------------------------------------------------------
    !
#endif
    !------------------------------------------------------------------------------------
END MODULE environ_pw_module
!----------------------------------------------------------------------------------------
