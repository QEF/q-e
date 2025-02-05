set(src_kcw
    wann2kcw.f90
    kcw_screen.f90
    kcw_ham.f90
    kcw_comm.f90
    print_clock_kcw.f90
    kcw_readin.f90
    kcw_setup.f90
    kcw_setup_screen.f90
    kcw_setup_ham.f90
    bcast_kcw_input.f90
    input_summary.f90
    read_wannier.f90
    rotate_ks.f90
    apply_u_matrix.f90
    ks_hamiltonian.f90
    screen_coeff.f90
    compute_map_ikq_single.f90
    calculate_phase.f90
    structure_factor.f90
    rho_of_q.f90
    kcw_prepare_q.f90
    bare_pot.f90
    solve_linter_koop_mod.f90
    kcw_initialize_ph.f90
    kcw_openfilq.f90
    kcw_q_setup.f90
    clean_pw_kcw.f90
    close_kcw.f90
    kcw_R_points.f90
    ham_R0_2nd.f90
    hamilt.f90
    read_alpha.f90
    kcw_io_new.f90
    bcast_wfc.f90
    interpolation.f90
    convert_kpts_names.f90
    write_hr_to_file.f90
    alpha_corr.f90
    koopmans_ham.f90
    full_ham.f90
    self_hartree.f90
    group_orbitals.f90
    kcw_run_nscf.f90
    kcw_init_q.f90
    kcw_allocate_q.f90
    kcw_deallocate_q.f90
    setup_coulomb_exx.f90
    coulomb.f90
    find_IBZ_q.f90 
    kcw_kpoint_grid.f90 
    kcw_deallocate_symmetry_arrays.f90 
    write_qlist_ibz.f90 
    read_qlist_ibz.f90 
    kcw_read_sym.f90 
    kcw_write_sym.f90 
    )

qe_add_library(qe_kcw ${src_kcw})
target_link_libraries(
    qe_kcw
    PRIVATE qe_pw
            qe_lr_modules
            qe_modules
            qe_fftx
            qe_upflib
            qe_xclib)

###########################################################
# kcw.x
###########################################################
set(src_kcw_x kcw.f90)
qe_add_executable(qe_kcw_exe ${src_kcw_x})
set_target_properties(qe_kcw_exe PROPERTIES OUTPUT_NAME kcw.x)
target_link_libraries(qe_kcw_exe PRIVATE qe_pw qe_lax qe_modules qe_kcw qe_xclib)
