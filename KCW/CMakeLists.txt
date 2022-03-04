set(src_kcw
    src/wann2kcw.f90
    src/kcw_screen.f90
    src/kcw_ham.f90
    src/kcw_comm.f90 
    src/print_clock_kcw.f90 
    src/kcw_readin.f90 
    src/kcw_setup.f90 
    src/kcw_setup_screen.f90 
    src/kcw_setup_ham.f90 
    src/bcast_kcw_input.f90 
    src/input_summary.f90 
    src/read_wannier.f90 
    src/rotate_ks.f90 
    src/apply_u_matrix.f90 
    src/ks_hamiltonian.f90 
    src/screen_coeff.f90 
    src/compute_map_ikq_single.f90 
    src/calculate_phase.f90 
    src/structure_factor.f90 
    src/rho_of_q.f90 
    src/kcw_prepare_q.f90 
    src/bare_pot.f90 
    src/solve_linter_koop_mod.f90 
    src/kcw_initialize_ph.f90 
    src/kcw_openfilq.f90 
    src/kcw_q_setup.f90 
    src/clean_pw_kcw.f90 
    src/close_kcw.f90 
    src/kcw_R_points.f90
    src/ham_R0_2nd.f90
    src/hamilt.f90
    src/read_alpha.f90
    src/kcw_io_new.f90
    src/bcast_wfc.f90
    src/interpolation.f90
    src/convert_kpts_names.f90
    src/write_hr_to_file.f90
    src/alpha_corr.f90 
    src/koopmans_ham.f90 
    src/full_ham.f90 
    src/self_hartree.f90
    src/group_orbitals.f90
    src/kcw_run_nscf.f90
    src/kcw_init_q.f90
    src/kcw_allocate_q.f90
    src/kcw_deallocate_q.f90
    src/setup_coulomb_exx.f90
    src/coulomb.f90)


qe_add_library(qe_kcw ${src_kcw})
target_link_libraries(qe_kcw
    PRIVATE
        qe_pw
        qe_lax
        qe_lr_modules
        qe_modules
        qe_fftx
        qe_upflib
        qe_xclib)

###########################################################
# kcw.x
###########################################################
set(src_kcw_x src/kcw.f90)
qe_enable_cuda_fortran("${src_kcw_x}")
qe_add_executable(qe_kcw_exe ${src_kcw_x})
set_target_properties(qe_kcw_exe 
    PROPERTIES 
        OUTPUT_NAME kcw.x)
target_link_libraries(qe_kcw_exe
    PRIVATE
        qe_pw
        qe_lax
        qe_modules
        qe_kcw
        qe_xclib)

###########################################################

qe_install_targets(
    qe_kcw
    qe_kcw_exe)

add_custom_target(kcw
    DEPENDS
        pw
        pp
#        qe_pp_pw2wannier90_exe # compile only pw2wannier
        qe_kcw_exe
    COMMENT
        "Koopmans spectral functionals")
