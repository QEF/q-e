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


add_subdirectory(src)

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

set (src_kcwpp 
     PP/bcast_kcw_pp_input.f90
     PP/input_pp_summary.f90
     PP/kcw_pp_readin.f90
     PP/read_hr.f90
     PP/kcw_bands.f90
     PP/sh_setup.f90)

add_subdirectory(PP)

qe_add_library(qe_kcwpp ${src_kcwpp})
target_link_libraries(qe_kcwpp
    PRIVATE
        qe_pw
        qe_lax
        qe_lr_modules
        qe_modules
        qe_fftx
        qe_upflib
        qe_xclib
        qe_kcw)


###########################################################

qe_install_targets(
    qe_kcw
    qe_kcw_exe
    qe_kcwpp_sh_exe
    qe_kcwpp_interp_exe)


add_custom_target(kcw
    DEPENDS
        pw
        pp
        qe_wannierprog_exe
        qe_kcw_exe
        qe_kcwpp_sh_exe
        qe_kcwpp_interp_exe
    COMMENT
        "Koopmans spectral functionals")
