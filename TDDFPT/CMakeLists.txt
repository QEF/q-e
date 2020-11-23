set(sources
    src/lr_variables.f90
    src/lr_charg_resp.f90
    src/bcast_lr_input.f90
    src/lr_readin.f90
    src/lr_alloc_init.f90
    src/lr_calc_dens.f90
    src/lr_dealloc.f90
    src/lr_ortho.f90
    src/lr_read_wf.f90
    src/lr_normalise.f90
    src/lr_lanczos.f90
    src/lr_apply_liouvillian.f90
    src/lr_dv_setup.f90
    src/lr_compute_intq.f90
    src/lr_set_intq_nc.f90
    src/lr_transform_intq_nc.f90
    src/lr_transform_intq_so.f90
    src/lr_solve_e.f90
    src/lr_dvpsi_e.f90
    src/stop_lr.f90
    src/lr_read_d0psi.f90
    src/lr_restart.f90
    src/lr_write_restart.f90
    src/print_clock_lr.f90
    src/sd0psi.f90
    src/lr_set_boxes_density.f90
    src/lr_init_nfo.f90
    src/lr_exx_kernel.f90
    src/lr_calc_dens_eels.f90
    src/lr_calc_dens_eels_nc.f90
    src/lr_apply_liouvillian_eels.f90
    src/lr_addus_dvpsi.f90
    src/lr_run_nscf.f90
    src/lr_setup_nscf.f90
    src/lr_sym_eels.f90
    src/lr_smallgq.f90
    src/lr_psym_eels.f90
    src/lr_summary.f90
    src/lr_dvpsi_eels.f90
    src/lr_dav_variables.f90
    src/lr_dav_routines.f90
    src/lr_us.f90
    src/lr_test_restart.f90
    src/lr_dav_debug.f90
    src/plugin_tddfpt_potential.f90
    src/lr_sternheimer.f90
    src/linear_solvers.f90
    src/paw_add_symmetry.f90
    src/dveqpsi_us.f90
    src/dveqpsi_us_only.f90
    src/orthogonalize_omega.f90
    src/mix_pot.f90
)

qe_add_library(qe_tddfpt ${sources})
add_library(QE::TDDFPT ALIAS qe_tddfpt)
target_link_libraries(qe_tddfpt
    PRIVATE
        QE::LR_Modules
        QE::PW
        QE::Modules)

###########################################################
# turbo_lanczos.x
###########################################################
set(sources src/lr_main.f90)
qe_add_executable(qe_tddfpt_turbolanczos_exe ${sources})
set_target_properties(qe_tddfpt_turbolanczos_exe PROPERTIES OUTPUT_NAME turbo_lanczos.x)
target_link_libraries(qe_tddfpt_turbolanczos_exe
    PRIVATE
        QE::PW
        QE::Modules
        QE::TDDFPT)

###########################################################
# turbo_davidson.x
###########################################################
set(sources src/lr_dav_main.f90)
qe_add_executable(qe_tddfpt_turbodavidson_exe ${sources})
set_target_properties(qe_tddfpt_turbodavidson_exe PROPERTIES OUTPUT_NAME turbo_davidson.x)
target_link_libraries(qe_tddfpt_turbodavidson_exe
    PRIVATE
        QE::PW
        QE::Modules
        QE::TDDFPT)

###########################################################
# turbo_eels.x
###########################################################
set(sources src/lr_eels_main.f90)
qe_add_executable(qe_tddfpt_turboeels_exe ${sources})
set_target_properties(qe_tddfpt_turboeels_exe PROPERTIES OUTPUT_NAME turbo_eels.x)
target_link_libraries(qe_tddfpt_turboeels_exe
    PRIVATE
        QE::PW
        QE::LAX
        QE::Modules
        QE::TDDFPT
        QE::LR_Modules)

###########################################################
# tddfpt_calculate_spectrum.x
###########################################################
set(sources tools/tddfpt_calculate_spectrum.f90)
qe_add_executable(qe_tddfpt_calculatespectrum_exe ${sources})
set_target_properties(qe_tddfpt_calculatespectrum_exe PROPERTIES OUTPUT_NAME tddfpt_calculate_spectrum.x)
target_link_libraries(qe_tddfpt_calculatespectrum_exe
    PRIVATE
        QE::LAPACK
        QE::UTILX
        QE::Modules)

###########################################################

qe_install_targets(
    # Libraries
    qe_tddfpt
    # Executables
    qe_tddfpt_turbolanczos_exe
    qe_tddfpt_turbodavidson_exe
    qe_tddfpt_turboeels_exe
    qe_tddfpt_calculatespectrum_exe
)