set(src_tddfpt
    src/environ_td_module.f90
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
    src/lr_mix_pot.f90
    src/lr_magnons_routines.f90
    src/lr_apply_liouvillian_magnons.f90
    src/lr_dvpsi_magnons.f90
    src/lr_Opsi_magnons.f90
    src/lr_calc_dens_magnons.f90)
qe_enable_cuda_fortran("${src_tddfpt}")


qe_add_library(qe_tddfpt ${src_tddfpt})
target_link_libraries(qe_tddfpt
    PRIVATE
        qe_hdf5_fortran
        qe_openmp_fortran
        qe_openacc_fortran
        qe_mpi_fortran
        qe_lr_modules
        qe_pw
        qe_modules
        qe_upflib
        qe_fftx
        qe_xclib)

###########################################################
# turbo_lanczos.x
###########################################################
set(src_turbo_lanczos_x src/lr_main.f90)
qe_enable_cuda_fortran("${src_turbo_lanczos_x}")
qe_add_executable(qe_tddfpt_turbolanczos_exe ${src_turbo_lanczos_x})
set_target_properties(qe_tddfpt_turbolanczos_exe 
    PROPERTIES 
        OUTPUT_NAME turbo_lanczos.x)
target_link_libraries(qe_tddfpt_turbolanczos_exe
    PRIVATE
        qe_openacc_fortran
        qe_openmp_fortran
        qe_mpi_fortran
        qe_pw
        qe_modules
        qe_tddfpt
        qe_upflib
        qe_fftx
        qe_xclib)

###########################################################
# turbo_davidson.x
###########################################################
set(src_turbo_davidson_x src/lr_dav_main.f90)
qe_enable_cuda_fortran("${src_turbo_davidson_x}")
qe_add_executable(qe_tddfpt_turbodavidson_exe ${src_turbo_davidson_x})
set_target_properties(qe_tddfpt_turbodavidson_exe 
    PROPERTIES 
        OUTPUT_NAME turbo_davidson.x)
target_link_libraries(qe_tddfpt_turbodavidson_exe
    PRIVATE
        qe_openacc_fortran
        qe_openmp_fortran
        qe_mpi_fortran
        qe_pw
        qe_modules
        qe_tddfpt
        qe_upflib
        qe_fftx
        qe_xclib)

###########################################################
# turbo_eels.x
###########################################################
set(src_turbo_eels_x src/lr_eels_main.f90)
qe_enable_cuda_fortran("${src_turbo_eels_x}")
qe_add_executable(qe_tddfpt_turboeels_exe ${src_turbo_eels_x})
set_target_properties(qe_tddfpt_turboeels_exe 
    PROPERTIES 
        OUTPUT_NAME turbo_eels.x)
target_link_libraries(qe_tddfpt_turboeels_exe
    PRIVATE
        qe_openacc_fortran
        qe_openmp_fortran
        qe_mpi_fortran
        qe_pw
        qe_lax
        qe_modules
        qe_tddfpt
        qe_lr_modules
        qe_upflib
        qe_fftx
        qe_xclib)

###########################################################
# turbo_magnon.x
###########################################################
set(src_turbomagnon_x src/lr_magnons_main.f90)
# Because this main file includes laxlib.fh which is preprocessed with __CUDA definition
# the main file needs to be compiled with consistent definitions.
qe_enable_cuda_fortran(src/lr_magnons_main.f90)
qe_add_executable(qe_tddfpt_turbomagnon_exe ${src_turbomagnon_x})
set_target_properties(qe_tddfpt_turbomagnon_exe
    PROPERTIES
        OUTPUT_NAME turbo_magnon.x)
target_link_libraries(qe_tddfpt_turbomagnon_exe
    PRIVATE
        qe_openacc_fortran
        qe_openmp_fortran
        qe_mpi_fortran
        qe_tddfpt
        qe_lr_modules
        qe_pw
        qe_lax
        qe_modules
        qe_fftx)

###########################################################
# turbo_spectrum.x
###########################################################
set(src_tddfpt_turbospectrum_x src/turbo_spectrum.f90)
qe_add_executable(qe_tddfpt_turbospectrum_exe ${src_tddfpt_turbospectrum_x})
set_target_properties(qe_tddfpt_turbospectrum_exe
    PROPERTIES
        OUTPUT_NAME turbo_spectrum.x)
target_link_libraries(qe_tddfpt_turbospectrum_exe
    PRIVATE
        qe_modules
        qe_xclib
        qe_lapack)

###########################################################

set(TDDFPT_EXE_TARGETS
    qe_tddfpt_turbolanczos_exe
    qe_tddfpt_turbodavidson_exe
    qe_tddfpt_turboeels_exe
    qe_tddfpt_turbomagnon_exe
    qe_tddfpt_turbospectrum_exe)

qe_install_targets(
    # Libraries
    qe_tddfpt
    # Executables
    ${TDDFPT_EXE_TARGETS})

add_custom_target(tddfpt
    DEPENDS
        ${TDDFPT_EXE_TARGETS}
    COMMENT
        "time dependent dft code")
