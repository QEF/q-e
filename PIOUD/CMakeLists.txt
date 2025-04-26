set(src_pioud
src/compute_scf_pioud.f90
src/engine_to_trpmd_pos.f90
src/engine_to_trpmd_nat_alat.f90
src/fcp_opt_routines.f90
src/fcp_variables.f90
src/trpmd_base.f90
src/ring_formats.f90
src/pw_gen_inputs.f90
src/ring_input_parameters_module.f90
src/trpmd_io_routines.f90
src/ring_io_units_module.f90
src/pimd_variables.f90
src/ring_variables.f90
src/pimd_subrout.f90
src/pimd_utils.f90
src/set_defaults.f90
src/stop_run_path.f90
)

qe_add_library(qe_pioud ${src_pioud})
target_link_libraries(qe_pioud
    PRIVATE
        qe_pw
        qe_lr_modules
        qe_modules
        qe_fftx
        qe_upflib
        qe_xclib)

###########################################################
# pioud.x
###########################################################
set(src_pioud_x src/pioud.f90)
qe_enable_cuda_fortran("${src_pioud_x}")
qe_add_executable(qe_pioud_exe ${src_pioud_x})
set_target_properties(qe_pioud_exe 
    PROPERTIES 
        OUTPUT_NAME pioud.x)
target_link_libraries(qe_pioud_exe
    PRIVATE
        qe_pw
        qe_lax
        qe_modules
        qe_pioud
        qe_xclib)

###########################################################

qe_install_targets(
    qe_pioud
    qe_pioud_exe)

add_custom_target(pioud
    DEPENDS
        qe_pioud_exe
    COMMENT
        "Path Integral Molecular Dynamics with PIOUD algorithm")
