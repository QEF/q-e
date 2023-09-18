set(src_neb
    src/bcast_file.f90
    src/compute_scf.f90
    src/engine_to_path_pos.f90
    src/engine_to_path_alat.f90
    src/engine_to_path_nat.f90
    src/engine_to_path_fix_atom_pos.f90
    src/engine_to_path_tot_charge.f90
    src/fcp_opt_routines.f90
    src/fcp_variables.f90
    src/gcscf_variables.f90
    src/neb_input.f90
    src/path_base.f90
    src/path_formats.f90
    src/path_gen_inputs.f90
    src/path_input_parameters_module.f90
    src/path_io_routines.f90
    src/path_io_tools.f90
    src/path_io_units_module.f90
    src/path_opt_routines.f90
    src/path_reparametrisation.f90
    src/path_read_cards_module.f90
    src/path_read_namelists_module.f90
    src/path_to_engine_fix_atom_pos.f90
    src/path_variables.f90
    src/set_defaults.f90
    src/stop_run_path.f90)

qe_add_library(qe_neb ${src_neb})
target_link_libraries(qe_neb
    PRIVATE
        qe_pw
        qe_modules
        qe_upflib
        qe_fftx
        qe_xclib)

###########################################################
# neb.x
###########################################################
set(src_neb_x src/neb.f90)
qe_enable_cuda_fortran("${src_neb_x}")
qe_add_executable(qe_neb_exe ${src_neb_x})
set_target_properties(qe_neb_exe 
    PROPERTIES 
        OUTPUT_NAME neb.x)
target_link_libraries(qe_neb_exe
    PRIVATE
        qe_pw
        qe_lax
        qe_modules
        qe_neb)

###########################################################
# path_interpolation.x
###########################################################
set(src_path_interpolation_x src/path_interpolation.f90)
qe_add_executable(qe_neb_pathinterpolation_exe ${src_path_interpolation_x})
set_target_properties(qe_neb_pathinterpolation_exe 
    PROPERTIES 
        OUTPUT_NAME path_interpolation.x)
target_link_libraries(qe_neb_pathinterpolation_exe
    PRIVATE
        qe_modules
        qe_neb
        qe_upflib)

###########################################################

qe_install_targets(
    # Libraries
    qe_neb
    # Executables
    qe_neb_exe
    qe_neb_pathinterpolation_exe)

add_custom_target(neb
    DEPENDS
        qe_neb_exe
        qe_neb_pathinterpolation_exe
    COMMENT
        "code for Nudged Elastic Band method")
