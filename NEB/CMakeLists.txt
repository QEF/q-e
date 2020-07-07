set(sources
    src/compute_scf.f90
    src/engine_to_path_pos.f90
    src/engine_to_path_alat.f90
    src/engine_to_path_nat.f90
    src/engine_to_path_fix_atom_pos.f90
    src/fcp_opt_routines.f90
    src/input.f90
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
    src/stop_run_path.f90
)

qe_add_library(qe_neb ${sources})
add_library(QE::NEB ALIAS qe_neb)
target_link_libraries(qe_neb
    PRIVATE
        QE::PW
        QE::Modules)

###########################################################
# neb.x
###########################################################
set(sources src/neb.f90)
qe_add_executable(qe_neb_exe ${sources})
set_target_properties(qe_neb_exe PROPERTIES OUTPUT_NAME neb.x)
target_link_libraries(qe_neb_exe
    PRIVATE
        QE::PW
        QE::LAX
        QE::UTILX
        QE::Modules
        QE::NEB)

###########################################################
# path_interpolation.x
###########################################################
set(sources src/path_interpolation.f90)
qe_add_executable(qe_neb_pathinterpolation_exe ${sources})
set_target_properties(qe_neb_pathinterpolation_exe PROPERTIES OUTPUT_NAME path_interpolation.x)
target_link_libraries(qe_neb_pathinterpolation_exe
    PRIVATE
        QE::Modules
        QE::NEB
        QE::UPF)

###########################################################

qe_install_targets(
    # Libraries
    qe_neb
    # Executables
    qe_neb_exe
    qe_neb_pathinterpolation_exe
)
