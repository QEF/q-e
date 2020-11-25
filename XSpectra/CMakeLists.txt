set(sources_x
    src/gaunt_mod.f90
    src/xspectra_mod.f90
    src/radin_mod.f90
    src/mygetK.f90
    src/ipoolscatter.f90
    src/lr_sm1_psi.f90
    src/orthoUatwfc_k.f90
    src/read_k_points.f90
    src/reset_k_points_and_reinit.f90
    src/select_nl_init.f90
    src/lanczos.f90
    src/plot_xanes_cross_sections.f90
    src/xanes_dipole_general_edge.f90
    src/banner_xspectra.f90
    src/set_xspectra_namelists_defaults.f90
    src/read_input_and_bcast.f90
    src/stdout_routines.f90
    src/xanes_dipole.f90
    src/check_orthogonality_k_epsilon.f90
    src/assign_paw_radii_to_species.f90
    src/io_routines.f90
    src/xanes_quadrupole.f90)

set(sources_gipaw
    src/paw_gipaw.f90
    src/gipaw_module.f90
    src/init_gipaw_1.f90
    src/init_gipaw_2.f90)

qe_add_library(qe_xspectra_gipaw ${sources_gipaw})
add_library(QE::XSpectra::GIPAW ALIAS qe_xspectra_gipaw)
target_link_libraries(qe_xspectra_gipaw
    PRIVATE
        QE::PW
        QE::UTILX
        QE::Modules
        QE::UPF)

qe_add_library(qe_xspectra ${sources_x})
add_library(QE::XSpectra ALIAS qe_xspectra)
target_link_libraries(qe_xspectra
    PRIVATE
        QE::PW
        QE::UTILX
        QE::XSpectra::GIPAW
        QE::Modules
        QE::UPF)

###########################################################
# xspectra.x
###########################################################
set(sources src/xspectra.f90)
qe_add_executable(qe_xspectra_exe ${sources})
set_target_properties(qe_xspectra_exe PROPERTIES OUTPUT_NAME xspectra.x)
target_link_libraries(qe_xspectra_exe
    PRIVATE
        QE::PW
        QE::Modules
        QE::XSpectra
        QE::XSpectra::GIPAW
        QE::UPF
        QE::UTILX)

###########################################################
# spectra_correction.x
###########################################################
set(sources 
    src/spectra_correction.f90
    src/gaunt_mod.f90)
qe_add_executable(qe_xspectra_spectracorrection_exe ${sources})
set_target_properties(qe_xspectra_spectracorrection_exe PROPERTIES OUTPUT_NAME spectra_correction.x)
target_link_libraries(qe_xspectra_spectracorrection_exe
    PRIVATE
        QE::Modules
        QE::XSpectra
        QE::XSpectra::GIPAW
        QE::UTILX)

###########################################################
# molecularnexafs.x
###########################################################
set(sources src/molecularnexafs.f90)
qe_add_executable(qe_xspectra_molecularnexafs_exe ${sources})
set_target_properties(qe_xspectra_molecularnexafs_exe PROPERTIES OUTPUT_NAME molecularnexafs.x)
target_link_libraries(qe_xspectra_molecularnexafs_exe
    PRIVATE
        QE::XSpectra
        QE::XSpectra::GIPAW)

###########################################################

qe_install_targets(
    # Libraries
    qe_xspectra_gipaw
    qe_xspectra
    # Executables
    qe_xspectra_exe
    qe_xspectra_spectracorrection_exe
    qe_xspectra_molecularnexafs_exe)
