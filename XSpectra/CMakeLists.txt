set(src_xspectra
    src/xspectra.f90
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

set(src_gipaw
    src/paw_gipaw.f90
    src/gipaw_module.f90
    src/init_gipaw_1.f90
    src/init_gipaw_2.f90)

qe_add_library(qe_xspectra_gipaw ${src_gipaw})
target_link_libraries(qe_xspectra_gipaw
    PRIVATE
        qe_pw
        qe_modules
        qe_upflib
        qe_fftx
        qe_xclib)

qe_add_library(qe_xspectra ${src_xspectra})
target_link_libraries(qe_xspectra
    PRIVATE
        qe_pw
        qe_xspectra_gipaw
        qe_modules
        qe_upflib
        qe_fftx
        qe_xclib)

###########################################################
# xspectra.x
###########################################################
set(sources src/xspectra.f90)
qe_add_executable(qe_xspectra_exe ${sources})
set_target_properties(qe_xspectra_exe 
    PROPERTIES 
        OUTPUT_NAME xspectra.x)
target_link_libraries(qe_xspectra_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_xspectra
        qe_xspectra_gipaw
        qe_upflib)

###########################################################
# spectra_correction.x
###########################################################
set(src_spectra_correction_x src/spectra_correction.f90)
qe_add_executable(qe_xspectra_spectracorrection_exe ${src_spectra_correction_x})
set_target_properties(qe_xspectra_spectracorrection_exe 
    PROPERTIES 
        OUTPUT_NAME spectra_correction.x)
target_link_libraries(qe_xspectra_spectracorrection_exe
    PRIVATE
        qe_modules
        qe_xspectra
        qe_xspectra_gipaw
        qe_pw
        qe_upflib)

###########################################################
# molecularnexafs.x
###########################################################
set(src_molecularnexafs_x src/molecularnexafs.f90)
qe_add_executable(qe_xspectra_molecularnexafs_exe ${src_molecularnexafs_x})
set_target_properties(qe_xspectra_molecularnexafs_exe 
    PROPERTIES 
        OUTPUT_NAME molecularnexafs.x)
target_link_libraries(qe_xspectra_molecularnexafs_exe
    PRIVATE
        qe_xspectra
        qe_xspectra_gipaw)

###########################################################

qe_install_targets(
    # Libraries
    qe_xspectra_gipaw
    qe_xspectra
    # Executables
    qe_xspectra_exe
    qe_xspectra_spectracorrection_exe
    qe_xspectra_molecularnexafs_exe)
