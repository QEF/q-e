set(sources
    src/epwcom.f90
    src/constants_epw.f90
    src/elph2.f90
    src/io_var.f90
    src/low_lvl.f90
    src/division.f90
    src/rigid_epw.f90
    src/io_epw.f90
    src/io_transport.f90
    src/io_selfen.f90
    src/wigner.f90
    src/wan2bloch.f90
    src/cum_mod.f90
    src/close_epw.f90
    src/poolgathering.f90
    src/printing.f90
    src/kfold.f90
    src/dynmat_asr.f90
    src/io_eliashberg.f90
    src/utilities.f90
    src/eliashbergcom.f90
    src/supercond.f90
    src/supercond_iso.f90
    src/supercond_aniso.f90
    src/grid.f90
    src/rotate.f90
    src/transport.f90
    src/transport_iter.f90
    src/pw2wan2epw.f90
    src/wannierization.f90
    src/indabs.f90
    src/bcast_epw_input.f90
    src/bloch2wan.f90
    src/dvqpsi.f90
    src/eliashberg.f90
    src/elphel2_shuffle.f90
    src/elphon_shuffle.f90
    src/elphon_shuffle_wrap.f90
    src/selfen.f90
    src/spectral_func.f90
    src/ephwann_shuffle.f90
    src/epw_init.f90
    src/epw_readin.f90
    src/epw_setup.f90
    src/epw_summary.f90
    src/loadumat.f90
    src/stop_epw.f90
    src/wannierEPW.f90
    src/loadbm.f90
    src/bfield.f90
    src/io_indabs.f90
    src/polaron.f90
    src/qdabs.f90
    src/wfpt_mod.f90
    src/rotate_wavefunction.f90
    src/selfen_ahc.f90
    src/read_ahc_files.f90
    )

qe_add_library(qe_epw ${sources})
target_link_libraries(
    qe_epw
    PRIVATE qe_wannier90
            qe_modules
            qe_lr_modules
            qe_phonon_ph
            qe_pw
            qe_mpi_fortran
            qe_upflib
            qe_fftx
            qe_elpa
            qe_xclib)

###########################################################
# epw.x
###########################################################
set(sources src/epw.f90)
qe_add_executable(qe_epw_exe ${sources})
set_target_properties(qe_epw_exe PROPERTIES OUTPUT_NAME epw.x)
target_link_libraries(qe_epw_exe PRIVATE qe_modules qe_pw qe_fftx qe_epw qe_upflib)

###########################################################
# ZG.x
###########################################################
set(sources ZG/src/ZG.f90)
qe_add_executable(qe_zg_exe ${sources})
set_target_properties(qe_zg_exe PROPERTIES OUTPUT_NAME ZG.x)
target_link_libraries(qe_zg_exe PRIVATE qe_modules qe_pw qe_phonon_ph qe_lr_modules qe_lapack)

###########################################################
# disca.x
###########################################################
set(sources ZG/src/disca.f90)
qe_add_executable(qe_zg_disca_exe ${sources})
set_target_properties(qe_zg_disca_exe PROPERTIES OUTPUT_NAME disca.x)
target_link_libraries(qe_zg_disca_exe PRIVATE qe_modules qe_pw qe_phonon_ph qe_lr_modules qe_lapack)

###########################################################
# pp_disca.x
###########################################################
set(sources ZG/src/pp_disca.f90)
qe_add_executable(qe_zg_pp_disca_exe ${sources})
set_target_properties(qe_zg_pp_disca_exe PROPERTIES OUTPUT_NAME pp_disca.x)
target_link_libraries(qe_zg_pp_disca_exe PRIVATE qe_modules qe_lr_modules qe_lapack)

###########################################################
# bands_unfold.x
###########################################################
set(src_bands_unfold_x ZG/src/bands_unfold.f90)
qe_add_executable(qe_zg_bands_unfold_exe ${sources})
set_target_properties(qe_zg_bands_unfold_exe PROPERTIES OUTPUT_NAME bands_unfold.x)
target_link_libraries(qe_zg_bands_unfold_exe PRIVATE qe_pw qe_modules qe_pp qe_upflib)

###########################################################
# pp_spctrlfn.x
###########################################################
set(sources ZG/src/pp_spctrlfn.f90)
qe_add_executable(qe_zg_pp_spctrlfn_exe ${sources})
set_target_properties(qe_zg_pp_spctrlfn_exe PROPERTIES OUTPUT_NAME pp_spctrlfn.x)
target_link_libraries(qe_zg_pp_spctrlfn_exe PRIVATE qe_modules qe_lr_modules qe_lapack)

###########################################################
# epsilon_Gaus.x
###########################################################
set(src_epsilon_Gaus_x ZG/src/epsilon_Gaus.f90)
qe_add_executable(qe_zg_epsilon_Gaus_exe ${sources})
set_target_properties(qe_zg_epsilon_Gaus_exe PROPERTIES OUTPUT_NAME epsilon_Gaus.x)
target_link_libraries(qe_zg_epsilon_Gaus_exe PRIVATE qe_pw qe_modules qe_pp qe_upflib)

###########################################################

qe_install_targets(
    # Libraries
    qe_epw
    # Exacutables
    qe_epw_exe qe_zg_exe qe_zg_disca_exe qe_zg_pp_disca_exe qe_zg_bands_unfold_exe qe_zg_pp_spctrlfn_exe qe_zg_epsilon_Gaus_exe)

install(PROGRAMS bin/pp.py TYPE BIN RENAME epw_pp.py)

add_custom_target(epw
    DEPENDS
        qe_epw_exe
    COMMENT
        "electron-Phonon Coupling with wannier functions")
