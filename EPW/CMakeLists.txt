set(sources
    src/input.f90
    src/ep_constants.f90
    src/global_var.f90
    src/io/io_var.f90
    src/bloch2wannier.f90
    src/utilities/low_lvl.f90
    src/utilities/parallelism.f90
    src/longrange.f90
    src/io/io.f90
    src/io/io_transport.f90
    src/io/io_selfen.f90
    src/wigner.f90
    src/wannier2bloch.f90
    src/wannier2bloch_opt.f90
    src/cumulant.f90
    src/close.f90
    src/printing.f90
    src/utilities/kfold.f90
    src/utilities/bzgrid.f90
    src/io/io_supercond.f90
    src/utilities/utilities.f90
    src/supercond_common.f90
    src/supercond_vertex.f90
    src/supercond.f90
    src/supercond_iso.f90
    src/supercond_aniso.f90
    src/supercond_coul.f90
    src/utilities/sparse_ir.f90
    src/io/io_sparse_ir.f90
    src/stop.f90
    src/transport_mag.f90
    src/utilities/symmetry.f90
    src/transport_legacy.f90
    src/transport.f90
    src/pw2wan.f90
    src/wannierization.f90
    src/indabs.f90
    src/io/io_indabs.f90
    src/utilities/bcast_input.f90
    src/io/io_ahc.f90
    src/dvqpsi.f90
    src/ep_coarse.f90
    src/ep_coarse_unfolding.f90
    src/supercond_driver.f90
    src/selfen.f90
    src/spectral.f90
    src/wannier.f90
    src/use_wannier.f90
    src/init.f90
    src/readin.f90
    src/setups.f90
    src/summaries.f90
    src/wann_common.f90
    src/polaron.f90
    src/qdabs.f90
    src/wfpt.f90
    src/exphonon.f90
    src/expolaron.f90
    src/utilities/screening.f90
    src/tdbe_driver.f90
    src/tdbe_mod.f90
    src/tdbe_common.f90
    src/phph.f90
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
# nscf2supercond.x
###########################################################
set(sources src/utilities/nscf2supercond.f90)
qe_add_executable(qe_epw_nscf2supercond_exe ${sources})
set_target_properties(qe_epw_nscf2supercond_exe PROPERTIES OUTPUT_NAME nscf2supercond.x)
target_link_libraries(qe_epw_nscf2supercond_exe PRIVATE qe_pw qe_modules qe_pp qe_upflib)

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
    qe_epw_exe
    qe_epw_nscf2supercond_exe
    qe_zg_exe
    qe_zg_disca_exe
    qe_zg_pp_disca_exe
    qe_zg_bands_unfold_exe
    qe_zg_pp_spctrlfn_exe
    qe_zg_epsilon_Gaus_exe)

install(PROGRAMS bin/pp.py TYPE BIN RENAME epw_pp.py)

add_custom_target(epw
    DEPENDS
        qe_epw_exe
        qe_epw_nscf2supercond_exe
        qe_zg_exe
        qe_zg_disca_exe
        qe_zg_pp_disca_exe
        qe_zg_bands_unfold_exe
        qe_zg_pp_spctrlfn_exe
        qe_zg_epsilon_Gaus_exe
    COMMENT
        "electron-Phonon Coupling with wannier functions")
