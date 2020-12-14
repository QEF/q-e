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
    src/ephwann_shuffle_mem.f90
    src/epw_init.f90
    src/epw_readin.f90
    src/epw_setup.f90
    src/epw_summary.f90
    src/loadumat.f90
    src/stop_epw.f90
    src/wannierEPW.f90
    src/loadbm.f90
    src/ephblochkq.f90
    src/wfc_elec.f90
    src/test_tools.f90)

qe_add_library(qe_epw ${sources})
add_library(QE::EPW ALIAS qe_epw)
target_link_libraries(qe_epw
    PRIVATE
        QE::WANNIER90
        QE::Modules
        QE::LR_Modules
        QE::PHonon::PH
        QE::PW
        QE::MPI_Fortran
        QE::UTILX
        QE::UPF
        QE::FFTX)

###########################################################
# epw.x
###########################################################
set(sources src/epw.f90)
qe_add_executable(qe_epw_exe ${sources})
set_target_properties(qe_epw_exe PROPERTIES OUTPUT_NAME epw.x)
target_link_libraries(qe_epw_exe
    PRIVATE
        QE::Modules
        QE::PP
	    QE::PW
	    QE::FFTX
        QE::UTILX
        QE::EPW)

###########################################################

qe_install_targets(
    # Libraries
    qe_epw
    # Exacutables
    qe_epw_exe)