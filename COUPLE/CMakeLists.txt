set(sources
    src/libpwscf.f90
    src/libcpv.f90
    src/libqemod.f90)

qe_add_library(qe_couple ${sources})
add_library(QE::COUPLE ALIAS qe_couple)
target_link_libraries(qe_couple
    PUBLIC
        QE::Modules
        QE::LAX
	    QE::FFTX
        QE::UTILX
        QE::CPV)

###########################################################

qe_install_targets(qe_couple)
