set(src_couple
    src/libpwscf.f90
    src/libcpv.f90
    src/libqemod.f90)
qe_enable_cuda_fortran("${src_couple}")

qe_add_library(qe_couple ${src_couple})
target_link_libraries(qe_couple
    PRIVATE
        qe_modules
        qe_lax
	    qe_fftx
        qe_cpv)

###########################################################

qe_install_targets(qe_couple)

add_custom_target(couple
    DEPENDS
        qe_couple
    COMMENT
        "library interface for coupling to external codes")
