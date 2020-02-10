set(f_sources
    scatter_mod.f90
    fft_ggen.f90
    fft_fwinv.f90
    fft_scalar.f90
    fftw_interfaces.f90
    fft_parallel.f90
    fft_interfaces.f90
    fft_interpolate.f90
    stick_base.f90
    fft_smallbox.f90
    fft_smallbox_type.f90
    fft_support.f90
    fft_error.f90
    fft_types.f90
    tg_gather.f90
    fft_helper_subroutines.f90
    fft_param.f90
    # FFT backends:
    # Builtin:
    fft_scalar.FFTW.f90
    # TODO: platform-specific stuff
    # fft_scalar.ARM_LIB.f90
    # fft_scalar.DFTI.f90
    # fft_scalar.ESSL.f90
    # fft_scalar.FFTW3.f90
    # fft_scalar.SX6.f90
)

set(c_sources
    fft_stick.c
    fftw.c
    fftw_sp.c
    fftw_dp.c)

# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${f_sources} PROPERTIES COMPILE_FLAGS -cpp)

qe_add_library(qe_fftx ${f_sources} ${c_sources})
add_library(QE::FFTX ALIAS qe_fftx)
target_link_libraries(qe_fftx
    PRIVATE
        QE::LAPACK)
qe_install_targets(qe_fftx)

###########################################################
# Tests
# TODO move all tests to a proper location
###########################################################
if(QE_ENABLE_TEST AND QE_ENABLE_MPI)
    # TODO test.f90 seems to work only with MPI
    # TODO rename all sources *.f90 -> *.F90
    set_source_files_properties(test.f90 PROPERTIES COMPILE_FLAGS -cpp)
    qe_add_executable(qe_fftx_test test.f90)
    target_link_libraries(qe_fftx_test
        PRIVATE
            QE::FFTX
            QE::MPI_Fortran)
    add_test(NAME test-fftx COMMAND qe_fftx_test)
endif(QE_ENABLE_TEST AND QE_ENABLE_MPI)
