###########################################################
# FFT
# The following targets will be defined:
add_library(qe_fft INTERFACE)
add_library(QE::FFT ALIAS qe_fft)
qe_install_targets(qe_fft)
###########################################################
if(VENDOR_FOUND AND NOT FFTW_ROOT)
    if(VENDOR MATCHES "Intel")
        qe_add_global_compile_definitions(__DFTI)
        set(qe_fft_wrappers fft_scalar.DFTI.f90)
    elseif(VENDOR MATCHES "Armpl")
        qe_add_global_compile_definitions(__ARM_LIB)
        set(qe_fft_wrappers fft_scalar.ARM_LIB.f90)
    elseif(VENDOR MATCHES "IBMESSL")
        qe_add_global_compile_definitions(__LINUX_ESSL)
        set(qe_fft_wrappers fft_scalar.ESSL.f90)
    # elseif(VENDOR MATCHES "ACML")
    #    qe_add_global_compile_definitions(__ACML)
    else()
        message(FATAL_ERROR "The vendor specific FFT found isn't yet supported in QE, please open an issue: https://gitlab.com/QEF/q-e/issues")
    endif()

    target_link_libraries(qe_fft 
        INTERFACE 
            ${FFTW_LIBRARIES}
            ${FFTW_LINKER_FLAGS})
    target_include_directories(qe_fft
        INTERFACE
            ${FFTW_INCLUDE_DIRS})
else()
    # No vendor-specific fft flavor has been found,
    # just look for a generic FFTW3 implementation
    if(QE_ENABLE_OPENMP)
        find_package(FFTW QUIET COMPONENTS DOUBLE_OPENMP_LIB)
        if(FFTW_FOUND)
            qe_add_global_compile_definitions(__FFTW3)
            set(qe_fft_wrappers fft_scalar.FFTW3.f90)
            target_link_libraries(qe_fft INTERFACE FFTW::DoubleOpenMP)
            message(STATUS "Found Double+OpenMP FFTW3 library")
            message(STATUS "Found FFTW3 libraries: ${FFTW_DOUBLE_OPENMP_LIB}")
            message(STATUS "Found FFTW3 include: ${FFTW_INCLUDE_DIRS}")
        endif()
    else()
        find_package(FFTW QUIET COMPONENTS DOUBLE_LIB)
        if(FFTW_FOUND)
            qe_add_global_compile_definitions(__FFTW3)
            set(qe_fft_wrappers fft_scalar.FFTW3.f90)
            target_link_libraries(qe_fft INTERFACE FFTW::Double)
            message(STATUS "Found Double FFTW3 library")
            message(STATUS "Found FFTW3 libraries: ${FFTW_DOUBLE_LIB}")
            message(STATUS "Found FFTW3 include: ${FFTW_INCLUDE_DIRS}")
        endif()
    endif()
endif()

# Cannot find anything useful, just fall back to the internal FFTW
if(NOT FFTW_FOUND)
    message(STATUS "FFTW is falling back to QE internal implementation (FFTXLib)")
    qe_add_global_compile_definitions(__FFTW)
    set(qe_fft_wrappers fft_scalar.FFTW.f90)
endif(NOT FFTW_FOUND)

set(f_sources
    fft_scatter.f90
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
)

set(c_sources
    fft_stick.c
    fftw.c
    fftw_sp.c
    fftw_dp.c)

qe_add_library(qe_fftx ${f_sources} ${c_sources} ${qe_fft_wrappers})
add_library(QE::FFTX ALIAS qe_fftx)
target_link_libraries(qe_fftx
    PRIVATE
        QE::FFT
        QE::OpenMP_Fortran
        QE::MPI_Fortran
        QE::LAPACK)
qe_install_targets(qe_fftx)

###########################################################
# Tests
# TODO move all tests to a proper location
###########################################################
if(QE_ENABLE_TEST AND QE_ENABLE_MPI) # TODO fix test without MPI
    # TODO test.f90 seems to work only with MPI
    set(sources fft_test.f90)
    qe_add_executable(qe_fftx_test ${sources})
    target_link_libraries(qe_fftx_test
        PRIVATE
            QE::OpenMP_Fortran
            QE::MPI_Fortran
            QE::LAPACK
            QE::FFTX)
    add_test(NAME test-fftx COMMAND qe_fftx_test)
endif()
