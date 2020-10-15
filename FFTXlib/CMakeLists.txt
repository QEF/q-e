###########################################################
# FFT
# The following targets will be defined:
add_library(qe_fft INTERFACE)
add_library(QE::FFT ALIAS qe_fft)
qe_install_targets(qe_fft)
###########################################################
set(VENDOR)
if(BLA_VENDOR MATCHES "Intel")
    set(VENDOR "MKL")
endif()
if(BLA_VENDOR MATCHES "Arm")
    set(VENDOR "ARMPL")
endif()
if(BLA_VENDOR MATCHES "Essl")
    set(VENDOR "ESSL")
endif()

if(QE_ENABLE_OPENMP)
    find_package(FFTW COMPONENTS OMP VENDOR)
    if(FFTW_FOUND)
        qe_add_global_compile_definitions(__FFTW3)
        set(qe_fft_wrappers fft_scalar.FFTW3.f90)
        target_link_libraries(qe_fft INTERFACE ${FFTW_LIBRARIES})
        message(STATUS "Found OpenMP FFTW3 library")
        message(STATUS "Found FFTW3 libraries: ${FFTW_LIBRARIES}")
        message(STATUS "Found FFTW3 include: ${FFTW_INCLUDE_DIRS}")
    endif()
else()
    find_package(FFTW COMPONENTS VENDOR)
    if(FFTW_FOUND)
        qe_add_global_compile_definitions(__FFTW3)
        set(qe_fft_wrappers fft_scalar.FFTW3.f90)
        target_link_libraries(qe_fft INTERFACE ${FFTW_LIBRARIES})
        message(STATUS "Found FFTW3 library")
        message(STATUS "Found FFTW3 libraries: ${FFTW_LIBRARIES}")
        message(STATUS "Found FFTW3 include: ${FFTW_INCLUDE_DIRS}")
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
    fft_scatter_2d.f90
    scatter_mod.f90 
    fft_ggen.f90
    fft_fwinv.f90
    fft_scalar.f90
    fftw_interfaces.f90
    fft_parallel.f90
    fft_parallel_2d.f90
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

if(QE_ENABLE_CUDA)
  set(f_sources ${f_sources}
      fft_scatter_gpu.f90
      fft_scatter_2d_gpu.f90
      fft_buffers.f90
      fft_scalar.cuFFT.f90
  )
endif()

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
if(QE_ENABLE_TEST)
    set(sources fft_test.f90)
    qe_add_executable(qe_fftx_test ${sources})
    set_target_properties(qe_fftx_test PROPERTIES OUTPUT_NAME qe_fftx_test.x RUNTIME_OUTPUT_DIRECTORY ${QE_TESTS_DIR})
    target_link_libraries(qe_fftx_test
        PRIVATE
            QE::OpenMP_Fortran
            QE::MPI_Fortran
            QE::LAPACK
            QE::FFTX)
    add_test(NAME test_qe_fftx COMMAND ${QE_TESTS_DIR}/qe_fftx_test.x)
endif()
