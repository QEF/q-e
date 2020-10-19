###########################################################
# FFT
# The following targets will be defined:
add_library(qe_fft INTERFACE)
add_library(QE::FFT ALIAS qe_fft)
qe_install_targets(qe_fft)
###########################################################
# Try to find a vendor FFTW
find_package(VendorFFTW)
if(VendorFFTW_FOUND)
    if(QE_FFTW_QUAD)
        message(FATAL_ERROR "FFTW quad precision is not supported in vendor FFTW")
    endif()

    if(VendorFFTW_ID STREQUAL "Intel")
        qe_add_global_compile_definitions(__DFTI)
        set(qe_fft_wrappers fft_scalar.DFTI.f90)
    elseif(VendorFFTW_ID STREQUAL "Arm")
        qe_add_global_compile_definitions(__ARM_LIB)
        set(qe_fft_wrappers fft_scalar.ARM_LIB.f90)
    elseif(VendorFFTW_ID STREQUAL "Essl")
        qe_add_global_compile_definitions(__LINUX_ESSL)
        set(qe_fft_wrappers fft_scalar.ESSL.f90)
    endif()
    
    target_link_libraries(qe_fft INTERFACE VendorFFTW)
else()
    # Try to find a vendor FFTW3
    message("-- No vendor FFTW found")

    set(FFTW_COMPS)
    if(QE_FFTW_SINGLE)
        if(QE_FFTW_CONF)
            list(APPEND FFTW_COMPS "FLOAT_${QE_FFTW_CONF}")
        else()
            list(APPEND FFTW_COMPS "FLOAT")
        endif()
    endif()
    if(QE_FFTW_DOUBLE)
        if(QE_FFTW_CONF)
            list(APPEND FFTW_COMPS "DOUBLE_${QE_FFTW_CONF}")
        else()
            list(APPEND FFTW_COMPS "DOUBLE")
        endif()
    endif()
    if(QE_FFTW_QUAD)
        if(QE_FFTW_CONF)
            list(APPEND FFTW_COMPS "LONGDOUBLE_${QE_FFTW_CONF}")
        else()
            list(APPEND FFTW_COMPS "LONGDOUBLE")
        endif()
    endif()

    find_package(FFTW COMPONENTS ${FFTW_COMPS})
    if(FFTW_FOUND)
        qe_add_global_compile_definitions(__FFTW3)
        set(qe_fft_wrappers fft_scalar.FFTW3.f90)
        target_link_libraries(qe_fft INTERFACE FFTW)
    else()
        message("-- No FFTW3 library found")
    endif()
endif()

# Cannot find anything useful, just fall back to the internal FFTW
if(NOT VendorFFTW_FOUND AND NOT FFTW_FOUND)
    message(STATUS "FFTW is falling back to QE internal implementation (FFTXLib)")
    qe_add_global_compile_definitions(__FFTW)
    set(qe_fft_wrappers fft_scalar.FFTW.f90)
endif()

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
