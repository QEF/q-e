###########################################################
# FFT
# The following targets will be defined:
add_library(qe_fft INTERFACE)
qe_install_targets(qe_fft)
###########################################################
SET(BLA_VENDOR_SAVED ${BLA_VENDOR})
if(QE_FFTW_VENDOR STREQUAL "AUTO")
    # Try to find a FFTW library installed in the system
    find_package(VendorFFTW)
    if(VendorFFTW_FOUND)
        if(VendorFFTW_ID STREQUAL "Intel")
            qe_add_global_compile_definitions(__DFTI)
            set(qe_fft_wrappers fft_scalar.DFTI.f90)
        elseif(VendorFFTW_ID STREQUAL "Arm")
            ## ARMPL implements the same interface of the standard FFTW
            ## no need of QE ARMPL wrapper
            #qe_add_global_compile_definitions(__ARM_LIB)
            #set(qe_fft_wrappers fft_scalar.ARM_LIB.f90)
            qe_add_global_compile_definitions(__FFTW3)
            set(qe_fft_wrappers fft_scalar.FFTW3.f90)
        elseif(VendorFFTW_ID STREQUAL "IBMESSL")
            qe_add_global_compile_definitions(__LINUX_ESSL)
            set(qe_fft_wrappers fft_scalar.ESSL.f90)
        endif()
        target_link_libraries(qe_fft INTERFACE VendorFFTW)
    else()
        # Try to find the official FFTW3
        message(STATUS "No vendor FFTW found")
        if(QE_ENABLE_OPENMP)
            find_package(FFTW3 COMPONENTS DOUBLE_OPENMP DOUBLE)
        else()
            find_package(FFTW3 COMPONENTS DOUBLE)
        endif()
        if(FFTW3_FOUND)
            qe_add_global_compile_definitions(__FFTW3)
            set(qe_fft_wrappers fft_scalar.FFTW3.f90)
            target_link_libraries(qe_fft INTERFACE FFTW3)
        else()
            message(STATUS "CMake variable FFTW3_ROOT may be used to hint the root directory of FFTW3 installation.")
            # Cannot find anything useful, just point out internal FFTW fallback.
            message(FATAL_ERROR "Failed to find an external FFTW library. "
                                "Alternatively, '-DQE_FFTW_VENDOR=Internal' may be used to enable reference FFTW "
                                "at a performance loss compared to optimized libraries.")

        endif()
    endif()
elseif(QE_FFTW_VENDOR MATCHES "Intel")
    if(QE_FFTW_VENDOR STREQUAL "Intel_DFTI")
        qe_add_global_compile_definitions(__DFTI)
        set(qe_fft_wrappers fft_scalar.DFTI.f90)
    elseif(QE_FFTW_VENDOR STREQUAL "Intel_FFTW3")
        qe_add_global_compile_definitions(__FFTW3)
        set(qe_fft_wrappers fft_scalar.FFTW3.f90)
    else()
        message(FATAL_ERROR "The unknown Intel FFTW library '${QE_FFTW_VENDOR}' is not supported!")
    endif()
    if(QE_ENABLE_OPENMP)
        SET(BLA_VENDOR Intel10_64lp)
    else()
        SET(BLA_VENDOR Intel10_64lp_seq)
    endif()
    find_package(LAPACK QUIET)
    find_package(VendorFFTW REQUIRED COMPONENTS MKL)
    if(VendorFFTW_FOUND)
        target_link_libraries(qe_fft INTERFACE VendorFFTW)
        message(STATUS "Found ${QE_FFTW_VENDOR} library")
    else()
        message(FATAL_ERROR "Failed to find ${QE_FFTW_VENDOR} library. "
            "Alternatively, remove -DQE_FFTW_VENDOR to automatic search a FFTW library in the system.")
    endif()
elseif(QE_FFTW_VENDOR STREQUAL "ArmPL")
    if(QE_ENABLE_OPENMP)
        SET(BLA_VENDOR Arm_mp)
    else()
        SET(BLA_VENDOR Arm)
    endif()
    find_package(LAPACK QUIET)
    find_package(VendorFFTW REQUIRED COMPONENTS ArmPL)
    if(VendorFFTW_FOUND)
        qe_add_global_compile_definitions(__FFTW3)
        set(qe_fft_wrappers fft_scalar.FFTW3.f90)
        target_link_libraries(qe_fft INTERFACE VendorFFTW)
        message(STATUS "Found ${QE_FFTW_VENDOR} library")
    else()
        message(FATAL_ERROR "Failed to find ${QE_FFTW_VENDOR} library. "
            "Alternatively, remove -DQE_FFTW_VENDOR to automatic search a FFTW library installed in the system.")
    endif()
elseif(QE_FFTW_VENDOR STREQUAL "IBMESSL")
    SET(BLA_VENDOR IBMESSL)
    find_package(BLAS QUIET)
    find_package(VendorFFTW REQUIRED COMPONENTS IBMESSL)
    if(VendorFFTW_FOUND)
        qe_add_global_compile_definitions(__LINUX_ESSL)
        set(qe_fft_wrappers fft_scalar.ESSL.f90)
        target_link_libraries(qe_fft INTERFACE VendorFFTW)
        message(STATUS "Found ${QE_FFTW_VENDOR} library")
    else()
        message(FATAL_ERROR "Failed to find ${QE_FFTW_VENDOR} library. "
            "Alternatively, remove -DQE_FFTW_VENDOR to automatic search a FFTW library installed in the system.")
    endif()
elseif(QE_FFTW_VENDOR STREQUAL "FFTW3")
    if(QE_ENABLE_OPENMP)
        find_package(FFTW3 REQUIRED COMPONENTS DOUBLE_OPENMP DOUBLE)
    else()
        find_package(FFTW3 REQUIRED COMPONENTS DOUBLE)
    endif()
    if(FFTW3_FOUND)
        qe_add_global_compile_definitions(__FFTW3)
        set(qe_fft_wrappers fft_scalar.FFTW3.f90)
        target_link_libraries(qe_fft INTERFACE FFTW)
        message(STATUS "Found FFTW3 library")
    else()
        message(FATAL_ERROR "Failed to find ${QE_FFTW_VENDOR} library. "
            "Alternatively, remove -DQE_FFTW_VENDOR to automatic search a FFTW library installed in the system.")
    endif()
elseif(QE_FFTW_VENDOR STREQUAL "Internal")
    message(STATUS "QE internal implementation of FFTW (FFTXLib)")
    qe_add_global_compile_definitions(__FFTW)
    set(qe_fft_wrappers fft_scalar.FFTW.f90)
else()
    message(FATAL_ERROR "The FFTW vendor library '${QE_FFTW_VENDOR}' is not supported!")
endif()
SET(BLA_VENDOR ${BLA_VENDOR_SAVED})

set(f_src_fftx
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
    fft_param.f90)
if(QE_ENABLE_CUDA)
    set(f_src_fftx
        ${f_src_fftx}
        fft_scalar.cuFFT.f90
        fft_buffers.f90
        fft_scatter_gpu.f90
        fft_scatter_2d_gpu.f90)
endif()
qe_enable_cuda_fortran("${f_src_fftx}")

set(c_src_fftx
    fft_stick.c
    fftw.c
    fftw_sp.c
    fftw_dp.c)

qe_add_library(qe_fftx 
    ${f_src_fftx} 
    ${c_src_fftx} 
    ${qe_fft_wrappers})

target_link_libraries(qe_fftx
    PRIVATE
        qe_fft
        qe_openmp_fortran
        qe_mpi_fortran)
if(QE_ENABLE_CUDA)
    target_link_libraries(qe_fftx
        PRIVATE
            CUDA::cufft)
endif()
        
qe_install_targets(qe_fftx)

###########################################################
# Tests
# TODO move all tests to a proper location
###########################################################
if(QE_ENABLE_TEST)
    set(src_fftx_test fft_test.f90)
    qe_add_executable(qe_fftx_test ${src_fftx_test})
    set_target_properties(qe_fftx_test 
        PROPERTIES 
            OUTPUT_NAME qe_fftx_test.x 
            RUNTIME_OUTPUT_DIRECTORY ${QE_TESTS_DIR})
    target_link_libraries(qe_fftx_test
        PRIVATE
            qe_openmp_fortran
            qe_mpi_fortran
            qe_fftx)
    if(QE_ENABLE_CUDA)
        set_property(TARGET qe_fftx_test APPEND PROPERTY COMPILE_DEFINITIONS "__FFTX_USE_EXTERNAL_CLOCK")
        target_link_libraries(qe_fftx_test PRIVATE qe_utilx)
    endif()
    add_test(NAME test_qe_fftx 
        COMMAND ${QE_TESTS_DIR}/qe_fftx_test.x)
endif()
