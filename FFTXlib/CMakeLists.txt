###########################################################
# FFT
# The following targets will be defined:
add_library(qe_fft INTERFACE)
add_library(QE::FFT ALIAS qe_fft)
qe_install_targets(qe_fft)
###########################################################
if(QE_FFTW_VENDOR)
    # A vendor FFTW has been specified
    if(QE_FFTW_VENDOR MATCHES "Intel")
        if(QE_ENABLE_OPENMP)
            SET(BLA_VENDOR Intel10_64lp)
        else()
            SET(BLA_VENDOR Intel10_64lp_seq)
        endif()
        find_package(LAPACK QUIET)
        find_package(VendorFFTW REQUIRED COMPONENTS MKL)
        if(VendorFFTW_FOUND)
            message(STATUS "Found ${QE_FFTW_VENDOR} library")
        else()
            message(FATAL_ERROR "Failed to find ${QE_FFTW_VENDOR} library. "
                "Alternatively, remove -DQE_FFTW_VENDOR to automatic search a FFTW library installed in the system.")
        endif()
        if(QE_FFTW_VENDOR MATCHES "Intel_DFTI")
            qe_add_global_compile_definitions(__DFTI)
            set(qe_fft_wrappers fft_scalar.DFTI.f90)
        elseif(QE_FFTW_VENDOR MATCHES "Intel_FFTW3")
            qe_add_global_compile_definitions(__FFTW3)
            set(qe_fft_wrappers fft_scalar.FFTW3.f90)
        endif()
        target_link_libraries(qe_fft INTERFACE VendorFFTW)
    elseif(QE_FFTW_VENDOR MATCHES "ArmPL")
        if(QE_ENABLE_OPENMP)
            SET(BLA_VENDOR Arm_mp)
        else()
            SET(BLA_VENDOR Arm)
        endif()
        find_package(LAPACK QUIET)
        find_package(VendorFFTW REQUIRED COMPONENTS ArmPL)
        if(VendorFFTW_FOUND)
            message(STATUS "Found ${QE_FFTW_VENDOR} library")
        else()
            message(FATAL_ERROR "Failed to find ${QE_FFTW_VENDOR} library. "
                "Alternatively, remove -DQE_FFTW_VENDOR to automatic search a FFTW library installed in the system.")
        endif()
        qe_add_global_compile_definitions(__FFTW)
        set(qe_fft_wrappers fft_scalar.FFTW.f90)
        target_link_libraries(qe_fft INTERFACE VendorFFTW)
    elseif(QE_FFTW_VENDOR MATCHES "IBMESSL")
        SET(BLA_VENDOR IBMESSL)
        find_package(BLAS QUIET)
        find_package(VendorFFTW REQUIRED COMPONENTS IBMESSL)
        if(VendorFFTW_FOUND)
            message(STATUS "Found ${QE_FFTW_VENDOR} library")
        else()
            message(FATAL_ERROR "Failed to find ${QE_FFTW_VENDOR} library. "
                "Alternatively, remove -DQE_FFTW_VENDOR to automatic search a FFTW library installed in the system.")
        endif()
        qe_add_global_compile_definitions(__LINUX_ESSL)
        set(qe_fft_wrappers fft_scalar.ESSL.f90)
        target_link_libraries(qe_fft INTERFACE VendorFFTW)
    elseif(QE_FFTW_VENDOR MATCHES "FFTW3")
        if(QE_ENABLE_OPENMP)
            find_package(FFTW3 REQUIRED COMPONENTS DOUBLE_OPENMP)
        else()
            find_package(FFTW3 REQUIRED COMPONENTS DOUBLE)
        endif()
        if(NOT FFTW3_FOUND)
            message(FATAL_ERROR "Failed to find ${QE_FFTW_VENDOR} library. "
                "Alternatively, remove -DQE_FFTW_VENDOR to automatic search a FFTW library installed in the system.")
        endif()
        qe_add_global_compile_definitions(__FFTW3)
        set(qe_fft_wrappers fft_scalar.FFTW3.f90)
        target_link_libraries(qe_fft INTERFACE FFTW)
    elseif(QE_FFTW_VENDOR MATCHES "Internal")
        message(STATUS "QE internal implementation of FFTW (FFTXLib)")
        qe_add_global_compile_definitions(__FFTW)
        set(qe_fft_wrappers fft_scalar.FFTW.f90)
    else()
        message(FATAL_ERROR "The FFTW vendor library '${QE_FFTW_VENDOR}' is not supported!")
    endif()
else()
    # Try to find a vendor FFTW installed in the system
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
            qe_add_global_compile_definitions(__FFTW)
            set(qe_fft_wrappers fft_scalar.FFTW.f90)
        elseif(VendorFFTW_ID STREQUAL "IBMESSL")
            qe_add_global_compile_definitions(__LINUX_ESSL)
            set(qe_fft_wrappers fft_scalar.ESSL.f90)
        endif()
        
        target_link_libraries(qe_fft INTERFACE VendorFFTW)
    else()
        # Try to find the official FFTW3
        message(STATUS "No vendor FFTW found")
        if(QE_ENABLE_OPENMP)
            find_package(FFTW3 COMPONENTS DOUBLE_OPENMP)
        else()
            find_package(FFTW3 COMPONENTS DOUBLE)
        endif()
        if(FFTW3_FOUND)
            qe_add_global_compile_definitions(__FFTW3)
            set(qe_fft_wrappers fft_scalar.FFTW3.f90)
            target_link_libraries(qe_fft INTERFACE FFTW3)
        else()
            # Cannot find anything useful, just fall back to the internal FFTW
            message(STATUS "No FFTW3 library found. FFTW is falling back to QE internal implementation (FFTXLib)")
            qe_add_global_compile_definitions(__FFTW)
            set(qe_fft_wrappers fft_scalar.FFTW.f90)
        endif()
    endif()
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
