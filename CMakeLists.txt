cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
# CMake < 3.13 cannot install external targets:
# https://gitlab.kitware.com/cmake/cmake/merge_requests/2152

# CMake < 3.14 generates incorrect dependency graphs with
# alias targets:
# https://gitlab.kitware.com/cmake/cmake/merge_requests/2521

# CMake < 3.20 does not support NVHPC compiler id
# https://cmake.org/cmake/help/latest/release/3.20.html#compilers

# Policy CMP0048: The project() command manages VERSION variables
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)

project(qe
    VERSION 7.3.1
    DESCRIPTION "ESPRESSO: opEn-Source Package for Research in Electronic Structure, Simulation, and Optimization"
    LANGUAGES Fortran C)

if(${qe_BINARY_DIR} STREQUAL ${qe_SOURCE_DIR})
    message(FATAL_ERROR "QE source folder cannot be safely used as a build folder!")
endif()

##########################################################
# Define the paths for static libraries and executables
##########################################################
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${qe_BINARY_DIR}/lib 
    CACHE
    PATH "Single output directory for building all libraries.")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${qe_BINARY_DIR}/bin 
    CACHE
    PATH "Single output directory for building all executables.")

###########################################################
# Build helpers
###########################################################
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
include(cmake/qeHelpers.cmake)

###########################################################
# Build Type
# Ensure that a specific, default build type is set when
# none has been explicitly set by the user
###########################################################
qe_ensure_build_type("Release")

###########################################################
# Modules
###########################################################
include(CheckFunctionExists)
# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
include(GNUInstallDirs)

###########################################################
# Build Options
###########################################################
include(CMakeDependentOption)
option(QE_ENABLE_CUDA
    "enable CUDA acceleration on NVIDIA GPUs" OFF)
if(QE_ENABLE_CUDA)
    option(QE_ENABLE_OPENACC "enable OpenACC acceleration" ON)
    # OpenMP enabled by default if CUDA is enable
    option(QE_ENABLE_OPENMP
        "enable distributed execution support via OpenMP" ON)
else()
    option(QE_ENABLE_OPENACC "enable OpenACC acceleration" OFF)
    option(QE_ENABLE_OPENMP
        "enable distributed execution support via OpenMP" OFF)
endif()
cmake_dependent_option(QE_ENABLE_OFFLOAD
    "enable OpenMP offload"
    OFF "QE_ENABLE_OPENMP" OFF)
option(QE_ENABLE_MPI
    "enable distributed execution support via MPI" ON)
option(QE_ENABLE_MPI_GPU_AWARE
    "enable GPU aware MPI operations" OFF)
option(QE_ENABLE_TEST
    "enable unit and system tests" ON)
cmake_dependent_option(QE_ENABLE_BENCHMARK
    "enable benchmark tests" OFF "QE_ENABLE_TEST" OFF)
option(QE_ENABLE_TRACE
    "enable execution tracing output" OFF)
option(QE_ENABLE_PROFILE_NVTX
        "enable execution of NVIDIA NVTX profiler plugin" OFF)
option(QE_ENABLE_MPI_INPLACE
    "enable inplace MPI calls (ignored when QE_ENABLE_MPI=OFF)" OFF)
option(QE_ENABLE_MPI_MODULE
    "use MPI via Fortran module instead of mpif.h header inclusion" OFF)
option(QE_ENABLE_BARRIER
    "enable global synchronization between execution units" OFF)
option(QE_LAPACK_INTERNAL
    "enable internal reference LAPACK" OFF)
option(QE_ENABLE_SCALAPACK
    "enable SCALAPACK execution units" OFF)
cmake_dependent_option(QE_ENABLE_SCALAPACK_QRCP
    "enable SCALAPACK QRCP in pw2wannier90 (requires SCALAPACK>=2.1.0 or Intel MKL>=2020)"
    OFF "QE_ENABLE_SCALAPACK" OFF)
option(QE_ENABLE_ELPA
    "enable ELPA execution units" OFF)
option(QE_ENABLE_LIBXC
    "enable LIBXC execution units" OFF)
option(QE_ENABLE_HDF5
    "enable HDF5 data collection" OFF)
option(QE_ENABLE_STATIC_BUILD
    "enable fully static build of executables" OFF)
option(QE_ENABLE_DOC
    "enable documentation building" OFF)
option(QE_CLOCK_SECONDS
    "print program time in seconds" OFF)
set(QE_FFTW_VENDOR "AUTO" CACHE
    STRING "select a specific FFTW library [Intel_DFTI, Intel_FFTW3, ArmPL, IBMESSL, FFTW3, Internal]")
set(QE_ENABLE_SANITIZER "none" CACHE STRING "none,asan,ubsan,tsan,msan")
set(QE_ENABLE_PLUGINS "" CACHE STRING "Semicolon-separated list of plugins")
set (QE_EXTRA_DEFINITIONS "" CACHE STRING "Semicolon-separated list of extra global definitions")
option(QE_ENABLE_FOX
    "enable XML I/O via Fox library" OFF)
if(QE_ENABLE_FOX)
    if(FOX_ROOT)
        set(QE_FOX_INTERNAL OFF)
    endif()
    option(QE_FOX_INTERNAL
        "enable FoX internal library" ON)
endif()

if(WANNIER90_ROOT)
    set(QE_WANNIER90_INTERNAL OFF)
endif()
option(QE_WANNIER90_INTERNAL
    "enable Wannier90 internal library" ON)
if(MBD_ROOT)
    set(QE_MBD_INTERNAL OFF)
endif()
option(QE_MBD_INTERNAL
    "enable LibMBD internal library" ON)
if(DEVICEXLIB_ROOT)
    set(QE_DEVICEXLIB_INTERNAL OFF)
endif()
option(QE_DEVICEXLIB_INTERNAL
    "enable DeviceXlib internal library" ON)
if(ENVIRON_ROOT)
    set(ENVIRON_DEFAULT "EXTERNAL")
else()
    set(ENVIRON_DEFAULT "NO")
endif()
set(QE_ENABLE_ENVIRON "${ENVIRON_DEFAULT}" CACHE
    STRING "select a specific Environ library [NO, EXTERNAL, INTERNAL]")

# OSCDFT
option(QE_ENABLE_OSCDFT "enable OS-CDFT 10.1021/acs.jctc.9b00281" OFF)


# TODO change all ifdefs throughout code base to match
# cmake options
# TODO symbols beginning with '__' followed by a capital
# character are reserved for standard library use (at
# least in C, not sure about Fortran), change all feature
# macros to avoid weird behaviors

if(QE_ENABLE_CUDA)
    qe_add_global_compile_definitions(__CUDA)
endif()
if(QE_ENABLE_TRACE)
    qe_add_global_compile_definitions(__TRACE)
endif()
if(QE_ENABLE_PROFILE_NVTX)
   qe_add_global_compile_definitions(__PROFILE_NVTX) 
endif() 
if(QE_ENABLE_MPI_INPLACE)
    qe_add_global_compile_definitions(__USE_INPLACE_MPI)
endif()
if(QE_ENABLE_MPI_MODULE)
    qe_add_global_compile_definitions(__MPI_MODULE)
endif()
if(QE_ENABLE_BARRIER)
    qe_add_global_compile_definitions(__USE_BARRIER)
endif()
foreach(DEF IN LISTS QE_EXTRA_DEFINITIONS)
    qe_add_global_compile_definitions(${DEF})
endforeach()
if(QE_ENABLE_MPI)
    # OMPI_SKIP_MPICXX: skip CXX APIs on openmpi, cause trouble to C APIs
    qe_add_global_compile_definitions(__MPI OMPI_SKIP_MPICXX)
    if(QE_ENABLE_MPI_GPU_AWARE)
        qe_add_global_compile_definitions(__GPU_MPI)
    endif()
endif()
if(QE_ENABLE_SCALAPACK)
    qe_add_global_compile_definitions(__SCALAPACK)
endif()
if(QE_ENABLE_HDF5)
    qe_add_global_compile_definitions(__HDF5)
endif()
if(QE_ENABLE_ENVIRON)
    qe_add_global_compile_definitions(__ENVIRON)
endif()
if(QE_CLOCK_SECONDS)
   qe_add_global_compile_definitions(__CLOCK_SECONDS)
endif()

# Feature checks
check_function_exists(mallinfo HAVE_MALLINFO)
if(HAVE_MALLINFO)
    qe_add_global_compile_definitions(HAVE_MALLINFO)
endif()

# OSCDFT
if(QE_ENABLE_OSCDFT)
    qe_add_global_compile_definitions(__OSCDFT)
endif()

# Check options consistency
if(QE_ENABLE_STATIC_BUILD AND BUILD_SHARED_LIBS)
    message(FATAL_ERROR "Full static build of QE executables requires static QE internal libraries. QE_ENABLE_STATIC_BUILD and BUILD_SHARED_LIBS cannot be both ON")
endif()
if(QE_ENABLE_ELPA AND NOT QE_ENABLE_SCALAPACK)
    message(FATAL_ERROR "ELPA requires SCALAPACK support, enable it with '-DQE_ENABLE_SCALAPACK=ON' or disable ELPA with '-DQE_ENABLE_ELPA=OFF'")
endif()
if(QE_ENABLE_SCALAPACK AND NOT QE_ENABLE_MPI)
    message(FATAL_ERROR "SCALAPACK requires MPI support, enable it with '-DQE_ENABLE_MPI=ON' or disable SCALAPACK with '-DQE_ENABLE_SCALAPACK=OFF'")
endif()
if(QE_ENABLE_CUDA AND NOT (CMAKE_Fortran_COMPILER_ID MATCHES "PGI" OR CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC"))
    message(FATAL_ERROR "NVHPC compiler is mandatory when CUDA is enabled due QE is based on CUDA Fortran language")
endif()
if(QE_ENABLE_OPENACC AND NOT (CMAKE_Fortran_COMPILER_ID MATCHES "PGI" OR CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC"))
    message(FATAL_ERROR "NVHPC compiler is mandatory when OpenACC is enabled")
endif()
if(QE_ENABLE_MPI_GPU_AWARE AND NOT (QE_ENABLE_CUDA AND QE_ENABLE_MPI))
    message(FATAL_ERROR "GPU aware MPI requires both MPI and CUDA features enabled")
endif()
# if(QE_ENABLE_HDF5 AND NOT QE_ENABLE_MPI)
#    message(FATAL_ERROR "HDF5 requires MPI support, enable it with '-DQE_ENABLE_MPI=ON' or disable HDF5 with '-DQE_ENABLE_HDF5=OFF'")
# endif()

# Add optional sanitizers ASAN, UBSAN, MSAN
set(VALID_SANITIZERS "none" "asan" "ubsan" "tsan" "msan")
# Perform sanitizer option check, only works in debug mode
if(NOT QE_ENABLE_SANITIZER IN_LIST VALID_SANITIZERS)
  message(FATAL_ERROR "Invalid -DQE_ENABLE_SANITIZER=${QE_ENABLE_SANITIZER}, value must be one of ${VALID_SANITIZERS}")
else()
  message(STATUS "Enable sanitizer QE_ENABLE_SANITIZER=${QE_ENABLE_SANITIZER}")
endif()
# only GNU works right now
if(NOT QE_ENABLE_SANITIZER STREQUAL "none" AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  message(FATAL_ERROR "-DQE_ENABLE_SANITIZER=${QE_ENABLE_SANITIZER} only works with the GNU compiler")
endif()

# valid plugins checks
set(VALID_QE_PLUGINS "d3q" "pw2qmcpack" "gipaw" "legacy")
# Perform sanitizer option check, only works in debug mode
foreach(PLUGIN IN LISTS QE_ENABLE_PLUGINS)
  if(NOT PLUGIN IN_LIST VALID_QE_PLUGINS)
    message(FATAL_ERROR "Invalid QE plugin ${PLUGIN}, value must be one of \"${VALID_QE_PLUGINS}\".")
  else()
    message(STATUS "Enable QE plugin ${PLUGIN}")
  endif()
endforeach()

############################################################
# C preprocessor
# Note: reply on the compiler preprocessor whenever possible
############################################################
if(DEFINED ENV{CPP})
    set(QE_CPP_DEFAULT $ENV{CPP})
else()
    set(QE_CPP_DEFAULT cpp)
endif()
# QE_CPP_DEFAULT is only effective when cached QE_CPP doesn't exist.
set(QE_CPP ${QE_CPP_DEFAULT} CACHE
    STRING "C preprocessor for qe_preprocess_source in qeHelpers.cmake")
find_program(QE_CPP_FULL_PATH NAMES ${QE_CPP} DOC "C preprocessor full path")
if(QE_CPP_FULL_PATH)
    message(STATUS "C preprocessor used by qe_preprocess_source in qeHelpers.cmake: ${QE_CPP_FULL_PATH}")
else()
    set(QE_CPP_SAVED ${QE_CPP})
    unset(QE_CPP CACHE)
    message(FATAL_ERROR "C preprocessor ${QE_CPP_SAVED} not found. Pass a working one to CMake via QE_CPP!")
endif()

###########################################################
# language standard requirements
###########################################################
# TODO need to require all compilers using the same one
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI" OR CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
    set(CMAKE_C_STANDARD 11)
    set(CMAKE_C_STANDARD_REQUIRED ON)
    set(CMAKE_C_EXTENSIONS OFF)
endif()

###########################################################
# check Fortran compiler -isystem option support
###########################################################
include(CheckFortranCompilerFlag)
check_fortran_compiler_flag("-isystem ." Fortran_ISYSTEM_SUPPORTED)
if(NOT Fortran_ISYSTEM_SUPPORTED AND NOT DEFINED CMAKE_NO_SYSTEM_FROM_IMPORTED)
  set(CMAKE_NO_SYSTEM_FROM_IMPORTED ON)
endif()

###########################################################
# OpenMP
# The following targets will be defined:
add_library(qe_openmp_fortran INTERFACE)
add_library(qe_openmp_c INTERFACE)
qe_install_targets(qe_openmp_fortran qe_openmp_c)
###########################################################
if(QE_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED Fortran C)
    target_link_libraries(qe_openmp_fortran INTERFACE OpenMP::OpenMP_Fortran)
    target_link_libraries(qe_openmp_c INTERFACE OpenMP::OpenMP_C)
endif(QE_ENABLE_OPENMP)
if(QE_ENABLE_OFFLOAD)
    target_compile_definitions(qe_openmp_fortran INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:__OPENMP_GPU>")
endif()

###########################################################
# OpenACC
# The following targets will be defined:
add_library(qe_openacc_fortran INTERFACE)
add_library(qe_openacc_c INTERFACE)
qe_install_targets(qe_openacc_fortran qe_openacc_c)
###########################################################
if(QE_ENABLE_OPENACC)
    find_package(OpenACC REQUIRED Fortran C)
    target_link_libraries(qe_openacc_fortran INTERFACE OpenACC::OpenACC_Fortran)
    target_link_libraries(qe_openacc_c INTERFACE OpenACC::OpenACC_C)
endif(QE_ENABLE_OPENACC)

############################################################
# Compiler vendor specific options
############################################################
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    include(GNUFortranCompiler)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI" OR CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
    include(NVFortranCompiler)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Cray")
    include(CrayFortranCompiler)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    include(IntelFortranCompiler)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "XL")
    include(IBMFortranCompiler)
endif()

if(QE_ENABLE_STATIC_BUILD)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
endif()

###########################################################
# CUDA
###########################################################
if(QE_ENABLE_CUDA OR QE_ENABLE_PROFILE_NVTX)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI" OR CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
        add_library(CUDA::cufft INTERFACE IMPORTED)
        set_target_properties(CUDA::cufft PROPERTIES INTERFACE_LINK_LIBRARIES "${CUDA_FLAG}lib=cufft")
        add_library(CUDA::cublas INTERFACE IMPORTED)
        set_target_properties(CUDA::cublas PROPERTIES INTERFACE_LINK_LIBRARIES "${CUDA_FLAG}lib=cublas")
        add_library(CUDA::cusolver INTERFACE IMPORTED)
        set_target_properties(CUDA::cusolver PROPERTIES INTERFACE_LINK_LIBRARIES "${CUDA_FLAG}lib=cusolver")
        add_library(CUDA::curand INTERFACE IMPORTED)
        set_target_properties(CUDA::curand PROPERTIES INTERFACE_LINK_LIBRARIES "${CUDA_FLAG}lib=curand")
        if(QE_ENABLE_PROFILE_NVTX)
            add_library(CUDA::nvToolsExt INTERFACE IMPORTED)
            set_target_properties(CUDA::nvToolsExt PROPERTIES INTERFACE_LINK_LIBRARIES "-cuda;libnvToolsExt.so")
            set(CMAKE_REQUIRED_LIBRARIES "-cuda;libnvToolsExt.so")
            check_function_exists(nvtxRangePushEx NVTX_FOUND)
            unset(CMAKE_REQUIRED_LIBRARIES)
            if(NOT NVTX_FOUND)
                message(FATAL_ERROR "Check nvtxRangePushEx in libnvToolsExt.so failed")
            endif()
        endif()
    else()
        find_package(CUDAToolkit REQUIRED)
    endif()
endif(QE_ENABLE_CUDA OR QE_ENABLE_PROFILE_NVTX)

###########################################################
# MPI
# The following targets will be defined:
add_library(qe_mpi_fortran INTERFACE)
qe_install_targets(qe_mpi_fortran)
###########################################################
if(QE_ENABLE_MPI)
    find_package(MPI REQUIRED Fortran)
    target_link_libraries(qe_mpi_fortran
        INTERFACE MPI::MPI_Fortran)
    message(STATUS "MPI settings used by CTest")
    message("     MPIEXEC_EXECUTABLE : ${MPIEXEC_EXECUTABLE}")
    message("     MPIEXEC_NUMPROC_FLAG : ${MPIEXEC_NUMPROC_FLAG}")
    message("     MPIEXEC_PREFLAGS : ${MPIEXEC_PREFLAGS}")
    string(REPLACE ";" " " MPIEXEC_PREFLAGS_PRINT "${MPIEXEC_PREFLAGS}")
    message("   Tests run as : ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} <NUM_PROCS> ${MPIEXEC_PREFLAGS_PRINT} <EXECUTABLE>")
endif(QE_ENABLE_MPI)

###########################################################
# Git
###########################################################
find_package(Git 2.13 REQUIRED)
if(EXISTS ${qe_SOURCE_DIR}/.git)
  message(STATUS "Source files are cloned from a git repository.")
  set(IS_GIT_PROJECT 1)
  include(GitInfo)
else()
  message(STATUS "Source files are not cloned from a git repository.")
endif()

###########################################################
# Lapack
# The following targets will be defined:
add_library(qe_lapack INTERFACE)
qe_install_targets(qe_lapack)
###########################################################
if(NOT QE_LAPACK_INTERNAL)
    if(NOT BLA_VENDOR)
        if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64.*")
            message(STATUS "Trying to find LAPACK from Intel MKL")
            if(QE_ENABLE_OPENMP)
                SET(BLA_VENDOR Intel10_64lp)
            else()
                SET(BLA_VENDOR Intel10_64lp_seq)
            endif()
            find_package(LAPACK)
        elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "i686.*|i386.*|x86.*")
            message(STATUS "Trying to find LAPACK from Intel MKL - 32bit")
            SET(BLA_VENDOR Intel10_32)
            find_package(LAPACK)
        elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64.*|AARCH64.*|arm64.*|ARM64.*)")
            message(STATUS "Trying to find LAPACK from ARM Performance Library")
            if(QE_ENABLE_OPENMP)
                SET(BLA_VENDOR Arm_mp)
            else()
                SET(BLA_VENDOR Arm)
            endif()
            find_package(LAPACK)
        endif()
        if(NOT LAPACK_FOUND)
            message(STATUS "Trying to find alternative LAPACK libraries")
            SET(BLA_VENDOR All)
            if(QE_ENABLE_OPENMP)
                set(CMAKE_REQUIRED_LINK_OPTIONS ${OpenMP_Fortran_FLAGS})
            endif()
            find_package(LAPACK)
            unset(CMAKE_REQUIRED_LINK_OPTIONS)
        endif()
    else()
        if(QE_ENABLE_OPENMP)
            set(CMAKE_REQUIRED_LINK_OPTIONS ${OpenMP_Fortran_FLAGS})
        endif()
        find_package(LAPACK)
        unset(CMAKE_REQUIRED_LINK_OPTIONS)
    endif()
    if(LAPACK_FOUND)
        list(APPEND _lapack_libs
            ${BLAS_LIBRARIES}
            ${BLAS_LINKER_FLAGS}
            ${LAPACK_LIBRARIES}
            ${LAPACK_LINKER_FLAGS})
        if(QE_ENABLE_OPENMP)
            list(APPEND _lapack_libs ${OpenMP_Fortran_LIBRARIES})
        endif()
        list(REMOVE_DUPLICATES "${_lapack_libs}")
        message(STATUS "Found LAPACK: ${_lapack_libs}")
        target_link_libraries(qe_lapack INTERFACE ${_lapack_libs})
        set(CMAKE_REQUIRED_LIBRARIES ${_lapack_libs})
        check_fortran_function_exists(zhpev ZHPEV_FOUND)
        unset(CMAKE_REQUIRED_LIBRARIES)
        if(NOT ZHPEV_FOUND)
          unset(ZHPEV_FOUND CACHE)
          message(FATAL_ERROR "Incomplete LAPACK! function zhpev not found!")
        endif()
    else()
        message(FATAL_ERROR "Failed to find a complete set of external BLAS/LAPACK library by FindLAPACK. "
                            "Variables controlling FindLAPACK can be found at CMake online documentation. "
                            "Alternatively, '-DQE_LAPACK_INTERNAL=ON' may be used to enable reference LAPACK "
                            "at a performance loss compared to optimized libraries.")
    endif()
else()
    message(WARNING "Internal reference LAPACK is enabled! It is less performant than vendor optimized libraries.")
    if(CMAKE_Fortran_COMPILER_ID MATCHES "XL")
        message(FATAL_ERROR "IBM XL compilers cannot build internal LAPACK with QE "
                            "due to the conflict in flags for free vs fixed format. "
                            "Please use an optimized LAPACK or build internal reference LAPACK separately.")
    endif()
    message(STATUS "Installing LAPACK via submodule")
    qe_git_submodule_update(external/lapack)
    add_subdirectory(external/lapack)
    target_link_libraries(qe_lapack INTERFACE lapack)
    # make lapack ready for other external libraries like mbd
    set(LAPACK_LIBRARIES lapack)
endif()

###########################################################
# SCALAPACK
# The following targets will be defined:
add_library(qe_scalapack INTERFACE)
qe_install_targets(qe_scalapack)
###########################################################
if(QE_ENABLE_SCALAPACK)
    find_package(SCALAPACK REQUIRED QUIET)
    message(STATUS "Found SCALAPACK: ${SCALAPACK_LIBRARIES};${SCALAPACK_LINKER_FLAGS}")
    target_link_libraries(qe_scalapack
        INTERFACE
            ${SCALAPACK_LIBRARIES}
            ${SCALAPACK_LINKER_FLAGS})
    if(QE_ENABLE_SCALAPACK_QRCP)
        include(CheckFortranFunctionExists)
        set(CMAKE_REQUIRED_LIBRARIES "${SCALAPACK_LIBRARIES}")
        check_fortran_function_exists("pzgeqpf" SCALAPACK_PZGEQPF_WORKS)
        unset(CMAKE_REQUIRED_LIBRARIES)
        if(SCALAPACK_PZGEQPF_WORKS)
            message(STATUS "Found pzgeqpf, add ScaLAPACK pzgeqpf macro")
            qe_add_global_compile_definitions(__SCALAPACK_QRCP)
        else()
            message(FATAL_ERROR "QE_ENABLE_SCALAPACK_QRCP requested but the current ScaLAPACK installation doesn't contain pzgeqpf!")
        endif()
    endif()
endif(QE_ENABLE_SCALAPACK)

###########################################################
# ELPA
# The following targets will be defined:
add_library(qe_elpa INTERFACE)
qe_install_targets(qe_elpa)
###########################################################
if(QE_ENABLE_ELPA)
    find_package(ELPA REQUIRED)

    # Check if ELPA version is compatible with QE
    if(ELPA_VERSION VERSION_GREATER_EQUAL "2018.11")
        set(QE_ELPA_DEFINITIONS __ELPA)
    elseif(ELPA_VERSION VERSION_GREATER_EQUAL "2016.11")
        set(QE_ELPA_DEFINITIONS __ELPA_2016)
    elseif(ELPA_VERSION VERSION_GREATER_EQUAL "2015")
        set(QE_ELPA_DEFINITIONS __ELPA_2015)
    else()
        message(FATAL_ERROR "ELPA verion ${ELPA_VERSION} is not supported.")
    endif()
    message(STATUS "Add ELPA flag : ${QE_ELPA_DEFINITIONS}")
    qe_add_global_compile_definitions(${QE_ELPA_DEFINITIONS})

    # Add link libraries and include directories
    target_link_libraries(qe_elpa
        INTERFACE
            ${ELPA_LIBRARIES}
            ${ELPA_LIBRARIES_DEP}
            ${ELPA_LINKER_FLAGS}
            qe_scalapack)
    target_include_directories(qe_elpa
        INTERFACE
            ${ELPA_Fortran_MODS_DIR}
            ${ELPA_INCLUDE_DIRS}
            ${ELPA_INCLUDE_DIRS_DEP})
endif(QE_ENABLE_ELPA)

###########################################################
# LIBXC
## The following targets will be defined:
add_library(qe_external_libxc INTERFACE)
qe_install_targets(qe_external_libxc)
###########################################################
if(QE_ENABLE_LIBXC)
    target_compile_definitions(qe_external_libxc INTERFACE "__LIBXC")

    find_package(Libxc 5.1.2 COMPONENTS Fortran)
    if (NOT Libxc_FOUND)
        message(STATUS "Libxc searching failed in CMake Module mode, trying Config mode")
        find_package(Libxc 5.1.2 COMPONENTS Fortran CONFIG)
    endif()

    if(${Libxc_FOUND})
        if (${Libxc_VERSION} VERSION_GREATER_EQUAL "5.1.2" )
            message(STATUS "Libxc version ${Libxc_VERSION} found.")
        else()
            message(FATAL_ERROR "Libxc version ${Libxc_VERSION} found. "
                                "CMake compilation of QE tested for libxc v5.1.2 or later only.")
        endif()
        target_link_libraries(qe_external_libxc INTERFACE Libxc::xcf03)
        target_include_directories(qe_external_libxc INTERFACE ${Libxc_INCLUDE_DIR})
        target_compile_definitions(qe_external_libxc INTERFACE ${Libxc_DEFINITIONS})
    else()
         message(FATAL_ERROR "Failed to find Libxc package (>=5.1.2) with Fortran enabled.")
    endif()
endif(QE_ENABLE_LIBXC)

###########################################################
# HDF5
# The following targets will be defined:
add_library(qe_hdf5_fortran INTERFACE)
add_library(qe_hdf5_c INTERFACE)
qe_install_targets(qe_hdf5_fortran qe_hdf5_c)
########################################################### 
if(QE_ENABLE_HDF5)
    if(QE_ENABLE_MPI)
        option(HDF5_PREFER_PARALLEL "Prefer parallel HDF5" ON)
    endif()
    if(QE_ENABLE_STATIC_BUILD)
        set(HDF5_USE_STATIC_LIBRARIES TRUE)
    endif()
    find_package(HDF5 REQUIRED Fortran C)
    if(NOT HDF5_FOUND)
        message(FATAL_ERROR "HDF5 Fortran interface has not been found!")
    endif()

    if (NOT HDF5_IS_PARALLEL OR NOT QE_ENABLE_MPI)
        message(STATUS "Serial HDF5 enabled!")
        qe_add_global_compile_definitions(__HDF5_SERIAL)
    else()
        message(STATUS "Parallel HDF5 enabled!")
    endif()

    target_link_libraries(qe_hdf5_fortran
        INTERFACE
            ${HDF5_Fortran_LIBRARIES})
    target_include_directories(qe_hdf5_fortran
        INTERFACE
            ${HDF5_Fortran_INCLUDE_DIRS})
    target_compile_definitions(qe_hdf5_fortran
        INTERFACE
            ${HDF5_Fortran_DEFINITIONS})

    target_link_libraries(qe_hdf5_c
        INTERFACE
            ${HDF5_C_LIBRARIES})
    target_include_directories(qe_hdf5_c
        INTERFACE
            ${HDF5_C_INCLUDE_DIRS})
    target_compile_definitions(qe_hdf5_c
        INTERFACE
            ${HDF5_C_DEFINITIONS})
endif(QE_ENABLE_HDF5)

###########################################################
# Tests
# Any executable target marked as test runner via
# 'add_test()' will be run by the 'test' make target
###########################################################
if(QE_ENABLE_TEST)
    include(cmake/unit_test.cmake)
    enable_testing()
endif(QE_ENABLE_TEST)

###########################################################
# PROFILERS LIBRARIES
# the target for profiler libray will be defined if
# some profiler is enabled
add_library(qe_ext_prof_tool INTERFACE)
qe_install_targets(qe_ext_prof_tool)
###########################################################
# this should work with nvfortran
if(QE_ENABLE_PROFILE_NVTX)
   target_link_libraries(qe_ext_prof_tool
	INTERFACE
        CUDA::nvToolsExt)
endif(QE_ENABLE_PROFILE_NVTX)

###########################################################
# Components
###########################################################
add_subdirectory(external)
add_subdirectory(FFTXlib)
add_subdirectory(UtilXlib)
add_subdirectory(Modules)
add_subdirectory(LAXlib)
add_subdirectory(XClib)
add_subdirectory(KS_Solvers)
add_subdirectory(dft-d3)
add_subdirectory(PW)
add_subdirectory(CPV)
add_subdirectory(atomic)
add_subdirectory(upflib)
add_subdirectory(COUPLE)
add_subdirectory(LR_Modules)
add_subdirectory(PHonon)
add_subdirectory(PP)
add_subdirectory(EPW)
add_subdirectory(GWW)
add_subdirectory(HP)
add_subdirectory(NEB)
add_subdirectory(PWCOND)
add_subdirectory(TDDFPT)
add_subdirectory(XSpectra)
add_subdirectory(QEHeat)
add_subdirectory(KCW)
add_subdirectory(GUI)
if(QE_ENABLE_DOC)
    add_subdirectory(Doc)
endif()

###########################################################
# Tests
###########################################################
if(QE_ENABLE_TEST)
  message(STATUS "Enabling tests in test-suite")
  add_subdirectory(test-suite)
endif()

###########################################################
# Pkg-config
###########################################################
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/quantum_espresso.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/quantum_espresso.pc
    @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/quantum_espresso.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

###########################################################
# Exports
###########################################################
install(EXPORT qeTargets
        FILE qeTargets.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/qe)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    qeConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion)

configure_file(cmake/qeConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/qeConfig.cmake @ONLY)

install(FILES
            ${CMAKE_CURRENT_BINARY_DIR}/qeConfigVersion.cmake
            ${CMAKE_CURRENT_BINARY_DIR}/qeConfig.cmake
        DESTINATION
            ${CMAKE_INSTALL_LIBDIR}/cmake/qe)

###########################################################
# Dependency graph generation
# Defines the custom target 'depgraph'
###########################################################
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/CMakeGraphVizOptions.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/CMakeGraphVizOptions.cmake COPYONLY)
add_custom_target(depgraph
    "${CMAKE_COMMAND}" "--graphviz=depgraph" .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")

###########################################################
# Custom make targets
# The collection 'pwall' is defined here
# Each of the followings is defined inside its subdirectory
#   pw ph pp hp pwcond neb cp tddfpt gwl ld1 upf epw
#   xspectra couple all_currents
###########################################################
add_custom_target(pwall
    DEPENDS
        pw
        ph
        pp
        pwcond
        neb
    COMMENT
        "same as \"make pw ph pp pwcond neb\"")
