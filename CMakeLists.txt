cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
# CMake < 3.13 cannot install external targets:
# https://gitlab.kitware.com/cmake/cmake/merge_requests/2152

# CMake < 3.14 generates incorrect dependency graphs with
# alias targets:
# https://gitlab.kitware.com/cmake/cmake/merge_requests/2521

# Policy CMP0048: The project() command manages VERSION variables
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)

project(qe
    VERSION 6.6.0
    DESCRIPTION "ESPRESSO: opEn-Source Package for Research in Electronic Structure, Simulation, and Optimization"
    LANGUAGES Fortran C)

if(${qe_BINARY_DIR} STREQUAL ${qe_SOURCE_DIR})
  message(FATAL_ERROR "QE source folder cannot be safely used as a build folder!")
endif()

# CMake < v3.18 cannot discover the ARM Performance Library
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64.*|AARCH64.*|arm64.*|ARM64.*)")
    if(CMAKE_VERSION VERSION_LESS 3.18.0)
        message("-- CMake versions less then 3.18 cannot automatically discover the ARM Performance Library!")
    endif()
endif()

######################################################################
# Define the paths for static libraries and executables
######################################################################
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${qe_BINARY_DIR}/lib CACHE PATH "Single output directory for building all libraries.")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${qe_BINARY_DIR}/bin CACHE PATH "Single output directory for building all executables.")
set(QE_TESTS_DIR ${qe_BINARY_DIR}/tests/bin)


###########################################################
# Build helpers
###########################################################
set(PROJECT_CMAKE ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
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
option(QE_ENABLE_TEST
    "enable unit tests" ON)
option(QE_ENABLE_MPI
    "enable distributed execution support via MPI" ON)
option(QE_ENABLE_OPENMP
    "enable distributed execution support via OpenMP" OFF)
option(QE_ENABLE_TRACE
    "enable execution tracing output" OFF)
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
option(QE_ENABLE_ELPA
    "enable ELPA execution units" OFF)
option(QE_ENABLE_HDF5
    "enable HDF5 data collection" OFF)
option(QE_ENABLE_CUDA
    "enable CUDA acceleration on NVIDIA GPUs" OFF)
option(QE_ENABLE_DOC
    "enable documentation building" OFF)

# TODO change all ifdefs throughout code base to match
# cmake options
# TODO symbols beginning with '__' followed by a capital
# character are reserved for standard library use (at
# least in C, not sure about Fortran), change all feature
# macros to avoid weird behaviors

# Disable all configuration headers used to be generated
# by configure (see <qe>/include/)
qe_add_global_compile_definitions(QE_NO_CONFIG_H)

if(QE_ENABLE_TRACE)
    qe_add_global_compile_definitions(__TRACE)
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
if(QE_ENABLE_MPI)
    # OMPI_SKIP_MPICXX: skip CXX APIs on openmpi, cause trouble to C APIs
    qe_add_global_compile_definitions(__MPI OMPI_SKIP_MPICXX)
endif()
if(QE_ENABLE_SCALAPACK)
    qe_add_global_compile_definitions(__SCALAPACK)
endif()
if(QE_ENABLE_HDF5)
    qe_add_global_compile_definitions(__HDF5)
endif()

# Feature checks
check_function_exists(mallinfo HAVE_MALLINFO)
if(HAVE_MALLINFO)
    qe_add_global_compile_definitions(HAVE_MALLINFO)
endif()

# Check options consistency
if(QE_ENABLE_ELPA AND NOT QE_ENABLE_SCALAPACK)
    message(FATAL_ERROR "ELPA requires SCALAPACK support, enable it with '-DQE_ENABLE_SCALAPACK=ON' or disable ELPA with '-DQE_ENABLE_ELPA=OFF'")
endif()
if(QE_ENABLE_SCALAPACK AND NOT QE_ENABLE_MPI)
    message(FATAL_ERROR "SCALAPACK requires MPI support, enable it with '-DQE_ENABLE_MPI=ON' or disable SCALAPACK with '-DQE_ENABLE_SCALAPACK=OFF'")
endif()
# if(QE_ENABLE_HDF5 AND NOT QE_ENABLE_MPI)
#    message(FATAL_ERROR "HDF5 requires MPI support, enable it with '-DQE_ENABLE_MPI=ON' or disable HDF5 with '-DQE_ENABLE_HDF5=OFF'")
# endif()


###########################################################
# language standard requirements
###########################################################
# TODO need to require all compilers using the same one
if(CMAKE_C_COMPILER_ID MATCHES "PGI")
  set(CMAKE_C_STANDARD 11)
  set(CMAKE_C_STANDARD_REQUIRED ON)
  set(CMAKE_C_EXTENSIONS OFF)
endif()

############################################################
## Compiler vendor specific options
############################################################
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    include(${PROJECT_CMAKE}/GNUFortranCompiler.cmake)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    include(${PROJECT_CMAKE}/NVFortranCompiler.cmake)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "XL")
    include(${PROJECT_CMAKE}/IBMFortranCompiler.cmake)
endif()

###########################################################
# OpenMP
# The following targets will be defined:
add_library(qe_openmp_fortran INTERFACE)
add_library(QE::OpenMP_Fortran ALIAS qe_openmp_fortran)
add_library(qe_openmp_c INTERFACE)
add_library(QE::OpenMP_C ALIAS qe_openmp_c)
qe_install_targets(qe_openmp_fortran qe_openmp_c)
###########################################################
if(QE_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED Fortran C)
    target_link_libraries(qe_openmp_fortran
        INTERFACE OpenMP::OpenMP_Fortran)
    target_link_libraries(qe_openmp_c
        INTERFACE OpenMP::OpenMP_C)
endif(QE_ENABLE_OPENMP)

###########################################################
# MPI
# The following targets will be defined:
add_library(qe_mpi_fortran INTERFACE)
add_library(QE::MPI_Fortran ALIAS qe_mpi_fortran)
qe_install_targets(qe_mpi_fortran)
###########################################################
if(QE_ENABLE_MPI)
    find_package(MPI REQUIRED Fortran)
    target_link_libraries(qe_mpi_fortran
        INTERFACE MPI::MPI_Fortran)
endif(QE_ENABLE_MPI)

###########################################################
# Lapack
# The following targets will be defined:
add_library(qe_lapack INTERFACE)
add_library(QE::LAPACK ALIAS qe_lapack)
qe_install_targets(qe_lapack)
#######################################################################
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
            find_package(LAPACK)
        endif()
    else()
        find_package(LAPACK)
    endif()
    if(LAPACK_FOUND)
        list(APPEND _lapack_libs
            ${BLAS_LIBRARIES}
            ${BLAS_LINKER_FLAGS}
            ${LAPACK_LIBRARIES}
            ${LAPACK_LINKER_FLAGS})
        list(REMOVE_DUPLICATES "${_lapack_libs}")
        message(STATUS "Found LAPACK: ${_lapack_libs}")
        target_link_libraries(qe_lapack INTERFACE ${_lapack_libs})
    else()
        message(FATAL_ERROR "Failed to find a complete set of external BLAS/LAPACK library by FindLAPACK. "
                            "Varaibles controlling FindLAPACK can be found at CMake online documentation. "
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
    message(STATUS "Installing QE::LAPACK via submodule")
    qe_git_submodule_update(external/lapack)
    add_subdirectory(external/lapack EXCLUDE_FROM_ALL)
    target_link_libraries(qe_lapack INTERFACE lapack)
    qe_fix_fortran_modules(lapack)
    qe_install_targets(lapack)
endif()

###########################################################
# SCALAPACK
# The following targets will be defined:
add_library(qe_scalapack INTERFACE)
add_library(QE::SCALAPACK ALIAS qe_scalapack)
qe_install_targets(qe_scalapack)
###########################################################
if(QE_ENABLE_SCALAPACK)
    find_package(SCALAPACK REQUIRED QUIET)
    message(STATUS "Found SCALAPACK: ${SCALAPACK_LIBRARIES};${SCALAPACK_LINKER_FLAGS}")
    target_link_libraries(qe_scalapack
        INTERFACE
            ${SCALAPACK_LIBRARIES}
            ${SCALAPACK_LINKER_FLAGS})
endif(QE_ENABLE_SCALAPACK)

###########################################################
# ELPA
# The following targets will be defined:
add_library(qe_elpa INTERFACE)
add_library(QE::ELPA ALIAS qe_elpa)
qe_install_targets(qe_elpa)
###########################################################
if(QE_ENABLE_ELPA)
    find_package(ELPA REQUIRED)

    # Check if ELPA version is compatible with QE
    if(ELPA_VERSION_STRING VERSION_GREATER_EQUAL "2018.11")
        set(QE_ELPA_DEFINITIONS __ELPA)
    elseif(ELPA_VERSION_STRING VERSION_GREATER_EQUAL "2016.11")
        set(QE_ELPA_DEFINITIONS __ELPA_2016)
    elseif(ELPA_VERSION_STRING VERSION_GREATER_EQUAL "2015")
        set(QE_ELPA_DEFINITIONS __ELPA_2015)
    else()
        message(FATAL_ERROR "ELPA verion ${ELPA_VERSION_STRING} is not supported.")
    endif()
    message(STATUS "Add ELPA flag : ${QE_ELPA_DEFINITIONS}")
    qe_add_global_compile_definitions(${QE_ELPA_DEFINITIONS})

    # Looking for module directory
    file(GLOB_RECURSE ELPA_MODS_FILES "${ELPA_INCLUDE_DIRS}/*.mod" "${ELPA_INCLUDE_DIRS}/../modules/*.mod")
    if(ELPA_MODS_FILES)
        list(GET ELPA_MODS_FILES 0 ELPA_MOD)
        get_filename_component(ELPA_MOD_DIR "${ELPA_MOD}" PATH)
        set(ELPA_INCLUDE_DIRS "${ELPA_MOD_DIR};${ELPA_INCLUDE_DIRS}")
    else()
        message(FATAL_ERROR "Fortran module directory of ELPA not found!")
    endif()

    # Add link libraries and include directories
    target_link_libraries(qe_elpa
        INTERFACE
            ${ELPA_LIBRARIES}
            ${ELPA_LIBRARIES_DEP}
            ${ELPA_LINKER_FLAGS})
    target_include_directories(qe_elpa
        INTERFACE
            ${ELPA_INCLUDE_DIRS}
            ${ELPA_INCLUDE_DIRS_DEP})
endif(QE_ENABLE_ELPA)

###########################################################
# HDF5
# The following targets will be defined:
add_library(qe_hdf5_fortran INTERFACE)
add_library(QE::HDF5_Fortran ALIAS qe_hdf5_fortran)
add_library(qe_hdf5_c INTERFACE)
add_library(QE::HDF5_C ALIAS qe_hdf5_c)
qe_install_targets(qe_hdf5_fortran qe_hdf5_c)
###########################################################
if(QE_ENABLE_HDF5)
    if(QE_ENABLE_MPI)
        option(HDF5_PREFER_PARALLEL "Prefer parallel HDF5" ON)
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
            ${HDF5_Fortran_LIBRARIES}
            ${HDF5_Fortran_HL_LIBRARIES})
    target_include_directories(qe_hdf5_fortran
        INTERFACE
            ${HDF5_Fortran_INCLUDE_DIRS})
    target_compile_definitions(qe_hdf5_fortran
        INTERFACE
            ${HDF5_Fortran_DEFINITIONS})

    target_link_libraries(qe_hdf5_c
        INTERFACE
            ${HDF5_C_LIBRARIES}
            ${HDF5_C_HL_LIBRARIES})
    target_include_directories(qe_hdf5_c
        INTERFACE
            ${HDF5_C_INCLUDE_DIRS})
    target_compile_definitions(qe_hdf5_c
        INTERFACE
            ${HDF5_C_DEFINITIONS})
endif(QE_ENABLE_HDF5)

###########################################################
# CUDA
# FIXME just a place holder for the moment
###########################################################
if(QE_ENABLE_CUDA)
    qe_add_global_compile_definitions(__CUDA)
endif(QE_ENABLE_CUDA)

###########################################################
# Tests
# Any executable target marked as test runner via
# 'add_test()' will be run by the 'test' make target
###########################################################
if(QE_ENABLE_TEST)
    enable_testing()
endif(QE_ENABLE_TEST)

###########################################################
# Components
###########################################################
add_subdirectory(external)
add_subdirectory(clib)
add_subdirectory(FFTXlib)
add_subdirectory(UtilXlib)
add_subdirectory(Modules)
add_subdirectory(LAXlib)
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
add_subdirectory(PlotPhon)
add_subdirectory(PWCOND)
add_subdirectory(QHA)
add_subdirectory(TDDFPT)
add_subdirectory(XSpectra)
if(QE_ENABLE_DOC)
    add_subdirectory(Doc)
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
        NAMESPACE qe::
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
###########################################################
add_custom_target(pw
    DEPENDS
        qe_pw_exe
        qe_pw_tools_ibrav2cell_exe
        qe_pw_tools_cell2ibrav_exe
        qe_pw_tools_ev_exe
        qe_pw_tools_kpoints_exe
        qe_pw_tools_pwi2xsf_exe
    COMMENT
        "basic code for scf, structure optimization, MD")

add_custom_target(ph
    DEPENDS
        qe_phonon_ph_exe
        qe_phonon_dynmat_exe
        qe_phonon_q2r_exe
        qe_phonon_matdyn_exe
        qe_phonon_q2qstar_exe
        qe_phonon_lambda_exe
        qe_phonon_alpha2f_exe
        qe_phonon_epa_exe
        qe_phonon_fqha_exe
        qe_plotphon_kforbands_exe
        qe_plotphon_bandstognuplot_exe
        qe_plotphon_eminmax_exe
    COMMENT
        "phonon code, Gamma-only and third-order derivatives")

add_custom_target(hp
    DEPENDS
        qe_hp_exe
    COMMENT
        "calculation of the Hubbard parameters from DFPT")

add_custom_target(pwcond
    DEPENDS
        qe_pwcond_exe
    COMMENT
        "ballistic conductance")

add_custom_target(neb
    DEPENDS
        qe_neb_exe
        qe_neb_pathinterpolation_exe
    COMMENT
        "code for Nudged Elastic Band method")

add_custom_target(pp
    DEPENDS
        qe_pp_exe
        qe_pp_opengrid_exe
        qe_pp_average_exe
        qe_pp_bands_exe
        qe_pp_dos_exe
        qe_pp_pawplot_exe
        qe_pp_planavg_exe
        qe_pp_plotband_exe
        qe_pp_plotproj_exe
        qe_pp_plotrho_exe
        qe_pp_pmw_exe
        qe_pp_xctest_exe
        qe_pp_projwfc_exe
        qe_pp_pw2wannier90_exe
        qe_pp_pw2critic_exe
        qe_pp_wfck2r_exe
        qe_pp_initial_state_exe
        qe_pp_pw2gw_exe
        qe_pp_sumpdos_exe
        qe_pp_epsilon_exe
        qe_pp_wannierham_exe
        qe_pp_wannierplot_exe
        qe_pp_molecularpdos_exe
        qe_pp_pw2bgw_exe
        qe_pp_fermivelocity_exe
        qe_pp_fermisurface_exe
        qe_pp_fermiproj_exe
        qe_pp_ppacf_exe
        qe_pp_benchmarklibxc_exe
    COMMENT
        "postprocessing programs")

add_custom_target(pwall
    DEPENDS
        pw
        ph
        pp
        pwcond
        neb
    COMMENT
        "same as \"make pw ph pp pwcond neb\"")

add_custom_target(cp
    DEPENDS
        qe_cpv_exe
        qe_cpv_manycp_exe
        qe_cpv_cppp_exe
        qe_cpv_wfdd_exe
    COMMENT
        "CP code: Car-Parrinello molecular dynamics")

add_custom_target(tddfpt
    DEPENDS
        qe_tddfpt_turbolanczos_exe
        qe_tddfpt_turbodavidson_exe
        qe_tddfpt_turboeels_exe
        qe_tddfpt_calculatespectrum_exe
    COMMENT
        "time dependent dft code")

add_custom_target(gwl
    DEPENDS
        qe_gww_util_grap_exe
        qe_gww_util_abcoefftoeps_exe
        qe_gww_util_memorypw4gww_exe
    COMMENT
        "GW with Lanczos chains")

add_custom_target(ld1
    DEPENDS
        qe_atomic_exe
    COMMENT
        "utilities for pseudopotential generation")

add_custom_target(upf
    DEPENDS
        #Library
        qe_upflib
        #Executables
        qe_upflib_virtual_v2_exe
        qe_upflib_upfconv_exe
    COMMENT
        "utilities for pseudopotential conversion")

add_custom_target(xspectra
    DEPENDS
        qe_xspectra_exe
        qe_xspectra_spectracorrection_exe
        qe_xspectra_molecularnexafs_exe
    COMMENT
        "X-ray core-hole spectroscopy calculations")

add_custom_target(couple
    DEPENDS
        qe_couple
    COMMENT
        "library interface for coupling to external codes")

add_custom_target(epw
    DEPENDS
        qe_epw_exe
    COMMENT
        "electron-Phonon Coupling with wannier functions")
