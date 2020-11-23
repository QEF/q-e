# This finder needs that BLAS and/or LAPACK have been founded!
if(BLAS_FOUND OR LAPACK_FOUND)
    # Trasform spaces in semicolons
    string(REPLACE " " ";" _blas_libs "${BLAS_LIBRARIES}")
    string(REPLACE " " ";" _lapack_libs "${LAPACK_LIBRARIES}")

    foreach(_lib ${_blas_libs} ${_lapack_libs})
        # Check only libraries and directories
        set(_search_path)
        if(EXISTS "${_lib}")
            # Get three parent directories of the library
            get_filename_component(_parent_l1 ${_lib} DIRECTORY)
            get_filename_component(_parent_l2 ${_parent_l1} DIRECTORY)
            get_filename_component(_parent_l3 ${_parent_l2} DIRECTORY)
            list(APPEND _search_path ${_lib} ${_parent_l1} ${_parent_l2} ${_parent_l3})
        else()
            # Remove first two characters '-L' (or '-l') and try again
            string(REGEX MATCH "/[^ ]*" _lib_clean ${_lib})
            if(EXISTS "${_lib_clean}")
                get_filename_component(_parent_l1 ${_lib_clean} DIRECTORY)
                get_filename_component(_parent_l2 ${_parent_l1} DIRECTORY)
                get_filename_component(_parent_l3 ${_parent_l2} DIRECTORY)
                list(APPEND _search_path ${_lib_clean} ${_parent_l1} ${_parent_l2} ${_parent_l3})
            else()
                # This is not a file or a directory
                continue()
            endif()
        endif()

        # Start to search vendor libraries and include directories
        # Try to find the Intel MKL
        if(NOT VendorFFTW_FIND_COMPONENTS OR "MKL" IN_LIST VendorFFTW_FIND_COMPONENTS)
            if(_search_path MATCHES "mkl")
                find_path(VendorFFTW_INCLUDE_MKL_DFTI
                    NAMES
                        "mkl_dfti.f90"
                    HINTS
                        ${_search_path}
                    PATH_SUFFIXES
                        "include"
                        "fftw"
                        "include/fftw"
                    NO_DEFAULT_PATH
                )
                find_path(VendorFFTW_INCLUDE_FFTW3
                    NAMES
                        "fftw3.f"
                        "fftw3.h"
                    HINTS
                        ${_search_path}
                    PATH_SUFFIXES
                        "include"
                        "fftw"
                        "include/fftw"
                    NO_DEFAULT_PATH
                )
                set(VendorFFTW_INCLUDE_DIRS ${VendorFFTW_INCLUDE_MKL_DFTI} ${VendorFFTW_INCLUDE_FFTW3})

                add_library(VendorFFTW INTERFACE IMPORTED)
                list(APPEND VendorFFTW_LIBRARIES 
                    ${BLAS_LIBRARIES}
                    ${BLAS_LINKER_FLAGS}
                    ${LAPACK_LIBRARIES}
                    ${LAPACK_LINKER_FLAGS})
                list(REMOVE_DUPLICATES "${VendorFFTW_LIBRARIES}")
                target_link_libraries(VendorFFTW INTERFACE ${VendorFFTW_LIBRARIES})
                target_include_directories(VendorFFTW INTERFACE ${VendorFFTW_INCLUDE_DIRS})
                set(VendorFFTW_ID "Intel")

                break()
            endif()
        endif()

        # Try to find the ARM Performance Library
        if(NOT VendorFFTW_FIND_COMPONENTS OR "ArmPL" IN_LIST VendorFFTW_FIND_COMPONENTS)
            if(_search_path MATCHES "armpl")
                find_path(VendorFFTW_INCLUDE_DIRS
                    NAMES
                        "fftw3.f"
                        "fftw3.h"
                    HINTS
                        ${_search_path}
                    PATH_SUFFIXES
                        "include"
                        "fftw"
                        "include/fftw"
                    NO_DEFAULT_PATH
                )

                add_library(VendorFFTW INTERFACE IMPORTED)
                list(APPEND VendorFFTW_LIBRARIES 
                    ${BLAS_LIBRARIES}
                    ${BLAS_LINKER_FLAGS}
                    ${LAPACK_LIBRARIES}
                    ${LAPACK_LINKER_FLAGS})
                list(REMOVE_DUPLICATES "${VendorFFTW_LIBRARIES}")
                target_link_libraries(VendorFFTW INTERFACE ${VendorFFTW_LIBRARIES})
                target_include_directories(VendorFFTW INTERFACE ${VendorFFTW_INCLUDE_DIRS})
                set(VendorFFTW_ID "Arm")

                break()
            endif()
        endif()

        # Try to find the IBM ESSL library
        if(NOT VendorFFTW_FIND_COMPONENTS OR "IBMESSL" IN_LIST VendorFFTW_FIND_COMPONENTS)
            if(_search_path MATCHES "essl")
                find_path(VendorFFTW_INCLUDE_DIRS
                    NAMES
                        "fftw3.f"
                        "fftw3.h"
                    HINTS
                        ${_search_path}
                    PATH_SUFFIXES
                        "include"
                        "fftw"
                        "include/fftw"
                    NO_DEFAULT_PATH
                )

                add_library(VendorFFTW INTERFACE IMPORTED)
                list(APPEND VendorFFTW_LIBRARIES 
                    ${BLAS_LIBRARIES}
                    ${BLAS_LINKER_FLAGS}
                    ${LAPACK_LIBRARIES}
                    ${LAPACK_LINKER_FLAGS})
                list(REMOVE_DUPLICATES "${VendorFFTW_LIBRARIES}")
                target_include_directories(VendorFFTW INTERFACE ${VendorFFTW_INCLUDE_DIRS})
                set(VendorFFTW_ID "IBMESSL")

                break()
            endif()
        endif()
    endforeach()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(VendorFFTW 
    REQUIRED_VARS 
        VendorFFTW_LIBRARIES 
        VendorFFTW_INCLUDE_DIRS 
        VendorFFTW_ID)

mark_as_advanced(VendorFFTW_LIBRARIES VendorFFTW_ID)