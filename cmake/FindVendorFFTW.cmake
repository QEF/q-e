if(LAPACK_FOUND)
    foreach(_lib ${LAPACK_LIBRARIES})
        if(_lib MATCHES "mkl")
            get_filename_component(_dir_l1 ${_lib} DIRECTORY)
            get_filename_component(_dir_l2 ${_dir_l1} DIRECTORY)
            get_filename_component(_dir_l3 ${_dir_l2} DIRECTORY)

            find_path(VendorFFTW_INCLUDE_MKL_DFTI
                NAMES
                    "mkl_dfti.f90"
                HINTS
                    ${_dir_l1} ${_dir_l2} ${_dir_l3}
                PATH_SUFFIXES
                    "include"
                    "fftw"
                    "include/fftw"
                NO_DEFAULT_PATH
            )
            find_path(VendorFFTW_INCLUDE_FFTW3
                NAMES
                    "fftw3.f"
                HINTS
                    ${_dir_l1} ${_dir_l2} ${_dir_l3}
                PATH_SUFFIXES
                    "include"
                    "fftw"
                    "include/fftw"
                NO_DEFAULT_PATH
            )
            set(VendorFFTW_INCLUDE_DIRS ${VendorFFTW_INCLUDE_MKL_DFTI} ${VendorFFTW_INCLUDE_FFTW3})

            add_library(VendorFFTW INTERFACE IMPORTED)
            set(VendorFFTW_LIBRARIES ${LAPACK_LIBRARIES})
            target_link_libraries(VendorFFTW INTERFACE ${LAPACK_LIBRARIES})
            target_include_directories(VendorFFTW INTERFACE ${VendorFFTW_INCLUDE_DIRS})
            set(VendorFFTW_ID "Intel")

            break()
# FIXME: undefined reference to `zfft1mx_'
#        elseif(_lib MATCHES "armpl")
#            get_filename_component(_dir_l1 ${_lib} DIRECTORY)
#            get_filename_component(_dir_l2 ${_dir_l1} DIRECTORY)
#            get_filename_component(_dir_l3 ${_dir_l2} DIRECTORY)
#
#            find_path(VendorFFTW_INCLUDE_DIRS
#                NAMES
#                    "fftw3.f"
#                HINTS
#                    ${_dir_l1} ${_dir_l2} ${_dir_l3}
#                PATH_SUFFIXES
#                    "include"
#                    "fftw"
#                    "include/fftw"
#                NO_DEFAULT_PATH
#            )
#
#            add_library(VendorFFTW INTERFACE IMPORTED)
#            set(VendorFFTW_LIBRARIES ${LAPACK_LIBRARIES})
#            target_link_libraries(VendorFFTW INTERFACE ${VendorFFTW_LIBRARIES})
#            set(VendorFFTW_ID "Arm")
#
#            break()
        elseif(_lib MATCHES "essl")
            get_filename_component(_dir_l1 ${_lib} DIRECTORY)
            get_filename_component(_dir_l2 ${_dir_l1} DIRECTORY)
            get_filename_component(_dir_l3 ${_dir_l2} DIRECTORY)
            find_path(VendorFFTW_INCLUDE_DIRS
                NAMES
                    "fftw3.f"
                HINTS
                    ${_dir_l1} ${_dir_l2} ${_dir_l3}
                PATH_SUFFIXES
                    "include"
                    "fftw"
                    "include/fftw"
                NO_DEFAULT_PATH
            )

            add_library(VendorFFTW INTERFACE IMPORTED)
            set(VendorFFTW_LIBRARIES ${LAPACK_LIBRARIES})
            target_link_libraries(VendorFFTW INTERFACE ${VendorFFTW_LIBRARIES})
            #target_include_directories(VendorFFTW INTERFACE ${VendorFFTW_INCLUDE_DIRS})
            set(VendorFFTW_ID "Essl")

            break()
        endif()
    endforeach()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(VendorFFTW REQUIRED_VARS VendorFFTW_LIBRARIES VendorFFTW_INCLUDE_DIRS VendorFFTW_ID)

mark_as_advanced(VendorFFTW_LIBRARIES VendorFFTW_ID)