if(BLAS_FOUND OR LAPACK_FOUND)
    foreach(_lib ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
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
                    "fftw3.h"
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
## FIXME: In the 'fft_scalar.ARM_LIB.f90' should be removed the call to the subroutine 'zfft1mx'
##        otherwise the ARMPL throw the following error: undefined reference to `zfft1mx_'
##        This subroutine is not more implemented in the ARMPL (ARMPL misses the symbol)
##
#        elseif(_lib MATCHES "armpl")
#            get_filename_component(_dir_l1 ${_lib} DIRECTORY)
#            get_filename_component(_dir_l2 ${_dir_l1} DIRECTORY)
#            get_filename_component(_dir_l3 ${_dir_l2} DIRECTORY)
#
#            find_path(VendorFFTW_INCLUDE_DIRS
#                NAMES
#                    "fftw3.f"
#                    "fftw3.h"
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
#            list(APPEND VendorFFTW_LIBRARIES 
#                ${BLAS_LIBRARIES}
#                ${BLAS_LINKER_FLAGS}
#                ${LAPACK_LIBRARIES}
#                ${LAPACK_LINKER_FLAGS})
#            list(REMOVE_DUPLICATES "${VendorFFTW_LIBRARIES}")
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
                    "fftw3.h"
                HINTS
                    ${_dir_l1} ${_dir_l2} ${_dir_l3}
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
    endforeach()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(VendorFFTW 
    REQUIRED_VARS 
        VendorFFTW_LIBRARIES 
        VendorFFTW_INCLUDE_DIRS 
        VendorFFTW_ID)

mark_as_advanced(VendorFFTW_LIBRARIES VendorFFTW_ID)