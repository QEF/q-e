# Add the directory of the first library file to RPATH/RUNPATH if that directory is not part of OS default locations.
# This is useful for locating dynamically linked libraries like OpenMP runtime libraries.
function(AddRPATH LIBRARY_NAME)
  if(ARGC GREATER 1)
    list(GET ARGV 1 LIBRARY_FILE)
    cmake_path(GET LIBRARY_FILE PARENT_PATH LIBRARY_DIR)
    list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES ${LIBRARY_DIR} isSystemDir)
    if("${isSystemDir}" STREQUAL "-1")
      message(STATUS "Append ${LIBRARY_NAME} runtime library path ${LIBRARY_DIR}")
      # Make sure that packages with depedency on this library actually find this one
      # For example, BLAS/LAPACK library may depend on OpenMP runtime libraries.
      # Spack copies compiler runtime libraries to create gcc-runtime/intel-oneapi-runtime packages
      # Under certain circumstances, the following bad mix can happen:
      #     FindOpenMP.cmake found libiomp5 provided by the compiler while
      #     FindBLAS.cmake found the copied libiomp5 due to intel-oneapi-runtime module
      # we have control of neither FindOpenMP.cmake nor FindBLAS.cmake
      # The following lines ensure FindBLAS.cmake follows what FindOpenMP.cmake found.
      list(APPEND CMAKE_LIBRARY_PATH ${LIBRARY_DIR})
      set(CMAKE_LIBRARY_PATH "${CMAKE_LIBRARY_PATH}" PARENT_SCOPE)
      # Build and install RPATH
      list(APPEND CMAKE_BUILD_RPATH ${LIBRARY_DIR})
      list(APPEND CMAKE_INSTALL_RPATH ${LIBRARY_DIR})
      set(CMAKE_BUILD_RPATH "${CMAKE_BUILD_RPATH}" PARENT_SCOPE)
      set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_RPATH}" PARENT_SCOPE)
    endif()
  endif()
endfunction()
