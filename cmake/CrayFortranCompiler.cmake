# Check compiler version
if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 14.0.0)
  message(FATAL_ERROR "Requires CCE 14.0.0 or higher ")
endif()

message(WARNING "The Cray Fortran compiler is not ready for QE production use.")

qe_add_global_compile_definitions(__CRAY)

# set preprocessor specific flag
set(Fortran_PREPROCESSOR_FLAGS "-eZ")

if(QE_ENABLE_OPENACC)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -hacc")
else()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -hnoacc")
endif()

if(QE_ENABLE_OPENMP)
  #
else()
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -hnoomp")
endif()
