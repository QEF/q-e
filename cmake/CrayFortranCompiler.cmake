# Check compiler version
if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 14.0.0)
  message(FATAL_ERROR "Requires CCE 14.0.0 or higher ")
endif()

message(WARNING "The Cray Fortran compiler is not ready for QE production use.")

qe_add_global_compile_definitions(__CRAY)

# set preprocessor specific flag
set(Fortran_PREPROCESSOR_FLAGS "-eZ")

if(NOT QE_ENABLE_OPENACC)
  target_compile_options(qe_openacc_fortran INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-hnoacc>")
endif()

if(QE_ENABLE_OFFLOAD)
  execute_process(
    COMMAND ${CMAKE_Fortran_COMPILER} -craype-verbose --version
    RESULT_VARIABLE VERSION_QUERY_RETURN
    OUTPUT_VARIABLE VERSION_QUERY_OUTPUT)

  if(VERSION_QUERY_RETURN)
    message(WARNING "Command `${CMAKE_Fortran_COMPILER} -craype-verbose --version` failed!")
  elseif(NOT VERSION_QUERY_OUTPUT MATCHES "accel=")
    message(FATAL_ERROR "Cannot find -haccel=<gpu_arc> option being used by the ftn compiler wrapper. "
                        "Make sure the GPU architecture module is loaded."
                        "Command `${CMAKE_Fortran_COMPILER} -craype-verbose --version` returns\n"
                        "${VERSION_QUERY_OUTPUT}")
  endif()
endif()

if(QE_ENABLE_OPENMP)
  target_link_options(qe_openmp_fortran INTERFACE "$<$<LINK_LANGUAGE:Fortran>:-homp>")
else()
  target_compile_options(qe_openmp_fortran INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-hnoomp>")
endif()
