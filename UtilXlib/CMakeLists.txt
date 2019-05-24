set(sources
    clocks_handler.f90
    cuda_util.f90
    divide.f90
    data_buffer.f90
    error_handler.f90
    find_free_unit.f90
    fletcher32_mod.f90
    mem_counter.f90
    mp.f90
    mp_base.f90
    mp_base_gpu.f90
    mp_bands_util.f90
    parallel_include.f90
    util_param.f90
    thread_util.f90
)

# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)

add_fortran_library(qe_utilxlib ${sources})
add_library(QE::UtilXLib ALIAS qe_utilxlib)
target_link_libraries(qe_utilxlib
    PUBLIC
        QE::Features
        QE::MPI_Fortran)
