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

qe_add_library(qe_utilx ${sources})
add_library(QE::UTILX ALIAS qe_utilx)
target_link_libraries(qe_utilx
    PUBLIC
        QE::MPI_Fortran)
qe_install_target(qe_utilx)
