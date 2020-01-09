set(sources
    Davidson/david_param.f90
    Davidson/cegterg.f90
    Davidson/regterg.f90
    Davidson_RCI/david_rci.f90
    CG/cg_param.f90
    CG/ccgdiagg.f90
    CG/rcgdiagg.f90
    CG/rotate_wfc_gamma.f90
    CG/rotate_wfc_k.f90
    PPCG/ppcg_param.f90
    PPCG/ppcg_gamma.f90
    PPCG/ppcg_k.f90
)

# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)

qe_add_library(qe_solvers ${sources})
add_library(QE::Solvers ALIAS qe_solvers)
target_link_libraries(qe_solvers
    PUBLIC
        QE::LAX
        QE::UTILX
        QE::MPI_Fortran)
qe_install_target(qe_solvers)
