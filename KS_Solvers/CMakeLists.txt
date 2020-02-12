set(sources_davidson_param
    Davidson/david_param.f90)

set(sources_davidson
    Davidson_RCI/david_rci.f90
    Davidson/cegterg.f90
    Davidson/regterg.f90
)

set(sources_davidson_rci
    Davidson_RCI/david_rci.f90
)

set(sources_cg
    CG/cg_param.f90
    CG/ccgdiagg.f90
    CG/rcgdiagg.f90
    CG/rotate_wfc_gamma.f90
    CG/rotate_wfc_k.f90
)

set(sources_ppcg
    PPCG/ppcg_param.f90
    PPCG/ppcg_gamma.f90
    PPCG/ppcg_k.f90
)

# Both sources_davidson and sources_davidson_rci depend on this
# module, we need to factor it out to a separate target to ensure
# we do not break parallel builds
qe_add_library(qe_kssolver_davidson_param ${sources_davidson_param})
add_library(QE::KS_Solvers::DavidsonParam ALIAS qe_kssolver_davidson_param)
target_link_libraries(qe_kssolver_davidson_param
    PRIVATE
        QE::UTILX)

qe_add_library(qe_kssolver_davidson ${sources_davidson})
add_library(QE::KS_Solvers::Davidson ALIAS qe_kssolver_davidson)
target_link_libraries(qe_kssolver_davidson
    PRIVATE
        QE::KS_Solvers::DavidsonParam
        QE::LAX
        QE::UTILX
        QE::MPI_Fortran)

qe_add_library(qe_kssolver_davidsonrci ${sources_davidson_rci})
add_library(QE::KS_Solvers::Davidson_RCI ALIAS qe_kssolver_davidsonrci)
target_link_libraries(qe_kssolver_davidsonrci
    PRIVATE
        QE::KS_Solvers::DavidsonParam
        QE::LAX
        QE::UTILX
        QE::MPI_Fortran)

qe_add_library(qe_kssolver_cg ${sources_cg})
add_library(QE::KS_Solvers::CG ALIAS qe_kssolver_cg)
target_link_libraries(qe_kssolver_cg
    PRIVATE
        QE::LAX
        QE::UTILX
        QE::MPI_Fortran)

qe_add_library(qe_kssolver_ppcg ${sources_ppcg})
add_library(QE::KS_Solvers::PPCG ALIAS qe_kssolver_ppcg)
target_link_libraries(qe_kssolver_ppcg
    PRIVATE
        QE::LAX
        QE::UTILX
        QE::MPI_Fortran)

###########################################################

qe_install_targets(
    qe_kssolver_davidson_param
    qe_kssolver_davidson
    qe_kssolver_davidsonrci
    qe_kssolver_cg
    qe_kssolver_ppcg
)
