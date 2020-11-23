set(sources_davidson
    Davidson_RCI/david_rci.f90
    Davidson/cegterg.f90
    Davidson/regterg.f90)

set(sources_davidson_rci
    Davidson_RCI/david_rci.f90)

set(sources_cg
    CG/ccgdiagg.f90
    CG/rcgdiagg.f90)

set(sources_ppcg
    PPCG/ppcg_gamma.f90
    PPCG/ppcg_k.f90)

set(sources_paro
    ParO/bpcg_gamma.f90
    ParO/bpcg_k.f90
    ParO/paro_gamma.f90
    ParO/pcg_gamma.f90
    ParO/paro_k.f90
    ParO/pcg_k.f90
    ParO/paro_gamma_new.f90
    ParO/paro_k_new.f90)

set(sources_dense
    DENSE/rotate_HSpsi_gamma.f90
    DENSE/rotate_HSpsi_k.f90
    DENSE/rotate_wfc_gamma.f90
    DENSE/rotate_wfc_k.f90)

qe_add_library(qe_kssolver_davidson ${sources_davidson})
add_library(QE::KS_Solvers::Davidson ALIAS qe_kssolver_davidson)
target_link_libraries(qe_kssolver_davidson
    PRIVATE
        QE::OpenMP_Fortran
        QE::LAX
        QE::UTILX
        QE::MPI_Fortran)

qe_add_library(qe_kssolver_davidsonrci ${sources_davidson_rci})
add_library(QE::KS_Solvers::Davidson_RCI ALIAS qe_kssolver_davidsonrci)
target_link_libraries(qe_kssolver_davidsonrci
    PRIVATE
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
        QE::SCALAPACK
        QE::OpenMP_Fortran
        QE::LAX
        QE::UTILX
        QE::MPI_Fortran)

qe_add_library(qe_kssolver_dense ${sources_dense})
add_library(QE::KS_Solvers::DENSE ALIAS qe_kssolver_dense)
target_link_libraries(qe_kssolver_dense
    PRIVATE
        QE::LAX
        QE::UTILX
        QE::MPI_Fortran)

qe_add_library(qe_kssolver_paro ${sources_paro})
add_library(QE::KS_Solvers::ParO ALIAS qe_kssolver_paro)
target_link_libraries(qe_kssolver_paro
    PRIVATE
        QE::KS_Solvers::DENSE
        QE::LAX
        QE::UTILX
        QE::MPI_Fortran)

###########################################################

qe_install_targets(
    qe_kssolver_davidson
    qe_kssolver_davidsonrci
    qe_kssolver_cg
    qe_kssolver_ppcg
    qe_kssolver_paro
    qe_kssolver_dense)
