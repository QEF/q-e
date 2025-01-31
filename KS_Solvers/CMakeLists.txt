set(ks_headers
  ks_solver_interfaces.h)
set (ks_include_dir ${CMAKE_CURRENT_BINARY_DIR}/include)
#FIXME maybe better to move this part to an helper function 
foreach (in_h ${ks_headers})
  get_filename_component(in_h_basename ${in_h} NAME_WE)
  set (out_h "${ks_include_dir}/${in_h_basename}.fh")
  qe_preprocess_source("${CMAKE_CURRENT_SOURCE_DIR}/${in_h}" ${out_h})
  list (APPEND ks_out_headers ${out_h})
endforeach()
add_custom_target(qe_ks_headers 
  DEPENDS ${ks_out_headers}
  VERBATIM )

set(src_davidson
    Davidson_RCI/david_rci.f90
    Davidson/cegterg.f90
    Davidson/regterg.f90
    # GPU
    Davidson/regterg_gpu.f90
    Davidson/cegterg_gpu.f90)
qe_enable_cuda_fortran("${src_davidson}")

set(src_davidson_rci
    Davidson_RCI/david_rci.f90)

set(src_cg
    CG/ccgdiagg.f90
    CG/rcgdiagg.f90)
qe_enable_cuda_fortran("${src_cg}")

set(src_dense
    DENSE/rotate_HSpsi_gamma.f90
    DENSE/rotate_HSpsi_k.f90
    DENSE/rotate_wfc_gamma.f90
    DENSE/rotate_wfc_k.f90 
    DENSE/rotate_xpsi_k.f90
    DENSE/rotate_xpsi_gamma.f90
    DENSE/gram_schmidt_k.f90
    DENSE/gram_schmidt_gamma.f90
    DENSE/rotate_driver.f90)
qe_enable_cuda_fortran("${src_dense}")

set(src_paro
    ParO/bpcg_gamma.f90
    ParO/bpcg_k.f90
    ParO/paro_gamma.f90
    ParO/pcg_gamma.f90
    ParO/paro_k.f90
    ParO/pcg_k.f90
    ParO/paro_gamma_new.f90
    ParO/paro_k_new.f90)
qe_enable_cuda_fortran("${src_paro}")

set(src_rmmdiis
    RMM/crmmdiagg.f90 
    RMM/rrmmdiagg.f90) 
    qe_enable_cuda_fortran("${src_rmmdiis}")

qe_add_library(qe_kssolver_davidson ${src_davidson})
target_link_libraries(qe_kssolver_davidson
    PRIVATE
        qe_device_lapack
        qe_openacc_fortran
        qe_openmp_fortran
        qe_lax
        qe_utilx
        qe_mpi_fortran
        qe_devxlib)

qe_add_library(qe_kssolver_davidsonrci ${src_davidson_rci})
target_link_libraries(qe_kssolver_davidsonrci
    PRIVATE
        qe_lax
        qe_utilx
        qe_mpi_fortran)

qe_add_library(qe_kssolver_cg ${src_cg})
target_link_libraries(qe_kssolver_cg
    PRIVATE
        qe_lax
        qe_openacc_fortran
        qe_utilx
        qe_mpi_fortran
        qe_devxlib)

qe_add_library(qe_kssolver_dense ${src_dense})
target_link_libraries(qe_kssolver_dense
    PRIVATE
        qe_openacc_fortran
        qe_lax
        qe_utilx
        qe_mpi_fortran
        qe_devxlib)
target_include_directories(qe_kssolver_dense
  PUBLIC
    $<BUILD_INTERFACE:${ks_include_dir}>
    $<INSTALL_INTERFACE:include/qe>
)
set_target_properties(qe_kssolver_dense
  PROPERTIES PUBLIC_HEADER ${ks_out_headers}
)
add_dependencies(qe_kssolver_dense qe_ks_headers)

qe_add_library(qe_kssolver_paro ${src_paro})
target_link_libraries(qe_kssolver_paro
    PRIVATE 
        qe_device_lapack
        qe_openacc_fortran
        qe_kssolver_dense
        qe_lax
        qe_utilx
        qe_mpi_fortran)

qe_add_library(qe_kssolver_rmmdiis ${src_rmmdiis})
target_link_libraries(qe_kssolver_rmmdiis
    PRIVATE
	qe_lax
        qe_openacc_fortran
        qe_utilx
	qe_devxlib
        qe_mpi_fortran )  

###########################################################

qe_install_targets(
    qe_kssolver_davidson
    qe_kssolver_davidsonrci
    qe_kssolver_cg
    qe_kssolver_paro
    qe_kssolver_dense
    qe_kssolver_rmmdiis)
