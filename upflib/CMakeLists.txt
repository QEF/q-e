set(src_upflib
    atom.f90 
    atomic_number.f90 
    upf_erf.f90 
    upf_utils.f90 
    gth.f90 
    pseudo_types.f90 
    radial_grids.f90 
    read_cpmd.f90 
    read_fhi.f90 
    read_ncpp.f90 
    read_upf_new.f90 
    read_upf_v1.f90 
    read_uspp.f90 
    splinelib.f90 
    simpsn.f90 
    upf.f90 
    upf_auxtools.f90 
    upf_const.f90 
    upf_error.f90 
    upf_free_unit.f90 
    upf_invmat.f90 
    upf_io.f90 
    upf_kinds.f90 
    upf_params.f90 
    upf_parallel_include.f90 
    upf_to_internal.f90 
    uspp.f90 
    write_upf_new.f90 
    xmltools.f90
    ylmr2.f90
    # GPU
    uspp_gpu.f90
    ylmr2_gpu.f90
    simpsn_gpu.f90
    splinelib_gpu.f90)
qe_enable_cuda_fortran("${src_upflib}")

qe_add_library(qe_upflib ${src_upflib})
target_link_libraries(qe_upflib
    PRIVATE
        qe_mpi_fortran
        qe_lapack
        qe_utilx
        qe_fox
        qe_devxlib)

##############################################################
# virtual_v2.x 
##############################################################
set(src_virtual_v2_x virtual_v2.f90) 
qe_add_executable(qe_upflib_virtual_v2_exe ${src_virtual_v2_x}) 
set_target_properties(qe_upflib_virtual_v2_exe 
    PROPERTIES 
        OUTPUT_NAME virtual_v2.x)
target_link_libraries(qe_upflib_virtual_v2_exe 
    PRIVATE 
        qe_upflib)

###########################################################
# upfconv.x
###########################################################
set(src_upfconv_x 
    upfconv.f90 
    casino_pp.f90 )
qe_add_executable(qe_upflib_upfconv_exe ${src_upfconv_x}) 
set_target_properties(qe_upflib_upfconv_exe 
    PROPERTIES 
        OUTPUT_NAME upfconv.x)
target_link_libraries(qe_upflib_upfconv_exe
    PRIVATE
        qe_upflib)

############################################################################
qe_install_targets(qe_upflib
    qe_upflib_virtual_v2_exe
    qe_upflib_upfconv_exe)
