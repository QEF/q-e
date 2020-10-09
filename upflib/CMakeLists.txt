set(sources
    atom.f90 
    atomic_number.f90 
    emend_upf.f90 
    upf_erf.f90 
    upf_utils.f90 
    gth.f90 
    pseudo_types.f90 
    radial_grids.f90 
    read_cpmd.f90 
    read_fhi.f90 
    read_ncpp.f90 
    read_upf_new.f90 
    read_upf_schema.f90 
    read_upf_v1.f90 
    read_upf_v2.f90 
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
    write_upf.f90 
    write_upf_new.f90 
    write_upf_schema.f90 
    write_upf_v2.f90 
    xmltools.f90
)

qe_add_library(qe_upflib ${sources})
add_library(QE::UPF ALIAS qe_upflib)
target_link_libraries(qe_upflib
    PRIVATE
        QE::MPI_Fortran
        QE::LAPACK
    PUBLIC
        QE::FOX)

##############################################################
# virtual_v2.x 
##############################################################
set(sources virtual_v2.f90) 
qe_add_executable(qe_upflib_virtual_v2_exe ${sources}) 
set_target_properties(qe_upflib_virtual_v2_exe PROPERTIES OUTPUT_NAME virtual_v2.x) 
target_link_libraries(qe_upflib_virtual_v2_exe 
      PRIVATE 
        QE::UPF
)
###########################################################
# upfconv.x
###########################################################
set(sources 
    upfconv.f90 
    casino_pp.f90 
)
qe_add_executable(qe_upflib_upfconv_exe ${sources}) 
set_target_properties(qe_upflib_upfconv_exe PROPERTIES OUTPUT_NAME upfconv.x)
target_link_libraries(qe_upflib_upfconv_exe
                      PRIVATE
                      QE::UPF
)

############################################################################
qe_install_targets(qe_upflib
                   qe_upflib_virtual_v2_exe
                   qe_upflib_upfconv_exe
)
