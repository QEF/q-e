set(src_upflib
    atom.f90
    atomic_number.f90
    atwfc_mod.f90
    dqvan2.f90
    gen_us_dj.f90
    gen_us_dy.f90
    init_us_0.f90
    init_us_b0.f90
    init_us_1.f90
    beta_mod.f90
    init_us_2_acc.f90
    gth.f90
    paw_variables.f90
    pseudo_types.f90
    qvan2.f90
    radial_grids.f90
    read_cpmd.f90
    read_fhi.f90
    read_ncpp.f90
    read_ps.f90
    read_psml.f90
    read_upf_new.f90
    read_upf_v1.f90
    read_uspp.f90
    spinor.f90
    sph_ind.f90
    sph_bes.f90
    splinelib.f90
    simpsn.f90
    upf_auxtools.f90
    upf_const.f90
    upf_error.f90
    upf_invmat.f90
    upf_io.f90
    upf_ions.f90
    upf_kinds.f90
    upf_spinorb.f90
    upf_params.f90
    upf_parallel_include.f90
    upf_to_internal.f90
    upf_utils.f90
    uspp.f90
    uspp_param.f90
    vloc_mod.f90
    rhoc_mod.f90
    rhoat_mod.f90
    qrad_mod.f90
    write_upf_new.f90
    ylmr2.f90
    dylmr2.f90
    # CUDA Fortran routines
    ylmr2_gpu.f90)
qe_enable_cuda_fortran("${src_upflib}")

qe_add_library(qe_xml xmltools.f90 dom.f90 wxml.f90)
target_compile_definitions(qe_xml PRIVATE __XML_STANDALONE)

qe_add_library(qe_upflib ${src_upflib})
target_link_libraries(qe_upflib
    PUBLIC
        qe_openacc_fortran
        qe_xml
    PRIVATE
        qe_mpi_fortran
        qe_lapack
        qe_utilx)

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

###########################################################
# casino2upf.x
###########################################################
set(src_casino2upf_x
    casino2upf.f90
    casino_pp.f90 )
qe_add_executable(qe_upflib_casino2upf_exe ${src_casino2upf_x})
set_target_properties(qe_upflib_casino2upf_exe
    PROPERTIES
        OUTPUT_NAME casino2upf.x)
target_link_libraries(qe_upflib_casino2upf_exe
    PRIVATE
        qe_upflib)

###########################################################
# hgh2qe.x
###########################################################
qe_add_executable(qe_upflib_hgh2qe_exe hgh2qe.f90)
set_target_properties(qe_upflib_hgh2qe_exe
    PROPERTIES
        OUTPUT_NAME hgh2qe.x)

############################################################################
qe_install_targets(qe_xml
    qe_upflib
    qe_upflib_virtual_v2_exe
    qe_upflib_upfconv_exe
    qe_upflib_casino2upf_exe
    qe_upflib_hgh2qe_exe)

add_custom_target(upf
    DEPENDS
        #Library
        qe_upflib
        #Executables
        qe_upflib_virtual_v2_exe
        qe_upflib_upfconv_exe
    COMMENT
        "utilities for pseudopotential conversion")
