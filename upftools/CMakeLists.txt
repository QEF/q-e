###########################################################
# extract_core.x
###########################################################
set(sources extract_core.f90)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_extractcore_exe ${sources})
set_target_properties(qe_upftools_extractcore_exe PROPERTIES OUTPUT_NAME extract_core.x)
target_link_libraries(qe_upftools_extractcore_exe
    PRIVATE
        QE::Modules)

###########################################################
# casino2upf.x
###########################################################
set(sources 
    casino2upf.f90 
    casino_pp.f90
)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_casino2upf_exe ${sources})
set_target_properties(qe_upftools_casino2upf_exe PROPERTIES OUTPUT_NAME casino2upf.x)
target_link_libraries(qe_upftools_casino2upf_exe
    PRIVATE
        QE::Modules)

###########################################################
# upf2upf2.x
###########################################################
set(sources upf2upf2.f90)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_upf2upf2_exe ${sources})
set_target_properties(qe_upftools_upf2upf2_exe PROPERTIES OUTPUT_NAME upf2upf2.x)
target_link_libraries(qe_upftools_upf2upf2_exe
    PRIVATE
        QE::Modules)

###########################################################
# cpmd2upf.x
###########################################################
set(sources cpmd2upf.f90)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_cpmd2upf_exe ${sources})
set_target_properties(qe_upftools_cpmd2upf_exe PROPERTIES OUTPUT_NAME cpmd2upf.x)
target_link_libraries(qe_upftools_cpmd2upf_exe
    PRIVATE
        QE::Modules)

###########################################################
# fhi2upf.x
###########################################################
set(sources fhi2upf.f90)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_fhi2upf_exe ${sources})
set_target_properties(qe_upftools_fhi2upf_exe PROPERTIES OUTPUT_NAME fhi2upf.x)
target_link_libraries(qe_upftools_fhi2upf_exe
    PRIVATE
        QE::Modules)

###########################################################
# fpmd2upf.x
###########################################################
set(sources 
    fpmd2upf.f90
    vanderbilt.f90
    write_upf.f90
    read_upf.f90
    casino_pp.f90
)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_fpmd2upf_exe ${sources})
set_target_properties(qe_upftools_fpmd2upf_exe PROPERTIES OUTPUT_NAME fpmd2upf.x)
target_link_libraries(qe_upftools_fpmd2upf_exe
    PRIVATE
        QE::Modules)

###########################################################
# ncpp2upf.x
###########################################################
set(sources 
    ncpp2upf.f90
    write_upf.f90
)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_ncpp2upf_exe ${sources})
set_target_properties(qe_upftools_ncpp2upf_exe PROPERTIES OUTPUT_NAME ncpp2upf.x)
target_link_libraries(qe_upftools_ncpp2upf_exe
    PRIVATE
        QE::Modules)

###########################################################
# oldcp2upf.x
###########################################################
set(sources 
    oldcp2upf.f90
    write_upf.f90
)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_oldcp2upf_exe ${sources})
set_target_properties(qe_upftools_oldcp2upf_exe PROPERTIES OUTPUT_NAME oldcp2upf.x)
target_link_libraries(qe_upftools_oldcp2upf_exe
    PRIVATE
        QE::Modules)

###########################################################
# read_upf.x
###########################################################
set(sources 
    read_ps.f90
    read_upf.f90
)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_readupf_exe ${sources})
set_target_properties(qe_upftools_readupf_exe PROPERTIES OUTPUT_NAME read_upf.x)
target_link_libraries(qe_upftools_readupf_exe
    PRIVATE
        QE::Modules)

###########################################################
# read_upf_tofile.x
###########################################################
set(sources read_upf_tofile.f90)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_readupftofile_exe ${sources})
set_target_properties(qe_upftools_readupftofile_exe PROPERTIES OUTPUT_NAME read_upf_tofile.x)
target_link_libraries(qe_upftools_readupftofile_exe
    PRIVATE
        QE::Modules)

###########################################################
# fix_upf.x
###########################################################
set(sources fix_upf.f90)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_fixupf_exe ${sources})
set_target_properties(qe_upftools_fixupf_exe PROPERTIES OUTPUT_NAME fix_upf.x)
target_link_libraries(qe_upftools_fixupf_exe
    PRIVATE
        QE::Modules)

###########################################################
# rrkj2upf.x
###########################################################
set(sources 
    rrkj2upf.f90
    write_upf.f90
)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_rrkj2upf_exe ${sources})
set_target_properties(qe_upftools_rrkj2upf_exe PROPERTIES OUTPUT_NAME rrkj2upf.x)
target_link_libraries(qe_upftools_rrkj2upf_exe
    PRIVATE
        QE::Modules)

###########################################################
# upf2casino.x
###########################################################
set(sources 
    upf2casino.f90
    vanderbilt.f90
    write_upf.f90
    casino_pp.f90
)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_upf2casino_exe ${sources})
set_target_properties(qe_upftools_upf2casino_exe PROPERTIES OUTPUT_NAME upf2casino.x)
target_link_libraries(qe_upftools_upf2casino_exe
    PRIVATE
        QE::Modules)

###########################################################
# uspp2upf.x
###########################################################
set(sources 
    uspp2upf.f90
    vanderbilt.f90
    write_upf.f90
)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_uspp2upf_exe ${sources})
set_target_properties(qe_upftools_uspp2upf_exe PROPERTIES OUTPUT_NAME uspp2upf.x)
target_link_libraries(qe_upftools_uspp2upf_exe
    PRIVATE
        QE::Modules)

###########################################################
# vdb2upf.x
###########################################################
set(sources 
    vdb2upf.f90
    vanderbilt.f90
    write_upf.f90
)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_vdb2upf_exe ${sources})
set_target_properties(qe_upftools_vdb2upf_exe PROPERTIES OUTPUT_NAME vdb2upf.x)
target_link_libraries(qe_upftools_vdb2upf_exe
    PRIVATE
        QE::Modules)

###########################################################
# virtual_v2.x
###########################################################
set(sources virtual_v2.f90)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_virtualv2_exe ${sources})
set_target_properties(qe_upftools_virtualv2_exe PROPERTIES OUTPUT_NAME virtual_v2.x)
target_link_libraries(qe_upftools_virtualv2_exe
    PRIVATE
        QE::Modules)

###########################################################
# interpolate.x
###########################################################
set(sources 
    interpolate.f90
    vanderbilt.f90
    write_upf.f90
)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_interpolate_exe ${sources})
set_target_properties(qe_upftools_interpolate_exe PROPERTIES OUTPUT_NAME interpolate.x)
target_link_libraries(qe_upftools_interpolate_exe
    PRIVATE
        QE::Modules)

###########################################################
# hgh2qe.x
###########################################################
set(sources 
    HGH2QE/hgh2qe.f90
    HGH2QE/lininterpol.f90
)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_upftools_hgh2qelininterpol_gh2qe_exe ${sources})
set_target_properties(qe_upftools_hgh2qelininterpol_gh2qe_exe PROPERTIES OUTPUT_NAME hgh2qe_lininterpol.x)
target_link_libraries(qe_upftools_hgh2qelininterpol_gh2qe_exe
    PRIVATE
        QE::Modules)

###########################################################

qe_install_targets(
    qe_upftools_extractcore_exe
    qe_upftools_casino2upf_exe
    qe_upftools_upf2upf2_exe
    qe_upftools_cpmd2upf_exe
    qe_upftools_fhi2upf_exe
    qe_upftools_fpmd2upf_exe
    qe_upftools_ncpp2upf_exe
    qe_upftools_oldcp2upf_exe
    qe_upftools_readupf_exe
    qe_upftools_readupftofile_exe
    qe_upftools_fixupf_exe
    qe_upftools_rrkj2upf_exe
    qe_upftools_upf2casino_exe
    qe_upftools_uspp2upf_exe
    qe_upftools_vdb2upf_exe
    qe_upftools_virtualv2_exe
    qe_upftools_interpolate_exe
    qe_upftools_hgh2qelininterpol_gh2qe_exe
)