set(src_kcwpp bcast_kcw_pp_input.f90 input_pp_summary.f90 kcw_pp_readin.f90 read_hr.f90 kcw_bands.f90 sh_setup.f90)

qe_add_library(qe_kcwpp ${src_kcwpp})
target_link_libraries(
    qe_kcwpp
    PRIVATE qe_pw
            qe_lax
            qe_lr_modules
            qe_modules
            qe_fftx
            qe_upflib
            qe_xclib
            qe_kcw)

###########################################################
# kcwpp_sh.x
###########################################################
set(src_kcwpp_sh_x compute_self_hartree.f90)
qe_enable_cuda_fortran("${src_kcwpp_sh_x}")
qe_add_executable(qe_kcwpp_sh_exe ${src_kcwpp_sh_x})
set_target_properties(qe_kcwpp_sh_exe PROPERTIES OUTPUT_NAME kcwpp_sh.x)
target_link_libraries(
    qe_kcwpp_sh_exe
    PRIVATE qe_pw
            qe_lax
            qe_lr_modules
            qe_modules
            qe_kcw
            qe_xclib
            qe_kcwpp)

###########################################################
# kcwpp_interp.x
###########################################################
set(src_kcwpp_interp_x interp_r_to_k.f90)
qe_enable_cuda_fortran("${src_kcwpp_interp_x}")
qe_add_executable(qe_kcwpp_interp_exe ${src_kcwpp_interp_x})
set_target_properties(qe_kcwpp_interp_exe PROPERTIES OUTPUT_NAME kcwpp_interp.x)
target_link_libraries(
    qe_kcwpp_interp_exe
    PRIVATE qe_pw
            qe_lax
            qe_lr_modules
            qe_modules
            qe_kcw
            qe_xclib
            qe_kcwpp)

###########################################################
# ups.x
###########################################################
set(src_ups_x ups.f90)
qe_add_executable(qe_kcwpp_ups_exe ${src_ups_x})
set_target_properties(qe_kcwpp_ups_exe PROPERTIES OUTPUT_NAME ups.x)
target_link_libraries(qe_kcwpp_ups_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_upflib)

###########################################################
# merge_wann.x
###########################################################
set(src_merge_wann_x merge_wann.f90)
qe_add_executable(qe_kcwpp_merge_wann_exe ${src_merge_wann_x})
set_target_properties(qe_kcwpp_merge_wann_exe PROPERTIES OUTPUT_NAME merge_wann.x)
target_link_libraries(qe_kcwpp_merge_wann_exe)
