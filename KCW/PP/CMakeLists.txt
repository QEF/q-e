###########################################################
# kcwpp_sh.x
###########################################################
set(src_kcwpp_sh_x compute_self_hartree.f90)
qe_enable_cuda_fortran("${src_kcwpp_sh_x}")
qe_add_executable(qe_kcwpp_sh_exe ${src_kcwpp_sh_x})
set_target_properties(qe_kcwpp_sh_exe 
    PROPERTIES 
        OUTPUT_NAME kcwpp_sh.x)
target_link_libraries(qe_kcwpp_sh_exe
    PRIVATE
        qe_pw
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
set_target_properties(qe_kcwpp_interp_exe 
    PROPERTIES 
        OUTPUT_NAME kcwpp_interp.x)
target_link_libraries(qe_kcwpp_interp_exe
    PRIVATE
        qe_pw
        qe_lax
        qe_lr_modules
        qe_modules
        qe_kcw
        qe_xclib
        qe_kcwpp)
