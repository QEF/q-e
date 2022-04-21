###########################################################
# kcw.x
###########################################################
set(src_kcw_x kcw.f90)
qe_enable_cuda_fortran("${src_kcw_x}")
qe_add_executable(qe_kcw_exe ${src_kcw_x})
set_target_properties(qe_kcw_exe 
    PROPERTIES 
        OUTPUT_NAME kcw.x)
target_link_libraries(qe_kcw_exe
    PRIVATE
        qe_pw
        qe_lax
        qe_modules
        qe_kcw
        qe_xclib)

