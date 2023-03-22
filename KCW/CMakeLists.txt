add_subdirectory(src)
add_subdirectory(PP)

###########################################################

qe_install_targets(qe_kcw qe_kcw_exe qe_kcwpp_sh_exe qe_kcwpp_interp_exe)

add_custom_target(
    kcw
    DEPENDS pw qe_pp_pw2wannier90_exe w90 qe_kcw_exe qe_kcwpp_sh_exe qe_kcwpp_interp_exe
    COMMENT "Koopmans spectral functionals")
