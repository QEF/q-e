set(src_qeheat
    src/traj_object.f90
    src/cpv_traj.f90
    src/test_h_psi_s_psi_commutator_Hx_psi.f90
    src/hartree_xc_mod.f90
    src/kohn_sham_mod.f90
    src/zero_mod.f90
    src/ionic_mod.f90
    src/project.f90
    src/init_us_3.f90
    src/ec_functionals.f90
    src/compute_charge.f90
    src/averages.f90
    src/scf_result.f90
)

qe_add_library(qe_qeheat ${src_qeheat})
target_link_libraries(qe_qeheat
    PRIVATE
    qe_phonon_ph
    qe_modules
    qe_pw
    qe_lr_modules
    qe_hdf5_fortran
    qe_openmp_fortran
    qe_upflib
    qe_fftx
    qe_fox
    qe_xclib
    qe_device_lapack
    PUBLIC
    qe_openacc_fortran
)

###########################################################
# all_currents.x
###########################################################
set(src_all_currents_x src/all_currents.f90)
qe_enable_cuda_fortran("${src_all_currents_x}")
qe_add_executable(qe_qeheat_exe ${src_all_currents_x})
set_target_properties(qe_qeheat_exe PROPERTIES OUTPUT_NAME all_currents.x)
target_link_libraries(qe_qeheat_exe
    PRIVATE
    qe_pw
    qe_lax
    qe_lr_modules
    qe_modules
    qe_qeheat
    qe_fftx
)

qe_install_targets(
    qe_qeheat_exe
)

add_custom_target(all_currents
    DEPENDS
	qe_qeheat_exe
    COMMENT
	"QEHeat code to compute energy and electronic density currents")
