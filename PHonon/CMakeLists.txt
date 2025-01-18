set(src_ph
    PH/acfdtest.f90
    PH/add_dkmds.f90
    PH/add_for_charges.f90
    PH/add_zstar_ue.f90
    PH/add_zstar_ue_us.f90
    PH/addcore.f90
    PH/dynmat_nlcc.f90
    PH/addnlcc_zstar_eu_us.f90
    PH/addusddens.f90
    PH/addusddens_cond.f90
    PH/addusddense.f90
    PH/addusdynmat.f90
    PH/allocate_part.f90
    PH/allocate_pert.f90
    PH/allocate_phq.f90
    PH/yambo.f90
    PH/bcast_ph_input.f90
    PH/cch_psi_all.f90
    PH/check_if_partial_dyn.f90
    PH/check_initial_status.f90
    PH/check_restart_recover.f90
    PH/clean_pw_ph.f90
    PH/clinear.f90
    PH/close_phq.f90
    PH/compute_alphasum.f90
    PH/compute_alphasum_cond.f90
    PH/compute_becalp.f90
    PH/compute_becsum_ph.f90
    PH/compute_becsum_ph_cond.f90
    PH/compute_drhous.f90
    PH/compute_drhous_nc.f90
    PH/compute_dvloc.f90
    PH/compute_nldyn.f90
    PH/compute_weight.f90
    PH/deallocate_part.f90
    PH/deallocate_phq.f90
    PH/d2ionq.f90
    PH/d2ionq_disp.f90
    PH/d2nsq_bare.f90
    PH/dnsq_bare.f90
    PH/dnsq_orth.f90
    PH/dnsq_orth_set_irr.f90
    PH/dnsq_scf.f90
    PH/davcio_drho.f90
    PH/dfile_autoname.f90
    PH/dfile_star.f90
    PH/dfpt_kernels.f90
    PH/dielec.f90
    PH/do_phonon.f90
    PH/do_q2r.f90
    PH/drho.f90
    PH/drhodv.f90
    PH/drhodvloc.f90
    PH/drhodvnl.f90
    PH/drhodvus.f90
    PH/dvanqq.f90
    PH/dvkb3.f90
    PH/dvpsi_e.f90
    PH/dvqpsi_us.f90
    PH/dvqpsi_us_only.f90
    PH/dyndia.f90
    PH/dynmat0.f90
    PH/dynmat_us.f90
    PH/dynmatcc.f90
    PH/dynmatrix.f90
    PH/dynmat_hub_bare.f90
    PH/dynmat_hub_scf.f90
    PH/dwfc.f90
    PH/dvqhub_barepsi_us.f90
    PH/dvqhub_barepsi_us2.f90
    PH/delta_sphi.f90
    PH/dpsi_orth.f90
    PH/dump_drho.f90
    PH/ep_matrix_element_wannier.f90
    PH/elph.f90
    PH/el_ph_collect.f90
    PH/elphon.f90
    PH/elph_scdft_mod.f90
    PH/elph_tetra_mod.f90
    PH/find_equiv_sites.f90
    PH/find_mode_sym.f90
    PH/find_irrep.f90
    PH/generate_dynamical_matrix_c.f90
    PH/generate_effective_charges_c.f90
    PH/gmressolve_all.f90
    PH/incdrhous.f90
    PH/incdrhous_nc.f90
    PH/initialize_ph.f90
    PH/init_representations.f90
    PH/io_dyn_mat.f90
    PH/io_dyn_mat_old.f90
    PH/io_pattern.f90
    PH/obsolete.f90
    PH/openfilq.f90
    PH/phcom.f90
    PH/ph_restart.f90
    PH/ph_set_upert.f90
    PH/phescf.f90
    PH/phq_init.f90
    PH/phq_readin.f90
    PH/phq_recover.f90
    PH/phq_setup.f90
    PH/phq_summary.f90
    PH/phqscf.f90
    PH/polariz.f90
    PH/print_clock_ph.f90
    PH/prepare_q.f90
    PH/prepare_sym_analysis.f90
    PH/psidspsi.f90
    PH/punch_plot_e.f90
    PH/q_points.f90
    PH/q_points_wannier.f90
    PH/q2qstar_ph.f90
    PH/random_matrix.f90
    PH/read_wfc_rspace_and_fwfft.f90
    PH/rigid.f90
    PH/rotate_dvscf_star.f90
    PH/rotate_and_add_dyn.f90
    PH/rotate_pattern_add.f90
    PH/run_nscf.f90
    PH/save_ph_input.f90
    PH/set_asr_c.f90
    PH/set_drhoc.f90
    PH/set_int12_nc.f90
    PH/set_irr.f90
    PH/set_irr_nosym.f90
    PH/set_irr_sym.f90
    PH/setlocq.f90
    PH/solve_e.f90
    PH/solve_e_fpol.f90
    PH/solve_linter.f90
    PH/stop_ph.f90
    PH/summarize.f90
    PH/sym_and_write_zue.f90
    PH/symdyn_munu.f90
    PH/symdynph_gq.f90
    PH/symm.f90
    PH/symmorphic_or_nzb.f90
    PH/swfc.f90
    PH/sym_dns.f90
    PH/sym_dns_wrapper.f90
    PH/transform_int_so.f90
    PH/transform_int_nc.f90
    PH/transform_alphasum_nc.f90
    PH/transform_alphasum_so.f90
    PH/tra_write_matrix.f90
    PH/trntnsc.f90
    PH/twochem_postproc_dfpt.f90
    PH/write_epsilon_and_zeu.f90
    PH/write_qplot_data.f90
    PH/write_matrix.f90
    PH/write_modes.f90
    PH/write_rec.f90
    PH/write_hub.f90
    PH/zstar_eu.f90
    PH/zstar_eu_us.f90
    # Raman
    PH/psyme2.f90
    PH/chi_test.f90
    PH/dielec_test.f90
    PH/dvpsi_e2.f90
    PH/el_opt.f90
    PH/hdiag.f90
    PH/pcgreen.f90
    PH/ramanm.f90
    PH/syme2.f90
    PH/solve_e2.f90
    PH/solve_e_nscf.f90
    PH/write_ramtns.f90
    PH/d2mxc.f90
    PH/raman_mat.f90
    PH/raman.f90
    PH/dhdrhopsi.f90
    PH/set_dvscf.f90
    PH/dvscf_interpolate.f90
    PH/ahc.f90)
qe_enable_cuda_fortran("${src_ph}")

set(src_phaux
    PH/write_eigenvectors.f90)

set(src_gamma
    Gamma/a_h.f90
    Gamma/cg_readin.f90
    Gamma/cg_setup.f90
    Gamma/cgcom.f90
    Gamma/cgsolve.f90
    Gamma/d2ion.f90
    Gamma/dgradcorr.f90
    Gamma/dielec.f90
    Gamma/drhodv.f90
    Gamma/dvb_cc.f90
    Gamma/dvpsi_e.f90
    Gamma/dvpsi_kb.f90
    Gamma/dyndiar.f90
    Gamma/dynmat_init.f90
    Gamma/dynmatcc.f90
    Gamma/find_equiv_sites.f90
    Gamma/generate_dynamical_matrix.f90
    Gamma/generate_effective_charges.f90
    Gamma/h_h.f90
    Gamma/macro.f90
    Gamma/rhod2vkb.f90
    Gamma/solve_e.f90
    Gamma/solve_ph.f90
    Gamma/writedyn.f90)
qe_enable_cuda_fortran("${src_gamma}")

set(src_fd
    FD/stop_pp.f90)

qe_add_library(qe_phonon_ph ${src_ph})
target_link_libraries(qe_phonon_ph
    PRIVATE
        qe_hdf5_fortran
        qe_openmp_fortran
	qe_openacc_fortran
        qe_mpi_fortran
        qe_pw
        qe_modules
        qe_lr_modules
        qe_upflib
        qe_fftx
        qe_dftd3
        qe_fox
        qe_xclib)

qe_add_library(qe_phonon_phaux ${src_phaux})
target_link_libraries(qe_phonon_phaux
    PRIVATE
        qe_modules
        qe_lr_modules)

qe_add_library(qe_phonon_gamma ${src_gamma})
target_link_libraries(qe_phonon_gamma
    PRIVATE
        qe_pw
        qe_modules
        qe_lr_modules
        qe_upflib
        qe_fftx
        qe_lax
        qe_xclib)

qe_add_library(qe_phonon_fd ${src_fd})
target_link_libraries(qe_phonon_fd
    PRIVATE
        qe_pw
        qe_modules
        qe_lr_modules)

if(QE_ENABLE_CUDA)
    target_link_libraries(qe_modules
        PRIVATE
            CUDA::cublas)
endif()

###########################################################
# ph.x
###########################################################
set(src_phonon_x PH/phonon.f90)
qe_enable_cuda_fortran("${src_phonon_x}")
qe_add_executable(qe_phonon_ph_exe ${src_phonon_x})
set_target_properties(qe_phonon_ph_exe PROPERTIES OUTPUT_NAME ph.x)
target_link_libraries(qe_phonon_ph_exe
    PRIVATE
        qe_openacc_fortran
	qe_openmp_fortran
	qe_mpi_fortran
        qe_pw
        qe_lax
        qe_modules
        qe_phonon_ph)

###########################################################
# dynmat.x
###########################################################
set(src_dynmat_x PH/dynmat.f90)
qe_add_executable(qe_phonon_dynmat_exe ${src_dynmat_x})
set_target_properties(qe_phonon_dynmat_exe PROPERTIES OUTPUT_NAME dynmat.x)
target_link_libraries(qe_phonon_dynmat_exe
    PRIVATE
        qe_lr_modules
        qe_modules
        qe_phonon_ph
        qe_phonon_phaux
        qe_fox)

###########################################################
# q2r.x
###########################################################
set(src_q2r_x PH/q2r.f90)
qe_add_executable(qe_phonon_q2r_exe ${src_q2r_x})
set_target_properties(qe_phonon_q2r_exe PROPERTIES OUTPUT_NAME q2r.x)
target_link_libraries(qe_phonon_q2r_exe
    PRIVATE
        qe_modules
        qe_phonon_ph)

###########################################################
# dvscf_q2r.x
###########################################################
set(src_dvscf_q2r_x PH/dvscf_q2r.f90)
qe_add_executable(qe_phonon_dvscf_q2r_exe ${src_dvscf_q2r_x})
set_target_properties(qe_phonon_dvscf_q2r_exe PROPERTIES OUTPUT_NAME dvscf_q2r.x)
target_link_libraries(qe_phonon_dvscf_q2r_exe
    PRIVATE
        qe_pw
        qe_lr_modules
        qe_modules
        qe_fftx
        qe_lapack
        qe_phonon_ph)

##################################################################
# matdyn.x 
##################################################################
set(src_matdyn_x PH/matdyn.f90)
qe_add_executable(qe_phonon_matdyn_exe ${src_matdyn_x})
set_target_properties(qe_phonon_matdyn_exe PROPERTIES OUTPUT_NAME matdyn.x) 
target_link_libraries(qe_phonon_matdyn_exe
    PRIVATE
        qe_modules
        qe_pw
        qe_phonon_ph
        qe_phonon_phaux
        qe_lr_modules)

###########################################################
# q2qstar.x
###########################################################
set(src_q2qstar_x PH/q2qstar.f90)
qe_add_executable(qe_phonon_q2qstar_exe ${src_q2qstar_x})
set_target_properties(qe_phonon_q2qstar_exe PROPERTIES OUTPUT_NAME q2qstar.x)
target_link_libraries(qe_phonon_q2qstar_exe
    PRIVATE
        qe_pw
        qe_lr_modules
        qe_modules
        qe_phonon_ph)

###########################################################
# lambda.x
###########################################################
set(src_lambda_x PH/lambda.f90)
qe_add_executable(qe_phonon_lambda_exe ${src_lambda_x})
set_target_properties(qe_phonon_lambda_exe PROPERTIES OUTPUT_NAME lambda.x)
target_link_libraries(qe_phonon_lambda_exe
    PRIVATE
        qe_modules
        qe_phonon_ph)

###########################################################
# alpha2f.x
###########################################################
set(src_alpha2f_x PH/alpha2f.f90)
qe_add_executable(qe_phonon_alpha2f_exe ${src_alpha2f_x})
set_target_properties(qe_phonon_alpha2f_exe PROPERTIES OUTPUT_NAME alpha2f.x)
target_link_libraries(qe_phonon_alpha2f_exe
    PRIVATE
        qe_openacc_fortran
        qe_pw
        qe_modules
        qe_phonon_ph)

###########################################################
# epa.x
###########################################################
set(src_epa_x PH/epa.f90)
qe_add_executable(qe_phonon_epa_exe ${src_epa_x})
set_target_properties(qe_phonon_epa_exe PROPERTIES OUTPUT_NAME epa.x)
target_link_libraries(qe_phonon_epa_exe
    PRIVATE
        qe_modules
        qe_phonon_ph)

###########################################################
# fqha.x
###########################################################
set(src_fqha_x PH/fqha.f90)
qe_add_executable(qe_phonon_fqha_exe ${src_fqha_x})
set_target_properties(qe_phonon_fqha_exe PROPERTIES OUTPUT_NAME fqha.x)
target_link_libraries(qe_phonon_fqha_exe
    PRIVATE
        qe_utilx
        qe_phonon_ph)

###########################################################
# postahc.x
###########################################################
set(src_postahc_x PH/postahc.f90)
qe_add_executable(qe_phonon_postahc_exe ${src_postahc_x})
set_target_properties(qe_phonon_postahc_exe PROPERTIES OUTPUT_NAME postahc.x)
target_link_libraries(qe_phonon_postahc_exe
    PRIVATE
        qe_modules
        qe_pw
        qe_lr_modules
        qe_phonon_ph)

###########################################################
# phcg.x
###########################################################
set(src_phonon_phcg_x Gamma/phcg.f90)
qe_enable_cuda_fortran("${src_phonon_phcg_x}")
qe_add_executable(qe_phonon_phcg_exe ${src_phonon_phcg_x})
set_target_properties(qe_phonon_phcg_exe PROPERTIES OUTPUT_NAME phcg.x)
target_link_libraries(qe_phonon_phcg_exe
    PRIVATE
        qe_pw
        qe_fftx
        qe_xclib
        qe_lr_modules
        qe_phonon_gamma
        qe_modules)

###########################################################
# fd.x
###########################################################
set(src_fd_x FD/fd.f90)
qe_add_executable(qe_phonon_fd_exe ${src_fd_x})
set_target_properties(qe_phonon_fd_exe PROPERTIES OUTPUT_NAME fd.x)
target_link_libraries(qe_phonon_fd_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_phonon_fd
        qe_phonon_ph
        qe_fftx)

###########################################################
# fd_ef.x
###########################################################
set(src_fd_ef_x FD/fd_ef.f90)
qe_add_executable(qe_phonon_fdef_exe ${src_fd_ef_x})
set_target_properties(qe_phonon_fdef_exe PROPERTIES OUTPUT_NAME fd_ef.x)
target_link_libraries(qe_phonon_fdef_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_phonon_fd
        qe_phonon_ph
        qe_fftx)

###########################################################
# fd_ifc.x
###########################################################
set(src_fd_ifc_x FD/fd_ifc.f90)
qe_add_executable(qe_phonon_fdifc_exe ${src_fd_ifc_x})
set_target_properties(qe_phonon_fdifc_exe PROPERTIES OUTPUT_NAME fd_ifc.x)
target_link_libraries(qe_phonon_fdifc_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_phonon_fd
        qe_phonon_ph
        qe_fftx)

###########################################################

set(PH_EXE_TARGETS
    qe_phonon_ph_exe
    qe_phonon_phcg_exe
    qe_phonon_dynmat_exe
    qe_phonon_q2r_exe
    qe_phonon_dvscf_q2r_exe
    qe_phonon_matdyn_exe
    qe_phonon_q2qstar_exe
    qe_phonon_lambda_exe
    qe_phonon_alpha2f_exe
    qe_phonon_epa_exe
    qe_phonon_fqha_exe
    qe_phonon_fd_exe
    qe_phonon_fdef_exe
    qe_phonon_fdifc_exe
    qe_phonon_postahc_exe)

qe_install_targets(
    # Libraries
    qe_phonon_ph
    qe_phonon_phaux
    qe_phonon_gamma
    qe_phonon_fd
    # Executables
    ${PH_EXE_TARGETS})

add_custom_target(ph
    DEPENDS
        ${PH_EXE_TARGETS}
    COMMENT
        "phonon code, Gamma-only and third-order derivatives")
