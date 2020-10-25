set(sources_ph
    PH/acfdtest.f90
    PH/add_dkmds.f90
    PH/add_for_charges.f90
    PH/add_zstar_ue.f90
    PH/add_zstar_ue_us.f90
    PH/addcore.f90
    PH/addnlcc.f90
    PH/addnlcc_zstar_eu_us.f90
    PH/addusddens.f90
    PH/addusddense.f90
    PH/addusdynmat.f90
    PH/adddvhubscf.f90
    PH/allocate_part.f90
    PH/allocate_pert.f90
    PH/allocate_phq.f90
    PH/yambo.f90
    PH/bcast_ph_input.f90
    PH/cch_psi_all.f90
    PH/ccg_psi.f90
    PH/check_if_partial_dyn.f90
    PH/check_initial_status.f90
    PH/check_restart_recover.f90
    PH/clean_pw_ph.f90
    PH/clinear.f90
    PH/close_phq.f90
    PH/compute_alphasum.f90
    PH/compute_becalp.f90
    PH/compute_becsum_ph.f90
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
    PH/dnsq_scf.f90
    PH/davcio_drho.f90
    PH/dfile_autoname.f90
    PH/dfile_star.f90
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
    PH/ef_shift.f90
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
    PH/psymdvscf.f90
    PH/psyme.f90
    PH/psym_dmag.f90
    PH/psym_dmage.f90
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
    PH/sym_def.f90
    PH/sym_dmag.f90
    PH/sym_dmage.f90
    PH/symdvscf.f90
    PH/symdyn_munu.f90
    PH/symdynph_gq.f90
    PH/syme.f90
    PH/symm.f90
    PH/symmorphic_or_nzb.f90
    PH/swfc.f90
    PH/syme_dns.f90
    PH/sym_dns.f90
    PH/sym_dns_wrapper.f90
    PH/transform_int_so.f90
    PH/transform_int_nc.f90
    PH/transform_alphasum_nc.f90
    PH/transform_alphasum_so.f90
    PH/tra_write_matrix.f90
    PH/trntnsc.f90
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
    PH/non_scf_ph.f90
    PH/adddvscf_ph_mag.f90
    PH/apply_trev.f90
    PH/c_bands_ph.f90
    PH/dvscf_interpolate.f90
    PH/ahc.f90
    PH/postahc.f90)

set(sources_phaux
    PH/write_eigenvectors.f90)

set(sources_gamma
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

set(sources_fd
    FD/stop_pp.f90
)
    

qe_add_library(qe_phonon_ph ${sources_ph})
add_library(QE::PHonon::PH ALIAS qe_phonon_ph)
target_link_libraries(qe_phonon_ph
    PRIVATE
        QE::HDF5_Fortran
        QE::OpenMP_Fortran
        QE::PW
        QE::Modules
        QE::LR_Modules)

qe_add_library(qe_phonon_phaux ${sources_phaux})
add_library(QE::PHonon::PHAUX ALIAS qe_phonon_phaux)
target_link_libraries(qe_phonon_phaux
    PRIVATE
        QE::Modules
        QE::LR_Modules)

qe_add_library(qe_phonon_gamma ${sources_gamma})
add_library(QE::PHonon::Gamma ALIAS qe_phonon_gamma)
target_link_libraries(qe_phonon_gamma
    PRIVATE
        QE::PW
        QE::Modules
        QE::LR_Modules)

qe_add_library(qe_phonon_fd ${sources_fd})
add_library(QE::PHonon::FD ALIAS qe_phonon_fd)
target_link_libraries(qe_phonon_fd
    PRIVATE
        QE::PW
        QE::UTILX
        QE::Modules
        QE::LR_Modules)
        

###########################################################
# ph.x
###########################################################
set(sources PH/phonon.f90)
qe_add_executable(qe_phonon_ph_exe ${sources})
set_target_properties(qe_phonon_ph_exe PROPERTIES OUTPUT_NAME ph.x)
target_link_libraries(qe_phonon_ph_exe
    PRIVATE
        QE::PW
        QE::LAX
        QE::Modules
        QE::PHonon::PH
)

###########################################################
# dynmat.x
###########################################################
set(sources PH/dynmat.f90)
qe_add_executable(qe_phonon_dynmat_exe ${sources})
set_target_properties(qe_phonon_dynmat_exe PROPERTIES OUTPUT_NAME dynmat.x)
target_link_libraries(qe_phonon_dynmat_exe
    PRIVATE
        QE::LR_Modules
        QE::Modules
        QE::UTILX
        QE::PHonon::PH
        QE::PHonon::PHAUX
)

###########################################################
# q2r.x
###########################################################
set(sources PH/q2r.f90)
qe_add_executable(qe_phonon_q2r_exe ${sources})
set_target_properties(qe_phonon_q2r_exe PROPERTIES OUTPUT_NAME q2r.x)
target_link_libraries(qe_phonon_q2r_exe
    PRIVATE
        QE::Modules
        QE::UTILX
        QE::PHonon::PH
)

##################################################################
# matdyn.x 
##################################################################
set(sources PH/matdyn.f90)
qe_add_executable(qe_phonon_matdyn_exe ${sources})
set_target_properties(qe_phonon_matdyn_exe PROPERTIES OUTPUT_NAME matdyn.x) 
target_link_libraries(qe_phonon_matdyn_exe
    PRIVATE
        QE::Modules
        QE::UTILX
        QE::PW
        QE::PHonon::PH
        QE::PHonon::PHAUX
) 

###########################################################
# q2qstar.x
###########################################################
set(sources PH/q2qstar.f90)
qe_add_executable(qe_phonon_q2qstar_exe ${sources})
set_target_properties(qe_phonon_q2qstar_exe PROPERTIES OUTPUT_NAME q2qstar.x)
target_link_libraries(qe_phonon_q2qstar_exe
    PRIVATE
        QE::PW
        QE::LR_Modules
        QE::Modules
        QE::UTILX
        QE::PHonon::PH
)

###########################################################
# lambda.x
###########################################################
set(sources PH/lambda.f90)
qe_add_executable(qe_phonon_lambda_exe ${sources})
set_target_properties(qe_phonon_lambda_exe PROPERTIES OUTPUT_NAME lambda.x)
target_link_libraries(qe_phonon_lambda_exe
    PRIVATE
        QE::Modules
        QE::UTILX
        QE::PHonon::PH
)

###########################################################
# alpha2f.x
###########################################################
set(sources PH/alpha2f.f90)
qe_add_executable(qe_phonon_alpha2f_exe ${sources})
set_target_properties(qe_phonon_alpha2f_exe PROPERTIES OUTPUT_NAME alpha2f.x)
target_link_libraries(qe_phonon_alpha2f_exe
    PRIVATE
        QE::PW
        QE::Modules
        QE::UTILX
        QE::PHonon::PH
)

###########################################################
# epa.x
###########################################################
set(sources PH/epa.f90)
qe_add_executable(qe_phonon_epa_exe ${sources})
set_target_properties(qe_phonon_epa_exe PROPERTIES OUTPUT_NAME epa.x)
target_link_libraries(qe_phonon_epa_exe
    PRIVATE
        QE::Modules
        QE::UTILX
        QE::PHonon::PH
)

###########################################################
# fqha.x
###########################################################
set(sources PH/fqha.f90)
qe_add_executable(qe_phonon_fqha_exe ${sources})
set_target_properties(qe_phonon_fqha_exe PROPERTIES OUTPUT_NAME fqha.x)
target_link_libraries(qe_phonon_fqha_exe
    PRIVATE
        QE::UTILX
        QE::PHonon::PH
)

###########################################################
# fd.x
###########################################################
set(sources FD/fd.f90)
qe_add_executable(qe_phonon_fd_exe ${sources})
set_target_properties(qe_phonon_fd_exe PROPERTIES OUTPUT_NAME fd.x)
target_link_libraries(qe_phonon_fd_exe
    PRIVATE
        QE::PW
        QE::UTILX
        QE::Modules
        QE::PHonon::FD
        QE::PHonon::PH
)

###########################################################
# fd_ef.x
###########################################################
set(sources FD/fd_ef.f90)
qe_add_executable(qe_phonon_fdef_exe ${sources})
set_target_properties(qe_phonon_fdef_exe PROPERTIES OUTPUT_NAME fd_ef.x)
target_link_libraries(qe_phonon_fdef_exe
    PRIVATE
        QE::PW
        QE::UTILX
        QE::Modules
        QE::PHonon::FD
        QE::PHonon::PH
)

###########################################################
# fd_ifc.x
###########################################################
set(sources FD/fd_ifc.f90)
qe_add_executable(qe_phonon_fdifc_exe ${sources})
set_target_properties(qe_phonon_fdifc_exe PROPERTIES OUTPUT_NAME fd_ifc.x)
target_link_libraries(qe_phonon_fdifc_exe
    PRIVATE
        QE::PW
        QE::UTILX
        QE::Modules
        QE::PHonon::FD
        QE::PHonon::PH
)

###########################################################

qe_install_targets(
    # Libraries
    qe_phonon_ph
    qe_phonon_phaux
    qe_phonon_gamma
    qe_phonon_fd
    # Executables
    qe_phonon_ph_exe
    qe_phonon_dynmat_exe
    qe_phonon_q2r_exe
    qe_phonon_q2qstar_exe
    qe_phonon_lambda_exe
    qe_phonon_alpha2f_exe
    qe_phonon_epa_exe
    qe_phonon_fqha_exe
    qe_phonon_fd_exe
    qe_phonon_fdef_exe
    qe_phonon_fdifc_exe
)
