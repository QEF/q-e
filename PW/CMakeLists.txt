set(sources
   src/a2fmod.f90 
   src/add_bfield.f90 
   src/add_efield.f90 
   src/add_vuspsi.f90 
   src/add_gatefield.f90 
   src/add_paw_to_deeq.f90 
   src/add_vhub_to_deeq.f90 
   src/addusdens.f90 
   src/addusforce.f90 
   src/addusstress.f90 
   src/allocate_fft.f90 
   src/allocate_locpot.f90 
   src/allocate_nlpot.f90 
   src/allocate_wfc.f90 
   src/atomic_rho.f90 
   src/atomic_wfc.f90 
   src/atomic_wfc_mod.f90 
   src/average_pp.f90 
   src/acfdt_in_pw.f90 
   src/newd.f90 
   src/beef.f90 
   src/bp_mod.f90 
   src/bp_c_phase.f90 
   src/bp_calc_btq.f90 
   src/bp_qvan3.f90 
   src/bp_strings.f90 
   src/buffers.f90 
   src/c_bands.f90 
   src/c_phase_field.f90 
   src/orbm_kubo.f90 
   src/cdiagh.f90 
   src/clean_pw.f90 
   src/close_files.f90 
   src/compute_becsum.f90 
   src/compute_deff.f90 
   src/compute_dip.f90 
   src/compute_rho.f90 
   src/compute_qdipol.f90 
   src/compute_qdipol_so.f90 
   src/compute_ux.f90 
   src/coset.f90 
   src/Coul_cut_2D.f90 
   src/d_matrix.f90 
   src/data_structure.f90 
   src/deriv_drhoc.f90 
   src/divide_class.f90 
   src/divide_class_so.f90 
   src/divide_et_impera.f90 
   src/dqvan2.f90 
   src/drhoc.f90 
   src/rotate_wfc.f90 
   src/run_driver.f90 
   src/dvloc_of_g.f90 
   src/dynamics_module.f90 
   src/efermig.f90 
   src/efermit.f90 
   src/electrons.f90 
   src/eqvect.f90 
   src/esm.f90 
   src/ewald.f90 
   src/ewald_dipole.f90 
   src/extfield.f90 
   src/exx_base.f90 
   src/exx_band.f90 
   src/exx.f90 
   src/fcp.f90 
   src/find_group.f90 
   src/forces_bp_efield.f90 
   src/force_cc.f90 
   src/force_corr.f90 
   src/force_ew.f90 
   src/force_hub.f90 
   src/force_lc.f90 
   src/force_us.f90 
   src/forces.f90 
   src/g_psi.f90 
   src/g_psi_mod.f90 
   src/gen_at_dj.f90 
   src/gen_at_dy.f90 
   src/gen_us_dj.f90 
   src/gen_us_dy.f90 
   src/get_locals.f90 
   src/gk_sort.f90 
   src/gradcorr.f90 
   src/gweights.f90 
   src/g2_kin.f90 
   src/hs_psi.f90 
   src/hs_1psi.f90 
   src/h_epsi_her_apply.f90 
   src/h_epsi_her_set.f90 
   src/h_psi.f90 
   src/h_psi_meta.f90 
   src/hinit0.f90 
   src/hinit1.f90 
   src/init_ns.f90 
   src/init_q_aeps.f90 
   src/init_run.f90 
   src/init_us_0.f90 
   src/init_us_b0.f90 
   src/init_us_1.f90 
   src/init_us_2.f90 
   src/init_at_1.f90 
   src/init_vloc.f90 
   src/input.f90 
   src/io_rho_xml.f90 
   src/irrek.f90 
   src/iweights.f90 
   src/intersite_V.f90 
   src/init_nsg.f90 
   src/nsg_adj.f90 
   src/start_k.f90 
   src/kpoint_grid.f90 
   src/lchk_tauxk.f90 
   src/ldaU.f90 
   src/make_pointlists.f90 
   src/makov_payne.f90 
   src/manypw.f90 
   src/martyna_tuckerman.f90 
   src/memory_report.f90 
   src/mix_rho.f90 
   src/move_ions.f90 
   src/multable.f90 
   src/n_plane_waves.f90 
   src/new_ns.f90 
   src/new_nsb.f90 
   src/new_nsg.f90 
   src/new_occ.f90 
   src/ns_adj.f90 
   src/non_scf.f90 
   src/offset_atom_wfc.f90 
   src/openfil.f90 
   src/orthoatwfc.f90 
   src/output_tau.f90 
   src/para.f90 
   src/paw_exx.f90 
   src/paw_init.f90 
   src/paw_onecenter.f90 
   src/paw_symmetry.f90 
   src/plugin_print_energies.f90 
   src/plugin_scf_energy.f90 
   src/plugin_scf_potential.f90 
   src/plugin_init_ions.f90 
   src/plugin_init_cell.f90 
   src/plugin_init_potential.f90 
   src/plugin_initbase.f90 
   src/plugin_clean.f90 
   src/plugin_check.f90 
   src/plugin_clock.f90 
   src/plugin_summary.f90 
   src/plugin_initialization.f90 
   src/plugin_ext_forces.f90 
   src/plugin_int_forces.f90 
   src/plugin_read_input.f90 
   src/plus_u_full.f90 
   src/potinit.f90 
   src/print_clock_pw.f90 
   src/print_ks_energies.f90 
   src/punch.f90 
   src/pw_restart_new.f90 
   src/add_qexsd_step.f90 
   src/pw_init_qexsd_input.f90 
   src/pwcom.f90 
   src/pw2blip.f90 
   src/pw2casino.f90 
   src/pw2casino_write.f90 
   src/qvan2.f90 
   src/rdiagh.f90 
   src/read_conf_from_file.f90 
   src/read_file_new.f90 
   src/realus.f90 
   src/remove_atomic_rho.f90 
   src/report_mag.f90 
   src/restart_in_electrons.f90 
   src/rho2zeta.f90 
   src/ruotaijk.f90 
   src/run_pwscf.f90 
   src/s_1psi.f90 
   src/s_psi.f90 
   src/save_in_cbands.f90 
   src/save_in_electrons.f90 
   src/scale_h.f90 
   src/loc_scdm.f90 
   src/loc_scdm_k.f90 
   src/scf_mod.f90 
   src/set_kplusq.f90 
   src/set_kup_and_kdw.f90 
   src/set_occupations.f90 
   src/set_rhoc.f90 
   src/set_spin_vars.f90 
   src/set_vrs.f90 
   src/setlocal.f90 
   src/setup.f90 
   src/spinor.f90 
   src/sph_ind.f90 
   src/stop_run.f90 
   src/stres_cc.f90 
   src/stres_ewa.f90 
   src/stres_gradcorr.f90 
   src/stres_har.f90 
   src/stres_hub.f90 
   src/stres_knl.f90 
   src/stres_loc.f90 
   src/stres_us.f90 
   src/stres_nonloc_dft.f90 
   src/stres_mgga.f90 
   src/stress.f90 
   src/struct_fact.f90 
   src/sum_band.f90 
   src/sumkg.f90 
   src/sumkt.f90 
   src/summary.f90 
   src/symme.f90 
   src/symm_base.f90 
   src/symmetrize_at.f90 
   src/tabd.f90 
   src/tetra.f90 
   src/transform_becsum_so.f90 
   src/transform_becsum_nc.f90 
   src/transform_qq_so.f90 
   src/trnvecc.f90 
   src/update_pot.f90 
   src/us_exx.f90 
   src/usnldiag.f90 
   src/v_of_rho.f90 
   src/vcsmd.f90 
   src/vcsubs.f90 
   src/vhpsi.f90 
   src/vloc_of_g.f90 
   src/vloc_psi.f90 
   src/utils.f90 
   src/xdm_dispersion.f90 
   src/wfcinit.f90 
   src/write_ns.f90 
   src/wsweight.f90 
   src/weights.f90 
   src/ortho_wfc.f90 
   src/wannier_init.f90 
   src/wannier_check.f90 
   src/wannier_clean.f90 
   src/wannier_occ.f90
)

qe_add_library(qe_pw ${sources})
add_library(QE::PW ALIAS qe_pw)
target_link_libraries(qe_pw
    PUBLIC
        QE::HDF5_Fortran
        QE::ELPA
        QE::SCALAPACK
        QE::LAX
        QE::UTILX  # TODO direct or transitive via QE::LAX?
        QE::KS_Solvers::Davidson
        QE::KS_Solvers::CG
        QE::KS_Solvers::PPCG
        QE::KS_Solvers::ParO
        QE::UPF # TODO turn to vendorized library from MaX gitrepo ? 
        QE::Modules
        QE::FFTX
        QE::DFTD3
        QE::OpenMP_Fortran
        QE::MPI_Fortran)

###########################################################
# pw.x
###########################################################
set(sources src/pwscf.f90)
qe_add_executable(qe_pw_exe ${sources})
set_target_properties(qe_pw_exe PROPERTIES OUTPUT_NAME pw.x)
target_link_libraries(qe_pw_exe
    PRIVATE
        QE::PW)

###########################################################
# ibrav2cell.x
###########################################################
set(sources tools/ibrav2cell.f90)
qe_add_executable(qe_pw_tools_ibrav2cell_exe ${sources})
set_target_properties(qe_pw_tools_ibrav2cell_exe PROPERTIES OUTPUT_NAME ibrav2cell.x)
target_link_libraries(qe_pw_tools_ibrav2cell_exe
    PRIVATE
        QE::PW)

###########################################################
# cell2ibrav.x
###########################################################
set(sources tools/cell2ibrav.f90)
qe_add_executable(qe_pw_tools_cell2ibrav_exe ${sources})
set_target_properties(qe_pw_tools_cell2ibrav_exe PROPERTIES OUTPUT_NAME cell2ibrav.x)
target_link_libraries(qe_pw_tools_cell2ibrav_exe
    PRIVATE
        QE::PW)

###########################################################
# ev.x
###########################################################
set(sources tools/ev.f90)
qe_add_executable(qe_pw_tools_ev_exe ${sources})
set_target_properties(qe_pw_tools_ev_exe PROPERTIES OUTPUT_NAME ev.x)
target_link_libraries(qe_pw_tools_ev_exe
    PRIVATE
        QE::PW)

###########################################################
# kpoints.x
###########################################################
set(sources tools/kpoints.f90)
qe_add_executable(qe_pw_tools_kpoints_exe ${sources})
set_target_properties(qe_pw_tools_kpoints_exe PROPERTIES OUTPUT_NAME kpoints.x)
target_link_libraries(qe_pw_tools_kpoints_exe
    PRIVATE
        QE::PW)

###########################################################
# pwi2xsf.x
###########################################################
set(sources tools/pwi2xsf.f90)
qe_add_executable(qe_pw_tools_pwi2xsf_exe ${sources})
set_target_properties(qe_pw_tools_pwi2xsf_exe PROPERTIES OUTPUT_NAME pwi2xsf.x)
target_link_libraries(qe_pw_tools_pwi2xsf_exe
    PRIVATE
        QE::PW)

###########################################################

qe_install_targets(
    # Libraries
    qe_pw
    # Executables
    qe_pw_exe
    qe_pw_tools_ibrav2cell_exe
    qe_pw_tools_cell2ibrav_exe
    qe_pw_tools_ev_exe
    qe_pw_tools_kpoints_exe
    qe_pw_tools_pwi2xsf_exe
)
