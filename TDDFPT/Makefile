# Makefile for PW
#Modified for v4.0 by obm


include ../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../Modules \
          $(MOD_FLAG)../EE $(MOD_FLAG). $(MOD_FLAG)../PW \
          $(MOD_FLAG)../PH

PWOBJS2 = ../PW/*.o

PWOBJS = \
../PW/a2fmod.o \
../PW/add_bfield.o \
../PW/add_efield.o \
../PW/add_vuspsi.o \
../PW/add_vuspsi_nc.o \
../PW/addusdens.o \
../PW/addusforce.o \
../PW/addusstress.o \
../PW/allocate_bp_efield.o \
../PW/allocate_fft.o \
../PW/allocate_locpot.o \
../PW/allocate_nlpot.o \
../PW/allocate_wfc.o \
../PW/atomic_rho.o \
../PW/atomic_wfc.o \
../PW/average_pp.o \
../PW/becmod.o \
../PW/bp_c_phase.o \
../PW/bp_calc_btq.o \
../PW/bp_qvan3.o \
../PW/bp_strings.o \
../PW/buffers.o \
../PW/c_bands.o \
../PW/c_phase_field.o \
../PW/ccgdiagg.o \
../PW/cdiagh.o \
../PW/cdiaghg.o \
../PW/cegterg.o \
../PW/cft3.o \
../PW/cft3s.o \
../PW/checkallsym.o \
../PW/checksym.o \
../PW/clean_pw.o \
../PW/close_files.o \
../PW/compute_becsum.o \
../PW/compute_deff.o \
../PW/compute_dip.o \
../PW/compute_rho.o \
../PW/compute_rho_new.o \
../PW/compute_fes_grads.o \
../PW/compute_scf.o \
../PW/compute_qdipol.o \
../PW/compute_qdipol_so.o \
../PW/compute_ux.o \
../PW/coset.o \
../PW/cryst_to_car.o \
../PW/cubicsym.o \
../PW/d_matrix.o \
../PW/data_structure.o \
../PW/davcio.o \
../PW/deriv_drhoc.o \
../PW/divide_class.o \
../PW/divide_class_so.o \
../PW/realus.o \
../PW/diropn.o \
../PW/divide.o \
../PW/divide_et_impera.o \
../PW/dqvan2.o \
../PW/drhoc.o \
../PW/dvloc_of_g.o \
../PW/dynamics_module.o \
../PW/efermig.o \
../PW/efermit.o \
../PW/electrons.o \
../PW/eqvect.o \
../PW/ewald.o \
../PW/ewald_dipole.o \
../PW/exx.o \
../PW/exx_loop.o \
../PW/find_group.o \
../PW/forces_bp_efield.o \
../PW/force_cc.o \
../PW/force_corr.o \
../PW/force_ew.o \
../PW/force_hub.o \
../PW/force_lc.o \
../PW/force_us.o \
../PW/forces.o \
../PW/g_psi.o \
../PW/g_psi_mod.o \
../PW/gen_at_dj.o \
../PW/gen_at_dy.o \
../PW/gen_us_dj.o \
../PW/gen_us_dy.o \
../PW/get_locals.o \
../PW/ggen.o \
../PW/gk_sort.o \
../PW/gradcorr.o \
../PW/gweights.o \
../PW/g2_kin.o \
../PW/h_epsi_her_apply.o \
../PW/h_epsi_her_set.o \
../PW/h_1psi.o \
../PW/h_psi.o \
../PW/h_psi_meta.o \
../PW/hexsym.o \
../PW/hinit0.o \
../PW/hinit1.o \
../PW/init_ns.o \
../PW/init_run.o \
../PW/init_us_1.o \
../PW/init_us_2.o \
../PW/init_at_1.o \
../PW/init_vloc.o \
../PW/input.o \
../PW/interpolate.o \
../PW/inverse_s.o \
../PW/io_rho_xml.o \
../PW/ions.o \
../PW/irrek.o \
../PW/iweights.o \
../PW/kpoint_grid.o \
../PW/lchk_tauxk.o \
../PW/make_pointlists.o \
../PW/makov_payne.o \
../PW/memory_report.o \
../PW/mix_pot.o \
../PW/mix_rho.o \
../PW/move_ions.o \
../PW/multable.o \
../PW/n_plane_waves.o \
../PW/new_ns.o \
../PW/new_occ.o \
../PW/ns_adj.o \
../PW/newd.o \
../PW/noncol.o \
../PW/non_scf.o \
../PW/openfil.o \
../PW/orthoatwfc.o \
../PW/output_tau.o \
../PW/para.o \
../PW/paw_init.o \
../PW/paw_onecenter.o \
../PW/potinit.o \
../PW/print_clock_pw.o \
../PW/print_ks_energies.o \
../PW/psymrho.o \
../PW/psymrho_mag.o \
../PW/punch.o \
../PW/pw_restart.o \
../PW/pwcom.o \
../PW/qvan2.o \
../PW/rcgdiagg.o \
../PW/rdiagh.o \
../PW/rdiaghg.o \
../PW/read_conf_from_file.o \
../PW/read_file.o \
../PW/read_pseudo.o \
../PW/regterg.o \
../PW/remove_atomic_rho.o \
../PW/report_mag.o \
../PW/reset_k_points.o \
../PW/restart_from_file.o \
../PW/restart_in_electrons.o \
../PW/restart_in_ions.o \
../PW/rgen.o \
../PW/rho2zeta.o \
../PW/rotate_wfc.o \
../PW/rotate_wfc_k.o \
../PW/rotate_wfc_gamma.o \
../PW/ruotaijk.o \
../PW/s_1psi.o \
../PW/s_axis_to_ca.o \
../PW/s_psi.o \
../PW/save_in_cbands.o \
../PW/save_in_electrons.o \
../PW/save_in_ions.o \
../PW/scale_h.o \
../PW/scf_mod.o \
../PW/seqopn.o \
../PW/set_fft_dim.o \
../PW/set_hubbard_l.o \
../PW/set_kplusq.o \
../PW/set_kup_and_kdw.o \
../PW/set_rhoc.o \
../PW/set_vrs.o \
../PW/setlocal.o \
../PW/setqf.o \
../PW/setup.o \
../PW/setupkpt.o \
../PW/sgam_at.o \
../PW/sgam_at_mag.o \
../PW/sgama.o \
../PW/spinor.o \
../PW/sph_ind.o \
../PW/startup.o \
../PW/stop_run.o \
../PW/stres_cc.o \
../PW/stres_ewa.o \
../PW/stres_gradcorr.o \
../PW/stres_har.o \
../PW/stres_hub.o \
../PW/stres_knl.o \
../PW/stres_loc.o \
../PW/stres_us.o \
../PW/stress.o \
../PW/struct_fact.o \
../PW/sum_band.o \
../PW/sumkg.o \
../PW/sumkt.o \
../PW/summary.o \
../PW/symme.o \
../PW/symmetrize_at.o \
../PW/symrho.o \
../PW/symrho_mag.o \
../PW/symscalar.o \
../PW/symtns.o \
../PW/symvect.o \
../PW/symz.o \
../PW/tabd.o \
../PW/transform_becsum_so.o \
../PW/transform_becsum_nc.o \
../PW/trntns.o \
../PW/trnvecc.o \
../PW/trnvect.o \
../PW/tweights.o \
../PW/update_pot.o \
../PW/usnldiag.o \
../PW/v_of_rho.o \
../PW/vcsmd.o \
../PW/vcsubs.o \
../PW/vhpsi.o \
../PW/vloc_of_g.o \
../PW/vloc_psi.o \
../PW/w0gauss.o \
../PW/w1gauss.o \
../PW/wfcinit.o \
../PW/wgauss.o \
../PW/write_ns.o \
../PW/wsweight.o \
../PW/weights.o \
../PW/ortho_wfc.o \
../PW/wannier_proj.o \
../PW/wannier_init.o \
../PW/wannier_check.o \
../PW/wannier_clean.o \
../PW/wannier_occ.o \
../PW/wannier_enrg.o
                                                         
EEOBJS = \
../EE/dvdr_tao.o                           \
../EE/add_dcc_field.o                      \
../EE/v_h_from_rho_r.o                     \
../EE/multiscale.o                         \
../EE/writetofile.o                        \
../EE/mg_pb_solver.o                       \
../EE/add_boundary.o                       \
../EE/add_dccdil_forces.o                  \
../EE/init_ee.o                            \
../EE/calc_ecomp.o                         \
../EE/setlocalcoul.o                       \
../EE/set_fft_dim_coarse.o                 \
../EE/set_mltgrid_dim.o                    \
../EE/data_structure_coarse.o              \
../EE/add_ele_corr.o                       \
../EE/write_ee_summary.o                   \
../EE/ggen_coarse.o

PHOBJS2=../PH/*.o

PHOBJS= \
../PH/add_dkmds.o \
../PH/add_for_charges.o \
../PH/add_zstar_ue.o \
../PH/add_zstar_ue_us.o \
../PH/addcore.o \
../PH/adddvepsi_us.o \
../PH/adddvscf.o \
../PH/addnlcc.o \
../PH/addnlcc_zstar_eu_us.o \
../PH/addusdbec.o \
../PH/addusdbec_nc.o \
../PH/addusddens.o \
../PH/addusddense.o \
../PH/addusdynmat.o \
../PH/addusldos.o \
../PH/allocate_part.o \
../PH/allocate_phq.o \
../PH/apply_dpot.o \
../PH/bcast_ph_input.o \
../PH/bcast_ph_input1.o \
../PH/cg_psi.o \
../PH/ccg_psi.o \
../PH/cgsolve_all.o \
../PH/cgsolve_all_imfreq.o \
../PH/ch_psi_all.o \
../PH/check_restart_recover.o \
../PH/cch_psi_all.o \
../PH/clinear.o \
../PH/close_phq.o \
../PH/compute_alphasum.o \
../PH/compute_becalp.o \
../PH/compute_becsum_ph.o \
../PH/compute_drhous.o \
../PH/compute_drhous_nc.o \
../PH/compute_dvloc.o \
../PH/compute_nldyn.o \
../PH/compute_vsgga.o \
../PH/compute_weight.o \
../PH/deallocate_part.o \
../PH/deallocate_phq.o \
../PH/d2ionq.o \
../PH/davcio_drho.o \
../PH/dielec.o \
../PH/drho.o \
../PH/drhodv.o \
../PH/drhodvloc.o \
../PH/drhodvnl.o \
../PH/drhodvus.o \
../PH/dv_of_drho.o \
../PH/dvanqq.o \
../PH/dvkb3.o \
../PH/dvpsi_e.o \
../PH/dvqpsi_us.o \
../PH/dvqpsi_us_only.o \
../PH/dyndia.o \
../PH/dynmat0.o \
../PH/dynmat_us.o \
../PH/dynmatcc.o \
../PH/dynmatrix.o \
../PH/ef_shift.o \
../PH/elph.o \
../PH/elphon.o \
../PH/find_equiv_sites.o \
../PH/find_mode_sym.o \
../PH/generate_dynamical_matrix_c.o \
../PH/generate_effective_charges_c.o \
../PH/gmressolve_all.o \
../PH/h_psiq.o \
../PH/incdrhoscf.o \
../PH/incdrhoscf_nc.o \
../PH/incdrhous.o \
../PH/incdrhous_nc.o \
../PH/io_pattern.o \
../PH/localdos.o \
../PH/mode_group.o \
../PH/newdq.o \
../PH/openfilq.o \
../PH/orthogonalize.o \
../PH/phcom.o \
../PH/ph_restart.o \
../PH/phq_init.o \
../PH/phq_readin.o \
../PH/phq_recover.o \
../PH/phq_setup.o \
../PH/phq_summary.o \
../PH/phqscf.o \
../PH/polariz.o \
../PH/print_clock_ph.o \
../PH/psidspsi.o \
../PH/psymdvscf.o \
../PH/psyme.o \
../PH/psym_dmag.o \
../PH/psym_dmage.o \
../PH/punch_plot_e.o \
../PH/punch_plot_ph.o \
../PH/q_points.o \
../PH/q2qstar_ph.o \
../PH/random_matrix.o \
../PH/rotate_and_add_dyn.o \
../PH/save_ph_input.o \
../PH/set_asr_c.o \
../PH/set_defaults_pw.o \
../PH/set_drhoc.o \
../PH/set_int12_nc.o \
../PH/set_irr.o \
../PH/set_irr_mode.o \
../PH/set_irr_nosym.o \
../PH/setlocq.o \
../PH/setqmod.o \
../PH/setup_dgc.o \
../PH/sgam_ph.o \
../PH/smallgq.o \
../PH/solve_e.o \
../PH/solve_e_fpol.o \
../PH/solve_linter.o \
../PH/star_q.o \
../PH/stop_ph.o \
../PH/summarize.o \
../PH/sym_and_write_zue.o \
../PH/sym_def.o \
../PH/sym_dmag.o \
../PH/sym_dmage.o \
../PH/symdvscf.o \
../PH/symdyn_munu.o \
../PH/symdynph_gq.o \
../PH/syme.o \
../PH/symm.o \
../PH/transform_int_so.o \
../PH/transform_int_nc.o \
../PH/transform_alphasum_nc.o \
../PH/transform_alphasum_so.o \
../PH/transform_dbecsum_so.o \
../PH/transform_dbecsum_nc.o \
../PH/tra_write_matrix.o \
../PH/trntnsc.o \
../PH/write_dyn_on_file.o \
../PH/write_epsilon_and_zeu.o \
../PH/write_matrix.o \
../PH/write_rec.o \
../PH/zstar_eu.o \
../PH/dgradcorr.o \
../PH/zstar_eu_us.o

RAMAN =    \
../PH/cft_wave.o \
../PH/psyme2.o   \
../PH/sym_elop.o \
../PH/trntnsr_3.o \
../PH/chi_test.o \
../PH/dielec_test.o \
../PH/dvpsi_e2.o \
../PH/el_opt.o   \
../PH/hdiag.o    \
../PH/pcgreen.o  \
../PH/ramanm.o   \
../PH/syme2.o    \
../PH/solve_e2.o \
../PH/solve_e_nscf.o \
../PH/write_ramtns.o \
../PH/symram.o   \
../PH/d2mxc.o    \
../PH/raman_mat.o\
../PH/raman.o    \
../PH/dhdrhopsi.o\
../PH/set_dvscf.o
                                      
EEMODS=../EE/ee_mod.o ../EE/gcoarse_mod.o

LROBJS= \
lr_newd.o \
lr_dgradcorr.o \
lr_variables.o \
bcast_lr_input.o \
lr_readin.o \
lr_addusdens.o \
lr_alloc_init.o \
lr_calc_dens.o \
lr_dot.o \
lr_dealloc.o \
lr_ortho.o \
lr_read_wf.o \
lr_normalise.o \
lr_lanczos.o \
lr_apply_liouvillian.o \
lr_main.o \
lr_dv_setup.o \
lr_setup_dgc.o \
lr_solve_e.o \
lr_dvpsi_e.o \
lr_ch_psi_all.o \
lr_cg_psi.o \
lr_cgsolve_all.o \
lr_h_psiq.o \
lr_sm1_psi.o \
lr_dv_of_drho.o \
stop_lr.o \
lr_read_d0psi.o \
lr_restart.o \
lr_write_restart.o \
print_clock_lr.o \
sd0psi.o \
lr_set_boxes_density.o \
lr_charg_resp.o \
lr_init_nfo.o 

MODULES = ../Modules/*.o

all :  tddfpt.x

#tddfpt.x : $(PWOBJS) $(EEOBJS) $(PHOBJS) $(LROBJS) 
#	$(MPIF90) $(LDFLAGS) -o tddfpt.x \
#	$(PWOBJS) $(EEOBJS) $(PHOBJS) $(LROBJS) $(MODULES) $(EEMODS) $(LIBOBJS) $(LIBS) 
#	- ( cd ../bin; ln -fs ../TDDFPT/tddfpt.x . )
#
tddfpt.x : $(LROBJS) $(PWOBJS) $(MODULES) $(EEOBJS) $(PHOBJS) $(RAMAN) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		$(PWOBJS) $(EEOBJS) $(PHOBJS) $(RAMAN) $(LROBJS) $(MODULES) $(EEMODS) $(LIBOBJS) $(LIBS) 
	- ( cd ../bin; ln -fs ../TDDFPT/tddfpt.x . )

dene.x : obm_tests.o 
	$(MPIF90) $(LDFLAGS) -o dene.x \
        $(PWOBJS) $(PHOBJS) obm_tests.o $(MODULES) $(LIBOBJS) $(LIBS)
clean :
	- /bin/rm -f tddfpt.x *.o *~ *.d *.mod *.i work.pc

include make.depend

# DO NOT DELETE
