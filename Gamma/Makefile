#
# temporary Makefile for Gamma-only code
#
include ../make.rules
include ../make.sys
#
# targets
#
PWOBJS=../PW/pwcom.o \
../PW/para.o \
../PW/pencils.o \
../PW/io.o \
../PW/funct.o \
../PW/g_psi_mod.o \
../PW/aainit.o \
../PW/addusforce.o \
../PW/addusstress.o \
../PW/allocate_locpot.o \
../PW/allocate_nlpot.o \
../PW/allowed.o \
../PW/atomic_wfc.o \
../PW/bachel.o \
../PW/bfgs.o \
../PW/cdiagh.o \
../PW/cft.o \
../PW/cft3.o \
../PW/cft3s.o \
../PW/cft_3.o \
../PW/cfts_3.o \
../PW/cft_fftw.o \
../PW/cft_sgi.o \
../PW/cft_sp.o \
../PW/checkallsym.o \
../PW/checksym.o \
../PW/clocks.o \
../PW/constrain.o \
../PW/conv_to_num.o \
../PW/coset.o \
../PW/cryst_to_car.o \
../PW/cubicsym.o \
../PW/d_matrix.o \
../PW/date_and_tim.o \
../PW/davcio.o \
../PW/delta_e.o \
../PW/deriv_drhoc.o \
../PW/diropn.o \
../PW/divide.o \
../PW/divide_et_impera.o \
../PW/dqvan2.o \
../PW/drhoc.o \
../PW/dsum.o \
../PW/dvloc_of_g.o \
../PW/dylmr2.o \
../PW/dynamics.o \
../PW/efermig.o \
../PW/efermit.o \
../PW/eqvect.o \
../PW/error.o \
../PW/estimate.o \
../PW/ewald.o \
../PW/force_cc.o \
../PW/force_corr.o \
../PW/force_ew.o \
../PW/force_lc.o \
../PW/forces.o \
../PW/functionals.o \
../PW/g_psi.o \
../PW/gen_us_dj.o \
../PW/gen_us_dy.o \
../PW/gk_sort.o \
../PW/good_fft_dimension.o \
../PW/gweights.o \
../PW/hexsym.o \
../PW/hinit0.o \
../PW/hinit1.o \
../PW/init_ns.o \
../PW/init_us_1.o \
../PW/init_us_2.o \
../PW/init_vloc.o \
../PW/invmat.o \
../PW/io_pot.o \
../PW/ions.o \
../PW/irrek.o \
../PW/iweights.o \
../PW/kpoint_grid.o \
../PW/latgen.o \
../PW/lchk_tauxk.o \
../PW/linmin.o \
../PW/lsda_functionals.o \
../PW/maximum.o \
../PW/mix_pot.o \
../PW/mode_group.o \
../PW/move_ions.o \
../PW/multable.o \
../PW/n_plane_waves.o \
../PW/new_ns.o \
../PW/openfil.o \
../PW/potinit.o \
../PW/print_clock_pw.o \
../PW/punch.o \
../PW/qvan2.o \
../PW/random.o \
../PW/read_ncpp.o \
../PW/rgen.o \
../PW/readin.o \
../PW/readnewvan.o \
../PW/readvan.o \
../PW/recips.o \
../PW/remove_atomic_rho.o \
../PW/restart_from_file.o \
../PW/restart_in_electrons.o \
../PW/restart_in_ions.o \
../PW/rho2zeta.o \
../PW/setv.o \
../PW/ruotaijk.o \
../PW/s_axis_to_ca.o \
../PW/save_in_cbands.o \
../PW/save_in_electrons.o \
../PW/save_in_ions.o \
../PW/saveall.o \
../PW/scale_h.o \
../PW/swap.o \
../PW/scnds.o \
../PW/scopy_t3e.o \
../PW/seqopn.o \
../PW/set_fft_dim.o \
../PW/set_kplusq.o \
../PW/set_kup_and_kdw.o \
../PW/set_pencils.o \
../PW/setqf.o \
../PW/set_vrs.o \
../PW/setup.o \
../PW/setupkpt.o \
../PW/sgama.o \
../PW/sgam_at.o \
../PW/sgam_ph.o \
../PW/show_memory.o \
../PW/smallg_q.o \
../PW/sph_bes.o \
../PW/stop_pw.o \
../PW/struct_fact.o \
../PW/sumkg.o \
../PW/read_pseudo.o \
../PW/stress.o \
../PW/stres_cc.o \
../PW/stres_ewa.o \
../PW/stres_gradcorr.o \
../PW/stres_har.o \
../PW/stres_knl.o \
../PW/stres_loc.o \
../PW/summary.o \
../PW/sumkt.o \
../PW/symrho.o \
../PW/symtns.o \
../PW/symvect.o \
../PW/tabd.o \
../PW/trntns.o \
../PW/trnvecc.o \
../PW/trnvect.o \
../PW/tweights.o \
../PW/update_pot.o \
../PW/updathes.o \
../PW/usnldiag.o \
../PW/vcsmd.o \
../PW/vcsubs.o \
../PW/vhpsi.o \
../PW/vloc_of_g.o \
../PW/volume.o \
../PW/vpack.o \
../PW/w0gauss.o \
../PW/w1gauss.o \
../PW/wgauss.o \
../PW/which_dft.o \
../PW/write_ns.o \
../PW/ylmr2.o \
../PW/input.o \
../PW/write_config_to_file.o \
../PW/read_conf_from_file.o  \
../PW/startup.o \
../PW/reduce.o \
../PW/poolreduce.o \
../PW/broadcast.o \
../PW/poolrecover.o \
../PW/poolscatter.o \
../PW/psymrho.o \
../PW/poolbcast.o \
../PW/fft_scatter.o \
../PW/poolextreme.o \
../PW/check.o   \
../PW/gather.o  \
../PW/scatter.o \
../PW/init_pool.o

MODULES = ../Modules/*.o ../PW/error_handler.o \
          ../PW/restart.o ../PW/upf_to_internal.o

GAMMA=rbecmod.o \
gamma.o \
fake.o \
add_vuspsi.o \
addusdens.o \
allocate_fft.o \
allocate_wfc.o \
atomic_rho.o \
c_bands.o \
rdiaghg.o \
regterg.o \
data_structure_para.o \
data_structure_scal.o \
electrons.o \
force_us.o \
ggen.o \
gradcorr.o \
h_psi.o \
init_run.o \
interpolate.o \
mix_rho.o \
newd.o \
pw_gemm.o \
read_file.o \
rotate_wfc.o \
s_psi.o \
set_rhoc.o \
setlocal.o \
stres_us.o \
sum_band.o \
v_of_rho.o \
vloc_psi.o \
wfcinit.o

CGOBJS = cgcom.o cg_readin.o a_h.o cgsolve.o d2ion.o dynmat_init.o \
         drhodv.o dvb_cc.o dvpsi_kb.o dyndiar.o dynmatcc.o  dielec.o \
         zvscal.o find_equiv_sites.o generate_dynamical_matrix.o   \
         generate_effective_charges.o h_h.o macro.o dvpsi_e.o \
         rhod2vkb.o solve_ph.o solve_e.o cg_summary.o cg_setup.o \
         cg_setupdgc.o  writedyn.o pw_dot.o \
         dmxc.o dgradcorr.o dgcxc_spin.o dgcxc.o raman.o
#
# targets
#
include .dependencies

all: pwg.x phcg.x

pwg.x:  $(GAMMA) $(PWOBJS) pwscf.o
	$(LD) -o pwg.x pwscf.o $(GAMMA) $(PWOBJS) $(MODULES) $(LDFLAGS)

phcg.x: $(GAMMA) $(PWOBJS) $(CGOBJS)
	$(LD) -o phcg.x $(CGOBJS) $(GAMMA) $(PWOBJS) $(MODULES) $(LDFLAGS)

clean_:
	-/bin/rm -f *.x *.o *.d *~ *.F90 *.mod work.pc

