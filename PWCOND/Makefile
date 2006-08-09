# Makefile for postprocessing (PP)

include ../make.sys

PPOBJS = \
../PP/start_postproc.o \
../PP/stop_pp.o

PWCONDOBJS = \
allocate_cond.o \
bessj.o \
compbs.o \
compbs_2.o \
condcom.o \
condmain.o \
cond_out.o \
do_cond.o \
eigenchnl.o \
form_zk.o \
four.o \
free_mem.o \
gep_x.o \
gramsh.o \
hev_ab.o \
init_cond.o \
init_orbitals.o \
init_gper.o \
integrals.o \
jbloch.o \
kbloch.o \
local.o \
local_set.o \
poten.o \
print_clock_pwcond.o \
rotate.o \
rotproc.o \
save_cond.o \
scatter_back.o \
scatter_forw.o \
summary_band.o \
summary_tran.o \
sunitary.o \
transmit.o

MODULES = \
../Modules/atom.o \
../Modules/basic_algebra_routines.o \
../Modules/berry_phase.o \
../Modules/cell_base.o \
../Modules/clocks.o \
../Modules/check_stop.o \
../Modules/constants.o \
../Modules/constraints_module.o \
../Modules/control_flags.o \
../Modules/electrons_base.o \
../Modules/fft_base.o \
../Modules/fft_scalar.o \
../Modules/fft_types.o \
../Modules/functionals.o \
../Modules/io_files.o \
../Modules/io_global.o \
../Modules/ions_base.o \
../Modules/input_parameters.o \
../Modules/kind.o \
../Modules/mp_global.o \
../Modules/mp_wave.o \
../Modules/mp.o \
../Modules/parallel_include.o \
../Modules/parameters.o \
../Modules/parser.o \
../Modules/printout_base.o \
../Modules/pseudo_types.o \
../Modules/ptoolkit.o \
../Modules/random_numbers.o \
../Modules/read_upf.o \
../Modules/read_uspp.o \
../Modules/recvec.o \
../Modules/shmem_include.o \
../Modules/task_groups.o \
../Modules/splinelib.o \
../Modules/stick_base.o \
../Modules/timestep.o \
../Modules/uspp.o \
../Modules/version.o \
../Modules/wavefunctions.o \
../Modules/xml_io_base.o


PWOBJS = \
../PW/add_bfield.o \
../PW/add_efield.o \
../PW/add_vuspsi.o \
../PW/add_vuspsi_nc.o \
../PW/addusdens.o \
../PW/allocate_fft.o \
../PW/allocate_locpot.o \
../PW/allocate_nlpot.o \
../PW/allocate_wfc.o \
../PW/atomic_rho.o \
../PW/atomic_wfc.o \
../PW/atomic_wfc_nc.o \
../PW/becmod.o \
../PW/bp_calc_btq.o \
../PW/bp_qvan3.o \
../PW/bp_radin.o \
../PW/bp_ylm_q.o \
../PW/c_gemm.o \
../PW/ccalbec.o \
../PW/ccalbec_nc.o \
../PW/ccgdiagg.o \
../PW/cdiagh.o \
../PW/cdiaghg.o \
../PW/cdiisg_nc.o \
../PW/cegterg.o \
../PW/cft3.o \
../PW/cft3s.o \
../PW/cgramg1.o \
../PW/cgramg1_nc.o \
../PW/checkallsym.o \
../PW/checksym.o \
../PW/cinitcgg.o \
../PW/clean_pw.o \
../PW/compute_dip.o \
../PW/compute_rho.o \
../PW/coset.o \
../PW/cryst_to_car.o \
../PW/cubicsym.o \
../PW/d_matrix.o \
../PW/data_structure.o \
../PW/davcio.o \
../PW/deriv_drhoc.o \
../PW/diropn.o \
../PW/divide.o \
../PW/divide_class.o \
../PW/divide_et_impera.o \
../PW/dndepsilon.o \
../PW/dprojdepsilon.o \
../PW/dprojdtau.o \
../PW/dqvan2.o \
../PW/drhoc.o \
../PW/dvloc_of_g.o \
../PW/dynamics_module.o \
../PW/efermig.o \
../PW/efermit.o \
../PW/eqvect.o \
../PW/error_handler.o \
../PW/estimate.o \
../PW/exx.o \
../PW/ewald.o \
../PW/find_group.o \
../PW/fftw.o \
../PW/g_psi.o \
../PW/g_psi_mod.o \
../PW/g_psi_nc.o \
../PW/gen_at_dj.o \
../PW/gen_at_dy.o \
../PW/gen_us_dj.o \
../PW/gen_us_dy.o \
../PW/get_locals.o \
../PW/ggen.o \
../PW/gk_sort.o \
../PW/gradcorr.o \
../PW/gweights.o \
../PW/h_1psi.o \
../PW/h_1psi_nc.o \
../PW/h_epsi_her.o \
../PW/h_psi.o \
../PW/h_psi_nc.o \
../PW/hexsym.o \
../PW/hinit0.o \
../PW/hinit1.o \
../PW/init_ns.o \
../PW/init_us_1.o \
../PW/init_us_2.o \
../PW/init_at_1.o \
../PW/init_vloc.o \
../PW/interpolate.o \
../PW/inverse_s.o \
../PW/io_rho_xml.o \
../PW/irrek.o \
../PW/iweights.o \
../PW/kpoint_grid.o \
../PW/lchk_tauxk.o \
../PW/linmin.o \
../PW/make_pointlists.o \
../PW/mode_group.o \
../PW/multable.o \
../PW/n_plane_waves.o \
../PW/new_ns.o \
../PW/ns_adj.o \
../PW/newd.o \
../PW/noncol.o \
../PW/openfil.o \
../PW/ortho.o \
../PW/orthoatwfc.o \
../PW/output_tau.o \
../PW/para.o \
../PW/potinit.o \
../PW/print_clock_pw.o \
../PW/psymrho.o \
../PW/psymrho_mag.o \
../PW/pwcom.o \
../PW/pw_gemm.o \
../PW/pw_restart.o \
../PW/qvan2.o \
../PW/rcgdiagg.o \
../PW/rdiagh.o \
../PW/rdiaghg.o \
../PW/read_file.o \
../PW/read_ncpp.o \
../PW/readin.o \
../PW/realus.o \
../PW/remove_atomic_rho.o \
../PW/regterg.o \
../PW/report_mag.o \
../PW/restart_from_file.o \
../PW/restart_in_electrons.o \
../PW/restart_in_ions.o \
../PW/rgen.o \
../PW/rho2zeta.o \
../PW/rinitcgg.o \
../PW/rotate_wfc.o \
../PW/rotate_wfc_gamma.o \
../PW/rotate_wfc_nc.o \
../PW/ruotaijk.o \
../PW/s_1psi.o \
../PW/s_axis_to_ca.o \
../PW/s_gemm.o \
../PW/s_psi.o \
../PW/s_psi_nc.o \
../PW/save_in_cbands.o \
../PW/save_in_electrons.o \
../PW/save_in_ions.o \
../PW/scale_h.o \
../PW/seqopn.o \
../PW/set_fft_dim.o \
../PW/set_hubbard_l.o \
../PW/set_kplusq.o \
../PW/set_kup_and_kdw.o \
../PW/set_rhoc.o \
../PW/set_vrs.o \
../PW/setlocal.o \
../PW/setqf.o \
../PW/setupkpt.o \
../PW/show_memory.o \
../PW/smallg_q.o \
../PW/spinor.o \
../PW/sph_ind.o \
../PW/startup.o \
../PW/struct_fact.o \
../PW/sum_band.o \
../PW/sumkg.o \
../PW/sumkt.o \
../PW/summary.o \
../PW/swap.o \
../PW/symrho.o \
../PW/symrho_mag.o \
../PW/symtns.o \
../PW/symvect.o \
../PW/tabd.o \
../PW/trntns.o \
../PW/trnvecc.o \
../PW/trnvect.o \
../PW/tweights.o \
../PW/update_pot.o \
../PW/updathes.o \
../PW/upf_to_internal.o \
../PW/usnldiag.o \
../PW/usnldiag_nc.o \
../PW/v_of_rho.o \
../PW/vhpsi.o \
../PW/vhpsi_nc.o \
../PW/vloc_of_g.o \
../PW/vloc_psi.o \
../PW/vpack.o \
../PW/w0gauss.o \
../PW/w1gauss.o \
../PW/wfcinit.o \
../PW/wgauss.o \
../PW/write_ns.o \
../PW/wsweight.o

TLDEPS= bindir mods libs pw pp

all : tldeps pwcond.x

pwcond.x : $(PWCONDOBJS) $(PWOBJS) $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o pwcond.x \
		$(PWCONDOBJS) $(PWOBJS) $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PWCOND/pwcond.x . )

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :
 
clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.a work.pc

include make.depend
