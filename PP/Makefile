# Makefile for postprocessing (PP)

include ../make.sys


PPOBJS = \
addusdens1d.o \
add_shift_cc.o \
add_shift_lc.o \
add_shift_us.o \
atomic_wfc_nc_proj.o \
cft.o \
cgracsc.o \
compute_sigma_avg.o \
cube.o \
cubicspinsym.o \
dosg.o \
do_initial_state.o \
do_shift_ew.o \
d_matrix_nc.o \
d_matrix_so.o \
elf.o \
ggen1d.o \
hexspinsym.o \
local_dos.o \
local_dos_mag.o \
local_dos1d.o \
openfil_pp.o\
plot_io.o \
plot_whole_cell.o \
punch_plot.o \
start_postproc.o \
stm.o \
stop_pp.o \
work_function.o \
xsf.o 

MODULES = \
../Modules/atom.o \
../Modules/basic_algebra_routines.o \
../Modules/cell_base.o \
../Modules/check_stop.o \
../Modules/clocks.o \
../Modules/constants.o \
../Modules/constraints_module.o \
../Modules/control_flags.o \
../Modules/fft_base.o \
../Modules/fft_scalar.o \
../Modules/fft_types.o \
../Modules/functionals.o \
../Modules/io_base.o \
../Modules/io_files.o \
../Modules/io_global.o \
../Modules/ions_base.o \
../Modules/input_parameters.o \
../Modules/kind.o \
../Modules/mp_buffers.o \
../Modules/mp_global.o \
../Modules/mp_wave.o \
../Modules/mp.o \
../Modules/parallel_include.o \
../Modules/parameters.o \
../Modules/parser.o \
../Modules/printout_base.o \
../Modules/pseudo_types.o \
../Modules/readpseudo.o \
../Modules/recvec.o \
../Modules/shmem_include.o \
../Modules/stick_base.o \
../Modules/splinelib.o \
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
../PW/bp_c_phase.o \
../PW/bp_calc_btq.o \
../PW/bp_qvan3.o \
../PW/bp_radin.o \
../PW/bp_strings.o \
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
../PW/cft_3.o \
../PW/cft_sgi.o \
../PW/cft_sp.o \
../PW/cft_sun.o \
../PW/cft_t3e.o \
../PW/cfts_3.o \
../PW/cgramg1.o \
../PW/cgramg1_nc.o \
../PW/checkallsym.o \
../PW/checksym.o \
../PW/cinitcgg.o \
../PW/clean_pw.o \
../PW/close_files.o \
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
../PW/ewald_dipole.o \
../PW/fftw.o \
../PW/g_psi.o \
../PW/g_psi_nc.o \
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
../PW/h_1psi.o \
../PW/h_psi.o \
../PW/h_epsi_her.o \
../PW/h_1psi_nc.o \
../PW/h_psi_nc.o \
../PW/hexsym.o \
../PW/hinit0.o \
../PW/hinit1.o \
../PW/init_ns.o \
../PW/init_paw_1.o \
../PW/init_paw_2.o \
../PW/init_pool.o \
../PW/init_us_1.o \
../PW/init_us_2.o \
../PW/init_at_1.o \
../PW/init_vloc.o \
../PW/interpolate.o \
../PW/inverse_s.o \
../PW/io_pot.o \
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
../PW/paw.o \
../PW/potinit.o \
../PW/print_clock_pw.o \
../PW/psymrho.o \
../PW/psymrho_mag.o \
../PW/pwcom.o \
../PW/pw_gemm.o \
../PW/pw_restart.o \
../PW/qvan2.o \
../PW/rcgdiagg.o \
../PW/rdiaghg.o \
../PW/read_file.o \
../PW/readin.o \
../PW/read_pseudo.o \
../PW/readvan.o \
../PW/readnewvan.o \
../PW/read_ncpp.o \
../PW/realus.o \
../PW/report_mag.o \
../PW/restart.o\
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
../PW/scala_cdiag.o \
../PW/scala_cdiaghg.o \
../PW/scala_utils.o \
../PW/scale_h.o \
../PW/scopy_t3e.o \
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
../PW/symscalar.o \
../PW/symtns.o \
../PW/symvect.o \
../PW/symz.o \
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

TLDEPS= bindir mods libs pw

all : tldeps average.x bands.x dos.x efg.x plotband.x plotrho.x \
      pmw.x pp.x projwfc.x pw2casino.x pw2wan.x pw2wannier90.x pw_export.x \
      voronoy.x initial_state.x dipole.x sumpdos.x plan_avg.x

average.x : average.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		average.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

bands.x : bands.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		bands.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

dipole.x : dipole.o $(PPOBJS) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		dipole.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

dos.x : dos.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		dos.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

efg.x : efg.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		efg.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

plan_avg.x : plan_avg.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		plan_avg.o $(PPOBJS) $(PWOBJS) $(MODULES) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

plotband.x : plotband.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		plotband.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

plotrho.x : plotrho.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		plotrho.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pmw.x : poormanwannier.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		poormanwannier.o $(PPOBJS) $(PWOBJS) $(MODULES) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

xctest.x : xctest.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		xctest.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pp.x : postproc.o chdens.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		postproc.o chdens.o $(PPOBJS) $(PWOBJS) $(MODULES) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

projwfc.x : projwfc.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		projwfc.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pw2casino.x : pw2casino.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		pw2casino.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pw2wan.x : pw2wan.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		pw2wan.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pw2wannier90.x : pw2wannier90.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		pw2wannier90.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

voronoy.x : voronoy.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		voronoy.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

initial_state.x : initial_state.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		initial_state.o $(PPOBJS) $(PWOBJS) $(MODULES) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pw2gw.x : pw2gw.o $(PPOBJS) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		pw2gw.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pw_export.x : pw_export.o $(PWOBJS) $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		pw_export.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

sumpdos.x : sumpdos.o 
	$(MPIF90) $(LDFLAGS) -o $@ sumpdos.o
	- ( cd ../bin ; ln -fs ../PP/$@ . )
 
tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i work.pc

include make.depend
