# Makefile for postprocessing (PP)

include ../make.rules
include ../make.sys

PPOBJS = \
addusdens1d.o \
add_shift_cc.o \
add_shift_lc.o \
add_shift_us.o \
cft.o \
cgracsc.o \
cube.o \
dosg.o \
do_initial_state.o \
do_shift_ew.o \
elf.o \
ggen1d.o \
local_dos.o \
local_dos1d.o \
openfil_pp.o\
plan_avg.o \
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
../Modules/clocks.o \
../Modules/check_stop.o \
../Modules/constants.o \
../Modules/control_flags.o \
../Modules/fft_base.o \
../Modules/fft_scalar.o \
../Modules/fft_types.o \
../Modules/functionals.o \
../Modules/io_base.o \
../Modules/io_files.o \
../Modules/io_global.o \
../Modules/ions_base.o \
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
../Modules/uspp.o \
../Modules/version.o \
../Modules/wavefunctions.o

PWOBJS = \
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
../PW/bachel.o \
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
../PW/ccgdiagg_nc.o \
../PW/cdiagh.o \
../PW/cdiaghg.o \
../PW/cdiisg_nc.o \
../PW/cegterg.o \
../PW/cegterg_nc.o \
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
../PW/cinitcgg_nc.o \
../PW/close_files.o \
../PW/compute_dip.o \
../PW/compute_rho.o \
../PW/constraints_module.o \
../PW/coset.o \
../PW/cryst_to_car.o \
../PW/cubicsym.o \
../PW/d_matrix.o \
../PW/data_structure.o \
../PW/davcio.o \
../PW/delta_e.o \
../PW/deriv_drhoc.o \
../PW/diropn.o \
../PW/divide.o \
../PW/divide_et_impera.o \
../PW/dndepsilon.o \
../PW/dndtau.o \
../PW/dprojdepsilon.o \
../PW/dprojdtau.o \
../PW/dqvan2.o \
../PW/drhoc.o \
../PW/dvloc_of_g.o \
../PW/efermig.o \
../PW/efermit.o \
../PW/eqvect.o \
../PW/error_handler.o \
../PW/estimate.o \
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
../PW/gradcorr_nc.o \
../PW/gweights.o \
../PW/h_1psi.o \
../PW/h_psi.o \
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
../PW/io_pot.o \
../PW/irrek.o \
../PW/iweights.o \
../PW/kpoint_grid.o \
../PW/lchk_tauxk.o \
../PW/linmin.o \
../PW/make_pointlists.o \
../PW/mix_pot.o \
../PW/mix_rho.o \
../PW/mix_rho_nc.o \
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
../PW/qvan2.o \
../PW/rdiaghg.o \
../PW/read_file.o \
../PW/report_mag.o \
../PW/restart.o\
../PW/rgen.o \
../PW/rho2zeta.o \
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
../PW/v_of_rho_nc.o \
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

EXOBJS = \
 ../Modules/iotk_base.o \
 ../Modules/iotk_interface.o \
 ../Modules/iotk_module.o \
 ../Modules/iotk_external.o \
 ../Modules/iotk_multitype.o \
 ../Modules/iotk_multitype_character.o \
 ../Modules/iotk_multitype_complex.o \
 ../Modules/iotk_multitype_integer.o \
 ../Modules/iotk_multitype_logical.o \
 ../Modules/iotk_multitype_real.o 

all : average.x bands.x chdens.x dos.x efg.x plotband.x plotrho.x \
      pmw.x pp.x projwfc.x pw2casino.x pw2wan.x voronoy.x initial_state.x

average.x : average.o $(PPOBJS)
	$(LD) -o $@ average.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

bands.x : bands.o $(PPOBJS)
	$(LD) -o $@ bands.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

chdens.x : chdens.o $(PPOBJS)
	$(LD) -o $@ chdens.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

dos.x : dos.o $(PPOBJS)
	$(LD) -o $@ dos.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

efg.x : efg.o $(PPOBJS)
	$(LD) -o $@ efg.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

plotband.x : plotband.o $(PPOBJS)
	$(LD) -o $@ plotband.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

plotrho.x : plotrho.o $(PPOBJS)
	$(LD) -o $@ plotrho.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pmw.x : poormanwannier.o $(PPOBJS)
	$(LD) -o $@ poormanwannier.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

xctest.x : xctest.o $(PPOBJS)
	$(LD) -o $@ xctest.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pp.x : postproc.o $(PPOBJS)
	$(LD) -o $@ postproc.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

projwfc.x : projwfc.o $(PPOBJS)
	$(LD) -o $@ projwfc.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pw2casino.x : pw2casino.o $(PPOBJS)
	$(LD) -o $@ pw2casino.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pw2wan.x : pw2wan.o $(PPOBJS)
	$(LD) -o $@ pw2wan.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

voronoy.x : voronoy.o $(PPOBJS)
	$(LD) -o $@ voronoy.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

initial_state.x : initial_state.o $(PPOBJS)
	$(LD) -o $@ initial_state.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

pw_export.x : pw_export.o $(EXOBJS) $(PPOBJS)
	$(LD) -o $@ pw_export.o $(EXOBJS) $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../PP/$@ . )

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i work.pc

include make.depend
