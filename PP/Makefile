#
# Makefile for postprocessing (PP)
#
include ../make.rules
include ../make.sys
#
PPOBJS = addusdens1d.o  \
average.o  \
cft.o \
cgracsc.o  \
chdens.o  \
xsf.o \
plot_whole_cell.o \
dosg.o  \
dost.o  \
elf.o  \
ggen1d.o  \
local_dos1d.o  \
local_dos.o  \
plan_avg.o  \
plot_io.o  \
projwave.o  \
punch_band.o  \
punch_plot.o  \
start_postproc.o  \
do_postproc.o  \
stm.o  \
stop_pp.o  \
voronoy.o  \
work_function.o \
wsweight.o

MODULES = ../Modules/*.o

PWOBJS = ../PW/pwcom.o  \
../PW/aainit.o  \
../PW/add_efield.o \
../PW/addusdens.o  \
../PW/addusforce.o  \
../PW/addusstress.o  \
../PW/add_vuspsi.o  \
../PW/allocate_fft.o  \
../PW/allocate_locpot.o  \
../PW/allocate_nlpot.o  \
../PW/allocate_wfc.o  \
../PW/atomic_rho.o  \
../PW/atomic_wfc.o  \
../PW/bachel.o  \
../PW/becmod.o  \
../PW/bfgs.o  \
../PW/broadcast.o  \
../PW/c_bands.o  \
../PW/ccalbec.o  \
../PW/ccgdiagg.o  \
../PW/cdiagh.o  \
../PW/cdiaghg.o  \
../PW/cdiisg.o  \
../PW/cegterg.o  \
../PW/cft_3.o  \
../PW/cft3.o  \
../PW/cft3s.o  \
../PW/cfts_3.o  \
../PW/cft_sgi.o  \
../PW/cft_sp.o  \
../PW/cft_sun.o  \
../PW/cft_t3e.o  \
../PW/cgather_sym.o  \
../PW/c_gemm.o  \
../PW/cgramg1.o  \
../PW/checkallsym.o  \
../PW/check.o  \
../PW/checksym.o  \
../PW/cinitcgg.o  \
../PW/constrain.o  \
../PW/conv_to_num.o  \
../PW/coset.o  \
../PW/cryst_to_car.o  \
../PW/cubicsym.o  \
../PW/data_structure.o  \
../PW/date_and_tim.o  \
../PW/davcio.o  \
../PW/delta_e.o  \
../PW/deriv_drhoc.o  \
../PW/diropn.o  \
../PW/divide_et_impera.o  \
../PW/divide.o  \
../PW/d_matrix.o  \
../PW/dndepsilon.o  \
../PW/dndtau.o  \
../PW/dprojdepsilon.o  \
../PW/dprojdtau.o  \
../PW/dqvan2.o  \
../PW/drhoc.o  \
../PW/dsum.o  \
../PW/dvloc_of_g.o  \
../PW/dylmr2.o  \
../PW/dynamics.o  \
../PW/efermig.o  \
../PW/efermit.o  \
../PW/electrons.o  \
../PW/eqvect.o  \
../PW/error.o  \
../PW/error_handler.o  \
../PW/estimate.o  \
../PW/ewald.o  \
../PW/fft_scatter.o  \
../PW/fftw.o  \
../PW/force_cc.o  \
../PW/force_corr.o  \
../PW/force_ew.o  \
../PW/force_hub.o  \
../PW/force_lc.o  \
../PW/forces.o  \
../PW/force_us.o  \
../PW/funct.o  \
../PW/functionals.o  \
../PW/gather.o  \
../PW/gen_at_dj.o  \
../PW/gen_at_dy.o  \
../PW/gen_us_dj.o  \
../PW/gen_us_dy.o  \
../PW/ggen.o  \
../PW/gk_sort.o  \
../PW/g_psi.o  \
../PW/g_psi_mod.o  \
../PW/gradcorr.o  \
../PW/gweights.o  \
../PW/h_1psi.o  \
../PW/hexsym.o  \
../PW/hinit0.o  \
../PW/hinit1.o  \
../PW/h_psi.o  \
../PW/init_ns.o  \
../PW/init_pool.o  \
../PW/init_run.o  \
../PW/init_us_1.o  \
../PW/init_us_2.o  \
../PW/init_vloc.o  \
../PW/input.o  \
../PW/interpolate.o  \
../PW/invmat.o  \
../PW/io.o  \
../PW/ions.o  \
../PW/io_pot.o  \
../PW/irrek.o  \
../PW/iweights.o  \
../PW/kpoint_grid.o  \
../PW/latgen.o  \
../PW/lchk_tauxk.o  \
../PW/linmin.o  \
../PW/lsda_functionals.o  \
../PW/maximum.o  \
../PW/mix_pot.o  \
../PW/mix_rho.o  \
../PW/mode_group.o  \
../PW/move_ions.o  \
../PW/multable.o  \
../PW/newd.o  \
../PW/new_ns.o  \
../PW/n_plane_waves.o  \
../PW/openfil.o  \
../PW/orthoatwfc.o  \
../PW/ortho.o  \
../PW/para.o  \
../PW/pencils.o  \
../PW/poolbcast.o  \
../PW/poolextreme.o  \
../PW/poolrecover.o  \
../PW/poolreduce.o  \
../PW/poolscatter.o  \
../PW/potinit.o  \
../PW/print_clock_pw.o  \
../PW/psymrho.o  \
../PW/punch.o  \
../PW/qvan2.o  \
../PW/random.o  \
../PW/read_conf_from_file.o  \
../PW/read_file.o  \
../PW/readin.o  \
../PW/read_ncpp.o  \
../PW/readnewvan.o  \
../PW/read_pseudo.o  \
../PW/readvan.o  \
../PW/recips.o  \
../PW/reduce.o  \
../PW/remove_atomic_rho.o  \
../PW/restart.o  \
../PW/restart_from_file.o  \
../PW/restart_in_electrons.o  \
../PW/restart_in_ions.o  \
../PW/rgen.o  \
../PW/rho2zeta.o  \
../PW/rotate_wfc.o  \
../PW/ruotaijk.o  \
../PW/s_1psi.o  \
../PW/saveall.o  \
../PW/save_in_cbands.o  \
../PW/save_in_electrons.o  \
../PW/save_in_ions.o  \
../PW/s_axis_to_ca.o  \
../PW/scala_cdiag.o  \
../PW/scala_cdiaghg.o  \
../PW/scala_utils.o  \
../PW/scale_h.o  \
../PW/scatter.o  \
../PW/scopy_t3e.o  \
../PW/seqopn.o  \
../PW/set_fft_dim.o  \
../PW/set_hubbard_l.o  \
../PW/set_kplusq.o  \
../PW/set_kup_and_kdw.o  \
../PW/setlocal.o  \
../PW/set_pencils.o  \
../PW/setqf.o  \
../PW/set_rhoc.o  \
../PW/setup.o  \
../PW/setupkpt.o  \
../PW/setv.o  \
../PW/set_vrs.o  \
../PW/sgama.o  \
../PW/sgam_at.o  \
../PW/sgam_ph.o  \
../PW/s_gemm.o  \
../PW/show_memory.o  \
../PW/smallg_q.o  \
../PW/s_psi.o  \
../PW/startup.o  \
../PW/stop_pw.o  \
../PW/stres_cc.o  \
../PW/stres_ewa.o  \
../PW/stres_gradcorr.o  \
../PW/stres_har.o  \
../PW/stres_hub.o  \
../PW/stres_knl.o  \
../PW/stres_loc.o  \
../PW/stress.o  \
../PW/stres_us.o  \
../PW/struct_fact.o  \
../PW/sum_band.o  \
../PW/sumkg.o  \
../PW/sumkt.o  \
../PW/summary.o  \
../PW/swap.o  \
../PW/symrho.o  \
../PW/symtns.o  \
../PW/symvect.o  \
../PW/tabd.o  \
../PW/trntns.o  \
../PW/trnvecc.o  \
../PW/trnvect.o  \
../PW/tweights.o  \
../PW/update_pot.o  \
../PW/updathes.o  \
../PW/upf_to_internal.o  \
../PW/usnldiag.o  \
../PW/vcsmd.o  \
../PW/vcsubs.o  \
../PW/vhpsi.o  \
../PW/vloc_of_g.o  \
../PW/v_of_rho.o  \
../PW/volume.o  \
../PW/vpack.o  \
../PW/w0gauss.o  \
../PW/w1gauss.o  \
../PW/wfcinit.o  \
../PW/wgauss.o  \
../PW/which_dft.o  \
../PW/write_config_to_file.o  \
../PW/write_ns.o  \
../PW/ylmr2.o   \
../PW/bp_bess.o \
../PW/bp_calc_btq.o \
../PW/bp_c_phase.o \
../PW/bp_dbess.o \
../PW/bp_qvan3.o \
../PW/bp_radin.o \
../PW/bp_strings.o \
../PW/bp_ylm_q.o \
../PW/bp_zgedi.o \
../PW/bp_zgefa.o 


#
# targets
#
include .dependencies

all: pp.x projwfc.x chdens.x average.x plotrho.x voronoy.x dos.x bands.x pw2wan.x

pp.x: $(PPOBJS) postproc.o
	$(LD) -o pp.x postproc.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)

projwfc.x:  $(PPOBJS) projwfc.o
	$(LD) -o projwfc.x projwfc.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)

chdens.x: $(PPOBJS) chmain.o
	$(LD) -o chdens.x chmain.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)

average.x: $(PPOBJS) avmain.o
	$(LD) -o average.x avmain.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)

voronoy.x:  $(PPOBJS) voromain.o
	$(LD) -o voronoy.x voromain.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)

dos.x:  $(PPOBJS) dos.o
	$(LD) -o dos.x dos.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)

bands.x: $(PPOBJS) bands.o
	$(LD) -o bands.x bands.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)

plotrho.x: $(PPOBJS) plotrho.o
	$(LD) -o plotrho.x plotrho.o ../Modules/parameters.o $(LDFLAGS)

pw2wan.x: $(PPOBJS) pw2wan.o
	$(LD) -o pw2wan.x pw2wan.o $(PPOBJS) $(PWOBJS) $(MODULES) $(LDFLAGS)


clean_:
	-/bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i work.pc
#
