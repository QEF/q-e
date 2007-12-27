# Makefile for phonon (PH)

include ../make.sys

PHOBJS = \
add_dkmds.o \
add_for_charges.o \
add_zstar_ue.o \
add_zstar_ue_us.o \
addcore.o \
adddvepsi_us.o \
adddvscf.o \
addnlcc.o \
addnlcc_zstar_eu_us.o \
addusdbec.o \
addusdbec_nc.o \
addusddens.o \
addusddense.o \
addusdynmat.o \
addusldos.o \
allocate_part.o \
allocate_phq.o \
bcast_ph_input.o \
bcast_ph_input1.o \
cg_psi.o \
ccg_psi.o \
cgsolve_all.o \
ch_psi_all.o \
cch_psi_all.o \
clinear.o \
close_phq.o \
compute_alphasum.o \
compute_becalp.o \
compute_becsum_ph.o \
compute_drhous.o \
compute_drhous_nc.o \
compute_dvloc.o \
compute_nldyn.o \
compute_vsgga.o \
compute_weight.o \
deallocate_part.o \
deallocate_phq.o \
d2ionq.o \
davcio_drho.o \
dgcxc.o \
dgcxc_spin.o \
dgradcorr.o \
dielec.o \
drho.o \
drhodv.o \
drhodvloc.o \
drhodvnl.o \
drhodvus.o \
dv_of_drho.o \
dvanqq.o \
dvkb3.o \
dvpsi_e.o \
dvqpsi_us.o \
dvqpsi_us_only.o \
dyndia.o \
dynmat0.o \
dynmat_us.o \
dynmatcc.o \
dynmatrix.o \
ef_shift.o \
elph.o \
elphon.o \
find_equiv_sites.o \
find_mode_sym.o \
generate_dynamical_matrix_c.o \
gmressolve_all.o \
h_psiq.o \
incdrhoscf.o \
incdrhoscf_nc.o \
incdrhous.o \
incdrhous_nc.o \
io_pattern.o \
localdos.o \
newdq.o \
openfilq.o \
phcom.o \
phq_init.o \
phq_readin.o \
phq_recover.o \
phq_setup.o \
phq_summary.o \
phqscf.o \
polariz.o \
print_clock_ph.o \
psidspsi.o \
psymdvscf.o \
psyme.o \
psym_dmag.o \
psym_dmage.o \
punch_plot_e.o \
punch_plot_ph.o \
q_points.o \
q2qstar_ph.o \
random_matrix.o \
rotate_and_add_dyn.o \
set_asr_c.o \
set_drhoc.o \
set_int12_nc.o \
set_irr.o \
set_irr_mode.o \
set_irr_nosym.o \
setlocq.o \
setqmod.o \
setup_dgc.o \
smallgq.o \
solve_e.o \
solve_e_fpol.o \
solve_linter.o \
star_q.o \
stop_ph.o \
sym_and_write_zue.o \
sym_def.o \
sym_dmag.o \
sym_dmage.o \
symdvscf.o \
symdyn_munu.o \
symdynph_gq.o \
syme.o \
symm.o \
transform_int_so.o \
transform_int_nc.o \
transform_alphasum_nc.o \
transform_alphasum_so.o \
transform_dbecsum_so.o \
transform_dbecsum_nc.o \
tra_write_matrix.o \
trntnsc.o \
xk_wk_collect.o \
write_dyn_on_file.o \
write_epsilon_and_zeu.o \
write_matrix.o \
zstar_eu.o \
zstar_eu_us.o

MODULES = \
../Modules/atom.o \
../Modules/autopilot.o \
../Modules/basic_algebra_routines.o \
../Modules/bfgs_module.o \
../Modules/berry_phase.o \
../Modules/cell_base.o \
../Modules/check_stop.o \
../Modules/clocks.o \
../Modules/constants.o \
../Modules/constraints_module.o \
../Modules/control_flags.o \
../Modules/descriptors.o \
../Modules/dspev_drv.o \
../Modules/electrons_base.o \
../Modules/exc_t.o \
../Modules/fft_base.o \
../Modules/fft_parallel.o \
../Modules/fft_scalar.o \
../Modules/fft_types.o \
../Modules/functionals.o \
../Modules/input_parameters.o \
../Modules/io_files.o \
../Modules/io_global.o \
../Modules/ions_base.o \
../Modules/ions_nose.o \
../Modules/kind.o \
../Modules/metadyn_base.o \
../Modules/metadyn_io.o \
../Modules/metadyn_vars.o \
../Modules/metagga.o \
../Modules/mp_global.o \
../Modules/mp_wave.o \
../Modules/mp.o \
../Modules/mp_base.o \
../Modules/parallel_include.o \
../Modules/parameters.o \
../Modules/parser.o \
../Modules/path_base.o \
../Modules/path_formats.o \
../Modules/path_variables.o \
../Modules/path_opt_routines.o \
../Modules/path_io_routines.o \
../Modules/path_reparametrisation.o \
../Modules/paw_variables.o \
../Modules/printout_base.o \
../Modules/pseudo_types.o \
../Modules/ptoolkit.o \
../Modules/radial_grids.o \
../Modules/random_numbers.o \
../Modules/read_cards.o \
../Modules/read_namelists.o \
../Modules/read_uspp.o \
../Modules/read_upf.o \
../Modules/recvec.o \
../Modules/splinelib.o \
../Modules/stick_base.o \
../Modules/shmem_include.o \
../Modules/task_groups.o \
../Modules/timestep.o \
../Modules/upf_to_internal.o \
../Modules/uspp.o \
../Modules/version.o \
../Modules/vxc_t.o \
../Modules/vxcgc.o \
../Modules/wavefunctions.o \
../Modules/wave_base.o \
../Modules/xml_io_base.o \
../Modules/zhpev_drv.o

RAMAN =    \
cft_wave.o \
psyme2.o   \
sym_elop.o \
trntnsr_3.o \
chi_test.o \
dielec_test.o \
dvpsi_e2.o \
el_opt.o   \
hdiag.o    \
pcgreen.o  \
ramanm.o   \
syme2.o    \
solve_e2.o \
solve_e_nscf.o \
write_ramtns.o \
symram.o   \
d2mxc.o    \
raman_mat.o\
raman.o    \
dhdrhopsi.o\
set_dvscf.o

PWOBJS = ../PW/libpw.a

TLDEPS= bindir mods libs pw

all : tldeps libph.a ph.x dynmat.x matdyn.x q2r.x 

ph.x : phonon.o libph.a $(PWOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o ph.x phonon.o \
		$(MODULES) libph.a $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PH/ph.x . )

dynmat.x : dynmat.o rigid.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
                dynmat.o rigid.o $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PH/$@ . )

matdyn.x : matdyn.o rigid.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		matdyn.o rigid.o $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PH/$@ . )

q2r.x : q2r.o rigid.o $(PWOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		q2r.o rigid.o $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PH/$@ . )

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

libph.a : $(PHOBJS) $(RAMAN)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L

include make.depend
