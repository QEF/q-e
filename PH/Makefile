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
compute_becsum.o \
compute_drhous.o \
compute_dvloc.o \
compute_nldyn.o \
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
incdrhous.o \
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
punch_plot_e.o \
punch_plot_ph.o \
q_points.o \
q2qstar_ph.o \
random_matrix.o \
rotate_and_add_dyn.o \
set_asr_c.o \
set_drhoc.o \
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
symdvscf.o \
symdyn_munu.o \
symdynph_gq.o \
syme.o \
symm.o \
tra_write_matrix.o \
trntnsc.o \
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
../Modules/electrons_base.o \
../Modules/fft_base.o \
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
../Modules/mp_global.o \
../Modules/mp_wave.o \
../Modules/mp.o \
../Modules/parallel_include.o \
../Modules/parameters.o \
../Modules/parser.o \
../Modules/path_base.o \
../Modules/path_formats.o \
../Modules/path_variables.o \
../Modules/path_opt_routines.o \
../Modules/path_io_routines.o \
../Modules/path_reparametrisation.o \
../Modules/printout_base.o \
../Modules/pseudo_types.o \
../Modules/ptoolkit.o \
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
../Modules/uspp.o \
../Modules/version.o \
../Modules/wavefunctions.o \
../Modules/wave_base.o \
../Modules/xml_io_base.o

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

all : tldeps libph.a ph.x

ph.x : phonon.o libph.a $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o ph.x phonon.o \
		$(MODULES) libph.a $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../PH/ph.x . )

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

libph.a : $(PHOBJS) $(RAMAN)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f ph.x *.o *.a *~ *.F90 *.d *.mod *.i work.pc

include make.depend
