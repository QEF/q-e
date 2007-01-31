# Makefile for CP/FPMD

include	../make.sys

FOBJS = \
atoms_type.o \
bessel.o \
berryion.o \
bforceion.o \
cg.o \
cg_sub.o \
cglib.o \
chargedensity.o \
chargemix.o \
chi2.o \
compute_fes_grads.o \
compute_scf.o \
cp_autopilot.o \
cp_emass.o \
cp_fpmd.o \
cp_interfaces.o \
cp_restart.o \
cplib.o \
cpr_mod.o \
cpr.o \
cprsub.o \
dealloc.o \
dforceb.o \
efermi.o \
efield.o \
eigs0.o \
electrons.o \
emptystates.o \
ensemble_dft.o \
environment.o \
exch_corr.o \
fftdrv.o \
fft.o \
forces.o \
fromscra.o \
gradrho.o \
greenf.o \
gtable.o \
guess.o \
init.o \
init_run.o \
inner_loop.o \
inner_loop_cold.o \
input.o \
ions.o \
ions_positions.o \
ksstates.o \
main.o \
mainvar.o \
main_loops.o \
cplib_meta.o \
metaxc.o \
modules.o \
move_electrons.o \
nl_base.o \
nlcc.o \
nl.o \
ortho_base.o \
ortho.o \
para.o \
path_routines.o \
periodic.o \
phasefactor.o \
polarization.o \
potentials.o \
pres_ai_mod.o \
print_out.o \
problem_size.o \
pseudo_base.o \
pseudopot.o \
pseudopot_sub.o \
qmatrixd.o \
qqberry.o \
read_pseudo.o \
redis.o \
restart.o \
restart_sub.o \
runcp.o \
smcp.o \
smd_modules.o \
smd.o  \
smlam.o \
spharmonic.o \
spline.o \
stop_run.o \
stress.o \
turbo.o \
util.o \
vanderwaals.o \
vol_clu.o \
cp_version.o \
wannier_base.o \
wannier.o \
waveinit.o \
wave.o \
wave_types.o \
wf.o \
$(WRAPPERS)

LOBJS = \
adjef.o \
entropy.o \
forceconv.o \
geninv.o \
indices.o

MODULES = \
../Modules/atom.o \
../Modules/autopilot.o \
../Modules/basic_algebra_routines.o \
../Modules/berry_phase.o \
../Modules/cell_base.o \
../Modules/check_stop.o \
../Modules/clocks.o \
../Modules/constants.o \
../Modules/constraints_module.o \
../Modules/control_flags.o \
../Modules/descriptors.o \
../Modules/electrons_base.o \
../Modules/energies.o \
../Modules/fft_base.o \
../Modules/fft_scalar.o \
../Modules/fft_types.o \
../Modules/functionals.o \
../Modules/griddim.o \
../Modules/input_parameters.o \
../Modules/io_files.o \
../Modules/io_global.o \
../Modules/ions_base.o \
../Modules/ions_nose.o \
../Modules/kind.o \
../Modules/mp_global.o \
../Modules/mp_wave.o \
../Modules/mp.o \
../Modules/metadyn_base.o \
../Modules/metadyn_io.o \
../Modules/metadyn_vars.o \
../Modules/parallel_types.o \
../Modules/path_base.o \
../Modules/path_formats.o \
../Modules/path_variables.o \
../Modules/path_opt_routines.o \
../Modules/path_io_routines.o \
../Modules/path_reparametrisation.o \
../Modules/parallel_include.o \
../Modules/parameters.o \
../Modules/parser.o \
../Modules/printout_base.o \
../Modules/pseudo_types.o \
../Modules/ptoolkit.o \
../Modules/random_numbers.o \
../Modules/read_cards.o \
../Modules/read_namelists.o \
../Modules/read_upf.o \
../Modules/read_uspp.o \
../Modules/recvec.o \
../Modules/shmem_include.o \
../Modules/sic.o \
../Modules/smallbox.o \
../Modules/splinelib.o \
../Modules/stick_base.o \
../Modules/task_groups.o \
../Modules/timestep.o \
../Modules/uspp.o \
../Modules/version.o \
../Modules/wavefunctions.o \
../Modules/wave_base.o \
../Modules/xml_io_base.o

WRAPPERS = wrapper.o
TLDEPS= bindir mods libs libiotk

all : cp
cp : tldeps libcp.a cp.x cppp.x

cp.x : cprstart.o libcp.a $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o cp.x \
		cprstart.o $(MODULES) libcp.a $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../CPV/cp.x . )

libcp.a : $(FOBJS) $(LOBJS) 
	 $(AR) $(ARFLAGS) $@ $?
	 $(RANLIB) $@

cp_version.o : cpver.h

cpver.h :
	echo "CHARACTER(LEN=70), PARAMETER :: version_date = '"`date`"'" \
		> cpver.h

PPOBJS = \
../Modules/xml_io_base.o \
../Modules/mp.o \
../Modules/mp_global.o \
../Modules/io_global.o \
../Modules/io_files.o \
../Modules/mp_wave.o \
../Modules/parser.o \
../Modules/kind.o \
../Modules/parameters.o ../Modules/parallel_include.o ../Modules/constants.o

cppp.x : fpmdpp.o $(PPOBJS) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o cppp.x fpmdpp.o $(PPOBJS) $(LIBOBJS) $(LIBS)
	- (cd ../bin ; ln -fs ../CPV/cppp.x . )

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f cppp.x *.o *.mod cpver.h *.i core* *.F90 fort* \
		*.cpp *.d *.L *.a *.s cp.x

include make.depend
