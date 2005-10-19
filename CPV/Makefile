# Makefile for CP/FPMD

include	../make.sys

FOBJS = \
atoms_type.o \
band_type.o \
bessel.o \
berryion.o \
bforceion.o \
blacs.o \
brillouin.o \
cg.o \
cg_sub.o \
cglib.o \
chargedensity.o \
chargemix.o \
charge_types.o \
chi2.o \
compute_fes_grads.o \
compute_scf.o \
cp_emass.o \
cp_fpmd.o \
cp_restart.o \
cplib.o \
cpr_mod.o \
cpr.o \
cprsub.o \
crayfft.o \
dealloc.o \
diis.o \
dipol_matrix.o \
dforceb.o \
efermi.o \
efield.o \
eigsp.o \
eigs0.o \
electrons.o \
emptystates.o \
ensemble_dft.o \
environment.o \
exch_corr.o \
fftdrv.o \
fft.o \
fields_type.o \
fnl.o \
forces.o \
fromscra.o \
greenf.o \
grid.o \
gsmesh.o \
gtable.o \
guess.o \
init.o \
init_run.o \
input.o \
interfaces_main.o \
interfaces.o \
ions.o \
ions_positions.o \
ksstates.o \
macdep.o \
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
optical.o \
ortho_base.o \
ortho.o \
para.o \
path_routines.o \
periodic.o \
phasefactor.o \
polarization.o \
potentials.o \
print_out.o \
problem_size.o \
pseudo_base.o \
pseudopot.o \
qmatrixd.o \
qqberry.o \
read_pseudo.o \
redis.o \
restart.o \
restart_sub.o \
rsmesh.o \
runcg_ion.o \
runcg.o \
runcp.o \
rundiis.o \
runsd.o \
scalapack.o \
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
cp_version.o \
vofrho2.o \
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
indices.o \
miller.o

WRAPPERS = wrapper.o
TLDEPS= bindir mods libs libiotk

all : tldeps cp.x cppp.x
cp : tldeps cp.x cppp.x

cp.x : cprstart.o $(FOBJS) $(LOBJS) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o cp.x \
		cprstart.o $(FOBJS) $(LOBJS) ../Modules/*.o \
		../flib/eispack.o $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../CPV/cp.x . )

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
../Modules/kind.o

cppp.x : fpmdpp.o $(PPOBJS) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o cppp.x fpmdpp.o $(PPOBJS) $(LIBOBJS) $(LIBS)
	- (cd ../bin ; ln -fs ../CPV/cppp.x . )

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f cppp.x *.o *.mod version.h *.i core* *.F90 fort* \
		*.cpp *.d work.pc  *.s cp.x

include make.depend
