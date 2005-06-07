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
cptypes.o \
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
input.o \
interfaces_main.o \
interfaces.o \
ions.o \
ions_nose.o \
ions_positions.o \
ksstates.o \
macdep.o \
main.o \
mainvar.o \
main_loops.o \
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
pseudotab_base.o \
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
stop_pw.o \
stress.o \
turbo.o \
util.o \
vanderwaals.o \
cp_version.o \
vofrho2.o \
wannier.o \
waveinit.o \
wave.o \
wave_types.o \
wf.o \
$(WRAPPERS)


LIBOBJ = \
adjef.o \
enthropy.o \
forceconv.o \
geninv.o \
indices.o \
miller.o

WRAPPERS = wrapper.o


all : cp.x fpmd.x fpmdpp.x
fpmd : fpmd.x fpmdpp.x
cp : cp.x

cp.x : cprstart.o $(FOBJS) $(LIBOBJ)
	$(LD) -o cp.x cprstart.o $(FOBJS) $(LIBOBJ) ../Modules/*.o ../flib/eispack.o $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../CPV/cp.x . )


cp_version.o : cpver.h

cpver.h :
	echo "CHARACTER(LEN=70), PARAMETER :: version_date = '"`date`"'" \
		> cpver.h

fpmd.x : start.o $(FOBJS) $(LIBOBJ)
	$(LD) -o fpmd.x start.o $(FOBJS) $(LIBOBJ) ../Modules/*.o \
		../flib/eispack.o ../flib/ptools.a ../flib/flib.a ../clib/clib.a \
		$(LDFLAGS)
	- (cd ../bin ; ln -fs ../CPV/fpmd.x . )

fpmdpp.x : fpmdpp.o
	$(LD) -o fpmdpp.x fpmdpp.o $(LDFLAGS)
	- (cd ../bin ; ln -fs ../CPV/fpmdpp.x . )

clean :
	- /bin/rm -f fpmd.x fpmdpp.x *.o *.mod version.h *.i core* *.F90 fort* \
		*.cpp *.d work.pc  *.s cp.x

include make.depend
