# Makefile for CP/FPMD

include	../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../Modules $(MOD_FLAG).

FOBJS = \
atoms_type.o \
berryion.o \
berry_phase.o \
bforceion.o \
cell_nose.o \
cg.o \
cg_sub.o \
cglib.o \
chargedensity.o \
cp_autopilot.o \
cp_emass.o \
cp_interfaces.o \
cp_restart.o \
cplib.o \
cpr_mod.o \
cpr.o \
dealloc.o \
dforceb.o \
efermi.o \
efield.o \
eigs0.o \
electrons.o \
electrons_nose.o \
ensemble_dft.o \
exch_corr.o \
exx_ggrid.o \
exx_gs.o \
exx_es.o \
exx_cg.o \
exx_pair.o \
exx_psi.o \
exx_setup.o \
exx_vofr.o \
fft.o \
forces.o \
fromscra.o \
gradrho.o \
gram.o \
gtable.o \
gvecw.o\
init.o \
init_run.o \
inner_loop_cold.o \
input.o \
ions_nose.o \
ions_positions.o \
ksstates.o \
ldaU.o \
ldaUpen.o \
mainvar.o \
cpr_loop.o \
cplib_meta.o \
metaxc.o \
modules.o \
move_electrons.o \
newd.o \
nl_base.o \
nlcc.o \
ortho_base.o \
ortho.o \
phasefactor.o \
plugin_initialization.o \
plugin_forces.o \
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
restart.o \
restart_sub.o \
runcp.o \
smallbox.o \
smallbox_grid.o \
smallbox_gvec.o \
smallbox_subs.o \
smallbox_lib.o \
spline.o \
stop_run.o \
stress.o \
vol_clu.o \
vofrho.o \
wannier_base.o \
wannier.o \
wave.o \
wave_types.o \
wf.o \
makov_payne.o

LOBJS = \
entropy.o 

QEMODS=../Modules/libqemod.a

TLDEPS= bindir mods libs libiotk

all : cp
cp : tldeps libcp.a cp.x cppp.x

cp.x : cprstart.o libcp.a $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o cp.x cprstart.o \
		libcp.a $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../CPV/cp.x . )

libcp.a : $(FOBJS) $(LOBJS) 
	 $(AR) $(ARFLAGS) $@ $?
	 $(RANLIB) $@

cp_version.o : cpver.h

cpver.h :
	echo "CHARACTER(LEN=70), PARAMETER :: version_date = '"`date`"'" \
		> cpver.h

cppp.x : fpmdpp.o $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o cppp.x fpmdpp.o $(QEMODS) $(LIBOBJS) $(LIBS)
	- (cd ../bin ; ln -fs ../CPV/cppp.x . )

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f cppp.x *.o *.mod cpver.h *.i core* *.F90 fort* \
		*.cpp *.d *.L *.a *.s cp.x

include make.depend
