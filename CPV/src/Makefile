# Makefile for CP/FPMD

include	../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../LAXlib $(MOD_FLAG)../../Modules $(MOD_FLAG).
IFLAGS=

FOBJS = \
berryion.o \
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
cp_restart_new.o \
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
energies.o \
ensemble_dft.o \
exch_corr.o \
exx_cg.o \
exx_es.o \
exx_gs.o \
exx_module.o \
exx_pair.o \
exx_psi.o \
exx_vofr.o \
fft.o \
forces.o \
fromscra.o \
gradrho.o \
gram.o \
gtable.o \
init.o \
init_run.o \
inner_loop_cold.o \
input.o \
ions_nose.o \
ions_positions.o \
ksstates.o \
ldaU_cp.o \
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
plugin_add_potential.o \
plugin_clean.o \
plugin_clock.o \
plugin_energy.o \
plugin_ext_forces.o \
plugin_get_potential.o \
plugin_init_base.o \
plugin_init_cell.o \
plugin_initialization.o \
plugin_init_ions.o \
plugin_int_forces.o \
plugin_print_energies.o \
plugin_print_info.o \
plugin_read_input.o \
plugin_utilities.o \
potentials.o \
pres_ai_mod.o \
print_out.o \
printout_base.o \
problem_size.o \
pseudo_base.o \
pseudopot.o \
pseudopot_sub.o \
qmatrixd.o \
qqberry.o \
restart.o \
restart_sub.o \
runcp.o \
sic.o \
smallbox.o \
smallbox_grid.o \
smallbox_gvec.o \
smallbox_subs.o \
smallbox_lib.o \
spline.o \
spinsq.o \
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

QEMODS=../../Modules/libqemod.a
QEFFT=../../FFTXlib/libqefft.a
QELA=../../LAXlib/libqela.a

TLDEPS= bindir mods libs libiotk

all : tldeps libcp.a manycp.x cp.x cppp.x wfdd.x

manycp.x : manycp.o libcp.a $(QEMODS) $(QEFFT) $(QELA) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o manycp.x manycp.o \
		libcp.a $(QEMODS) $(QEFFT) $(QELA) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../CPV/src/manycp.x . )

cp.x : cprstart.o libcp.a $(QEMODS) $(QEFFT) $(QELA) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o cp.x cprstart.o \
		libcp.a $(QEMODS) $(QEFFT) $(QELA) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../CPV/src/cp.x . )

libcp.a : $(FOBJS) $(LOBJS) 
	 $(AR) $(ARFLAGS) $@ $?
	 $(RANLIB) $@

cp_version.o : cpver.h

cpver.h :
	echo "CHARACTER(LEN=70), PARAMETER :: version_date = '"`date`"'" \
		> cpver.h

cppp.x : fpmdpp.o $(QEMODS) $(QEFFT) $(QELA) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o cppp.x fpmdpp.o $(QEMODS) $(QEFFT) $(QELA) $(LIBOBJS) $(LIBS)
	- (cd ../../bin ; ln -fs ../CPV/src/cppp.x . )

wfdd.x : wfdd.o
	$(LD) $(LDFLAGS) -o $@ wfdd.o $(QEMODS) $(QEFFT) $(QELA) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../CPV/src/$@ . )

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f cppp.x *.o *.mod cpver.h *.i core* *_tmp.f90 fort* \
		*.cpp *.d *.L *.a *.s cp.x manycp.x wfdd.x
	-/bin/rm -rf ../../bin/manycp.x ../../bin/cppp.x ../../bin/cp.x ../../bin/wfdd.x

include make.depend
