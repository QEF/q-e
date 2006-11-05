#
# Makefile for atomic code
#

include ../make.sys

LD1OBJS = \
add_exchange.o \
all_electron.o \
ascheq.o \
ascheqlocps.o \
ascheqps.o \
c6_dft.o \
c6_tfvw.o \
cfdsol.o \
chargeps.o \
compute_chi.o \
compute_det.o \
compute_phi.o \
compute_phipot.o \
compute_phius.o \
compute_solution.o \
descreening.o \
dfx_new.o \
dir_outward.o \
dirsol.o \
dmixp.o \
drhoofv.o \
drho0ofvx.o \
do_mesh.o \
dvex.o \
el_config.o \
elsd.o \
elsdps.o \
esic.o \
exc_t.o \
find_qi.o \
gener_pseudo.o \
green.o \
hartree.o \
int_0_inf_dr.o \
integrate_inward.o \
integrate_outward.o \
inward.o \
outward.o \
intref.o \
ld1_readin.o \
ld1_setup.o \
ld1_writeout.o \
lderiv.o \
lderivps.o \
lschps.o \
newd_at.o \
new_potential.o \
nodenum.o \
normalize.o \
occ_spin.o \
pseudovloc.o \
read_newpseudo.o \
read_pseudo.o \
read_pseudoupf.o \
run_pseudo.o \
run_test.o \
scf.o \
seriebes.o \
series.o \
set_rho_core.o \
set_sl3.o \
sic_correction.o \
starting_potential.o \
start_potps.o \
test_bessel.o \
trou.o \
vdpack.o \
vext.o \
vpack.o \
vxcgc.o \
vxc_t.o \
write_cpmd.o \
write_paw_recon.o \
write_pseudo.o \
write_results.o \
write_resultsps.o \
write_upf.o

LD1MODS = \
atomic_paw.o \
ld1inc.o \
parameters.o 

MODULES = \
../Modules/constants.o \
../Modules/functionals.o \
../Modules/io_global.o \
../Modules/kind.o \
../Modules/io_files.o \
../Modules/pseudo_types.o \
../Modules/random_numbers.o \
../Modules/read_upf.o \
../Modules/version.o \
../Modules/parameters.o ../Modules/parallel_include.o

TLDEPS= bindir mods libs 

all: tldeps ld1.x

ld1.x: $(LD1OBJS) ld1.o $(LD1MODS) $(MODULES) $(PHOBJS) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o ld1.x \
		ld1.o $(LD1OBJS) $(LD1MODS) $(MODULES) $(PHOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../atomic/ld1.x . )

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean:
	- /bin/rm -f *.o *.F90 *.x $(LIB) *~ *.d *.mod *.tmp.* work.pc

include make.depend
