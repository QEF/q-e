#
# Makefile for atomic code
#

include ../make.sys

LD1OBJS = \
add_exchange.o \
all_electron.o \
ascheq.o \
ascheqps.o \
ascheqps_old.o \
ascheqps_drv.o \
c6_dft.o \
c6_tfvw.o \
cfdsol.o \
chargeps.o \
compute_chi.o \
compute_chi_tm.o \
compute_det.o \
compute_phi.o \
compute_phi_tm.o \
compute_phius.o \
compute_potps.o \
compute_solution.o \
descreening.o \
dfx_new.o \
dir_outward.o \
dirsol.o \
dmixp.o \
drhoofv.o \
drho0ofvx.o \
dvex.o \
el_config.o \
elsd.o \
elsd_highv.o \
elsdps.o \
elsdps_paw.o \
esic.o \
find_qi.o \
gener_pseudo.o \
green.o \
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
set_rho_core.o \
set_psi_in.o \
set_rc_rv.o \
set_sl3.o \
sic_correction.o \
starting_potential.o \
start_potps.o \
test_bessel.o \
trou.o \
vdpack.o \
vext.o \
vpack.o \
v_of_rho_at.o \
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
../Modules/atom.o \
../Modules/cell_base.o \
../Modules/constants.o \
../Modules/control_flags.o \
../Modules/exc_t.o \
../Modules/functionals.o \
../Modules/io_global.o \
../Modules/kind.o \
../Modules/ions_base.o \
../Modules/io_files.o \
../Modules/paw_variables.o \
../Modules/pseudo_types.o \
../Modules/radial_grids.o \
../Modules/random_numbers.o \
../Modules/read_upf.o \
../Modules/timestep.o \
../Modules/version.o \
../Modules/parameters.o \
../Modules/parallel_include.o \
../Modules/read_paw.o \
../Modules/mp_global.o \
../Modules/mp.o \
../Modules/uspp.o \
../Modules/vxc_t.o \
../Modules/vxcgc.o \
../PW/startup.o

TLDEPS= bindir mods libs 

all: tldeps ld1.x

ld1.x: $(LD1OBJS) ld1.o $(LD1MODS) $(MODULES) $(PHOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o ld1.x \
		ld1.o $(LD1OBJS) $(LD1MODS) $(MODULES) $(PHOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../atomic/ld1.x . )

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean:
	- /bin/rm -f *.o *.F90 *.x $(LIB) *~ *.d *.mod *.tmp.* *.L

include make.depend
