#
# Makefile for atomic code
#

include ../make.rules
include ../make.sys

LD1OBJS = \
all_electron.o \
ascheq.o \
ascheqlocps.o \
ascheqps.o \
cfdsol.o \
charge.o \
chargeps.o \
compute_chi.o \
compute_det.o \
compute_phi.o \
compute_phipot.o \
compute_phius.o \
compute_solution.o \
descreening.o \
dir_outward.o \
dirsol.o \
dmixp.o \
do_mesh.o \
el_config.o \
elsd.o \
elsdps.o \
esic.o \
exc_t.o \
find_qi.o \
gener_pseudo.o \
hartree.o \
int_0_inf_dr.o \
integrate_inward.o \
integrate_outward.o \
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
sic_correction.o \
starting_potential.o \
start_potps.o \
trou.o \
vext.o \
vpack.o \
vxcgc.o \
vxc_t.o \
write_paw_recon.o \
write_pseudo.o \
write_results.o \
write_resultsps.o \
write_upf.o

# sph_bes.o erf.o random.o bachel.o
# funct.o which_dft.o capital.o matches.o functionals.o lsda_functionals.o 

LD1MODS = \
ld1inc.o \
parameters.o 

MODULES = \
../Modules/constants.o \
../Modules/functionals.o \
../Modules/io_global.o \
../Modules/kind.o \
../Modules/io_files.o \
../Modules/pseudo_types.o \
../Modules/readpseudo.o \
../PW/bachel.o 

all: ld1.x

ld1.x: $(LD1OBJS) ld1.o
	$(LD) -o ld1.x ld1.o $(LD1OBJS) $(LD1MODS) $(MODULES) $(LDFLAGS)
	- ( cd ../bin ; ln -fs ../atomic/ld1.x . )

clean:
	- /bin/rm -f *.o *.F90 *.x $(LIB) *~ *.d *.mod *.tmp.* work.pc

include make.depend
