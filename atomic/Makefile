#
# Makefile for atomic code
#

include ../make.rules
include ../make.sys

PWMODS = \
../Modules/kind.o \
../Modules/constants.o \
../Modules/io_global.o \
../Modules/functionals.o \
../Modules/parser.o ../Modules/mp.o ../Modules/mp_global.o \
../PW/functionals.o \
../PW/lsda_functionals.o \
../PW/bachel.o

LD1MOD = \
parameters.o \
ld1inc.o \
pseudo_types.o \
read_pseudo_module.o

LD1OBJ = all_electron.o \
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
date_and_tim.o \
descreening.o \
dir_outward.o \
dirsol.o \
dmixp.o \
do_mesh.o \
el_config.o \
elsd.o \
elsdps.o \
errore.o \
esic.o \
exc_t.o \
find_qi.o \
gener_pseudo.o \
hartree.o \
int_0_inf_dr.o \
integrate_inward.o \
integrate_outward.o \
intref.o \
invmat1.o \
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
openfil.o \
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
sph_besr.o \
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
write_upf.o capital.o matches.o

# sph_bes.o erf.o random.o bachel.o
# funct.o which_dft.o capital.o matches.o functionals.o lsda_functionals.o 


all: $(LD1OBJ) ld1.o
	$(LD) -o ld1.x ld1.o $(LD1OBJ) $(PWMODS) $(LD1MOD) $(LDFLAGS)

clean:
	-/bin/rm -f *.o *.x $(LIB) *~ *.tmp.* make.* work.pc

include .dependencies
