#
# Makefile for atomic code
#

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../LAXlib $(MOD_FLAG).
IFLAGS=

LD1OBJS = \
add_exchange.o \
all_electron.o \
ascheq.o \
ascheqps.o \
ascheqps_drv.o \
c6_dft.o \
c6_tfvw.o \
calculate_gipaw_orbitals.o \
cfdsol.o \
chargeps.o \
compute_chi.o \
compute_chi_tm.o \
compute_phi.o \
compute_phi_tm.o \
compute_phius.o \
compute_potps.o \
compute_potps_new.o \
compute_q_3bess.o \
compute_relpert.o \
compute_solution.o \
descreening.o \
default_conf.o \
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
export_upf.o \
find_qi.o \
gener_pseudo.o \
grad_log.o \
green.o \
import_upf.o \
int_0_inf_dr.o \
integrate_inward.o \
integrate_outward.o \
inward.o \
outward.o \
intref.o \
kin_e_density.o \
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
paw_type.o \
partial_wave_expansion.o \
pseudovloc.o \
pseudo_q.o \
read_pseudo_ncpp.o \
read_pseudo_rrkj3.o \
run_lda_half.o \
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
start_scheq.o \
test_bessel.o \
trou.o \
vdpack.o \
vext.o \
vpack.o \
v_of_rho_at.o \
vxcgc.o \
write_cpmd.o \
write_ae_pseudo.o \
write_files.o \
write_paw_recon.o \
write_pseudo.o \
write_results.o \
write_resultsps.o

LD1MODS = \
atomic_paw.o \
ld1inc.o \
parameters.o 

QEMODS=../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a

TLDEPS= bindir mods libs 

all: tldeps ld1.x

ld1.x: $(LD1OBJS) ld1.o $(LD1MODS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o ld1.x \
		ld1.o $(LD1OBJS) $(LD1MODS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../atomic/src/ld1.x . )

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean:
	- /bin/rm -f *.o *_tmp.f90 *.x *~ *.d *.mod *.tmp.* *.L
	- /bin/rm -f ../../bin/ld1.x

include make.depend
