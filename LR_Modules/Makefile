# Makefile for LR_Modules

include ../make.inc

VPATH = @srcdir@

# location of needed modules and included files (if any)
MODFLAGS=$(BASEMOD_FLAGS) $(MOD_FLAG)$(BUILDDIR)/PW/src

LR_MODULES = \
adddvscf.o \
adddvhubscf.o \
addusdbec.o \
addusdbec_nc.o \
apply_dpot.o \
apply_dpot_mod.o \
apply_trev.o \
addusldos.o \
adddvepsi_us.o \
c_bands_ph.o \
cft_wave.o \
cgsolve_all.o \
ccgsolve_all.o \
cg_psi.o \
ccg_psi.o \
c_bands_ph.o \
ch_psi_all.o \
ch_psi_all_complex.o \
Coul_cut_2D_ph.o \
check_q_points_sym.o \
dnsq_scf.o \
dfpt_tetra_mod.o \
dfpt_type.o \
dynmat_sub.o \
h_prec.o \
h_prec_freq.o \
incdrhoscf_mod.o \
dv_of_drho.o \
dgradcorr.o \
dnonloccorr.o \
dv_vdW_DF.o \
dv_rVV10.o \
newdq.o \
non_scf_ph.o \
set_intq_nc.o \
lr_sym_mod.o \
sym_def.o \
psymdvscf.o \
sym_dns.o \
symdvscf.o \
sym_dmag.o \
orthogonalize.o \
orthogonalize_omega.o \
setup_nscf.o \
set_dbecsum_nc.o \
set_int3_nc.o \
set_kplusq_nc.o \
smallgq.o \
lr_sm1_psi.o \
lr_addusddens.o \
localdos.o \
setup_dgc.o \
compute_vsgga.o \
setup_dmuxc.o \
setup_nbnd_occ.o \
setup_alpha_pv.o \
set_small_group_of_q.o \
sgam_lr.o \
star_q.o \
setqmod.o \
lr_orthoUwfc.o \
qdipol_cryst.o \
mode_group.o \
mix_pot.o \
lrcom.o \
qpoint.o \
lr_dot.o \
lr_dot_magnons.o \
lr_nc_mag.o \
lr_two_chem.o \
twochem_postproc_dfpt.o \
lr_sym_paw.o \
lr_restart.o \
lanczos_nonhermitian.o \
lanczos_pseudohermitian.o \
lanczos_pseudohermitian_c.o \
lanczos_nonhermitian_c.o \
response_kernels.o \
non_scf_ph.o \
efermi_shift.o \
dfpt_kernels.o

# see comments in file tldeps for explanations
TLDEPS=mods pwlibs
include ../install/tldeps

all : liblrmod.a

liblrmod.a: $(LR_MODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L


include make.depend
