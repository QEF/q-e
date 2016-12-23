# Makefile for LR_Modules

include ../make.inc

# location of needed modules and included files (if any)
MODFLAGS=$(MOD_FLAG)../FFTXlib $(MOD_FLAG)../Modules $(MOD_FLAG)../PW/src $(MOD_FLAG).

IFLAGS=

LR_MODULES = \
adddvscf.o \
addusdbec.o \
addusdbec_nc.o \
apply_dpot.o \
cft_wave.o \
cgsolve_all.o \
cg_psi.o \
ch_psi_all.o \
commutator_Hx_psi.o \
dfpt_tetra_mod.o \
dynmat_sub.o \
h_prec.o \
incdrhoscf.o \
incdrhoscf_nc.o \
dv_of_drho.o \
dgradcorr.o \
dnonloccorr.o \
dv_vdW_DF.o \
dv_rVV10.o \
newdq.o \
orthogonalize.o \
setup_nscf.o \
set_dbecsum_nc.o \
set_int3_nc.o \
smallgq.o \
lr_sm1_psi.o \
setup_dgc.o \
compute_vsgga.o \
setup_dmuxc.o \
setup_nbnd_occ.o \
setup_alpha_pv.o \
qdipol_cryst.o \
adddvepsi_us.o \
setqmod.o \
lrcom.o

TLDEPS=libfft mods

all : liblrmod.a

liblrmod.a: $(LR_MODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L

include make.depend
