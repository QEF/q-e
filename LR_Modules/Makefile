# Makefile for LR_Modules

include ../make.sys

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
h_psiq.o \
incdrhoscf.o \
incdrhoscf_nc.o \
dv_of_drho.o \
dgradcorr.o \
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
	- /bin/rm -f *.o *.a *.d *.i *~ *.F90 *.mod *.L

include make.depend
