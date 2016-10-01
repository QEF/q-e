# Makefile for Gamma-only phonon code

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../LAXlib $(MOD_FLAG)../../FFTXlib \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG).
IFLAGS=

CGOBJS = \
a_h.o \
cg_readin.o \
cg_setup.o \
cg_setupdgc.o \
cgcom.o \
cgsolve.o \
d2ion.o \
dgradcorr.o \
dielec.o \
drhodv.o \
dvb_cc.o \
dvpsi_e.o \
dvpsi_kb.o \
dyndiar.o \
dynmat_init.o \
dynmatcc.o \
find_equiv_sites.o \
generate_dynamical_matrix.o \
generate_effective_charges.o \
h_h.o \
macro.o \
phcg.o \
rhod2vkb.o \
solve_e.o \
solve_ph.o \
writedyn.o 

PWOBJS = ../../PW/src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a

TLDEPS= bindir mods libs pw

all : tldeps phcg.x libcg.a

phcg.x : $(PWOBJS) $(CGOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o phcg.x \
		$(CGOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PHonon/Gamma/phcg.x . )
tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.d *.i *~ *_tmp.f90 *.mod *.L *.a
	- /bin/rm -f ../../bin/phcg.x

libcg.a : $(CGOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

include make.depend
