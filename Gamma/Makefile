# Makefile for Gamma-only phonon code

include ../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../Modules \
          $(MOD_FLAG)../PW $(MOD_FLAG).

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
pw_dot.o \
phcg.o \
rhod2vkb.o \
solve_e.o \
solve_ph.o \
writedyn.o \
zvscal.o

PWOBJS = ../PW/libpw.a
QEMODS = ../Modules/libqemod.a

TLDEPS= bindir mods libs pw

all : tldeps phcg.x

phcg.x : $(PWOBJS) $(CGOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o phcg.x \
		$(CGOBJS) $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../Gamma/phcg.x . )
tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *.d *.i *~ *.F90 *.mod *.L

include make.depend
