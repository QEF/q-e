# Makefile for Gamma-only phonon code

include ../../make.inc

VPATH = @srcdir@

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../dft-d3 \
          $(MOD_FLAG)../../LR_Modules

CGOBJS = \
a_h.o \
cg_readin.o \
cg_setup.o \
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
rhod2vkb.o \
solve_e.o \
solve_ph.o \
writedyn.o 

PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a 
QEMODS = ../../LR_Modules/liblrmod.a $(BASEMODS) 

# see comments in file tldeps for explanations
TLDEPS=pwlibs
include ../../install/tldeps

all : phcg.x

phcg.x : phcg.o libphcg.a $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
		phcg.o libphcg.a $(PWOBJS) $(QEMODS) $(QELIBS)
	- ( cd ../../bin ; ln -fs ../PHonon/Gamma/phcg.x . )

libphcg.a : $(CGOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.d *.i *~ *_tmp.f90 *.mod *.L *.a
	- /bin/rm -f ../../bin/phcg.x

include make.depend


