# Makefile for TF-VdW calculation.

include ../../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../FFTXlib \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG)../../PHonon/PH $(MOD_FLAG).

VdWOBJS = \
allocate_vdw.o \
cegterg.o \
check_v_eff.o \
ch_psi_all.o \
dv_of_drho.o \
dvpsi_e.o \
eff_pot.o \
h_psiq.o \
incdrhoscf.o \
pbcg_psi.o \
polariz.o \
print_clock_vdw.o \
solve_e.o \
stop_vdw.o \
vdwcom.o

PPOBJS = \
../../PP/src/openfil_pp.o

QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a
PWOBJS = ../../PW/src/libpw.a
PHOBJS = ../../PHonon/PH/libph.a

LIBOBJS = ../../flib/ptools.a ../../flib/flib.a ../../clib/clib.a ../../iotk/src/libiotk.a

TLDEPS= bindir mods libs pw ph pp

all : tldeps vdw.x

vdw.x : vdw.o $(VdWOBJS) $(PPOBJS) $(PWOBJS) $(PHOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ vdw.o \
		$(VdWOBJS) $(PPOBJS) $(PHOBJS) $(PWOBJS) $(QEMODS) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../VdW/src/$@ . )

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L

include make.depend
