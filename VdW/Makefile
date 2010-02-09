# Makefile for TF-VdW calculation.

include ../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../Modules \
          $(MOD_FLAG)../EE $(MOD_FLAG)../PW $(MOD_FLAG)../PH $(MOD_FLAG).

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
../PP/openfil_pp.o

QEMODS = ../Modules/libqemod.a
EEOBJS=../EE/libee.a
PWOBJS = ../PW/libpw.a
PHOBJS = ../PH/libph.a

TLDEPS= bindir mods libs pw ph pp

all : tldeps vdw.x
 
vdw.x : vdw.o $(VdWOBJS) $(PPOBJS) $(PWOBJS) $(PHOBJS) $(QEMODS) $(EEOBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ vdw.o \
		$(VdWOBJS) $(PPOBJS) $(PHOBJS) $(PWOBJS) $(EEOBJS) $(QEMODS) \
		$(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../VdW/$@ . )

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L

include make.depend
