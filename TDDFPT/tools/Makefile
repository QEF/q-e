# Makefile for TDDFPT tools

include ../../make.sys
IFLAGS         = -I../../include
LIBOBJS        = ../../flib/ptools.a ../../flib/flib.a ../../clib/clib.a ../../iotk/src/libiotk.a ../../Multigrid/mglib.a


# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG).

QEMODS = ../../Modules/libqemod.a
PWOBJS = ../../PW/libpw.a
EEOBJS = ../../EE/libee.a

TLDEPS=bindir libs pw

all : tldeps tddfpt_pp.x
	
tddfpt_pp.x : tddfpt_calculate_spectrum.o $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	tddfpt_calculate_spectrum.o $(QEMODS) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../TDDFPT/tools/$@ . )
	- if [ -d ../bin ] ; then  ( cd ../bin ; ln -fs ../tools/$@ . ); fi

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	rm -f ../bin/tddfpt_pp.x *.x *.o *.mod

