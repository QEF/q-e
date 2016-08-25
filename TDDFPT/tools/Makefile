# Makefile for TDDFPT tools

include ../../make.inc
IFLAGS         = -I../../include

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../LAXlib $(MOD_FLAG).

QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a
PWOBJS = ../../PW/src/libpw.a

TLDEPS=bindir libs pw

all : tldeps turbo_spectrum.x

turbo_spectrum.x : tddfpt_calculate_spectrum.o $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	tddfpt_calculate_spectrum.o $(QEMODS) $(PWOBJS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../TDDFPT/tools/$@ . )
	- if [ -d ../bin ] ; then  ( cd ../bin ; ln -fs ../tools/$@ . ); fi

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	rm -f ../bin/turbo_spectrum.x *.x *.o *.mod *_tmp.f90

