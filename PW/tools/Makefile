# Makefile for tools

include ../../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
	  $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../src $(MOD_FLAG). 

PWOBJS = ../src/libpw.a
QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a

TLDEPS= bindir mods libs pw libfft

LIBOBJS        = ../../flib/ptools.a ../../flib/flib.a ../../clib/clib.a ../../iotk/src/libiotk.a

all : tldeps dist.x ev.x kpoints.x pwi2xsf.x

dist.x : dist.o $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		dist.o $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PW/tools/$@ . )

ev.x : ev.o ev_xml.o $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		ev.o ev_xml.o $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PW/tools/$@ . )

kpoints.x : kpoints.o $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		kpoints.o $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PW/tools/$@ . )

pwi2xsf.x : pwi2xsf.o $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		pwi2xsf.o $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	- ( cd ../../bin ; ln -fs ../PW/tools/$@ . )

tldeps:
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f pwi2xsf pwi2xsf_old *.x *.o *~ *.F90 *.mod *.d *.i *.L

include make.depend
