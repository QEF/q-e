# Makefile for converters to UPF format

include ../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../Modules $(MOD_FLAG)../FFTXlib $(MOD_FLAG)../LAXlib $(MOD_FLAG).

OBJS = write_upf.o
QEMODS = ../Modules/libqemod.a ../FFTXlib/libqefft.a

TLDEPS = mods libs libiotk

all : tldeps casino2upf.x cpmd2upf.x fhi2upf.x fpmd2upf.x \
     ncpp2upf.x oldcp2upf.x read_upf_tofile.x rrkj2upf.x upf2casino.x \
     uspp2upf.x vdb2upf.x virtual.x interpolate.x

casino2upf.x : casino2upf.o casino_pp.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ casino2upf.o casino_pp.o $(OBJS) $(QEMODS) $(LIBOBJS) $(LIBS)

upf2upf2.x : upf2upf2.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ upf2upf2.o $(OBJS) $(QEMODS) $(LIBOBJS) $(LIBS)

cpmd2upf.x : cpmd2upf.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ cpmd2upf.o $(OBJS) $(QEMODS) $(LIBOBJS) $(LIBS)

fhi2upf.x : fhi2upf.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ fhi2upf.o $(OBJS) $(QEMODS) $(LIBOBJS) $(LIBS)

fpmd2upf.x : fpmd2upf.o $(OBJS) $(LIBOBJS) 
	$(LD) $(LDFLAGS) -o $@ fpmd2upf.o $(OBJS) $(QEMODS) $(LIBOBJS) $(LIBS)

ncpp2upf.x : ncpp2upf.o $(OBJS) $(QEMODS) $(LIBOBJS) 
	$(LD) $(LDFLAGS) -o $@ ncpp2upf.o $(OBJS) $(QEMODS) $(LIBOBJS) $(LIBS)

oldcp2upf.x : oldcp2upf.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ oldcp2upf.o $(OBJS) $(QEMODS) $(LIBOBJS) $(LIBS)

read_upf.x : read_upf.o
	$(LD) $(LDFLAGS) -o $@ read_upf.o $(LIBS)

read_upf_tofile.x : read_upf_tofile.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ read_upf_tofile.o $(OBJS) $(QEMODS) $(LIBOBJS) $(LIBS)

rrkj2upf.x : rrkj2upf.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ rrkj2upf.o $(OBJS) $(QEMODS) $(LIBOBJS) $(LIBS)

upf2casino.x : upf2casino.o casino_pp.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ upf2casino.o casino_pp.o $(OBJS) $(QEMODS) $(LIBOBJS) $(LIBS)

uspp2upf.x : uspp2upf.o vanderbilt.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		uspp2upf.o vanderbilt.o $(OBJS) $(QEMODS) $(LIBOBJS) $(LIBS)

vdb2upf.x : vdb2upf.o vanderbilt.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		vdb2upf.o vanderbilt.o $(OBJS) $(QEMODS) $(LIBOBJS) $(LIBS)

virtual.x : virtual.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		virtual.o $(OBJS) $(QEMODS) $(LIBOBJS) $(LIBS)

interpolate.x : interpolate.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		interpolate.o $(OBJS) $(QEMODS) $(LIBOBJS) $(LIBS)

tldeps:
	if test -n "$(TLDEPS)" ; then \
	( cd ../ ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f  *.x *.o *~ *_tmp.f90 *.mod *.d *.i *.L

include make.depend
