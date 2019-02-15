# Makefile for converters to UPF format

include ../make.inc

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS)

OBJS = write_upf.o read_upf.o
QEMODS = ../Modules/libqemod.a ../FFTXlib/libqefft.a  ../UtilXlib/libutil.a 

TLDEPS = libs mods

all : tldeps casino2upf.x cpmd2upf.x fhi2upf.x fpmd2upf.x \
     ncpp2upf.x oldcp2upf.x read_upf_tofile.x rrkj2upf.x upf2casino.x \
     uspp2upf.x vdb2upf.x virtual_v2.x interpolate.x fix_upf.x

casino2upf.x : casino2upf.o casino_pp.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ casino2upf.o casino_pp.o $(OBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

upf2upf2.x : upf2upf2.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ upf2upf2.o $(OBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

cpmd2upf.x : cpmd2upf.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ cpmd2upf.o $(OBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

fhi2upf.x : fhi2upf.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ fhi2upf.o $(OBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

fpmd2upf.x : fpmd2upf.o $(OBJS) $(LIBOBJS) 
	$(LD) $(LDFLAGS) -o $@ fpmd2upf.o $(OBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

ncpp2upf.x : ncpp2upf.o $(OBJS) $(QEMODS) $(LIBOBJS) 
	$(LD) $(LDFLAGS) -o $@ ncpp2upf.o $(OBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

oldcp2upf.x : oldcp2upf.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ oldcp2upf.o $(OBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

read_upf.x : read_ps.o
	$(LD) $(LDFLAGS) -o $@ read_ps.o $(OBJS)

read_upf_tofile.x : read_upf_tofile.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ read_upf_tofile.o $(OBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)


fix_upf.x: fix_upf.o $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ fix_upf.o $(QEMODS) $(LIBOBJS) $(QELIBS)

rrkj2upf.x : rrkj2upf.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ rrkj2upf.o $(OBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

upf2casino.x : upf2casino.o casino_pp.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ upf2casino.o casino_pp.o $(OBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

uspp2upf.x : uspp2upf.o vanderbilt.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		uspp2upf.o vanderbilt.o $(OBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

vdb2upf.x : vdb2upf.o vanderbilt.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		vdb2upf.o vanderbilt.o $(OBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)


virtual_v2.x: virtual_v2.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		virtual_v2.o $(OBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

interpolate.x : interpolate.o $(OBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		interpolate.o $(OBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)

tldeps:
	if test -n "$(TLDEPS)" ; then \
	( cd ../ ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f  *.x *.o *~ *_tmp.f90 *.mod *.d *.i *.L

include make.depend
