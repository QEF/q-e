# Makefile for converters to UPF format

include ../make.sys

OBJS = write_upf.o 
MODS = ../Modules/kind.o ../Modules/parameters.o ../Modules/pseudo_types.o \
       ../Modules/parser.o ../Modules/io_global.o ../Modules/mp_global.o \
       ../Modules/io_files.o \
       ../Modules/mp.o ../Modules/parallel_include.o ../Modules/constants.o \
       ../flib/gridsetup.o
SPLINE = ../Modules/splinelib.o
FUNCT = ../Modules/functionals.o
UPF_MOD =../Modules/read_upf.o ../flib/flush_unit.o ../flib/matches.o ../flib/error.o
TLDEPS = mods libs libiotk

all : tldeps cpmd2upf.x fhi2upf.x fpmd2upf.x ncpp2upf.x oldcp2upf.x \
      read_upf_tofile.x rrkj2upf.x uspp2upf.x vdb2upf.x virtual.x 

cpmd2upf.x : cpmd2upf.o $(OBJS) $(MODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ cpmd2upf.o $(OBJS) $(MODS) $(LIBOBJS) $(LIBS)

fhi2upf.x : fhi2upf.o $(OBJS) $(MODS) $(FUNCT) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ fhi2upf.o $(OBJS) $(MODS) $(FUNCT) $(LIBOBJS) $(LIBS)

fpmd2upf.x : fpmd2upf.o $(OBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ fpmd2upf.o $(OBJS) $(LIBOBJS) $(MODS) $(LIBS)

ncpp2upf.x : ncpp2upf.o $(OBJS) $(MODS) $(FUNCT) $(LIBOBJS) 
	$(LD) $(LDFLAGS) -o $@ ncpp2upf.o $(OBJS) $(MODS) $(FUNCT) $(LIBOBJS) $(LIBS)

oldcp2upf.x : oldcp2upf.o $(OBJS) $(MODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ oldcp2upf.o $(OBJS) $(MODS) $(LIBOBJS) $(LIBS)

read_upf.x : read_upf.o
	$(LD) $(LDFLAGS) -o $@ read_upf.o $(LIBS)

read_upf_tofile.x : read_upf_tofile.o $(MODS) $(UPF_MOD) 
	$(LD) $(LDFLAGS) -o $@ read_upf_tofile.o $(UPF_MOD) $(MODS) $(LIBS)

rrkj2upf.x : rrkj2upf.o $(OBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ rrkj2upf.o $(OBJS) $(LIBOBJS) $(LIBS)

uspp2upf.x : uspp2upf.o vanderbilt.o $(OBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		uspp2upf.o vanderbilt.o $(OBJS) $(LIBOBJS) $(LIBS)

vdb2upf.x : vdb2upf.o vanderbilt.o $(OBJS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		vdb2upf.o vanderbilt.o $(OBJS) $(LIBOBJS) $(LIBS)

virtual.x : virtual.o $(OBJS) $(MODS) $(SPLINE) $(FUNCT) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		virtual.o $(OBJS) $(MODS) $(SPLINE) $(FUNCT) $(LIBOBJS) $(LIBS)

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f  *.x *.o *~ *.F90 *.mod *.d *.i *.L

include make.depend
