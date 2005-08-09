# Makefile for converters to UPF format

include ../make.sys

OBJS = write_upf.o errore.o
MODS = ../Modules/kind.o ../Modules/parameters.o ../Modules/pseudo_types.o \
       ../Modules/parser.o ../Modules/io_global.o ../Modules/mp_global.o \
       ../Modules/mp.o ../Modules/parallel_include.o 
SPLINE = ../Modules/splinelib.o
FUNCT = ../Modules/functionals.o
TLDEPS = mods libs libiotk

all : tldeps cpmd2upf.x fhi2upf.x fpmd2upf.x ncpp2upf.x oldcp2upf.x \
      rrkj2upf.x uspp2upf.x vdb2upf.x virtual.x

cpmd2upf.x : cpmd2upf.o $(OBJS) $(LIBOBJS)
	$(LD) -o $@ cpmd2upf.o $(OBJS) $(LIBOBJS) $(LDFLAGS)

fhi2upf.x : fhi2upf.o $(OBJS) $(FUNCT) $(LIBOBJS)
	$(LD) -o $@ fhi2upf.o $(OBJS) $(FUNCT) $(LIBOBJS) $(LDFLAGS)

fpmd2upf.x : fpmd2upf.o $(OBJS) $(LIBOBJS)
	$(LD) -o $@ fpmd2upf.o $(OBJS) $(LIBOBJS) $(MODS) $(LDFLAGS)

ncpp2upf.x : ncpp2upf.o $(OBJS) $(FUNCT) $(LIBOBJS) 
	$(LD) -o $@ ncpp2upf.o $(OBJS) $(FUNCT) $(LIBOBJS) $(LDFLAGS)

oldcp2upf.x : oldcp2upf.o $(OBJS) $(LIBOBJS)
	$(LD) -o $@ oldcp2upf.o $(OBJS) $(LIBOBJS) $(LDFLAGS)

read_upf.x : read_upf.o
	$(LD) -o $@ read_upf.o $(LDFLAGS)

rrkj2upf.x : rrkj2upf.o $(OBJS) $(LIBOBJS)
	$(LD) -o $@ rrkj2upf.o $(OBJS) $(LIBOBJS) $(LDFLAGS)

uspp2upf.x : uspp2upf.o vanderbilt.o $(OBJS) $(LIBOBJS)
	$(LD) -o $@ uspp2upf.o vanderbilt.o $(OBJS) $(LIBOBJS) $(LDFLAGS)

vdb2upf.x : vdb2upf.o vanderbilt.o $(OBJS) $(LIBOBJS)
	$(LD) -o $@ vdb2upf.o vanderbilt.o $(OBJS) $(LIBOBJS) $(LDFLAGS)

virtual.x : virtual.o $(OBJS) $(SPLINE) $(FUNCT) $(LIBOBJS)
	$(LD) -o $@ virtual.o $(OBJS) $(SPLINE) $(FUNCT) $(LIBOBJS) $(LDFLAGS)

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f  *.x *.o *~ *.F90 *.mod *.d *.i work.pc

include make.depend
