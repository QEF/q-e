# Makefile for converters to UPF format

include ../make.sys

OBJS = write_upf.o 
MODS = ../Modules/kind.o ../Modules/parameters.o ../Modules/pseudo_types.o \
       ../Modules/parser.o ../Modules/io_global.o ../Modules/mp_global.o \
       ../Modules/mp.o ../Modules/parallel_include.o ../Modules/constants.o
SPLINE = ../Modules/splinelib.o
FUNCT = ../Modules/functionals.o
TLDEPS = mods libs libiotk

all : tldeps cpmd2upf.x fhi2upf.x fpmd2upf.x ncpp2upf.x oldcp2upf.x \
      rrkj2upf.x uspp2upf.x vdb2upf.x virtual.x

cpmd2upf.x : cpmd2upf.o $(OBJS) $(MODS) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ cpmd2upf.o $(OBJS) $(MODS) $(LIBOBJS) $(LIBS)

fhi2upf.x : fhi2upf.o $(OBJS) $(MODS) $(FUNCT) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ fhi2upf.o $(OBJS) $(MODS) $(FUNCT) $(LIBOBJS) $(LIBS)

fpmd2upf.x : fpmd2upf.o $(OBJS) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ fpmd2upf.o $(OBJS) $(LIBOBJS) $(MODS) $(LIBS)

ncpp2upf.x : ncpp2upf.o $(OBJS) $(MODS) $(FUNCT) $(LIBOBJS) 
	$(MPIF90) $(LDFLAGS) -o $@ ncpp2upf.o $(OBJS) $(MODS) $(FUNCT) $(LIBOBJS) $(LIBS)

oldcp2upf.x : oldcp2upf.o $(OBJS) $(MODS) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ oldcp2upf.o $(OBJS) $(MODS) $(LIBOBJS) $(LIBS)

read_upf.x : read_upf.o
	$(MPIF90) $(LDFLAGS) -o $@ read_upf.o $(LIBS)

rrkj2upf.x : rrkj2upf.o $(OBJS) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ rrkj2upf.o $(OBJS) $(LIBOBJS) $(LIBS)

uspp2upf.x : uspp2upf.o vanderbilt.o $(OBJS) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		uspp2upf.o vanderbilt.o $(OBJS) $(LIBOBJS) $(LIBS)

vdb2upf.x : vdb2upf.o vanderbilt.o $(OBJS) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		vdb2upf.o vanderbilt.o $(OBJS) $(LIBOBJS) $(LIBS)

virtual.x : virtual.o $(OBJS) $(MODS) $(SPLINE) $(FUNCT) $(LIBOBJS)
	$(MPIF90) $(LDFLAGS) -o $@ \
		virtual.o $(OBJS) $(MODS) $(SPLINE) $(FUNCT) $(LIBOBJS) $(LIBS)

tldeps:
	test -n "$(TLDEPS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f  *.x *.o *~ *.F90 *.mod *.d *.i work.pc

include make.depend
