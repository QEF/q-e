# Makefile for converters to UPF format

include ../make.rules
include ../make.sys

OBJS = write_upf.o errore.o
MODS = ../Modules/kind.o ../Modules/parameters.o ../Modules/pseudo_types.o \
       ../Modules/parser.o ../Modules/io_global.o ../Modules/mp_global.o \
       ../Modules/mp.o  ../Modules/parallel_include.o
SPLINE = ../Modules/splinelib.o

all : cpmd2upf.x fhi2upf.x fpmd2upf.x ncpp2upf.x oldcp2upf.x \
      rrkj2upf.x uspp2upf.x vdb2upf.x virtual.x

cpmd2upf.x : cpmd2upf.o nclib.o $(OBJS)
	$(LD) -o $@ cpmd2upf.o nclib.o $(OBJS) $(LDFLAGS)

fhi2upf.x : fhi2upf.o nclib.o $(OBJS)
	$(LD) -o $@ fhi2upf.o nclib.o $(OBJS) $(LDFLAGS)

fpmd2upf.x : fpmd2upf.o $(OBJS)
	$(LD) -o $@ fpmd2upf.o $(OBJS) $(MODS) $(LDFLAGS)

ncpp2upf.x : ncpp2upf.o nclib.o $(OBJS)
	$(LD) -o $@ ncpp2upf.o nclib.o $(OBJS) $(LDFLAGS)

oldcp2upf.x : oldcp2upf.o nclib.o $(OBJS)
	$(LD) -o $@ oldcp2upf.o nclib.o $(OBJS) $(LDFLAGS)

read_upf.x : read_upf.o
	$(LD) -o $@ read_upf.o $(LDFLAGS)

rrkj2upf.x : rrkj2upf.o $(OBJS)
	$(LD) -o $@ rrkj2upf.o $(OBJS) $(LDFLAGS)

uspp2upf.x : uspp2upf.o vanderbilt.o $(OBJS)
	$(LD) -o $@ uspp2upf.o vanderbilt.o $(OBJS) $(LDFLAGS)

vdb2upf.x : vdb2upf.o vanderbilt.o $(OBJS)
	$(LD) -o $@ vdb2upf.o vanderbilt.o $(OBJS) $(LDFLAGS)

virtual.x : virtual.o nclib.o $(OBJS) $(SPLINE)
	$(LD) -o $@ virtual.o nclib.o $(OBJS) $(SPLINE) $(LDFLAGS)

clean :
	- /bin/rm -f  *.x *.o *~ *.F90 *.mod *.d *.i work.pc

include make.depend
