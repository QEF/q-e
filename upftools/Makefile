# Makefile for converters to UPF format

include ../make.rules
include ../make.sys

OBJS = write_upf.o
MODS = ../Modules/kind.o ../Modules/parameters.o ../Modules/pseudo_types.o \
       ../Modules/parser.o ../Modules/io_global.o ../Modules/mp_global.o \
       ../Modules/mp.o 

all : cpmd2upf.x fhi2upf.x fpmd2upf.x ncpp2upf.x oldcp2upf.x \
      rrkj2upf.x uspp2upf.x vdb2upf.x

# dependencies
cpmd2upf.x : cpmd2upf.o nclib.o $(OBJS)
fhi2upf.x : fhi2upf.o nclib.o $(OBJS)
ncpp2upf.x : ncpp2upf.o nclib.o $(OBJS)
oldcp2upf.x : oldcp2upf.o nclib.o $(OBJS)

uspp2upf.x : uspp2upf.o vanderbilt.o $(OBJS)
vdb2upf.x : vdb2upf.o vanderbilt.o $(OBJS)

# compilation rules
cpmd2upf.x fhi2upf.x ncpp2upf.x oldcp2upf.x :
	$(LD) -o $@ $< nclib.o $(OBJS) $(LDFLAGS)

uspp2upf.x vdb2upf.x :
	$(LD) -o $@ $< vanderbilt.o $(OBJS) $(LDFLAGS)

# special cases
read_upf.x : read_upf.o
	$(LD) -o $@ read_upf.o $(LDFLAGS)

rrkj2upf.x : rrkj2upf.o $(OBJS)
	$(LD) -o $@ rrkj2upf.o $(OBJS) $(LDFLAGS)

fpmd2upf.x : fpmd2upf.o $(OBJS)
	$(LD) -o $@ fpmd2upf.o $(OBJS) $(MODS) $(LDFLAGS)

clean :
	- /bin/rm -f  *.x *.o *~ *.F90 *.mod *.d *.i work.pc

include .dependencies
