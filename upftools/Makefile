# Makefile for converters to UPF format

include ../make.rules
include ../make.sys

OBJS = write_upf.o
MODS = ../Modules/parser.o ../Modules/io_global.o ../Modules/mp_global.o \
       ../Modules/mp.o

all : rrkj2upf.x ncpp2upf.x vdb2upf.x uspp2upf.x fhi2upf.x cpmd2upf.x oldcp2upf.x fpmd2upf.x

rrkj2upf.x :  rrkj2upf.o $(OBJS)
	$(LD) -o $@  rrkj2upf.o $(OBJS) $(LDFLAGS)

ncpp2upf.x :  ncpp2upf.o $(OBJS) nclib.o
	$(LD) -o $@  ncpp2upf.o nclib.o $(OBJS) $(LDFLAGS)

vdb2upf.x :  vdb2upf.o $(OBJS) vanderbilt.o
	$(LD) -o $@  vdb2upf.o  vanderbilt.o write_upf.o $(LDFLAGS)

uspp2upf.x :  uspp2upf.o $(OBJS) vanderbilt.o
	$(LD) -o $@  uspp2upf.o vanderbilt.o write_upf.o $(LDFLAGS)

fhi2upf.x :  fhi2upf.o $(OBJS) nclib.o
	$(LD) -o $@  fhi2upf.o nclib.o write_upf.o $(LDFLAGS)

cpmd2upf.x :  cpmd2upf.o $(OBJS) nclib.o
	$(LD) -o $@  cpmd2upf.o nclib.o write_upf.o $(LDFLAGS)

oldcp2upf.x :  oldcp2upf.o $(OBJS) nclib.o
	$(LD) -o $@  oldcp2upf.o nclib.o write_upf.o $(LDFLAGS)

fpmd2upf.x : fpmd2upf.o
	$(LD) -o $@ fpmd2upf.o $(MODS)  $(LDFLAGS)

read_upf.x :
	$(LD) -o $@ read_upf.o $(LDFLAGS)

clean :
	- /bin/rm -f  *.x *.o *~ *.F90 *.mod *.d *.i work.pc

any2upf.o: $(OBJS)
cpmd2upf.o: $(OBJS)
fhi2upf.o: $(OBJS)
ncpp2upf.o: $(OBJS)
oldcp2upf.o: $(OBJS)
rrkj2upf.o: $(OBJS)
uspp2upf.o: $(OBJS)
vanderbilt.o: $(OBJS)
vdb2upf.o: $(OBJS)
