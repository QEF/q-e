# Makefile for converters to UPF format

include ../make.rules
include ../make.sys

OBJS = write_upf.o

all : rrkj2upf ncpp2upf vdb2upf uspp2upf fhi2upf cpmd2upf oldcp2upf

rrkj2upf :  rrkj2upf.o $(OBJS)
	$(LD) -o $@  rrkj2upf.o $(OBJS) $(LDFLAGS)

ncpp2upf :  ncpp2upf.o $(OBJS) nclib.o
	$(LD) -o $@  ncpp2upf.o nclib.o $(OBJS) $(LDFLAGS)

vdb2upf :  vdb2upf.o $(OBJS) vanderbilt.o
	$(LD) -o $@  vdb2upf.o  vanderbilt.o write_upf.o $(LDFLAGS)

uspp2upf :  uspp2upf.o $(OBJS) vanderbilt.o
	$(LD) -o $@  uspp2upf.o vanderbilt.o write_upf.o $(LDFLAGS)

fhi2upf :  fhi2upf.o $(OBJS) nclib.o
	$(LD) -o $@  fhi2upf.o nclib.o write_upf.o $(LDFLAGS)

cpmd2upf :  cpmd2upf.o $(OBJS) nclib.o
	$(LD) -o $@  cpmd2upf.o nclib.o write_upf.o $(LDFLAGS)

oldcp2upf :  oldcp2upf.o $(OBJS) nclib.o
	$(LD) -o $@  oldcp2upf.o nclib.o write_upf.o $(LDFLAGS)

fpmd2upf.x : fpmd2upf.o
	$(LD) -o $@ $(LDFLAGS) fpmd2upf.o \
        ../Modules/*.o ../flib/ptools.a ../flib/flib.a ../clib/clib.a $(LIBS)

read_upf :
	$(LD) -o $@  read_upf.o $(LDFLAGS)

clean :
	- /bin/rm -f  *2upf *.o *~ *.F90 *.mod *.d *.i work.pc

any2upf.o: $(OBJS)
cpmd2upf.o: $(OBJS)
fhi2upf.o: $(OBJS)
ncpp2upf.o: $(OBJS)
oldcp2upf.o: $(OBJS)
rrkj2upf.o: $(OBJS)
uspp2upf.o: $(OBJS)
vanderbilt.o: $(OBJS)
vdb2upf.o: $(OBJS)
