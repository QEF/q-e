#
# Makefile for converters to UPF format
#
include ../make.rules
include ../make.sys

OBJS =	write_upf.o
all: rrkj2upf ncpp2upf vdb2upf uspp2upf fhi2upf cpmd2upf oldcp2upf

rrkj2upf:  rrkj2upf.o $(OBJS)
	$(LD) -o $@  rrkj2upf.o $(OBJS) $(LFLAGS)

ncpp2upf:  ncpp2upf.o $(OBJS) nclib.o
	$(LD) -o $@  ncpp2upf.o nclib.o $(OBJS) $(LFLAGS)

vdb2upf:  vdb2upf.o $(OBJS) vanderbilt.o
	$(LD) -o $@  vdb2upf.o  vanderbilt.o write_upf.o $(LFLAGS)

uspp2upf:  uspp2upf.o $(OBJS) vanderbilt.o
	$(LD) -o $@  uspp2upf.o vanderbilt.o write_upf.o $(LFLAGS)

fhi2upf:  fhi2upf.o $(OBJS) nclib.o
	$(LD) -o $@  fhi2upf.o nclib.o write_upf.o $(LFLAGS)

cpmd2upf:  cpmd2upf.o $(OBJS) nclib.o
	$(LD) -o $@  cpmd2upf.o nclib.o write_upf.o $(LFLAGS)

oldcp2upf:  oldcp2upf.o $(OBJS) nclib.o
	$(LD) -o $@  oldcp2upf.o nclib.o write_upf.o $(LFLAGS)

fpmd2upf.x : fpmd2upf.o
	$(LNK) -o $@ $(LNKFLAGS) fpmd2upf.o \
        ../Modules/*.o ../flib/ptools.a ../flib/flib.a ../clib/clib.a $(LIB)


read_upf:
	$(LD) -o $@  read_upf.o $(LFLAGS)

clean_:
	rm -f  *2upf *.o *.mod *.d *.F90

any2upf.o:  $(OBJS)
cpmd2upf.o: $(OBJS)
oldcp2upf.o:$(OBJS)
rrkj2upf.o: $(OBJS)
ncpp2upf.o: $(OBJS)
vdb2upf.o:  $(OBJS)
uspp2upf.o: $(OBJS)
fhi2upf.o:  $(OBJS)
vanderbilt.o: $(OBJS)

