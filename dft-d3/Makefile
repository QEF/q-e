# Makefile for dft-d3

include ../make.inc

VPATH = @srcdir@

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS)

# list of modules

TESTOBJ = \
test_code.o

QELIBS= ../Modules/libqemod.a ../UtilXlib/libutil.a ../LAXlib/libqela.a

TESTMODS = \
libdftd3qe.a

DFTD3OBJ = \
common.o \
sizes.o \
pars.o \
core.o \
api.o \
dftd3_qe.o 

# TLDEPS should contain the same list of dependencies as in the main Makefile
TLDEPS= mods
# Ensure that "make" run in this subdirectory builds the needed dependencies
include ../install/tldeps

all : libdftd3qe.a 

libdftd3qe.a: $(DFTD3OBJ)
	$(AR) $(ARFLAGS) $@ $(DFTD3OBJ)        
	$(RANLIB) $@    

test.x : tldeps $(TESTOBJ) $(TESTMODS)
	$(LD) $(LDFLAGS) $(MODFLAGS) -o $@ \
        $(TESTOBJ) $(TESTMODS) $(QELIBS) 

clean :
	( /bin/rm -f *.o *.a *.d *.i *.x *~ *_tmp.f90 *.mod *.L ); \


include make.depend


