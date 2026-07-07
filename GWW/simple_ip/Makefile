# Makefile for simple_ip

include ../../make.inc

VPATH = @srcdir@

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../pw4gww $(MOD_FLAG)../../PW/src 

SIMPLEIPOBJS = \
simple_ip_objects.o \
start_end.o \
input_simple_ip.o \
tetra_mod1.o \
interpolation.o \
diagonalization.o \
dielectric.o  


QEMODS = $(BASEMODS)
PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a

LIBMIN= 

# TLDEPS should contain the same list of dependencies as in the main Makefile
TLDEPS=pwlibs
# Ensure that "make" run in this subdirectory builds the needed dependencies
include ../../install/tldeps

all : simple_ip.x  libsimple_ip.a

simple_ip.x : simple_ip.o $(SIMPLEIPOBJS) $(PWOBJS) $(QEMODS)  $(LIBMIN) 
	$(MPIF90) $(LDFLAGS) -o $@ \
		simple_ip.o $(SIMPLEIPOBJS) $(PWOBJS) $(QEMODS) $(LIBMIN) $(QELIBS) 
	- ( cd ../../bin; ln -fs ../GWW/simple_ip/$@ . )

libsimple_ip.a : $(SIMPLEIPOBJS)
	$(AR) $(ARFLAGS) $@ $(SIMPLEIPOBJS)

clean :
	- /bin/rm -f -v simple_ip.x  *.o *~ *.F90 *.d *.mod *.i work.pc

include make.depend
