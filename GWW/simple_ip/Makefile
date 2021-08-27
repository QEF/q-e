# Makefile for simple_ip

include ../../make.inc

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

TLDEPS=pwlibs

all : tldeps simple_ip.x  libsimple_ip.a

simple_ip.x : simple_ip.o $(SIMPLEIPOBJS) $(PWOBJS) $(QEMODS)  $(LIBMIN) 
	$(MPIF90) $(LDFLAGS) -o $@ \
		simple_ip.o $(SIMPLEIPOBJS) $(PWOBJS) $(QEMODS) $(LIBMIN) $(QELIBS) 
	- ( cd ../../bin; ln -fs ../GWW/simple_ip/$@ . )

libsimple_ip.a : $(SIMPLEIPOBJS)
	$(AR) $(ARFLAGS) $@ $(SIMPLEIPOBJS)

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f -v simple_ip.x  *.o *~ *.F90 *.d *.mod *.i work.pc

include make.depend

# DO NOT DELETE
