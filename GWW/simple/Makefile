# Makefile for  simple

include ../../make.inc

VPATH = @srcdir@

LIBMIN = ../minpack/minpacklib.a

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG)../pw4gww $(MOD_FLAG)../gww 

SIMPLEOBJS = \
  stop_pp.o \
  read_export.o \
  openfile_simple.o \
  input_simple.o \
  wfc_basis.o \
  product_basis.o \
  v_product.o \
  epe.o  \
  gk_sort_limit.o \
  khamiltonian.o \
  commutator.o 


QEMODS = $(BASEMODS)
PWOBJS = ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a 
GWWOBJ = ../gww/libgww.a

# TLDEPS should contain the same list of dependencies as in the main Makefile
TLDEPS= pwlibs gwwlib
# Ensure that "make" run in this subdirectory builds the needed dependencies
include ../../install/tldeps

all : simple.x

simple.x : simple.o libsimple.a  $(SIMPLEOBJS) $(PWOBJS) $(QEMODS) $(GWWOBJ)
	$(LD) $(LDFLAGS) -o $@ \
		simple.o libsimple.a $(PWOBJS) $(GWWOBJ) $(QEMODS) $(QELIBS) $(LIBMIN)
	- ( cd ../../bin ; ln -fs ../GWW/simple/$@ . )

libsimple.a : $(SIMPLEOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L libsimple.a

include make.depend
