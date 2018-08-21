# Makefile for  simple

include ../../make.inc
LIBMIN = ../minpack/minpacklib.a

# location of needed modules
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src $(MOD_FLAG)../pw4gww $(MOD_FLAG)../gww 

#location of needed libraries
LIBOBJS= ../../iotk/src/libiotk.a \
         ../../clib/clib.a

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
  init_us_2_max.o \
  commutator.o 


QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a \
        ../../KS_Solvers/libks_solvers.a \
	../../LAXlib/libqela.a ../../UtilXlib/libutil.a

PWOBJS = ../../PW/src/libpw.a
PW4GWWOBJS = ../pw4gww/libpw4gww.a
GWWOBJ = ../gww/libgww.a

TLDEPS= pwlibs gwwlib

all : tldeps simple.x

simple.x : simple.o libsimple.a  $(SIMPLEOBJS) $(PWOBJS) $(QEMODS) $(PW4GWWOBJS) $(GWWOBJ)
	$(LD) $(LDFLAGS) -o $@ \
		simple.o libsimple.a $(PWOBJS) $(PW4GWWOBJS) $(GWWOBJ) $(EEOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS) $(LIBMIN)
	- ( cd ../../bin ; ln -fs ../GWW/simple/$@ . )

tldeps :
	if test -n "$(TLDEPS)" ; then   ( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi


libsimple.a : $(SIMPLEOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L libsimple.a

include make.depend
# DO NOT DELETE

