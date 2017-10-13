# Makefile for dft-d3

include ../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../ELPA/src $(MOD_FLAG)../FFTXlib $(MOD_FLAG)../LAXlib $(MOD_FLAG)../UtilXlib \
          $(MOD_FLAG)../KS_Solvers/CG $(IFLAGS) $(MOD_FLAG)../Modules $(MOD_FLAG).
# list of modules


TESTOBJ = \
test_code.o

QEMODS= ../Modules/libqemod.a ../UtilXlib/libutil.a


TESTMODS = \
libdftd3qe.a


DFTD3OBJ = \
dftd3-lib-0.9/lib/common.o \
dftd3-lib-0.9/lib/sizes.o \
dftd3-lib-0.9/lib/pars.o \
dftd3-lib-0.9/lib/core.o \
dftd3-lib-0.9/lib/api.o \
dftd3_qe.o


TESTD3 = \
common.o \
sizes.o \
pars.o \
core.o \
api.o \
dftd3_qe.o 


TLDEPS=mods

all : linkstodftd3sources libdftd3qe.a test.x


linkstodftd3sources:
	if ((test -d dftd3-lib-0.9/lib) && (test ! -f linkstodftd3sources)); then \
	ln -fs dftd3-lib-0.9/lib/common.f90 . ; \
	ln -fs dftd3-lib-0.9/lib/core.f90 . ; \
	ln -fs dftd3-lib-0.9/lib/sizes.f90 . ; \
	ln -fs dftd3-lib-0.9/lib/pars.f90 . ; \
	ln -fs dftd3-lib-0.9/lib/api.f90 . ; \
	touch ./linkstodftd3sources ; fi


libdftd3qe.a: $(DFTD3OBJ)
	$(AR) $(ARFLAGS) $@ $(TESTD3)        
	$(RANLIB) $@    

test.x : $(TESTOBJ) $(TESTMODS) $(QEMODS)
	$(LD) $(LDFLAGS) $(MODFLAGS) -o $@ \
        $(TESTOBJ) $(TESTMODS) $(QEMODS) $(LIBS)	
      

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../; $(MAKE) $(TLDEPS) || exit 1 ) ; fi


clean :
	( /bin/rm -f *.o *.a *.d *.i *.x *~ *_tmp.f90 *.mod *.L ); \
	rm -f common.f90
	rm -f core.f90
	rm -f sizes.f90
	rm -f pars.f90
	rm -f api.f90
	rm -f linkstodftd3sources

touch-foo:
	$(dummy-variable)

include make.depend.dftd3
