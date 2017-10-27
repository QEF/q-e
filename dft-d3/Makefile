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
common.o \
sizes.o \
pars.o \
core.o \
api.o \
dftd3_qe.o 


TLDEPS=mods

all : tldeps libdftd3qe.a test.x


libdftd3qe.a: $(DFTD3OBJ)
	$(AR) $(ARFLAGS) $@ $(DFTD3OBJ)        
	$(RANLIB) $@    

test.x : $(TESTOBJ) $(TESTMODS) $(QEMODS)
	$(LD) $(LDFLAGS) $(MODFLAGS) -o $@ \
        $(TESTOBJ) $(TESTMODS) $(QEMODS) $(LIBS)	
      

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../; $(MAKE) $(TLDEPS) || exit 1 ) ; fi


clean :
	( /bin/rm -f *.o *.a *.d *.i *.x *~ *_tmp.f90 *.mod *.L ); \

include make.depend.dftd3
