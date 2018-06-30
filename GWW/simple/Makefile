# Makefile for  simple

include ../../make.inc

# location of include files
IFLAGS=-I../../../include

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../EE $(MOD_FLAG)../../PW/src $(MOD_FLAG)../pw4gww $(MOD_FLAG)../gww \
	$(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../LAXlib $(MOD_FLAG)../../UtilXlib $(MOD_FLAG).

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


QEMODS = ../../Modules/libqemod.a ../../FFTXlib/libqefft.a ../../KS_Solvers/CG/libcg.a ../../KS_Solvers/Davidson/libdavid.a \
	../../LAXlib/libqela.a ../../UtilXlib/libutil.a

PWOBJS = ../../PW/src/libpw.a
PW4GWWOBJS = ../pw4gww/libpw4gww.a
GWWOBJ = ../gww/libgww.a

TLDEPS= bindir libs pw

all : tldeps simple.x

simple.x : simple.o libsimple.a  $(SIMPLEOBJS) $(PWOBJS) $(QEMODS) $(PW4GWWOBJS) $(GWWOBJ)
	$(LD) $(LDFLAGS) -o $@ \
		simple.o libsimple.a $(PWOBJS) $(PW4GWWOBJS) $(GWWOBJ) $(EEOBJS) $(QEMODS) $(LIBOBJS) $(LIBS) $(LIBMIN)
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

