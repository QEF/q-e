# Makefile for ParO

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG) ../../LAXlib $(MOD_FLAG) ../../UtilXlib $(MOD_FLAG).

PARO = \
paro_gamma.o paro_gamma_new.o \
paro_k.o paro_k_new.o \
pcg_gamma.o \
pcg_k.o \
bpcg_gamma.o \
bpcg_k.o 

all : libparo.a

libparo.a:  $(PARO)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

include make.depend
