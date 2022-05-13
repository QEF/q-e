# Makefile for PPCG

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG) ../../LAXlib $(MOD_FLAG) ../../UtilXlib $(MOD_FLAG).

PPCG = \
generic_cublas.o \
ppcg_gamma.o \
ppcg_gamma_gpu.o \
ppcg_k.o \
ppcg_k_gpu.o \

all : libppcg.a


libppcg.a:  $(PPCG)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

include make.depend
