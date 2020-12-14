# Makefile for DENSE diagonalization

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG) ../../ELPA/src $(MOD_FLAG) ../../LAXlib $(MOD_FLAG) ../../UtilXlib $(MOD_FLAG).

DENSE = \
rotate_HSpsi_gamma.o \
rotate_HSpsi_k.o \
rotate_wfc_gamma.o \
rotate_wfc_k.o 

-include make.gpu

all : libdense.a


libdense.a:  $(DENSE)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

include make.depend
