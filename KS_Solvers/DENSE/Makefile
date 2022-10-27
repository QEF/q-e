# Makefile for DENSE diagonalization

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG) ../../LAXlib $(MOD_FLAG) ../../UtilXlib $(MOD_FLAG) ../../Modules $(MOD_FLAG).

DENSE = \
rotate_HSpsi_gamma.o \
rotate_HSpsi_k.o \
rotate_wfc_gamma.o \
rotate_wfc_k.o \
rotate_driver.o \
rotate_wfc_gamma_gpu.o \
rotate_wfc_k_gpu.o \
rotate_driver_cuf.o \
rotate_xpsi_k.o \
rotate_xpsi_k_gpu.o \
rotate_xpsi_gamma.o \
rotate_xpsi_gamma_gpu.o \
gram_schmidt_gamma.o \
gram_schmidt_k.o \
gram_schmidt_gamma_gpu.o \
gram_schmidt_k_gpu.o  

all : libdense.a


libdense.a:  $(DENSE)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

include make.depend
