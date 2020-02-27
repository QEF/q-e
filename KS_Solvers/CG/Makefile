# Makefile for CG

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG) ../../ELPA/src $(MOD_FLAG) ../../LAXlib $(MOD_FLAG) ../../UtilXlib $(MOD_FLAG).

CG = \
ccgdiagg.o \
rcgdiagg.o 

-include make.gpu

all : libcg.a


libcg.a:  $(CG)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

include make.depend
