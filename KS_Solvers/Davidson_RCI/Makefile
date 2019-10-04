# Makefile for DAVID_RCI

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG) ../../ELPA/src $(MOD_FLAG) ../../LAXlib $(MOD_FLAG) ../../UtilXlib $(MOD_FLAG).

DAVID_RCI = david_rci.o \
david_param.o

all : libdavid_rci.a


libdavid_rci.a:  $(DAVID_RCI)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

include make.depend
