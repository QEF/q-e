# Makefile for Direct diagonalization

include ../../make.inc

VPATH = @srcdir@

# location of needed modules and included files (if any)
MODFLAGS= $(MOD_FLAG) ../../LAXlib $(MOD_FLAG) ../../UtilXlib $(MOD_FLAG) ../../Modules $(MOD_FLAG).

DIRECT = \
direct_diag_k.o

all : libdirect.a


libdirect.a:  $(DIRECT)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x


include make.depend


