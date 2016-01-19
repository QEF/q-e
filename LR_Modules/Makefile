# Makefile for LR_Modules

include ../make.sys

# location of needed modules and included files (if any)
MODFLAGS=$(MOD_FLAG)../FFTXlib $(MOD_FLAG)../Modules $(MOD_FLAG)../PW/src $(MOD_FLAG).

IFLAGS=

LR_MODULES = \
apply_dpot.o \
cft_wave.o \
h_psiq.o \
lrcom.o

TLDEPS=libfft mods

all : liblrmod.a

liblrmod.a: $(LR_MODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *.F90 *.mod *.L

include make.depend
