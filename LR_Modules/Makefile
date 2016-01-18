# Makefile for LR_Modules

include ../make.sys

# location of needed modules
MODFLAGS=$(MOD_FLAG)../FFTXlib $(MOD_FLAG)../Modules $(MOD_FLAG)../PW/src $(MOD_FLAG).

LRMODULES = cft_wave.o

TLDEPS=libfft mods

all : liblrmod.a

liblrmod.a: $(LRMODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *.F90 *.mod *.L

include make.depend
