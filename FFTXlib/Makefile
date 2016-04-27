# Makefile for FFTXlib

include ../make.sys

# location of needed modules
#MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG).

FFTX = \
scatter_mod.o  \
fft_scalar.o  \
fft_parallel.o  \
fft_interfaces.o  \
stick_base.o  \
stick_set.o  \
fft_smallbox.o  \
fft_support.o  \
fft_error.o  \
fft_stick.o  \
fft_types.o 


all : libqefft.a

libqefft.a: 	$(FFTX)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

fft_scalar.o : fft_scalar.f90  fft_scalar.FFTW3.f90  fft_scalar.FFTW.f90  fft_scalar.SX6.f90 fft_scalar.DFTI.f90  fft_scalar.ESSL.f90


fft_stick.o : fft_stick.c fftw.c fftw.h konst.h

TEST : test.o libqefft.a
	$(LD) $(LDFLAGS) -o fft_test.x test.o libqefft.a $(LIBS)

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L 

# .PHONY forces execution of a rule irrespective of the presence of an
# updated file with the same name of the rule. In this way, the script 
# that generates version.f90 always runs, updating the version if you 
# execute "svn update". The update_version script takes care of not
# changing the file if the svn version did not change

.PHONY: all clean

include make.depend
