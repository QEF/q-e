# Makefile for FFTXlib testing

include ../../make.inc

MODFLAGS = $(MOD_FLAG).. $(MOD_FLAG).
LDFLAGS += -Mcudalib=cufft

SRCS = test_fft_scalar_gpu.f90 \
       test_fft_scatter_mod_gpu.f90 \
       test_fwinv_gpu.f90

EXECS = $(SRCS:.f90=.x)

all: common $(EXECS)

common: sort.o tester.o utils.o recips.o

%.x: %.o
	$(LD) $(LDFLAGS) $< utils.o tester.o sort.o recips.o -o $@ ../libqefft.a $(FFT_LIBS) $(BLAS_LIBS) $(MPI_LIBS) $(LD_LIBS)

clean:
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x rnd_seed*

.PHONY: all clean

