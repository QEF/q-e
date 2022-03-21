# Makefile for FFTXlib testing

include ../../make.inc

MODFLAGS = $(MOD_FLAG)../src $(MOD_FLAG).

SRCS = test_fft_scalar_gpu.f90 \
       test_fft_scatter_mod_gpu.f90 \
       test_fwinv_gpu.f90

EXECS = $(SRCS:.f90=.x)

all: common $(EXECS)

common: sort.o tester.o utils.o recips.o

%.x: %.o
	$(LD) $(LDFLAGS) -Mcudalib=cufft $< utils.o tester.o sort.o recips.o -o $@ ../src/libqefft.a $(FFT_LIBS) $(BLAS_LIBS) $(MPI_LIBS) $(LD_LIBS)

test : fft_test.o ../src/libqefft.a
	$(LD) $(LDFLAGS) -o fft_test.x fft_test.o ../src/libqefft.a $(FFT_LIBS) $(BLAS_LIBS) $(MPI_LIBS) $(LD_LIBS)

TEST : test # compatibility

clean:
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x rnd_seed*

.PHONY: all clean

