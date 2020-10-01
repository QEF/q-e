# Makefile for FFTXlib

include ../make.inc

FFTX = \
scatter_mod.o  \
fft_scatter_2d.o  \
fft_scatter_gpu.o  \
fft_scatter_2d_gpu.o  \
fft_ggen.o  \
fft_fwinv.o  \
fft_scalar.o  \
fft_scalar.ARM_LIB.o  \
fft_scalar.DFTI.o  \
fft_scalar.ESSL.o  \
fft_scalar.FFTW.o  \
fft_scatter.o  \
fft_scalar.cuFFT.o  \
fftw_interfaces.o  \
fft_scalar.FFTW3.o  \
fft_scalar.SX6.o  \
fft_parallel.o  \
fft_parallel_2d.o  \
fft_interfaces.o  \
fft_interpolate.o \
stick_base.o  \
fftw.o  \
fftw_sp.o  \
fftw_dp.o  \
fft_smallbox.o  \
fft_smallbox_type.o  \
fft_support.o  \
fft_error.o  \
fft_stick.o  \
fft_types.o \
tg_gather.o \
fft_helper_subroutines.o \
fft_param.o \
fft_buffers.o 


all : libqefft.a

libqefft.a: 	$(FFTX)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

fft_scalar.o : fft_scalar.f90  fft_scalar.FFTW3.f90  fft_scalar.FFTW.f90  fft_scalar.SX6.f90 fft_scalar.DFTI.f90  fft_scalar.ESSL.f90


fft_stick.o : fft_stick.c fftw_sp.h fftw_dp.h

fftw.o : fftw.c fftw.h

fftw_sp.o : fftw_sp.c fftw_sp.h fftw.h konst.h

fftw_dp.o : fftw_dp.c fftw_dp.h fftw.h konst.h

test : fft_test.o libqefft.a
	$(LD) $(LDFLAGS) -o fft_test.x fft_test.o libqefft.a $(FFT_LIBS) $(BLAS_LIBS) $(MPI_LIBS) $(LD_LIBS)

TEST : test # compatibility

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L 

include make.depend
