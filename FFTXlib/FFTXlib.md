project: FFTXlib
project_dir: .
output_dir: ./Doc
predocmark: >
docmark_alt: #
predocmark_alt: <
display: public
         private
exclude: fft_scalar.DFTI.f90
         fft_scalar.ESSL.f90
         fft_scalar.FFTW3.f90
         fft_scalar.FFTW.f90
         fft_scalar.SX6.f90
         mpif.h
include: /cineca/prod/compilers/intel/cs-xe-2015/binary/impi_5.0.2/include64/
graph: true

A self-contained library for handling FFT in QE
