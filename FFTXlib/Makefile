# Makefile for FFTXLIB
sinclude ../make.inc

default: all

all: fftx

fftx:
	( cd src ; $(MAKE) all || exit 1 )

fftx_clean:
	( cd src ; $(MAKE) clean )

clean: fftx_clean

distclean: clean
