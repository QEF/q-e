#
# Makefile for Modules
#
MODULES= constants.o \
kind.o \
parallel_include.o \
shmem_include.o \
parameters.o \
io_base.o  \
io_global.o  \
mp.o  \
mp_global.o  \
mp_wave.o  \
mp_buffers.o  \
parser.o \
readpseudo.o  \
pseudodata.o  \
cell_base.o  \
pseudo_types.o \
ptoolkit.o \
stick_base.o \
fft_types.o \
fft_base.o \
fft_scalar.o
#
include ../make.rules
#
include ../make.sys
#
# targets
#
all: $(MODULES)
#
clean_:
	-/bin/rm -f *.o *.d *~ *.F90 *.mod work.pc
#
include .dependencies

