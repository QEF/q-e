# Makefile for Modules

include ../make.rules
include ../make.sys

MODULES = \
berry_phase.o \
bfgs_module.o \
cell_base.o  \
clocks.o \
constants.o \
control_flags.o \
descriptors.o \
electrons_base.o \
energies.o \
fft_base.o \
fft_scalar.o \
fft_types.o \
formats.o \
griddim.o \
input_parameters.o \
io_base.o \
io_files.o \
io_global.o  \
ions_base.o \
kind.o \
mp.o \
mp_global.o \
mp_wave.o \
mp_buffers.o \
parallel_include.o \
parallel_types.o \
parameters.o \
parser.o \
pseudodata.o \
pseudo_types.o \
ptoolkit.o \
read_cards.o \
read_namelists.o \
readpseudo.o \
recvec.o \
shmem_include.o \
smallbox.o \
stick_base.o \
version.o \
wavefunctions.o

all : $(MODULES)

clean :
	- /bin/rm -f *.o *.d *~ *.F90 *.mod work.pc

include .dependencies
