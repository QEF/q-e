# Makefile for Modules

include ../make.rules
include ../make.sys

MODULES = \
atom.o \
berry_phase.o \
bfgs_module.o \
cell_base.o  \
check_stop.o  \
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
functionals.o \
griddim.o \
input_parameters.o \
io_base.o \
io_files.o \
io_global.o  \
io_routines.o  \
ions_base.o \
kind.o \
minimization_routines.o \
mp.o \
mp_global.o \
mp_wave.o \
mp_buffers.o \
neb_base.o \
neb_variables.o \
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
splinelib.o \
stick_base.o \
supercell.o \
version.o \
uspp.o \
wave_base.o \
wavefunctions.o

all : $(MODULES)

clean :
	- /bin/rm -f *.o *.d *~ *.F90 *.mod work.pc

include .dependencies
