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
cp_emass.o \
cp_restart.o \
cp_fpmd.o \
descriptors.o \
electrons_base.o \
energies.o \
fft_base.o \
fft_scalar.o \
fft_types.o \
functionals.o \
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
path_base.o \
path_formats.o \
path_io_routines.o  \
path_opt_routines.o \
path_variables.o \
parallel_include.o \
parallel_types.o \
parameters.o \
parser.o \
pseudodata.o \
pseudo_types.o \
ptoolkit.o \
printout_base.o \
read_cards.o \
read_namelists.o \
readpseudo.o \
recvec.o \
shmem_include.o \
smallbox.o \
splinelib.o \
stick_base.o \
supercell.o \
timestep.o \
version.o \
uspp.o \
wave_base.o \
wavefunctions.o

IOTK=iotk_module.o

all : $(MODULES)

$(IOTK): ../include/iotk_config.h

iotk : $(IOTK)

clean :
	- /bin/rm -f *.o *.d *.i *~ *.F90 *.mod work.pc

include make.depend
