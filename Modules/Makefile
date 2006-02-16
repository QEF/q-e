# Makefile for Modules

include ../make.sys

MODULES = \
atom.o \
autopilot.o \
berry_phase.o \
bfgs_module.o \
cell_base.o  \
check_stop.o  \
clocks.o \
constants.o \
constraints_module.o \
control_flags.o \
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
ions_nose.o \
kind.o \
mp.o \
mp_global.o \
mp_wave.o \
metadyn_base.o \
metadyn_io.o \
metadyn_vars.o \
path_base.o \
path_formats.o \
path_io_routines.o  \
path_opt_routines.o \
path_reparametrisation.o \
path_variables.o \
parallel_include.o \
parallel_types.o \
parameters.o \
parser.o \
pseudodata.o \
pseudo_types.o \
ptoolkit.o \
printout_base.o \
random_numbers.o \
read_cards.o \
read_namelists.o \
readpseudo.o \
recvec.o \
shmem_include.o \
sic.o \
smallbox.o \
splinelib.o \
stick_base.o \
timestep.o \
version.o \
uspp.o \
wave_base.o \
wavefunctions.o \
xml_io_base.o

all : $(MODULES)

clean :
	- /bin/rm -f *.o *.d *.i *~ *.F90 *.mod work.pc

include make.depend
