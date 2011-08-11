# Makefile for Modules

include ../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG).

MODULES = \
atom.o \
autopilot.o \
basic_algebra_routines.o \
bfgs_module.o \
cell_base.o  \
check_stop.o  \
clocks.o \
constants.o \
constraints_module.o \
control_flags.o \
coulomb_vcut.o \
descriptors.o \
dspev_drv.o \
electrons_base.o \
environment.o \
error_handler.o \
fft_base.o \
fft_interfaces.o \
fft_parallel.o \
fft_scalar.o \
fft_types.o \
funct.o \
griddim.o \
image_io_routines.o \
input_parameters.o \
io_files.o \
io_global.o  \
ions_base.o \
kernel_table.o \
kind.o \
mm_dispersion.o \
mp.o \
mp_base.o \
mp_global.o \
mp_image_global_module.o \
mp_wave.o \
open_close_input_file_interf.o \
open_close_input_file.o \
parallel_include.o \
parameters.o \
parser.o \
paw_variables.o \
plugin_flags.o \
plugin_arguments.o \
pseudo_types.o \
ptoolkit.o \
radial_grids.o \
random_numbers.o \
read_input.o \
read_namelists.o \
read_ncpp.o \
read_pseudo.o \
read_upf_v1.o \
read_upf_v2.o \
read_uspp.o \
read_xml.o \
read_xml_cards.o \
read_xml_fields.o \
recvec.o \
recvec_subs.o \
run_info.o \
set_signal.o \
sic.o \
splinelib.o \
stick_base.o \
stick_set.o \
timestep.o\
version.o \
upf.o \
upf_to_internal.o \
uspp.o \
wave_base.o \
wavefunctions.o \
write_upf_v2.o \
xc_vdW_DF.o \
xml_input.o \
xml_io_base.o \
zhpev_drv.o \
wannier_new.o \
wrappers.o\
ws_base.o \
read_cards.o

all : libqemod.a

libqemod.a: 	$(MODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *.F90 *.mod *.L

include make.depend
