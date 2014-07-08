# Makefile for Modules

include ../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../ELPA/src $(MOD_FLAG).

MODULES = \
atom.o \
autopilot.o \
basic_algebra_routines.o \
becmod.o \
bfgs_module.o \
bspline.o \
bz_form.o \
cell_base.o  \
check_stop.o  \
clocks.o \
command_line_options.o \
compute_dipole.o \
constants.o \
constraints_module.o \
control_flags.o \
coulomb_vcut.o \
descriptors.o \
dspev_drv.o \
electrons_base.o \
environment.o \
error_handler.o \
fd_gradient.o \
fft_base.o \
fft_custom.o \
fft_interfaces.o \
fft_parallel.o \
fft_scalar.o \
fft_types.o \
funct.o \
generate_function.o \
griddim.o \
input_parameters.o \
io_files.o \
io_global.o  \
ions_base.o \
kernel_table.o \
kind.o \
libxc.o \
libxc_funcs.o \
mm_dispersion.o \
mp.o \
mp_bands.o \
mp_base.o \
mp_diag.o \
mp_global.o \
mp_images.o \
mp_pots.o \
mp_pools.o \
mp_wave.o \
mp_world.o \
noncol.o \
open_close_input_file.o \
parallel_include.o \
parameters.o \
parser.o \
paw_variables.o \
plugin_flags.o \
plugin_arguments.o \
plugin_variables.o \
pseudo_types.o \
ptoolkit.o \
qexml.o \
qmmm.o \
radial_grids.o \
random_numbers.o \
read_cards.o \
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
tsvdw.o\
upf.o \
upf_to_internal.o \
uspp.o \
version.o \
wannier_gw.o\
wannier_new.o \
wave_base.o \
wavefunctions.o \
wrappers.o\
write_upf_v2.o \
ws_base.o \
xc_vdW_DF.o \
xc_rVV10.o \
xml_input.o \
xml_io_base.o \
zdotc_wrapper.o \
zhpev_drv.o 


all : version version.o libqemod.a

## If no update_version is found, or it is not executable,
## just copy version.f90.in (that contains 'version=unknown').
## We use cat and not cp so that version.f90 has a new timestamp,
## forcing recompilation of version.f90.
## We use the inner if to avoid recompilation at each make.
version:
	- ( if test -x ../install/update_version ; then \
	../install/update_version; \
	else if test ! -f version.f90 ; then \
	cat version.f90.in > version.f90 ; fi ; fi )	

version.f90: version

libqemod.a: 	$(MODULES)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *.F90 *.mod *.L version.f90

# .PHONY forces execution of a rule irrespective of the presence of an
# updated file with the same name of the rule. In this way, the script 
# that generates version.f90 always runs, updating the version if you 
# execute "svn update". The update_version script takes care of not
# changing the file if the svn version did not change

.PHONY: version all clean

include make.depend
