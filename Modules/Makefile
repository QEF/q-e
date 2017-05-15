# Makefile for Modules

include ../make.inc

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../ELPA/src $(MOD_FLAG)../FFTXlib $(MOD_FLAG)../LAXlib $(IFLAGS) $(MOD_FLAG).

# list of modules

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
dist.o \
electrons_base.o \
environment.o \
error_handler.o \
fcp_variables.o \
fd_gradient.o \
fft_base.o \
fft_custom.o \
fft_rho.o \
fsockets.o \
funct.o \
generate_function.o \
gth.o \
gvecw.o \
input_parameters.o \
invmat.o \
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
mp_exx.o \
mp_global.o \
mp_images.o \
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
pw_dot.o \
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
space_group.o \
set_signal.o \
setqf.o \
splinelib.o \
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
io_base.o \
qes_module.o  \
qes_types.o \
qes_libs.o  \
qexsd.o   \
qexsd_input.o \
qexsd_reader_module.o \
hdf5_qe.o\
xsf.o \
wyckoff.o \
wypos.o \
zdotc_wrapper.o  \
zvscal.o 

# list of subrooutines and functions (not modules) previously found in flib/

OBJS = \
atomic_number.o \
capital.o \
cryst_to_car.o \
erf.o \
expint.o \
find_free_unit.o \
functionals.o \
generate_k_along_lines.o \
lsda_functionals.o \
more_functionals.o \
has_xml.o \
inpfile.o \
int_to_char.o \
latgen.o \
linpack.o \
metagga.o \
matches.o \
plot_io.o \
radial_gradients.o \
rgen.o \
recips.o \
remove_tot_torque.o \
set_hubbard_l.o \
set_hubbard_n.o \
simpsn.o \
sort.o \
sph_bes.o \
sph_dbes.o \
trimcheck.o \
test_input_file.o \
date_and_tim.o \
volume.o \
dylmr2.o \
ylmr2.o \
wgauss.o \
w0gauss.o \
w1gauss.o \
deviatoric.o

TLDEPS=libfft

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

## The following is needed only for lapack compiled from sources

dlamch.o : dlamch.f
	$(F77) $(FFLAGS_NOOPT) -c $<

libqemod.a: 	$(MODULES) $(OBJS)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi


clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L version.f90

# .PHONY forces execution of a rule irrespective of the presence of an
# updated file with the same name of the rule. In this way, the script 
# that generates version.f90 always runs, updating the version if you 
# execute "svn update". The update_version script takes care of not
# changing the file if the svn version did not change

.PHONY: version all clean

include make.depend
