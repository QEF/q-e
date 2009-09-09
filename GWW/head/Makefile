# Makefile for head (head.x)
# 
# Author G. Stenuit
#
include ../../make.sys

# location of include files
IFLAGS=-I../../include

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../EE $(MOD_FLAG)../../PH $(MOD_FLAG)../../PW \
	  $(MOD_FLAG)../pw4gww $(MOD_FLAG).

#location of needed libraries
LIBOBJS= ../../iotk/src/libiotk.a ../../flib/flib.a \
         ../../clib/clib.a ../../flib/ptools.a \
         ../../Multigrid/mglib.a

HEADOBJS = \
bcast_ph_input.o \
close_phq.o \
openfilq.o \
phq_init.o \
phq_readin.o \
solve_head.o 

MODULES = \
../../Modules/atom.o \
../../Modules/autopilot.o \
../../Modules/basic_algebra_routines.o \
../../Modules/bfgs_module.o \
../../Modules/berry_phase.o \
../../Modules/cell_base.o \
../../Modules/check_stop.o \
../../Modules/clocks.o \
../../Modules/constants.o \
../../Modules/constraints_module.o \
../../Modules/control_flags.o \
../../Modules/descriptors.o \
../../Modules/dspev_drv.o \
../../Modules/electrons_base.o \
../../Modules/error_handler.o \
../../Modules/fft_base.o \
../../Modules/fft_parallel.o \
../../Modules/fft_scalar.o \
../../Modules/fft_types.o \
../../Modules/functionals.o \
../../Modules/input_parameters.o \
../../Modules/io_files.o \
../../Modules/io_global.o \
../../Modules/ions_base.o \
../../Modules/ions_nose.o \
../../Modules/kind.o \
../../Modules/metadyn_base.o \
../../Modules/metadyn_io.o \
../../Modules/metadyn_vars.o \
../../Modules/mm_dispersion.o \
../../Modules/mp_global.o \
../../Modules/mp_wave.o \
../../Modules/mp.o \
../../Modules/mp_base.o \
../../Modules/parallel_include.o \
../../Modules/parameters.o \
../../Modules/parser.o \
../../Modules/path_base.o \
../../Modules/path_formats.o \
../../Modules/path_variables.o \
../../Modules/path_opt_routines.o \
../../Modules/path_io_routines.o \
../../Modules/path_reparametrisation.o \
../../Modules/paw_variables.o \
../../Modules/printout_base.o \
../../Modules/pseudo_types.o \
../../Modules/ptoolkit.o \
../../Modules/radial_grids.o \
../../Modules/random_numbers.o \
../../Modules/read_cards.o \
../../Modules/read_namelists.o \
../../Modules/read_ncpp.o \
../../Modules/read_uspp.o \
../../Modules/read_upf_v1.o \
../../Modules/read_upf_v2.o \
../../Modules/recvec.o \
../../Modules/splinelib.o \
../../Modules/stick_base.o \
../../Modules/task_groups.o \
../../Modules/timestep.o \
../../Modules/upf_to_internal.o \
../../Modules/upf.o \
../../Modules/uspp.o \
../../Modules/version.o \
../../Modules/wavefunctions.o \
../../Modules/wave_base.o \
../../Modules/write_upf_v2.o \
../../Modules/xml_io_base.o \
../../Modules/zhpev_drv.o \
../../Modules/wannier_new.o \
../../Modules/wrappers.o \
../../Modules/ws_base.o 

EEMODS= ../../EE/ee_mod.o ../../EE/gcoarse_mod.o

LIBPWPH = ../pw4gww/libpw4gww.a ../../PH/libph.a ../../PW/libpw.a

LIBMIN= ../minpack/minpacklib.a

TLDEPS= bindir mods libs pw ph 

all : tldeps head.x  

head.x : head.o $(LIBOBJS) $(HEADOBJS) $(LIBPWPH) $(LIBMIN)
	$(LD) $(LDFLAGS) -o head.x head.o \
		 $(MODULES) $(EEMODS) $(HEADOBJS) $(LIBPWPH) $(LIBOBJS) $(LIBS) $(LIBMIN)
	- ( cd ../../bin ; ln -fs ../GWW/head/head.x . )

tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :
	- /bin/rm -f *.x *.o *.a *~ *.F90 *.d *.mod *.i *.L

include make.depend
