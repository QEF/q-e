# Makefile for TDDFPT tools

include ../../make.sys
IFLAGS         = -I../../include
LIBOBJS        = ../../flib/ptools.a ../../flib/flib.a ../../clib/clib.a ../../iotk/src/libiotk.a ../../Multigrid/mglib.a


# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG).

PWOBJS = ../../PW/libpw.a

MODULES = \
../../Modules/basic_algebra_routines.o \
../../Modules/cell_base.o \
../../Modules/constants.o \
../../Modules/control_flags.o \
../../Modules/error_handler.o \
../../Modules/io_global.o \
../../Modules/io_files.o \
../../Modules/ions_base.o \
../../Modules/fft_types.o \
../../Modules/fft_scalar.o \
../../Modules/kind.o \
../../Modules/mp.o \
../../Modules/mp_base.o \
../../Modules/mp_global.o \
../../Modules/parallel_include.o \
../../Modules/parameters.o \
../../Modules/parser.o \
../../Modules/path_formats.o \
../../Modules/printout_base.o \
../../Modules/pseudo_types.o \
../../Modules/random_numbers.o \
../../Modules/radial_grids.o \
../../Modules/recvec.o \
../../Modules/splinelib.o \
../../Modules/timestep.o \
../../Modules/uspp.o




all : tddfpt_pp.x
	
tddfpt_pp.x : tddfpt_calculate_spectrum.o $(MODULES)
	$(LD) $(LDFLAGS) -o $@ \
	tddfpt_calculate_spectrum.o $(MODULES) $(PWOBJS) $(LIBOBJS) $(LIBS)
clean :
	rm -f tddfpt_pp.x *.o *.mod

