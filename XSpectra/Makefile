# Makefile for postprocessing (PP)

#include ../make.rules
include ../make.sys

XOBJS = \
./xspectra_mod.o \
./radin_mod.o \
./mygetK.o \
./read_file_xspectra.o \
./ipoolscatter.o \
./init_xspectra_ldau.o

PPOBJS = \
../PP/openfil_pp.o\
../PP/start_postproc.o \
../PP/stop_pp.o

MODULES = \
../Modules/atom.o \
../Modules/cell_base.o \
../Modules/check_stop.o \
../Modules/clocks.o \
../Modules/constants.o \
../Modules/control_flags.o \
../Modules/fft_base.o \
../Modules/fft_scalar.o \
../Modules/fft_types.o \
../Modules/functionals.o \
../Modules/input_parameters.o \
../Modules/io_files.o \
../Modules/io_global.o \
../Modules/ions_base.o \
../Modules/kind.o \
../Modules/mp_global.o \
../Modules/mp_wave.o \
../Modules/mp.o \
../Modules/parallel_include.o \
../Modules/parameters.o \
../Modules/parser.o \
../Modules/ptoolkit.o \
../Modules/printout_base.o \
../Modules/pseudo_types.o \
../Modules/random_numbers.o \
../Modules/read_upf_v1.o \
../Modules/read_upf_v2.o \
../Modules/read_uspp.o \
../Modules/recvec.o \
../Modules/splinelib.o \
../Modules/shmem_include.o \
../Modules/stick_base.o \
../Modules/task_groups.o \
../Modules/timestep.o \
../Modules/uspp.o \
../Modules/version.o \
../Modules/wavefunctions.o \
../Modules/xml_io_base.o \
../Modules/error_handler.o \
../Modules/paw_variables.o \
../Modules/mp_base.o \
../Modules/basic_algebra_routines.o \
../Modules/fft_parallel.o\
../Modules/constraints_module.o \
../Modules/vxcgc.o \
../Modules/radial_grids.o \
../Modules/upf_to_internal.o \
../Modules/zhpev_drv.o \
../Modules/vxc_t.o \
../Modules/exc_t.o \
../Modules/electrons_base.o \
../Modules/metagga.o \
../Modules/descriptors.o \
../Modules/dspev_drv.o \
../Modules/upf.o \
../Modules/read_ncpp.o \
../Modules/write_upf_v2.o

PWOBJS = ../PW/libpw.a \
../GIPAW/paw_gipaw.o \
../GIPAW/init_gipaw_1.o \
../GIPAW/init_gipaw_2.o \
../PW/startup.o

EEMODS=../EE/ee_mod.o ../EE/gcoarse_mod.o

all : xspectra.x 

xspectra.x :   xspectra.o $(XOBJS) 
	$(MPIF90) -o $@  xspectra.o $(XOBJS) $(PWOBJS) $(PPOBJS) $(MODULES) $(EEMODS) $(LDFLAGS) $(LIBOBJS) $(LIBS)
	- ( cd ../bin ; ln -fs ../XSpectra/$@ . )

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i work.pc

include make.depend
# DO NOT DELETE
