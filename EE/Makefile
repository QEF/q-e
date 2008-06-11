# Makefile for Modules

include ../make.sys

EE = \
ee_mod.o                             \
gcoarse_mod.o                        \
dvdr_tao.o                           \
add_dccdil_forces.o                  \
v_h_from_rho_r.o                     \
init_ee.o                            \
calc_ecomp.o                         \
setlocalcoul.o                       \
set_fft_dim_coarse.o                 \
set_mltgrid_dim.o                    \
add_boundary.o                       \
add_dcc_field.o                      \
mg_pb_solver.o                       \
data_structure_coarse.o              \
multiscale.o                         \
writetofile.o                        \
add_ele_corr.o                       \
write_ee_summary.o                   \
ggen_coarse.o                        

MODULES = \
../Modules/fft_base.o \
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
../Modules/printout_base.o \
../Modules/pseudo_types.o \
../Modules/recvec.o \
../Modules/shmem_include.o \
../Modules/stick_base.o \
../Modules/timestep.o \
../Modules/uspp.o \
../Modules/version.o \
../Modules/wavefunctions.o \
../Modules/xml_io_base.o

all : $(MODULES) $(PWOBJS) $(EE) 

clean :
	- /bin/rm -f *.o *.d *.i *~ *.o *.mod work.pc

include make.depend

# DO NOT DELETE
