# Makefile for EE (electric potential with free boundary conditions)

include ../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../Modules $(MOD_FLAG).

EEOBJS = \
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

all : libee.a

libee.a: $(EEOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@
clean :
	- /bin/rm -f *.o *~ *.a *.d *.i *.L *.mod 

include make.depend

# DO NOT DELETE
