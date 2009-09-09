# Makefile for pw4gww (PW4GWW)
# Author G. Stenuit and L. Martin-Samos

include ../../make.sys

# location of include files
IFLAGS=-I../../include

# location of needed modules
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules \
          $(MOD_FLAG)../../EE $(MOD_FLAG)../../PW $(MOD_FLAG).

#location of needed libraries
LIBOBJS= ../../iotk/src/libiotk.a ../../flib/flib.a \
         ../../clib/clib.a ../../flib/ptools.a \
         ../../Multigrid/mglib.a

PW4GWWOBJS = \
allocate_wannier.o \
calculate_wing.o \
cprim.o \
dft_exchange.o \
dirdel.o \
distance_gw.o \
energies_u.o \
energies_xc.o \
go_wannier_product.o \
gram_schmidt_pwannier.o \
local_wannier.o \
matrix_wannier_gamma.o \
memory_limit.o \
produce_wannier_gamma.o \
product_wannier.o \
product_wannier_para.o \
product_wannier_para_c.o \
openfil_pw4gww.o \
stop_pp.o \
real_wfc.o \
rotate_wannier.o \
set_wannier_P.o \
start_pw4gww.o \
ultra_external.o \
ultralocalization.o \
ultralocalization_para.o \
wannier.o \
wannier_gw.o \
wannier_pmat_terms_ggrid.o \
wannier_uterms.o \
wannier_uterms_c.o \
wannier_valence_terms.o \
wfc_real.o \
write_vpot_matrix.o \
write_wannier_matrix.o \
write_wannier_matrix_c.o \
write_wannier_products.o 

MODULES = \
../../Modules/atom.o \
../../Modules/autopilot.o \
../../Modules/basic_algebra_routines.o \
../../Modules/berry_phase.o \
../../Modules/bfgs_module.o \
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
../../Modules/mm_dispersion.o \
../../Modules/mp_global.o \
../../Modules/mp_wave.o \
../../Modules/mp.o \
../../Modules/mp_base.o \
../../Modules/metadyn_base.o \
../../Modules/metadyn_io.o \
../../Modules/metadyn_vars.o \
../../Modules/path_base.o \
../../Modules/path_formats.o \
../../Modules/path_variables.o \
../../Modules/path_opt_routines.o \
../../Modules/path_io_routines.o \
../../Modules/path_reparametrisation.o \
../../Modules/parallel_include.o \
../../Modules/parameters.o \
../../Modules/parser.o \
../../Modules/paw_variables.o \
../../Modules/printout_base.o \
../../Modules/pseudo_types.o \
../../Modules/ptoolkit.o \
../../Modules/radial_grids.o \
../../Modules/random_numbers.o \
../../Modules/read_cards.o \
../../Modules/read_namelists.o \
../../Modules/read_ncpp.o \
../../Modules/read_upf_v1.o \
../../Modules/read_upf_v2.o \
../../Modules/read_uspp.o \
../../Modules/recvec.o \
../../Modules/splinelib.o \
../../Modules/stick_base.o \
../../Modules/task_groups.o \
../../Modules/timestep.o \
../../Modules/upf_to_internal.o \
../../Modules/uspp.o \
../../Modules/upf.o \
../../Modules/version.o \
../../Modules/wannier_new.o \
../../Modules/wavefunctions.o \
../../Modules/wave_base.o \
../../Modules/wrappers.o \
../../Modules/write_upf_v2.o \
../../Modules/xml_io_base.o \
../../Modules/zhpev_drv.o \
../../Modules/ws_base.o 

EEMODS=../../EE/ee_mod.o ../../EE/gcoarse_mod.o

PWOBJS = ../../PW/libpw.a

LIBMIN= ../minpack/minpacklib.a

TLDEPS= bindir libs pw

all : tldeps pw4gww.x

pw4gww.x : pw4gww.o libpw4gww.a $(PW4GWWOBJS) $(PWOBJS) $(MODULES) $(EEMODS) $(LIBMIN) 
	$(LD) $(LDFLAGS) -o $@ \
		pw4gww.o $(MODULES) $(EEMODS) libpw4gww.a  $(PWOBJS) $(LIBOBJS) $(LIBS) $(LIBMIN)
	- ( cd ../../bin ; ln -fs ../GWW/pw4gww/$@ . )
 
tldeps:
	test -n "$(TLDEPS)" && ( cd ../.. ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

libpw4gww.a : $(PW4GWWOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *~ *.F90 *.d *.mod *.i *.L libpw4gww.a

include make.depend
