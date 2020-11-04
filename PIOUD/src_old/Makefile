# Makefile for PIOUD

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../dft-d3/

PIOUDOBJS = \
pioud.o \

PIOUDLIBS = \
compute_scf.o \
engine_to_path_pos.o \
engine_to_path_alat.o \
engine_to_path_nat.o \
engine_to_path_fix_atom_pos.o \
fcp_opt_routines.o \
input.o \
path_base.o \
path_formats.o \
path_gen_inputs.o \
path_input_parameters_module.o \
path_io_routines.o \
path_io_tools.o \
path_io_units_module.o \
path_opt_routines.o \
path_reparametrisation.o \
path_read_cards_module.o \
path_read_namelists_module.o \
path_to_engine_fix_atom_pos.o \
path_variables.o \
pimd_variables.o \
pimd_subrout.o \
pimd_utils.o \
set_defaults.o \
stop_run_path.o 

QEMODS=../../Modules/libqemod.a ../../KS_Solvers/libks_solvers.a \
       ../../FFTXlib/libqefft.a ../../LAXlib/libqela.a ../../UtilXlib/libutil.a ../../dft-d3/libdftd3qe.a
PWOBJS= ../../PW/src/libpw.a

TLDEPS=pwlibs

all : tldeps pioud.x

pioud.x : $(PIOUDOBJS) libpioud.a $(LIBOBJS) $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	$(PIOUDOBJS) libpioud.a $(PWOBJS) $(QEMODS) $(LIBOBJS) $(QELIBS)
	- ( cd ../../bin; ln -fs ../PIOUD/src/$@ . )

libpioud.a : $(PIOUDLIBS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L
	- /bin/rm -f ../../bin/pioud.x

include make.depend
