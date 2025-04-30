# Makefile for PIOUD

include ../../make.inc

# location of needed modules and included files (if any)
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../PW/src \
          $(MOD_FLAG)../../dft-d3/

PIOUDOBJS = \
pioud.o \


PIOUDLIBS = \
compute_scf_pioud.o \
engine_to_trpmd_pos.o \
engine_to_trpmd_nat_alat.o \
fcp_opt_routines.o \
fcp_variables.o \
trpmd_base.o \
ring_formats.o \
pw_gen_inputs.o \
ring_input_parameters_module.o \
trpmd_io_routines.o \
ring_io_units_module.o \
pimd_variables.o \
ring_variables.o \
pimd_subrout.o \
pimd_utils.o \
set_defaults.o \
stop_run_path.o 

PWOBJS= ../../PW/src/libpw.a ../../KS_Solvers/libks_solvers.a ../../dft-d3/libdftd3qe.a
QEMODS = $(BASEMODS)

F90FLAGS += -pg
LDFLAGS += -pg

TLDEPS=pwlibs

TLDEPS=pwlibs

all : tldeps pioud.x

pioud.x : $(PIOUDOBJS) libpioud.a $(PWOBJS) $(QEMODS)
	$(LD) $(LDFLAGS) -o $@ \
	$(PIOUDOBJS) libpioud.a $(PWOBJS) $(QEMODS) $(QELIBS)
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
