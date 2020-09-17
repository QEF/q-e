# Makefile for XClib

include ../make.inc

MODFLAGS= $(BASEMOD_FLAGS)

LDAL = \
qe_constants.o \
qe_drivers_d_gga.o \
qe_drivers_d_lda_lsda.o \
qe_drivers_gga.o \
qe_drivers_lda_lsda.o \
qe_drivers_mgga.o \
qe_external_calls.o \
qe_funct_corr_gga.o \
qe_funct_corr_lda_lsda.o \
qe_funct_exch_gga.o \
qe_funct_exch_lda_lsda.o \
qe_funct_mgga.o \
qe_kind.o \
xc_beef_interface.o \
xc_get_dft.o \
xc_input_params_mod.o \
xc_interfaces.o \
xc_wrapper_d_gga.o \
xc_wrapper_d_lda_lsda.o \
xc_wrapper_gga.o \
xc_wrapper_lda_lsda.o \
xc_wrapper_mgga.o

all: xc_lib.a

xc_lib.a: 	$(LDAL)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L 

# .PHONY forces execution of a rule irrespective of the presence of an
# updated file with the same name of the rule. In this way, the script 
# that generates version.f90 always runs, updating the version if you 
# execute "svn update". The update_version script takes care of not
# changing the file if the svn version did not change

.PHONY: all clean

include make.depend
