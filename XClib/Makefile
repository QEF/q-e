# Makefile for XClib

include ../make.inc

MODFLAGS= $(BASEMOD_FLAGS)

LDAL = \
kind_l.o \
constants_l.o \
xc_beef_interface.o \
correlation_lda_lsda.o \
correlation_gga.o \
dft_par_mod.o \
drivers_lda_lsda.o \
drivers_d_lda_lsda.o \
drivers_gga.o \
drivers_d_gga.o \
drivers_mgga.o \
exchange_lda_lsda.o \
exchange_gga.o \
metagga.o \
external.o \
get_xc.o \
xc_lda_lsda_wrapper.o \
xc_gga_wrapper.o \
xc_mgga_wrapper.o \
xc_interfaces.o

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
