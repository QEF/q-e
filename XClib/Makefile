# Makefile for XClib

include ../make.inc

MODFLAGS= $(BASEMOD_FLAGS)

LDAL = \
kind_l.o \
constants_l.o \
correlation_lda_lsda.o \
correlation_gga.o \
dft_par_mod.o \
drivers_lda_lsda.o \
drivers_d_lda_lsda.o \
drivers_gga.o \
exchange_lda_lsda.o \
exchange_gga.o \
external.o \
get_xc.o \
xc_interfaces.o

all: libqeldaxc.a

libqeldaxc.a: 	$(LDAL)
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
