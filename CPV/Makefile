# Makefile for CP

include	../make.sys
include ../make.rules

# CODEDIR = CPV

WRAPPERS = wrapper.o

FOBJS = \
cplib.o \
cpr.o \
cpr_mod.o \
cprstart.o \
cprsub.o \
crayfft.o \
errore.o \
fft_cp.o \
input.o \
macdep.o \
modules.o \
para.o \
read_pseudo.o \
restart.o \
which_dft.o \
neb_routines.o \
$(WRAPPERS)

all : cp.x

cp.x : $(FOBJS) cpflush.o
	$(LD) -o cp.x $(FOBJS) ../Modules/*.o cpflush.o \
          ../flib/flib.a ../flib/eispack.o $(LDFLAGS)
	chmod a+rx cp.x

clean :
	- /bin/rm -f *.o *.i core* fort* *.mod *.s *.d work.pc cp.x

include .dependencies
