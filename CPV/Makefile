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
compute_scf.o \
crayfft.o \
fft_cp.o \
input.o \
macdep.o \
modules.o \
para.o \
read_pseudo.o \
restart.o \
neb_routines.o \
stop_pw.o \
smd.o  \
redis.o \
smcp.o \
smlam.o \
smd_modules.o \
wf.o \
wfcpr.o \
$(WRAPPERS)

all : cp.x

cp.x : $(FOBJS)
	$(LD) -o cp.x $(FOBJS) ../Modules/*.o \
        ../flib/eispack.o $(LDFLAGS)
	chmod a+rx cp.x

clean :
	- /bin/rm -f *.o *.i core* fort* *.mod *.s *.d work.pc cp.x

include .dependencies
