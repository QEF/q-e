# Makefile for flib

include ../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../Modules $(MOD_FLAG).

OBJS = \
avrec.o \
atomic_number.o \
bachel.o \
capital.o \
cryst_to_car.o \
dost.o \
erf.o \
expint.o \
functionals.o \
lsda_functionals.o \
more_functionals.o \
iglocal.o \
inpfile.o \
int_to_char.o \
invmat.o \
invmat_complex.o \
latgen.o \
linpack.o \
metagga.o \
matches.o \
plot_io.o \
radial_gradients.o \
recips.o \
remove_tot_torque.o \
simpsn.o \
sort.o \
sph_bes.o \
sph_dbes.o \
transto.o \
date_and_tim.o \
sort_gvec.o \
volume.o \
dylmr2.o \
ylmr2.o

POBJS = \
flush_unit.o \
distools.o 

all: flib.a ptools.a 

flib_only: flib.a ptools.a

flib.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

dlamch.o : dlamch.f
	$(F77) $(FFLAGS_NOOPT) -c $<

ptools.a : $(POBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.a *.o *.mod *.i *.F90 core* *.L

include make.depend
