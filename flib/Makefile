# Makefile for flib

include ../make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)../iotk/src $(MOD_FLAG)../Modules $(MOD_FLAG).

OBJS = \
avrec.o \
atomic_number.o \
bachel.o \
capital.o \
dost.o \
erf.o \
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

all: blas.a flib.a lapack_atlas.a lapack.a ptools.a

flib_only: flib.a ptools.a

flib.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

blas.a : blas.o
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@ 

lapack_atlas.a : lapack_atlas.o dlamch.o
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@ 

lapack.a : lapack_atlas.o lapack_all.o dlamch.o
	$(AR) $(ARFLAGS) $@ $?

dlamch.o : dlamch.f
	$(F77) $(FFLAGS_NOOPT) -c $<

ptools.a : $(POBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.a *.o *.mod *.i *.F90 core* *.L

include make.depend
