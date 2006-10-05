# Makefile for flib

include ../make.sys

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
iceil.o \
iglocal.o \
inpfile.o \
int_to_char.o \
invmat.o \
invmat_complex.o \
latgen.o \
linpack.o \
matches.o \
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
blockset.o \
flush_unit.o \
gridsetup.o \
hangup.o \
ilcm.o \
indxg2l.o \
indxg2p.o \
indxl2g.o \
infog1l.o \
infog2l.o \
localdim.o \
localindex.o \
npreroc.o \
numroc.o \
ownerof.o \
error.o

all: blas.a flib.a ptools.a $(MYLIB)

flib.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

blas.a : blas.o dlamch.o
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@ 

lapack : $(OBJS) lapack.o dlamch.o
	$(AR) $(ARFLAGS) flib.a $?
	$(RANLIB) flib.a
	touch $@

lapack_ibm : $(OBJS) lapack_ibm.o dlamch.o
	$(AR) $(ARFLAGS) flib.a $?
	$(RANLIB) flib.a
	touch $@

lapack_mkl : $(OBJS) lapack_mkl.o 
	$(AR) $(ARFLAGS) flib.a $?
	touch $@

lapack_atlas : $(OBJS) lapack_atlas.o dlamch.o
	$(AR) $(ARFLAGS) flib.a $?
	touch $@

dlamch.o : dlamch.f
	$(F77) $(FFLAGS_NOOPT) -c $<

ptools.a : $(POBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f ptools.a *.o *.mod *.i *.F90 core* work.pc flib.a \
		blas_and_lapack blas lapack lapack_ibm lapack_t3e \
		lapack_mkl lapack_atlas

include make.depend
