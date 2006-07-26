# Makefile for flib

include ../make.sys

OBJS = \
avrec.o \
atomic_number.o \
bachel.o \
capital.o \
dost.o \
erf.o \
flush_unit.o \
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
scnds.o \
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

all: flib.a ptools.a $(MYLIB)

flib.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

blas_and_lapack : blas lapack
	touch $@

blas : blas.o dlamch.o
	$(AR) $(ARFLAGS) flib.a $?
	$(RANLIB) flib.a
	touch $@

lapack : lapack.o dlamch.o
	$(AR) $(ARFLAGS) flib.a $?
	$(RANLIB) flib.a
	touch $@

lapack_ibm : lapack_ibm.o dlamch.o
	$(AR) $(ARFLAGS) flib.a $?
	$(RANLIB) flib.a
	touch $@

lapack_t3e : lapack_t3e.o 
	$(AR) $(ARFLAGS) flib.a $?
	touch $@

lapack_mkl : lapack_mkl.o 
	$(AR) $(ARFLAGS) flib.a $?
	touch $@

lapack_atlas : lapack_atlas.o dlamch.o
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
