# Makefile for flib

include ../make.sys
include ../make.rules

OBJS = \
avrec.o \
atomic_number.o \
capital.o \
dost.o \
erf.o \
iceil.o \
iglocal.o \
invmat.o \
latgen.o \
linpack.o \
matches.o \
ngnr_set.o \
recips.o \
rranf.o \
scnds.o \
simpsn.o \
sort.o \
sph_bes.o \
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
matmul.o \
npreroc.o \
numroc.o \
ownerof.o \
error.o \
reduce.o

all: flib.a ptools.a eispack $(MYLIB)

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

eispack : eispack.o 
	touch $@

dlamch.o : dlamch.f
	$(F77) $(F77FLAGS_NOOPT) -c $<

ptools.a : $(POBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f ptools.a *.o *.mod *.i *.F90 core* work.pc flib.a \
		blas_and_lapack blas lapack lapack_ibm lapack_t3e \
		lapack_mkl eispack lapack_atlas

include make.depend
