# Makefile for flib

include ../make.sys
include ../make.rules

OBJS = \
avrec.o \
erf.o \
iceil.o \
iglocal.o \
inv3.o \
itoa.o \
latgen.o \
ngnr_set.o \
recips.o \
scnds.o \
simpsn.o \
sort.o \
sph_bes.o \
transto.o \
volume.o

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
reduce.o

all: flib.a ptools.a eispack $(MYLIB)

flib.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $?

blas_and_lapack : blas lapack

blas : blas.o dlamch.o
	$(AR) $(ARFLAGS) flib.a $?

lapack : lapack.o dlamch.o
	$(AR) $(ARFLAGS) flib.a $?

lapack_ibm : lapack_ibm.o dlamch.o
	$(AR) $(ARFLAGS) flib.a $?

lapack_t3e : lapack_t3e.o 
	$(AR) $(ARFLAGS) flib.a $?

lapack_mkl : lapack_mkl.o 
	$(AR) $(ARFLAGS) flib.a $?

eispack : eispack.o 
	echo eispack

dlamch.o : dlamch.f
	$(F77) $(F77FLAGS_NOOPT) -c $<

ptools.a : $(POBJS)
	$(AR) $(ARFLAGS) $@ $?

clean :
	- /bin/rm -f ptools.a *.o *.mod *.i core* work.pc flib.a

include .dependencies
