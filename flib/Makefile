
include ../make.sys

all: flib.a ptools.a eispack $(MYLIB)

OBJS = inv3.o itoa.o sort.o simpsn.o \
       transto.o iceil.o erf.o

POBJS = gridsetup.o blockset.o hangup.o reduce.o \
       matmul.o localindex.o localdim.o ownerof.o \
       ilcm.o indxg2l.o indxg2p.o indxl2g.o infog1l.o infog2l.o \
       numroc.o npreroc.o


flib.a : $(OBJS)
	ar $(ARFLAGS) $@ $?
	$(RANLIB) $@

blas_and_lapack : blas lapack

blas : blas.o dlamch.o
	ar $(ARFLAGS) flib.a $?
	$(RANLIB) flib.a

lapack : lapack.o dlamch.o
	ar $(ARFLAGS) flib.a $?
	$(RANLIB) flib.a

lapack_ibm : lapack_ibm.o dlamch.o
	ar $(ARFLAGS) flib.a $?
	$(RANLIB) flib.a

lapack_t3e : lapack_t3e.o 
	ar $(ARFLAGS) flib.a $?
	$(RANLIB) flib.a

eispack : eispack.o
	echo eispack

dlamch.o : dlamch.f
	$(F77) $(F77FLAGS_NOOPT) -c $<

ptools.a : $(POBJS)
	ar $(ARFLAGS) $@ $?

include .dependencies
include ../make.rules

clean_:
	- rm -f ptools.a *.o *.mod *.i core* work.pc flib.a
