# Makefile for upflib

include ../make.inc

QEMODS = libupf.a

# list of modules

OBJS= \
atom.o \
atomic_number.o \
upf_erf.o \
upf_utils.o \
gth.o \
pseudo_types.o \
radial_grids.o \
read_cpmd.o \
read_fhi.o \
read_ncpp.o \
read_upf_new.o \
read_upf_v1.o \
read_uspp.o \
splinelib.o \
simpsn.o \
upf.o \
upf_auxtools.o \
upf_const.o \
upf_error.o \
upf_free_unit.o \
upf_invmat.o \
upf_io.o \
upf_kinds.o \
upf_params.o \
upf_parallel_include.o \
upf_to_internal.o \
uspp.o \
write_upf_new.o \
xmltools.o \
ylmr2.o

-include make.gpu

all : libupf.a virtual_v2.x upfconv.x

libupf.a: $(MODULES) $(OBJS)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

virtual_v2.x : virtual_v2.o $(QEMODS) 
	$(LD) $(LDFLAGS) -o $@ virtual_v2.o $(QEMODS) $(LAPACK_LIBS) $(BLAS_LIBS)

upfconv.x : upfconv.o casino_pp.o $(QEMODS) 
	$(LD) $(LDFLAGS) -o $@ upfconv.o casino_pp.o $(QEMODS) $(LAPACK_LIBS) $(BLAS_LIBS)

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L

include make.depend

