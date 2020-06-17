# Makefile for upflib

include ../make.inc

# location of needed modules
MODFLAGS=$(MOD_FLAG)$(TOPDIR)/FoX/finclude

QEMODS = libupf.a

# list of modules

OBJS= \
atom.o \
atomic_number.o \
emend_upf.o \
upf_erf.o \
upf_utils.o \
gth.o \
pseudo_types.o \
radial_grids.o \
read_cpmd.o \
read_fhi.o \
read_ncpp.o \
read_upf_new.o \
read_upf_schema.o \
read_upf_v1.o \
read_upf_v2.o \
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
write_upf.o \
write_upf_new.o \
write_upf_schema.o \
write_upf_v2.o \
xmltools.o

TLDEPS=libfox

all : libupf.a virtual_v2.x upfconv.x

libupf.a: $(MODULES) $(OBJS)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

virtual_v2.x : virtual_v2.o $(QEMODS) 
	$(LD) $(LDFLAGS) -o $@ virtual_v2.o $(QEMODS) $(FOX_LIB) $(LAPACK_LIBS) $(BLAS_LIBS) 

upfconv.x : upfconv.o casino_pp.o $(QEMODS) 
	$(LD) $(LDFLAGS) -o $@ upfconv.o casino_pp.o $(QEMODS) $(FOX_LIB) $(LAPACK_LIBS) $(BLAS_LIBS) 

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L

include make.depend

