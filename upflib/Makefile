# Makefile for upflib

include ../make.inc

QEMODS = libupf.a
MODFLAGS= $(MOD_FLAG)../UtilXlib

# list of modules

OBJS= \
atom.o \
atomic_number.o \
dqvan2.o \
gth.o \
init_us_0.o \
init_us_b0.o \
init_us_1.o \
init_at_1.o \
init_us_2_base.o \
interp_at_wfc.o \
paw_variables.o \
pseudo_types.o \
qvan2.o \
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
upf_invmat.o \
upf_io.o \
upf_ions.o \
upf_kinds.o \
upf_params.o \
upf_parallel_include.o \
upf_spinorb.o \
upf_to_internal.o \
upf_utils.o \
uspp.o \
uspp_data.o \
write_upf_new.o \
xmltools.o \
ylmr2.o

# GPU versions of routines
OBJS += \
  interp_at_wfc_gpu.o \
  simpsn_gpu.o \
  qvan2_gpu.o \
  init_us_2_base_gpu.o \
  ylmr2_gpu.o

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

