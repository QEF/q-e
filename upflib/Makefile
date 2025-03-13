# Makefile for upflib

include ../make.inc

MODFLAGS= $(MOD_FLAG)../UtilXlib $(MOD_FLAG).

# list of modules:
# OBJS_NODEP do not depend upon any external modules
# OBJS_DEP   depend upon modules in UtilXlib/
# OBJS_GPU   are GPU-specific

OBJS_DEP= \
atwfc_mod.o \
beta_mod.o \
qrad_mod.o  \
rhoat_mod.o \
rhoc_mod.o \
vloc_mod.o \
qvan2.o \
dqvan2.o \
gen_us_dj.o \
gen_us_dy.o \
init_us_0.o \
init_us_b0.o \
init_us_1.o \
init_us_2_acc.o

OBJS_NODEP= \
atom.o \
atomic_number.o \
dylmr2.o \
gth.o \
paw_variables.o \
pseudo_types.o \
radial_grids.o \
read_cpmd.o \
read_fhi.o \
read_ncpp.o \
read_ps.o \
read_psml.o \
read_upf_new.o \
read_upf_v1.o \
read_uspp.o \
spinor.o \
sph_ind.o \
sph_bes.o \
splinelib.o \
simpsn.o \
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
uspp_param.o \
write_upf_new.o \
xmltools.o \
ylmr2.o

# CUDA Fortran routines
OBJS_GPU= \
  ylmr2_gpu.o

OBJS_NODEP+= $(OBJS_GPU) dom.o wxml.o

all : libupf.a virtual_v2.x upfconv.x casino2upf.x

libupf.a: $(OBJS_DEP) $(OBJS_NODEP) $(OBJS_GPU)
	$(AR) $(ARFLAGS) $@ $?       
	$(RANLIB) $@    

# The following utilities are linked with $(OBJS_NODEP) instead of libupf.a
# because the latter has external dependencies that produce "missing symbols"
# errors for some linkers that try to resolve symbols even if never used

virtual_v2.x : virtual_v2.o $(OBJS_NODEP)
	$(LD) $(LDFLAGS) -o $@ virtual_v2.o $(OBJS_NODEP) $(LAPACK_LIBS) $(BLAS_LIBS)

casino2upf.x : casino2upf.o casino_pp.o $(OBJS_NODEP)
	$(LD) $(LDFLAGS) -o $@ casino2upf.o casino_pp.o $(OBJS_NODEP) $(LAPACK_LIBS) $(BLAS_LIBS)
upfconv.x : upfconv.o casino_pp.o $(OBJS_NODEP) 
	$(LD) $(LDFLAGS) -o $@ upfconv.o casino_pp.o $(OBJS_NODEP) $(LAPACK_LIBS) $(BLAS_LIBS)

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L

include make.depend

