# Makefile for upflib

include ../make.inc

MODFLAGS= $(MOD_FLAG)../UtilXlib $(MOD_FLAG)../external/devxlib/src

# list of modules:
# OBJS_NODEP do not depend upon any external modules
# OBJS_DEP   depend upon modules in UtilXlib/
# OBJS_GPU   depend upon modules in devxlib

OBJS_DEP= \
gen_us_dj.o \
gen_us_dy.o \
init_us_0.o \
init_us_b0.o \
init_us_1.o \
init_us_2_base.o \
init_tab_atwfc.o \
init_tab_beta.o \
init_tab_qrad.o

OBJS_NODEP= \
atom.o \
atomic_number.o \
dqvan2.o \
dylmr2.o \
gth.o \
interp_atwfc.o \
paw_variables.o \
pseudo_types.o \
qvan2.o \
radial_grids.o \
read_cpmd.o \
read_fhi.o \
read_ncpp.o \
read_ps.o \
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
uspp_data.o \
uspp_param.o \
write_upf_new.o \
xmltools.o \
ylmr2.o

# GPU versions of routines
OBJS_GPU= \
  dylmr2_gpu.o \
  init_us_2_base_gpu.o \
  interp_atwfc_gpu.o \
  qvan2_gpu.o \
  simpsn_gpu.o \
  sph_bes_gpu.o \
  ylmr2_gpu.o \
  upf_acc_interfaces.o

OBJS_NODEP+= dom.o wxml.o

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

