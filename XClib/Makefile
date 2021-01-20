# Makefile for XClib

include ../make.inc

XCL = \
xclib_para_m.o \
xclib_error.o \
qe_kind.o \
qe_constants.o \
xc_input_params_mod.o \
dft_mod.o \
qe_funct_corr_lda_lsda.o \
qe_funct_exch_lda_lsda.o \
qe_funct_corr_gga.o \
qe_funct_exch_gga.o \
qe_funct_mgga.o \
xc_beef_interface.o \
qe_drivers_lda_lsda.o \
qe_drivers_d_lda_lsda.o \
qe_drivers_gga.o \
qe_drivers_d_gga.o \
qe_drivers_mgga.o \
xc_wrapper_d_gga.o \
xc_wrapper_d_lda_lsda.o \
xc_wrapper_gga.o \
xc_wrapper_lda_lsda.o \
xc_wrapper_mgga.o \
xc_lib.o

BEEF = \
pbecor.o \
beefun.o

all: xc_lib.a

xc_lib.a: $(XCL) $(BEEF)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

xclib_test.x : xclib_test.o xc_lib.a 
	$(LD) $(LD_LIBS) $(LDFLAGS) $(DFLAGS) -o $@ \
		xclib_test.o xc_lib.a $(BLAS_LIBS) $(LD_LIBS) $(LIBXC_LIBS)
	- ( cd ../bin ; ln -fs ../XClib/$@ . )

test : xclib_test.x

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

include make.depend

# DO NOT DELETE
