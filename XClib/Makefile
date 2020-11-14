# Makefile for XClib

include ../make.inc

MODFLAGS = $(MOD_FLAG)../UtilXlib

XCL = \
dft_mod.o \
qe_constants.o \
qe_drivers_d_gga.o \
qe_drivers_d_lda_lsda.o \
qe_drivers_gga.o \
qe_drivers_lda_lsda.o \
qe_drivers_mgga.o \
qe_funct_corr_gga.o \
qe_funct_corr_lda_lsda.o \
qe_funct_exch_gga.o \
qe_funct_exch_lda_lsda.o \
qe_funct_mgga.o \
qe_kind.o \
xc_beef_interface.o \
xc_input_params_mod.o \
xc_lib.o \
xclib_error.o \
xclib_para_m.o \
xc_wrapper_d_gga.o \
xc_wrapper_d_lda_lsda.o \
xc_wrapper_gga.o \
xc_wrapper_lda_lsda.o \
xc_wrapper_mgga.o \
xclib_test.o

BEEF = \
pbecor.o \
beefun.o

# PGI requires -c11 to understand __thread
# See https://www.pgroup.com/userforum/viewtopic.php?t=5504
ifneq (,$(findstring -D__PGI,$(DFLAGS)))
	CFLAGS := -c11 $(CFLAGS)
endif

all: xc_lib.a xclib_test.x

xc_lib.a: $(XCL) $(BEEF)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

xclib_test.x : xclib_test.o xc_lib.a 
	$(LD) $(LD_LIBS) $(LDFLAGS) $(DFLAGS) -o $@ \
		xclib_test.o xc_lib.a $(QELIBS)
	- ( cd ../bin ; ln -fs ../XClib/$@ . )

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

