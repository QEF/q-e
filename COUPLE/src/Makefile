# Makefile for library coupling interface
include ../../make.inc

# location of required modules and include files
MODFLAGS= $(MOD_FLAG)../../iotk/src $(MOD_FLAG)../../Modules $(MOD_FLAG)../../FFTXlib $(MOD_FLAG)../../CPV/src
IFLAGS=#-DDEBUG_QECOUPLE
TLDEPS=pw cp

COUPLEOBJS= \
libpwscf.o \
libcpv.o \
libqemod.o

########################################################################

all : tldeps libqecouple.a

libqecouple.a : $(COUPLEOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

tldeps :
	if test -n "$(TLDEPS)" ; then \
	( cd ../.. ; $(MAKE) $(TLDEPS) || exit 1 ) ; fi

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L

# we currently need to specify explicit dependencies here
libpwscf.o: libpwscf.f90 ../../Modules/environment.o ../../Modules/mp_global.o \
	../../Modules/read_input.o ../../Modules/command_line_options.o \
	../../Modules/parallel_include.o

libcpv.o: libcpv.f90 ../../Modules/environment.o ../../Modules/mp_global.o \
	../../Modules/read_input.o ../../Modules/command_line_options.o \
	../../Modules/mp_images.o ../../Modules/check_stop.o \
	../../Modules/io_global.o ../../Modules/parallel_include.o \
	../../CPV/src/input.o

libqemod.o: libqemod.f90 ../../Modules/qmmm.o

# DO NOT DELETE
