# Makefile for library coupling interface
include ../../make.inc

VPATH = @srcdir@

# location of required modules and include files
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)../../CPV/src

COUPLEOBJS= \
libpwscf.o \
libcpv.o \
libqemod.o

########################################################################

# see comments in file tldeps for explanations
TLDEPS=pw cp
include ../../install/tldeps

all : libqecouple.a

libqecouple.a : $(COUPLEOBJS)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

clean :
	- /bin/rm -f *.x *.o *.a *~ *_tmp.f90 *.d *.mod *.i *.L

# we currently need to specify explicit dependencies here
libpwscf.o: libpwscf.f90 ../../Modules/environment.o ../../Modules/mp_global.o \
	../../Modules/read_input.o ../../Modules/command_line_options.o \
	../../UtilXlib/parallel_include.o

libcpv.o: libcpv.f90 ../../Modules/environment.o ../../Modules/mp_global.o \
	../../Modules/read_input.o ../../Modules/command_line_options.o \
	../../Modules/mp_images.o ../../Modules/check_stop.o \
	../../Modules/io_global.o ../../UtilXlib/parallel_include.o \
	../../CPV/src/input.o

libqemod.o: libqemod.f90 ../../Modules/qmmm.o

