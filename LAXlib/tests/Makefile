# Makefile for LAXlib

include ../../make.inc

MODFLAGS= $(MOD_FLAG)../../UtilXlib $(MOD_FLAG)../../ELPA/src  $(MOD_FLAG) ..

SRCS = test_diaghg_1.f90 \
       test_diaghg_2.f90 \
       test_diaghg_gpu_1.f90 \
       test_diaghg_gpu_2.f90

EXECS = $(SRCS:.f90=.x)

all: common $(EXECS)

common: tester.o mp_world.o utils.o

%.x: %.o
	$(LD) $(LDFLAGS) $< tester.o mp_world.o utils.o -o $@ ../libqela.a \
	../../UtilXlib/libutil.a $(LIBOBJS) $(LIBS)

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

# .PHONY forces execution of a rule irrespective of the presence of an
# updated file with the same name of the rule. In this way, the script 
# that generates version.f90 always runs, updating the version if you 
# execute "svn update". The update_version script takes care of not
# changing the file if the svn version did not change

.PHONY: all clean

