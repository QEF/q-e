# Makefile for FFTXlib testing

include ../../make.inc

MODFLAGS= $(MOD_FLAG).. $(MOD_FLAG).

SRCS = test_fft_scalar_gpu.f90

EXECS = $(SRCS:.f90=.x)

all: common $(EXECS)

common: tester.o utils.o

%.x: %.o
	$(LD) $(LDFLAGS) $< utils.o tester.o -o $@ ../libqefft.a $(LIBS)

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

# .PHONY forces execution of a rule irrespective of the presence of an
# updated file with the same name of the rule. In this way, the script 
# that generates version.f90 always runs, updating the version if you 
# execute "svn update". The update_version script takes care of not
# changing the file if the svn version did not change

.PHONY: all clean

