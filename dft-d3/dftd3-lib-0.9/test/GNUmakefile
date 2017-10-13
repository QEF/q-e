################################################################################
#
# Program makefile
#
# Needs following variables
#
# FC  Fortran compiler
# FCFLAGS  Fortran compiler flags
# LN  Linker
# LNFLAGS  Linker flags
#
################################################################################
.SUFFIX:
.SUFFIX: .f90 .o

.PHONY: clean distclean

LIBDIR = ../lib

all: testapi

OBJS_TESTAPI = testapi.o

testapi: $(OBJS_TESTAPI)
	$(LN) $(LNFLAGS) -o $@ $^ -L$(LIBDIR) -ldftd3

clean:
	rm -f *.o *.mod

distclean: clean
	rm -f testapi


%.o: %.f90
	$(FC) $(FCFLAGS) -I$(LIBDIR) -c $< -o $@


# Dependencies
testapi.o:
