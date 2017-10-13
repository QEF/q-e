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

all: dftd3

OBJS_DFTD3  = extras.o main.o

dftd3: $(OBJS_DFTD3) $(LIBDIR)/libdftd3.a
	$(LN) $(LNFLAGS) -o $@ $^

clean:
	rm -f *.o *.mod

distclean: clean
	rm -f dftd3


%.o: %.f90
	$(FC) $(FCFLAGS) -I$(LIBDIR) -c $< -o $@


# Dependencies
dftd3.o: extras.o
