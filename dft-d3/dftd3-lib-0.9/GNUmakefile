.SUFFIX:
.SUFFIX: .f90 .o

.PHONY: dftd3 lib testapi

all: lib dftd3 testapi

include make.arch

lib:
	$(MAKE) -C lib FC="$(FC)" FCFLAGS="$(FCFLAGS)" LN="$(LN)" \
            LNFLAGS="$(LNFLAGS)" SRCDIR="."

dftd3: lib
	$(MAKE) -C prg FC="$(FC)" FCFLAGS="$(FCFLAGS)" LN="$(LN)" \
	    LNFLAGS="$(LNFLAGS)" dftd3

testapi: lib
	$(MAKE) -C test FC="$(FC)" FCFLAGS="$(FCFLAGS)" LN="$(LN)" \
	    LNFLAGS="$(LNFLAGS)" testapi


.PHONY: clean distclean
clean:
	$(MAKE) -C lib clean
	$(MAKE) -C prg clean
	$(MAKE) -C test clean

distclean:
	$(MAKE) -C lib distclean
	$(MAKE) -C prg distclean
	$(MAKE) -C test distclean
