# Makefile for PW
sinclude ../make.inc

# see comments in file tldeps for explanations
TLDEPS=bindir mods libks_solvers dftd3
include ../install/tldeps

all: pw pwtools

pw: pwlibs
	( cd src ; $(MAKE) all || exit 1 )

pwtools: pwlibs
	( cd tools ; $(MAKE) all || exit 1 )

pwlibs:
	( cd src ; $(MAKE) libpw.a || exit 1 )

doc:
	( cd Doc ; $(MAKE) all )

doc_clean:
	( cd Doc ; $(MAKE) clean )

clean : examples_clean
	( cd src ; $(MAKE) clean )
	( cd tools ; $(MAKE) clean )

examples_clean:
	if test -d examples ; then \
	( cd examples ; ./clean_all ) ; fi 

distclean: clean doc_clean
