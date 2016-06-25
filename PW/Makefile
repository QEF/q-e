# Makefile for PW
sinclude ../make.inc

default: all

all: pw pwtools

pw:	
	( cd src ; $(MAKE) all || exit 1 )

pw-lib:
	( cd src ; $(MAKE) libpw.a || exit 1 )

pwtools: pw
	( cd tools ; $(MAKE) all || exit 1 )

doc:
	( cd Doc ; $(MAKE) all || exit 1 )

doc_clean:
	( cd Doc ; $(MAKE) clean )

clean : examples_clean
	( cd src ; $(MAKE) clean )
	( cd tools ; $(MAKE) clean )

examples_clean:
	if test -d examples ; then \
	( cd examples ; ./clean_all ) ; fi 

distclean: clean doc_clean
