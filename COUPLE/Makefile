# Makefile for COUPLE

sinclude ../make.inc

default: all

all: libcouple

libcouple:	
	( cd src ; $(MAKE) all || exit 1 ) 

doc:
	if test -d Doc ; then \
	( cd Doc ; $(MAKE) all || exit 1) ; fi

doc_clean:
	if test -d Doc ; then \
	( cd Doc ; $(MAKE) clean ; fi ) ; fi

clean : examples_clean tests_clean
	if test -d src ; then \
	( cd src ; $(MAKE) clean ) ; fi ;\

examples_clean:
	if test -d examples ; then \
	( cd examples ; $(MAKE) -f Makefile.gfortran clean ) ; fi 

tests_clean:
	if test -d tests ; then \
	( cd tests ; ./clean_all ) ; fi 


distclean: clean doc_clean
