# Makefile for PIOUD
sinclude ../make.inc

default: all

all:
	if test -d src ; then \
	( cd src ; $(MAKE) || exit 1 ) ; fi

clean : examples_clean
	if test -d src ; then \
	( cd src ; $(MAKE) clean ) ; fi

examples_clean:
	if test -d examples ; then \
	( cd examples ; ./clean_all ) ; fi 

distclean: clean doc_clean
