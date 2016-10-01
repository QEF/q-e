# Makefile for atomic

sinclude ../make.inc

default: all

all:
	( cd src ; $(MAKE) all || exit 1 )

clean : examples_clean
	( cd src ; $(MAKE) clean )

examples_clean:
	if test -d examples ; then \
	( cd examples ; ./make_clean ) ; fi
	if test -d pseudo_library ; then \
	( cd pseudo_library ; ./make_clean) ; fi
	
doc:
	( cd Doc ; $(MAKE) all || exit 1 )

doc_clean:
	( cd Doc ; $(MAKE) clean )

distclean: clean doc_clean



