# Makefile for PostProc
sinclude ../make.inc

default: all

all:
	if test -d src ; then \
	( cd src ; $(MAKE) || exit 1 ) ; fi

doc:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) all || exit 1 ) ; fi

doc_clean:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) clean ) ; fi

clean : examples_clean
	if test -d src ; then \
	( cd src ; $(MAKE) clean ) ; fi

examples_clean:
	if test -d examples ; then \
	( cd examples ; ./clean_all ) ; fi

distclean: clean doc_clean
