# Makefile for HP
sinclude ../make.inc

default: all

all: hp

hp:
	( cd src ; $(MAKE) all || exit 1 )

hp-lib:
	( cd src ; $(MAKE) libs-hp || exit 1 )

clean: hp_clean

hp_clean:
	( cd src ; $(MAKE) clean )

examples_clean:
	if test -d examples ; then \
	( cd examples ; ./clean_all ) ; fi

doc:
	if test -d Doc ; then \
	( cd Doc ; $(MAKE) all || exit 1 ) ; fi

doc_clean:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) clean ) ; fi

distclean: clean doc_clean
