# Makefile for HP
sinclude ../make.inc

# see comments in file tldeps for explanations
TLDEPS=pwlibs lrmods
include ../install/tldeps

all: hp

hp: hplibs
	( cd src ; $(MAKE) all || exit 1 )

hplibs:
	( cd src ; $(MAKE) hplibs || exit 1 )

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
