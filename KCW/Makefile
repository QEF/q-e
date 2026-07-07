# Makefile for KCW
sinclude ../make.inc

# see comments in file tldeps for explanations
TLDEPS=pwlibs lrmods
include ../install/tldeps

all: kcw pp

kcw: kcwlib
	( cd src ; $(MAKE) all || exit 1 )

kcwlib:
	( cd src ; $(MAKE) kcwlib || exit 1 )

pp: kcwlib
	( cd PP ; $(MAKE) all || exit 1 )

clean: kcw_clean pp_clean

kcw_clean:
	( cd src ; $(MAKE) clean )

pp_clean:
	( cd PP ; $(MAKE) clean )

doc:
	if test -d Doc ; then \
	( cd Doc ; $(MAKE) all || exit 1 ) ; fi

doc_clean:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) clean ) ; fi

distclean: clean doc_clean
