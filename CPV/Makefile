# Makefile for CPV
sinclude ../make.inc

# see comments in file tldeps for explanations
TLDEPS=bindir mods
include ../install/tldeps

all:
	( cd src ; $(MAKE) all || exit 1 )

clean :
	( cd src ; $(MAKE) clean )

doc:
	(cd Doc ; $(MAKE) all )

doc_clean:
	(cd Doc ; $(MAKE) clean )

distclean: clean doc_clean



