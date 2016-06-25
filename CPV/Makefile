# Makefile for CPV
sinclude ../make.inc

default: all

all:
	( cd src ; $(MAKE) all || exit 1 )

clean :
	( cd src ; $(MAKE) clean )

doc:
	(cd Doc ; $(MAKE) all || exit 1 )

doc_clean:
	(cd Doc ; $(MAKE) clean )

distclean: clean doc_clean



