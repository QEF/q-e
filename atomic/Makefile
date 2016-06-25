# Makefile for atomic

sinclude ../make.inc

default: all

all:
	( cd src ; $(MAKE) all || exit 1 )

clean : examples_clean
	( cd src ; $(MAKE) clean )

examples_clean:
	( cd examples ; ./make_clean )
	( cd pseudo_library ; ./make_clean )

doc:
	( cd Doc ; $(MAKE) all || exit 1 )
doc_clean:
	( cd Doc ; $(MAKE) clean )

distclean: clean doc_clean



