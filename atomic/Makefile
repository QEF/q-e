# Makefile for atomic

default: all

all:
	( cd src ; if test "$(MAKE)" = "" ; then make $@ || exit 1; \
	else $(MAKE) $@ || exit 1; fi )

clean : examples_clean
	( cd src ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi )

examples_clean:
	( cd examples ; ./make_clean )
	( cd pseudo_library ; ./make_clean )

doc:
	( cd Doc ; if test "$(MAKE)" = "" ; then make all || exit 1 ; \
	else $(MAKE) all || exit 1 ; fi )

doc_clean:
	( cd Doc ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi )

distclean: clean doc_clean
