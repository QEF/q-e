# Makefile for XSpectra

default: all

all:
	( cd src ; if test "$(MAKE)" = "" ; then make $@ || exit 1 ; \
	else $(MAKE) $@ || exit 1 ; fi )

clean :
	( cd src ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi )

distclean: clean



