# Makefile for TDDFPT
#

default: all

all: bindir
	( cd src ; if test "$(MAKE)" = "" ; then make $@ || exit 1; \
	else $(MAKE) $@ || exit 1; fi )
	( cd tools ; if test "$(MAKE)" = "" ; then make $@; \
	else $(MAKE) $@  || exit 1; fi )

bindir:
	test -d bin || mkdir bin

clean :
	( cd src ; $(MAKE) clean )
	( cd tools ; $(MAKE) clean )

distclean: clean



