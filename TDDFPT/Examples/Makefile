# Makefile for TDDFPT Examples
#Modified for v4.0 by obm

default: all

all:	
	for dir in \
		CH4\
	; do \
	    if test -d $$dir ; then \
		( cd $$dir ; \
		if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
		else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) \
	    fi \
	done

clean :
	for dir in \
		CH4\
	; do \
	    if test -d $$dir ; then \
		( cd $$dir ; \
		if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= clean ; \
		else $(MAKE) $(MFLAGS) TLDEPS= clean ; fi ) \
	    fi \
	done

distclean: clean



