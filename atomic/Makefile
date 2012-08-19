# Makefile for atomic
# Adapted from TDDFPT main Makefile

default: all

all:
	if test -d src ; then \
	( cd src ; if test "$(MAKE)" = "" ; then make $(MFLAGS) $@; \
	else $(MAKE) $(MFLAGS) $@ ; fi ) ; fi ; \

clean : examples_clean
	if test -d src ; then \
	( cd src ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi ;\

examples_clean:
	if test -d examples ; then \
	( cd examples ; make_clean ) ; fi 
	if test -d pseudo_library ; then \
	( cd pseudo_library ; make_clean ) ; fi 

doc:
	if test -d Doc ; then \
	(cd Doc ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all ; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi
doc_clean:
	if test -d Doc ; then \
	(cd Doc ; if test "$(MAKE)" = "" ; then make $(MFLAGS) clean ; \
	else $(MAKE) $(MFLAGS) clean ; fi ) ; fi

distclean: clean doc_clean



