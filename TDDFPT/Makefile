# Makefile for TDDFPT
#Modified for v4.0 by obm

default: all

all:
	if test -d src ; then \
	( cd src ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi ; \
	mv src/tddfpt.x bin/;\
	( cd ../bin; ln -fs ../TDDFPT/bin/tddfpt.x . ) ;\
	if test -d tools ; then \
	( cd tools ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi ; \
	mv tools/tddfpt_pp.x bin/;\
	( cd ../bin; ln -fs ../TDDFPT/bin/tddfpt_pp.x . ) ;\

clean :
	if test -d src ; then \
	( cd src ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi ;\
	if test -d tools ; then \
	( cd tools ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi ;\
	/bin/rm -f bin/tddfpt.x bin/tddfpt_pp.x

distclean: clean



