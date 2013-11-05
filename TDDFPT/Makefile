# Makefile for TDDFPT
#Modified for v4.0 by obm
#

default: all

all:
	if [ ! -d bin ] ; then mkdir bin; fi;
	if test -d src ; then \
	( cd src ; $(MAKE) $(MFLAGS) $@ ; ) ; fi ; \
	if test -d tools ; then \
	( cd tools ; $(MAKE) $(MFLAGS) $@ ; ) ; fi ; \

clean :
	if test -d src ; then \
	( cd src ; $(MAKE) clean ; ) ; fi ;\
	if test -d tools ; then \
	( cd tools ; $(MAKE) clean ; ) ; fi

distclean: clean



