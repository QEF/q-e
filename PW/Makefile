# Makefile for PW
# Adapted from TDDFPT main Makefile

default: all

all: pw pwtools

pw:	
	if test -d src ; then \
	( cd src ; if test "$(MAKE)" = "" ; then make $(MFLAGS) $@; \
	else $(MAKE) $(MFLAGS) ; fi ) ; fi ; \

pwtools:
	if test -d tools ; then \
	( cd tools ; if test "$(MAKE)" = "" ; then make $(MFLAGS) $@; \
	else $(MAKE) $(MFLAGS) ; fi ) ; fi ; \


clean :
	if test -d src ; then \
	( cd src ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi ;\
	if test -d tools ; then \
	( cd tools ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi ;\


distclean: clean
