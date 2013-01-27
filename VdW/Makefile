# Makefile for VdW
# Adapted from TDDFPT main Makefile

default: all

all: vdw

vdw:	
	@echo "VdW no longer maintained"
#	if test -d src ; then \
#	( cd src ; if test "$(MAKE)" = "" ; then make $(MFLAGS) $@; \
#	else $(MAKE) $(MFLAGS) ; fi ) ; fi ; \

clean :
	if test -d src ; then \
	( cd src ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi ;\

distclean: clean
