# Makefile for NEB

default: all

all:
	if test -d src ; then \
	( cd src ; if test "$(MAKE)" = "" ; then make $(MFLAGS) $@; \
	else $(MAKE) $(MFLAGS) $@ ; fi ) ; fi ; \

doc:
	if test -d Doc ; then \
	( cd Doc ; if test "$(MAKE)" = "" ; then make $(MFLAGS) $@; \
	else $(MAKE) $(MFLAGS) $@ ; fi) ; fi ; \

clean : examples_clean
	if test -d src ; then \
	( cd src ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi ;\

examples_clean:
	if test -d examples ; then \
	( cd examples ; ./clean_all ) ; fi 

doc_clean :
	if test -d Doc ; then \
	( cd Doc ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi ;\


distclean: clean doc_clean
