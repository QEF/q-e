# Makefile for PW

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

doc:
	if test -d Doc ; then \
	(cd Doc ; if test "$(MAKE)" = "" ; then make $(MFLAGS) all ; \
	else $(MAKE) $(MFLAGS) all ; fi ) ; fi

doc_clean:
	if test -d Doc ; then \
	(cd Doc ; if test "$(MAKE)" = "" ; then make $(MFLAGS) clean ; \
	else $(MAKE) $(MFLAGS) clean ; fi ) ; fi

clean : examples_clean
	if test -d src ; then \
	( cd src ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi ;\
	if test -d tools ; then \
	( cd tools ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi ) ; fi ;\

examples_clean:
	if test -d examples ; then \
	( cd examples ; ./clean_all ) ; fi 


distclean: clean doc_clean
