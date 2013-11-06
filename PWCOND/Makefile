# Makefile for PWCOND

default: all

all:
	( cd src ; if test "$(MAKE)" = "" ; then make all || exit 1; \
	else $(MAKE) all || exit 1 ; fi )

clean : clean_examples
	( cd src ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi )

clean_examples :
	( cd examples ; ./clean_all)

doc:
	(cd Doc ; if test "$(MAKE)" = "" ; then make all || exit 1 ; \
	else $(MAKE) all || exit 1 ; fi )

doc_clean:
	(cd Doc ; if test "$(MAKE)" = "" ; then make clean ; \
	else $(MAKE) clean ; fi )

distclean: clean doc_clean



