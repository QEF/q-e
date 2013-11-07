# Makefile for PWCOND
# Adapted from TDDFPT main Makefile

default: all

all:
	( cd src ; $(MAKE) all || exit 1 )

doc:
	( cd Doc ; $(MAKE) all || exit 1 )

doc_clean:
	( cd Doc ; $(MAKE) clean )

clean : examples_clean
	( cd src ; $(MAKE) clean )

examples_clean:
	if test -d examples ; then \
	( cd examples ; ./clean_all ) ; fi 

distclean: clean doc_clean
