# Makefile for COUPLE

default: all

all: libcouple

libcouple:	
	if test -d src ; then \
	( cd src ; $(MAKE) $(MFLAGS)  ) ; fi 

doc:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) $(MFLAGS) all  ) ; fi

doc_clean:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) $(MFLAGS) clean ; fi ) ; fi

clean : examples_clean tests_clean
	if test -d src ; then \
	( cd src ; $(MAKE) clean ) ; fi ;\

examples_clean:
	if test -d examples ; then \
	( cd examples ; $(MAKE) -f Makefile.gfortran clean ) ; fi 

tests_clean:
	if test -d tests ; then \
	( cd tests ; ./clean_all ) ; fi 


distclean: clean doc_clean
