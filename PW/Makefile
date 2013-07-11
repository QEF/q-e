# Makefile for PW

default: all

all: pw pwtools

pw:	
	if test -d src ; then \
	( cd src ; $(MAKE) $(MFLAGS)  ) ; fi 

pwtools:
	if test -d tools ; then \
	( cd tools ; $(MAKE) $(MFLAGS)  ) ; fi

doc:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) $(MFLAGS) all  ) ; fi

doc_clean:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) $(MFLAGS) clean ; fi ) ; fi

clean : examples_clean
	if test -d src ; then \
	( cd src ; $(MAKE) clean ) ; fi ;\
	if test -d tools ; then \
	( cd tools ; $(MAKE) clean ) ; fi

examples_clean:
	if test -d examples ; then \
	( cd examples ; ./clean_all ) ; fi 


distclean: clean doc_clean
