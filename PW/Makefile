# Makefile for PW

default: all

all: pw pwtools

pw:	
	( cd src ; $(MAKE) all || exit 1 )

pwtools: pw
	( cd tools ; $(MAKE) all || exit 1 )

doc:
	(cd Doc ; $(MAKE) $(MFLAGS) all || exit 1 )

doc_clean:
	(cd Doc ; $(MAKE) $(MFLAGS) clean ; fi )

clean : examples_clean
	( cd src ; $(MAKE) clean )
	( cd tools ; $(MAKE) clean )

examples_clean:
	( cd examples ; ./clean_all )


distclean: clean doc_clean
