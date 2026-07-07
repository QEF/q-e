# Makefile for TDDFPT
include ../make.inc

# see comments in file tldeps for explanations
include ../install/tldeps
TLDEPS=lrmods

all: bindir tddfpt

bindir:
	if [ ! -d bin ] ; then mkdir bin; fi;

tddfpt:
	( cd src ; $(MAKE) all || exit 1 ) 

clean :
	( cd src ; $(MAKE) clean ; )

doc:
	if test -d Doc ; then \
        ( cd Doc ; $(MAKE) all || exit 1 ) ; fi

doc_clean:
	if test -d Doc ; then \
        (cd Doc ; $(MAKE) clean ) ; fi

distclean: clean
	if [ -d bin ] ; then rm -rf bin; fi;

