# Makefile for KCWANN
sinclude ../make.inc

default: all

all: kc kc_pp

kc:
	( cd src ; $(MAKE) all || exit 1 )

kc_pp:
	( cd PostProc; $(MAKE) all || exit 1 )

clean: kc_clean kc_pp_clean

kc_clean:
	( cd src ; $(MAKE) clean )

kc_pp_clean:
	( cd PostProc ; $(MAKE) clean )

doc:
	if test -d Doc ; then \
	( cd Doc ; $(MAKE) all || exit 1 ) ; fi

doc_clean:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) clean ) ; fi

distclean: clean doc_clean
