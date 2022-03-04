# Makefile for KCWANN
sinclude ../make.inc

default: all

all: kcw pp

kcw:
	( cd src ; $(MAKE) all || exit 1 )

pp:
	( cd PP; $(MAKE) all || exit 1 )

clean: kcw_clean pp_clean

kcw_clean:
	( cd src ; $(MAKE) clean )

pp_clean:
	( cd PP ; $(MAKE) clean )

doc:
	if test -d Doc ; then \
	( cd Doc ; $(MAKE) all || exit 1 ) ; fi

doc_clean:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) clean ) ; fi

distclean: clean doc_clean
