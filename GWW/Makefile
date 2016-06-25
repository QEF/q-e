include ../make.inc

all: libminpack pw4gwwa heada gwwa  bsea

libminpack:
	( cd minpack ; $(MAKE) all || exit 1 )

pw4gwwa:
	( cd pw4gww ; $(MAKE) all || exit 1 )

heada: libminpack
	( cd head ; $(MAKE) all || exit 1 )

gwwa:
	( cd gww ; $(MAKE) all || exit 1 )

bsea:
	( cd bse ; $(MAKE) all || exit 1 )

clean: 
	( cd pw4gww  ; $(MAKE) clean )
	( cd head    ; $(MAKE) clean )
	( cd gww     ; $(MAKE) clean )
	( cd minpack ; $(MAKE) clean )
	( cd bse     ; $(MAKE) clean )

doc:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) all || exit 1 ) ; fi

doc_clean:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) clean ) ; fi

distclean: clean doc_clean
