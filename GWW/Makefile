include ../make.inc

all: libminpack pw4gwwa heada gwwa  bsea simplea simple_bsea simple_ipa

libminpack:
	( cd minpack ; $(MAKE) all || exit 1 )

pw4gwwa:
	( cd pw4gww ; $(MAKE) all || exit 1 )

heada: libminpack
	( cd head ; $(MAKE) all || exit 1 )

gwwa:
	( cd gww ; $(MAKE) all || exit 1 )

bsea: gwwa
	( cd bse ; $(MAKE) all || exit 1 )

simplea: 
	( cd simple ; $(MAKE) all || exit 1 )

simple_bsea:
	( cd simple_bse ; $(MAKE) all || exit 1 )

simple_ipa:
	( cd simple_ip ; $(MAKE) all || exit 1 )
clean: 
	( cd pw4gww  ; $(MAKE) clean )
	( cd head    ; $(MAKE) clean )
	( cd gww     ; $(MAKE) clean )
	( cd minpack ; $(MAKE) clean )
	( cd bse     ; $(MAKE) clean )
	( cd simple     ; $(MAKE) clean )
	( cd simple_bse     ; $(MAKE) clean )
	( cd simple_ip     ; $(MAKE) clean )

doc:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) all || exit 1 ) ; fi

doc_clean:
	if test -d Doc ; then \
	(cd Doc ; $(MAKE) clean ) ; fi

distclean: clean doc_clean
