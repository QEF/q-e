# Makefile for TDDFPT Examples
#Modified for v4.0 by obm

default: normal_test


normal_test:	
	for dir in \
		CH4 CH4-TDHF CH4-DAV EELS-Si-NCPP-LDA \
	; do \
	    if test -d $$dir ; then \
		( cd $$dir ; \
		if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
		else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) \
	    fi \
	done

exx:	
	for dir in \
		CH4-TDHF CH4-PBE0 CH4-ECUTFOCK\
	; do \
	    if test -d $$dir ; then \
		( cd $$dir ; \
		if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
		else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) \
	    fi \
	done

small_exx:
	@for dir in \
		CH4-TDHF CH4-PBE0 CH4-ECUTFOCK\
	; do \
	    if test -d $$dir ; then \
		( cd $$dir ; \
		if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= small_test ; \
		else $(MAKE) $(MFLAGS) TLDEPS= small_test ; fi ) \
	    fi \
	done

clean :
	@for dir in \
		CH4 SiH4 CH4-PR Benzene CH4-TDHF CH4-PBE0 CH4-ECUTFOCK\
	; do \
	    if test -d $$dir ; then \
		( cd $$dir ; \
		if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= clean ; \
		else $(MAKE) $(MFLAGS) TLDEPS= clean ; fi ) \
	    fi \
	done

distclean: clean



