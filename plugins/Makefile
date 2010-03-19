
#===============================
# espresso
#===============================
#
# Makefile for plugins
#

include ../make.sys
W90=wannier90-1.2
WANT=want-2.1.1
#
# MAIN target
#
all: w90 want

w90:
	if test -e archive/$(W90).tar.gz ; then \
	(if test ! -d ../$(W90); then \
	 ( gzip -dc archive/$(W90).tar.gz | (cd ../; tar -xvf -)); fi) ; fi
	if test -e install/make_wannier90.sys; then \
	(cp install/make_wannier90.sys ../$(W90)/make.sys); fi 
	cd ../$(W90); $(MAKE) all
	- (cd ../bin; ln -fs ../$(W90)/wannier90.x .) 

want:
	if test -e archive/$(WANT).tar.gz ; then \
	(if test ! -d ../$(WANT); then \
	( gzip -dc archive/$(WANT).tar.gz | (cd ../; tar -xvf -)); fi) ; fi
	cd ../$(WANT) ; \
	./conf/configure CC="$(CC)" \
#	CFLAGS="$(CFLAGS) $(DFLAGS) $(IFLAGS)" \
#	CPP="$(CPP)" \
#	CPPFLAGS="$(CPPFLAGS) $(DFLAGS) $(IFLAGS)" \
	F90="$(F90)" \
	MPIF90="$(MPIF90)" \
	F90FLAGS="$(F90FLAGS) $($PRE_FDFLAGS) $(FDFLAGS) $(IFLAGS) $(MODFLAGS)" \
	F77="$(F77)" \
	FFLAGS="$(FFLGAS)" \
	FFLAGS_NOOPT="$(FFLAGS_NOOPT)" \
	LD="$(LD)" \
	LDFLAGS="$(LDFLAGS)" \
	AR="$(AR)" \
	ARFLAGS="$(ARFLAGS)" \
	RANLIB="$(RANLIB)" \
	DFLAGS="$(DFLAGS)" \
	FDFLAGS="$(DFLAGS)" \
	IFLAGS="-I../include $(IFLAGS)"
	cd ../$(WANT); $(MAKE) all

###################################
# cleaning
###################################
# each lib independently
w90_clean:
	if test -d ../$(W90); then (cd ../$(W90); \
	$(MAKE) clean); fi
	if test -e ../$(W90)/wannier90.x; then \
	rm -f ../$(W90)/wannier90.x; fi
	if test -e ../$(W90)/libwannier.a; then \
	rm -f ../$(W90)/libwannier.a; fi
	rm -f ../bin/wannier90.x
w90_veryclean:
	if test -d ../$(W90); then (rm -R -f ../$(W90)); fi
	rm -f ../bin/wannier90.x
want_clean:
	if test -d ../$(WANT); then (cd ../$(WANT); \
	$(MAKE) clean); fi
#	if test -e ../$(W90)/wannier90.x; then \
#	rm -f ../$(W90)/wannier90.x; fi
#	if test -e ../$(W90)/libwannier.a; then \
#	rm -f ../$(W90)/libwannier.a; fi
#	rm -f ../bin/wannier90.x

want_veryclean:
	if test -d ../$(WANT); then (rm -R -f ../$(WANT)); fi

# general cleaning
clean: w90_clean want_clean
veryclean: w90_veryclean want_veryclean
