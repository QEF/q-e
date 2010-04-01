
#===============================
# espresso
#===============================
#
# Makefile for plugins
#

include ../make.sys
W90=wannier90-1.2
WANT=want-2.1.1
YAMBO=yambo-3.2.1-r.448
#
# MAIN target
#

all: w90 want

EXTLIBS=liblapack libblas libiotk

###########################################################
# wannier90
###########################################################

w90:
	if test -e archive/$(W90).tar.gz ; then \
	(if test ! -d ../$(W90); then \
	 ( gzip -dc archive/$(W90).tar.gz | (cd ../; tar -xvf -)); fi) ; fi
	if test -e install/make_wannier90.sys; then \
	(cp install/make_wannier90.sys ../$(W90)/make.sys); fi 
	cd ../$(W90); $(MAKE) all
	- (cd ../bin; ln -fs ../$(W90)/wannier90.x .) 

###########################################################
# WANT
###########################################################

want: extlibs uncompress-want configure-want
	if test -d ../$(WANT); then \
	cd ../$(WANT); $(MAKE) all; fi
#	- ( cd ../bin ; ln -sf ../$(WANT)/bin/* . ) 
	touch ./make-want
uncompress-want:
	if test -e archive/$(WANT).tar.gz ; then \
	(if test ! -d ../$(WANT); then \
	( gzip -dc archive/$(WANT).tar.gz | (cd ../; tar -xvf -)); fi) ; fi
	touch ./uncompress-want
configure-want:
	cd ../$(WANT) ; \
	./conf/configure \
	CC="$(CC)" \
	F90="$(F90)" \
	MPIF90="$(MPIF90)" \
	F90FLAGS="$(F90FLAGS) $($PRE_FDFLAGS) $(FDFLAGS) $(IFLAGS) $(MODFLAGS)" \
	F77="$(F77)" \
	FFLAGS="$(FFLAGS)" \
	FFLAGS_NOOPT="$(FFLAGS_NOOPT)" \
	LD="$(LD)" \
	LDFLAGS="$(LDFLAGS)" \
	AR="$(AR)" \
	ARFLAGS="$(ARFLAGS)" \
	RANLIB="$(RANLIB)" \
	DFLAGS="$(DFLAGS)" \
	FDFLAGS="$(DFLAGS)" \
	IFLAGS="-I../include $(IFLAGS)" \
	LAPACK_LIBS="$(LAPACK_LIBS)" \
	BLAS_LIBS="$(BLAS_LIBS)" \
	FFT_LIBS="$(FFT_LIBS)" 
	touch ./configure-want

###################################################
# YAMBO
###################################################

yambo: extlibs uncompress-yambo configure-yambo
	if test -d ../$(YAMBO); then \
	cd ../$(YAMBO); $(MAKE) all; fi
#	- ( cd ../bin ; ln -sf ../$(YAMBO)/bin/* . ) 
	touch ./make-yambo
uncompress-yambo:
	if test -e archive/$(YAMBO).tar.gz ; then \
	(if test ! -d ../$(YAMBO); then \
	( gzip -dc archive/$(YAMBO).tar.gz | (cd ../; tar -xvf -)); fi) ; fi
	touch ./uncompress-yambo
configure-yambo:
	cd ../$(YAMBO) ; \
	./configure --with-blas="$(BLAS_LIBS)" \
	--with-lapack="$(LAPACK_LIBS)" \
	--with-iotk="$(TOPDIR)/iotk" \
	PFC="$(MPIF90)" \
	FC="$(LD)" \
	FCFLAGS="$(F90FLAGS) $(FDFLAGS) $(IFLAGS) $(MODFLAGS)"
	F77="$(F77)" \
	FFLAGS="$(FFLAGS)"
	AR="$(AR)" \
	ARFLAGS="$(ARFLAGS)" 
	touch ./configure-yambo
##################################
# extlibs 	FCFLAGS="$(F90FLAGS) $(FDFLAGS) $(IFLAGS) $(MODFLAGS)" \
###################################
extlibs:
	test -n "$(EXTLIBS)" && ( cd .. ; $(MAKE) $(MFLAGS) $(EXTLIBS) || exit 1) || :
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
	- rm -rf ./make-want ./configure-want

want_veryclean:
	if test -d ../$(WANT); then (rm -R -f ../$(WANT)); fi
	- rm -rf ./make-want ./configure-want ./uncompress-want
yambo_clean:
	if test -d ../$(YAMBO); then (cd ../$(YAMBO); \
	$(MAKE) clean); fi
	- rm -rf ./make-yambo ./configure-yambo

yambo_veryclean:
	if test -d ../$(YAMBO); then (rm -R -f ../$(YAMBO)); fi
	- rm -rf ./make-yambo ./configure-yambo ./uncompress-yambo

##########################################################
# general cleaning
##########################################################
clean: w90_clean want_clean yambo_clean
veryclean: w90_veryclean want_veryclean yambo_veryclean
