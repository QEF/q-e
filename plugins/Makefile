
#===============================
# espresso
#===============================
#
# Makefile for plugins
#

include ../make.sys
W90=wannier90-1.2
#
# MAIN target
#
all: 
	$(MAKE) w90

w90:
	if test -e archive/$(W90).tar.gz ; then \
	(if test ! -d ../$(W90); then \
	 ( gzip -dc archive/$(W90).tar.gz | (cd ../; tar -xvf -)); fi) ; fi
	if test -e install/make_wannier90.sys; then \
	(cp install/make_wannier90.sys ../$(W90)/make.sys); fi 
	cd ../$(W90); $(MAKE) all
	- (cd ../bin; ln -fs ../$(W90)/wannier90.x .) 
 
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

# general cleaning
clean: w90_clean
veryclean: w90_veryclean
