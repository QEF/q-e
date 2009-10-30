
#===============================
# espresso
#===============================
#
# Makefile for plugins
#

include ../make.sys

#
# MAIN target
#
all: 
	$(MAKE) w90

w90:
	if test ! -d ../wannier90-1.1; then \
	( cat archive/wannier90-1.1.tar | (cd ../; tar -xvf -)); fi
	if test -e install/make_wannier90.sys; then \
	(cp install/make_wannier90.sys ../wannier90-1.1/make.sys); fi 
	cd ../wannier90-1.1; $(MAKE) all
	- (cd ../bin; ln -fs ../wannier90-1.1/wannier90.x .) 
 
###################################
# cleaning
###################################
# each lib independently
w90_clean:
	if test -d ../wannier90-1.1; then (cd ../wannier90-1.1; \
	$(MAKE) clean); fi
	if test -e ../wannier90-1.1/wannier90.x; then \
	rm -f ../wannier90-1.1/wannier90.x; fi
	if test -e ../wannier90-1.1/libwannier.a; then \
	rm -f ../wannier90-1.1/libwannier.a; fi
	rm -f ../bin/wannier90.x
w90_veryclean:
	if test -d ../wannier90-1.1; then (rm -R -f ../wannier90-1.1); fi
	rm -f ../bin/wannier90.x

# general cleaning
clean: w90_clean
veryclean: w90_veryclean
