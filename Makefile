# Copyright (C) 2001-2016 Quantum ESPRESSO group
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

include make.inc

default :
	@echo 'to install Quantum ESPRESSO, type at the shell prompt:'
	@echo '  ./configure [--prefix=]'
	@echo '  make [-j] target'
	@echo ' '
	@echo 'where target identifies one or multiple CORE PACKAGES:'
	@echo '  pw           basic code for scf, structure optimization, MD'
	@echo '  ph           phonon code, Gamma-only and third-order derivatives'
	@echo '  pwcond       ballistic conductance'
	@echo '  neb          code for Nudged Elastic Band method'
	@echo '  pp           postprocessing programs'
	@echo '  pwall        same as "make pw ph pp pwcond neb"'
	@echo '  cp           CP code: CP MD with ultrasoft pseudopotentials'
	@echo '  tddfpt       time dependent dft code'
	@echo '  gwl          GW with Lanczos chains'
	@echo '  ld1          utilities for pseudopotential generation'
	@echo '  upf          utilities for pseudopotential conversion'
	@echo '  xspectra     X-ray core-hole spectroscopy calculations'
	@echo '  couple       Library interface for coupling to external codes'
	@echo '  epw          Electron-Phonon Coupling with wannier functions'
	@echo '  gui          Graphical User Interface'
	@echo '  examples     fetch from web examples for all core packages'
	@echo '  test-suite   run semi-automated test-suite for regression testing'
	@echo '  all          same as "make pwall cp ld1 upf tddfpt"'
	@echo ' '
	@echo 'where target identifies one or multiple THIRD-PARTIES PACKAGES:'
	@echo '  gipaw        NMR and EPR spectra'
	@echo '  w90          Maximally localised Wannier Functions'
	@echo '  want         Quantum Transport with Wannier functions'
	@echo '  west         Many-body perturbation corrections Without Empty STates'
#	@echo '  SaX          Standard GW-BSE with plane waves'
	@echo '  yambo        electronic excitations with plane waves'
	@echo '  yambo-devel  yambo devel version'
	@echo '  SternheimerGW calculate GW using Sternheimer equations'
	@echo '  plumed       Metadynamics plugin for pw or cp'
	@echo '  d3q          general third-order code and thermal transport codes'
	@echo ' '
	@echo 'where target is one of the following suite operation:'
	@echo '  doc          build documentation'
	@echo '  links        create links to all executables in bin/'
	@echo '  tar          create a tarball of the source tree'
	@if test -d GUI/; then \
		echo '  tar-gui      create a standalone PWgui tarball from the GUI sources'; \
		echo '  tar-qe-modes create a tarball for QE-modes (Emacs major modes for Quantum ESPRESSO)'; fi
	@echo '  clean        remove executables and objects'
	@echo '  veryclean    remove files produced by "configure" as well'
	@echo '  distclean    revert distribution to the original status'

###########################################################
# Main targets
###########################################################

# The syntax "( cd PW ; $(MAKE) TLDEPS= all || exit 1)" below
# guarantees that error code 1 is returned in case of error and make stops
# If "|| exit 1" is not present, the error code from make in subdirectories
# is not returned and make goes on even if compilation has failed

pw : bindir libs mods libdavid libcg dftd3
	if test -d PW ; then \
	( cd PW ; $(MAKE) TLDEPS= all || exit 1) ; fi

cp : bindir libs mods libdavid libcg
	if test -d CPV ; then \
	( cd CPV ; $(MAKE) TLDEPS= all || exit 1) ; fi

ph : pw lrmods
	if test -d PHonon; then \
	(cd PHonon; $(MAKE) all || exit 1) ; fi

neb : pw
	if test -d NEB; then \
  (cd NEB; $(MAKE) all || exit 1) ; fi

tddfpt : pw
	if test -d TDDFPT; then \
	(cd TDDFPT; $(MAKE) all || exit 1) ; fi

pp : pw
	if test -d PP ; then \
	( cd PP ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

pwcond : pp
	if test -d PWCOND ; then \
	( cd PWCOND ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

acfdt : ph
	if test -d ACFDT ; then \
	( cd ACFDT ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

# target still present for backward compatibility
gww:
	@echo '"make gww" is obsolete, use "make gwl" instead '

gwl : ph
	if test -d GWW ; then \
	( cd GWW ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

gipaw : pw
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

d3q : ph
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

ld1 : bindir libs mods
	if test -d atomic ; then \
	( cd atomic ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

upf : libs mods
	if test -d upftools ; then \
	( cd upftools ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

pw_export : pw
	if test -d PP ; then \
	( cd PP ; $(MAKE) TLDEPS= pw_export.x || exit 1 ) ; fi

xspectra : pw
	if test -d XSpectra ; then \
	( cd XSpectra ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

couple : pw cp
	if test -d COUPLE ; then \
	( cd COUPLE ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

epw: ph
	if test -d EPW ; then \
	( cd EPW ; $(MAKE) all || exit 1; \
		cd ../bin; ln -fs ../EPW/bin/epw.x . ); fi

travis : pwall epw
	if test -d test-suite ; then \
	( cd test-suite ; make run-travis || exit 1 ) ; fi

gui :
	@echo 'Check "GUI/README" how to access the Graphical User Interface'
#@echo 'Check "PWgui-X.Y/README" how to access the Graphical User Interface'

examples :
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

pwall : pw neb ph pp pwcond acfdt

all   : pwall cp ld1 upf tddfpt xspectra gwl 

###########################################################
# Auxiliary targets used by main targets:
# compile modules, libraries, directory for binaries, etc
###########################################################

mods : libiotk libfox libutil libla libfft
	( cd Modules ; $(MAKE) TLDEPS= all || exit 1 )

libdavid_rci : libs libutil libla
	( cd KS_Solvers/Davidson_RCI ; $(MAKE) TLDEPS= all || exit 1 )

libdavid : libs libutil libla
	( cd KS_Solvers/Davidson ; $(MAKE) TLDEPS= all || exit 1 )

libcg : libs libutil libla
	( cd KS_Solvers/CG ; $(MAKE) TLDEPS= all || exit 1 )

libla : liblapack libutil
	( cd LAXlib ; $(MAKE) TLDEPS= all || exit 1 )

libfft : 
	( cd FFTXlib ; $(MAKE) TLDEPS= all || exit 1 )

libutil : 
	( cd UtilXlib ; $(MAKE) TLDEPS= all || exit 1 )

libs :
	( cd clib ; $(MAKE) TLDEPS= all || exit 1 )

lrmods : libs libutil libla libfft
	( cd LR_Modules ; $(MAKE) TLDEPS= all || exit 1 )

dftd3 : mods
	( cd dft-d3 ; $(MAKE) TLDEPS= all || exit 1 )

bindir :
	test -d bin || mkdir bin

#############################################################
# Targets for external libraries
############################################################

libblas : 
	cd install ; $(MAKE) -f extlibs_makefile $@

liblapack: 
	cd install ; $(MAKE) -f extlibs_makefile $@

libiotk: 
	cd install ; $(MAKE) -f extlibs_makefile $@
libfox: 
	cd install ; $(MAKE) -f extlibs_makefile $@

# In case of trouble with iotk and compilers, add
# FFLAGS="$(FFLAGS_NOOPT)" after $(MFLAGS)

#########################################################
# plugins
#########################################################

w90: bindir liblapack
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

want : 
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

SaX : 
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

yambo: 
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

yambo-devel: 
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

plumed: 
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

west: pw
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

SternheimerGW: pw lrmods 
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

#########################################################
# "make links" produces links to all executables in bin/
#########################################################

# Contains workaround for name conflicts (dos.x and bands.x) with WANT
links : bindir
	( cd bin/ ; \
	rm -f *.x ; \
	for exe in ../*/*/*.x ../*/bin/* ; do \
	    if test ! -L $$exe ; then ln -fs $$exe . ; fi \
	done ; \
	[ -f ../WANT/wannier/dos.x ] && \
		ln -fs ../WANT/wannier/dos.x ../bin/dos_want.x ; \
	[ -f ../PP/src/dos.x ] &&  \
		ln -fs ../PP/src/dos.x ../bin/dos.x ; \
	[ -f ../WANT/wannier/bands.x ] && \
		ln -fs ../WANT/wannier/bands.x ../bin/bands_want.x ; \
	[ -f ../PP/src/dos.x ] &&  ln -fs ../PP/src/bands.x ../bin/bands.x ; \
	[ -f ../W90/wannier90.x ] &&  ln -fs ../W90/wannier90.x ../bin/wannier90.x ;\
	)

#########################################################
# 'make install' works based on --with-prefix
# - If the final directory does not exists it creates it
#########################################################

install : 
	@if test -d bin ; then mkdir -p $(PREFIX)/bin ; \
	for x in `find * ! -path "test-suite/*" -name *.x -type f` ; do \
		cp $$x $(PREFIX)/bin/ ; done ; \
	fi
	@echo 'Quantum ESPRESSO binaries installed in $(PREFIX)/bin'

#########################################################
# Run test-suite for numerical regression testing
# NB: it is assumed that reference outputs have been 
#     already computed once (usualy during release)
#########################################################

test-suite: pw cp 
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

#########################################################
# Other targets: clean up
#########################################################

# remove object files and executables
clean : 
	touch make.inc 
	for dir in \
		CPV LAXlib FFTXlib UtilXlib Modules PP PW EPW \
                KS_Solvers/CG KS_Solvers/Davidson KS_Solvers/Davidson_RCI \
		NEB ACFDT COUPLE GWW XSpectra PWCOND dft-d3 \
		atomic clib LR_Modules pwtools upftools \
		dev-tools extlibs Environ TDDFPT PHonon GWW \
	; do \
	    if test -d $$dir ; then \
		( cd $$dir ; \
		$(MAKE) TLDEPS= clean ) \
	    fi \
	done
	- @(cd install ; $(MAKE) -f plugins_makefile clean)
	- @(cd install ; $(MAKE) -f extlibs_makefile clean)
	- /bin/rm -rf bin/*.x tempdir

# remove files produced by "configure" as well
veryclean : clean
	- @(cd install ; $(MAKE) -f plugins_makefile veryclean)
	- @(cd install ; $(MAKE) -f extlibs_makefile veryclean)
	- rm -rf install/patch-plumed
	- cd install ; rm -f config.log configure.msg config.status \
		CPV/version.h ChangeLog* intel.pcl */intel.pcl
	- rm -rf include/configure.h install/make_wannier90.inc
	- cd install ; rm -fr autom4te.cache
	- cd pseudo; ./clean_ps ; cd -
	- cd install; ./clean.sh ; cd -
	- cd include; ./clean.sh ; cd -
	- rm -f espresso.tar.gz -
	- rm -rf make.inc -
	- rm -rf FoX
# remove everything not in the original distribution
distclean : veryclean
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

tar :
	@if test -f espresso.tar.gz ; then /bin/rm espresso.tar.gz ; fi
	# do not include unneeded stuff 
	find ./ -type f | grep -v -e /.svn/ -e'/\.' -e'\.o$$' -e'\.mod$$'\
		-e /.git/ -e'\.a$$' -e'\.d$$' -e'\.i$$' -e'_tmp\.f90$$' -e'\.x$$' \
		-e'~$$' -e'\./GUI' -e '\./tempdir' | xargs tar rvf espresso.tar
	gzip espresso.tar

#########################################################
# Tools for the developers
#########################################################
tar-gui :
	@if test -d GUI/PWgui ; then \
	    cd GUI/PWgui ; \
	    $(MAKE) TLDEPS= clean svninit pwgui-source; \
	    mv PWgui-*.tgz ../.. ; \
	else \
	    echo ; \
	    echo "  Sorry, tar-gui works only for svn sources !!!" ; \
	    echo ; \
	fi

tar-qe-modes :
	@if test -d GUI/QE-modes ; then \
	    cd GUI/QE-modes ; \
	    $(MAKE) TLDEPS= veryclean tar; \
	    mv QE-modes-*.tar.gz ../.. ; \
	else \
	    echo ; \
	    echo "  Sorry, tar-qe-modes works only for svn sources !!!" ; \
	    echo ; \
	fi

# NOTICE about "make doc": in order to build the .html and .txt
# documentation in Doc, "tcl", "tcllib", "xsltproc" are needed;
# in order to build the .pdf files in Doc, "pdflatex" is needed;
# in order to build html files for user guide and developer manual,
# "latex2html" and "convert" (from Image-Magick) are needed.
doc : 
	if test -d Doc ; then \
	( cd Doc ; $(MAKE) TLDEPS= all ) ; fi
	for dir in */Doc; do \
	( if test -f $$dir/Makefile ; then \
	( cd $$dir; $(MAKE) TLDEPS= all ) ; fi ) ;  done

doc_clean :
	if test -d Doc ; then \
	( cd Doc ; $(MAKE) TLDEPS= clean ) ; fi
	for dir in */Doc; do \
	( if test -f $$dir/Makefile ; then \
	( cd $$dir; $(MAKE) TLDEPS= clean ) ; fi ) ;  done

depend: libiotk
	@echo 'Checking dependencies...'
	- ( if test -x install/makedeps.sh ; then install/makedeps.sh ; fi)
