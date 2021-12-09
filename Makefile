# Copyright (C) 2001-2020 Quantum ESPRESSO Foundation
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

include make.inc

# execute a target irrespective of the presence of a file or directory 
# with the same name
.PHONY: install

default :
	@echo 'to install Quantum ESPRESSO, type at the shell prompt:'
	@echo '  ./configure [--prefix=]'
	@echo '  make [-j] target'
	@echo ' '
	@echo 'where target identifies one or multiple CORE PACKAGES:'
	@echo '  pw           basic code for scf, structure optimization, MD'
	@echo '  ph           phonon code, Gamma-only and third-order derivatives'
	@echo '  hp           calculation of the Hubbard parameters from DFPT'
	@echo '  pwcond       ballistic conductance'
	@echo '  neb          code for Nudged Elastic Band method'
	@echo '  pp           postprocessing programs'
	@echo '  pwall        same as "make pw ph pp pwcond neb"'
	@echo '  cp           CP code: Car-Parrinello molecular dynamics'
	@echo '  all_currents QEHeat code: energy flux and charge current'
	@echo '  tddfpt       time dependent dft code'
	@echo '  gwl          GW with Lanczos chains'
	@echo '  ld1          utilities for pseudopotential generation'
	@echo '  xspectra     X-ray core-hole spectroscopy calculations'
	@echo '  couple       Library interface for coupling to external codes'
	@echo '  epw          Electron-Phonon Coupling with Wannier functions'
	@echo '  gui          Graphical User Interface'
	@echo '  all          same as "make pwall cp ld1 tddfpt xspectra hp"'
	@echo ' '
	@echo 'where target identifies one or multiple THIRD-PARTIES PACKAGES:'
	@echo '  gipaw        NMR and EPR spectra'
	@echo '  w90          Maximally localised Wannier Functions'
	@echo '  want         Quantum Transport with Wannier functions'
	@echo '  yambo        electronic excitations with plane waves'
	@echo '  d3q          general third-order code and thermal transport codes'
	@echo ' '
	@echo 'where target is one of the following suite operation:'
	@echo '  doc          build documentation'
	@echo '  links        create links to all executables in bin/'
	@echo '  install      copy all executables to PREFIX/bin/'
	@echo '               (works with "configure --prefix=PREFIX)"'
	@echo '  tar          create a tarball of the source tree'
	@echo '  depend       generate dependencies (make.depend files)'
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

pw : pwlibs
	if test -d PW ; then \
	( cd PW ; $(MAKE) TLDEPS= all || exit 1) ; fi

cp : bindir mods
	if test -d CPV ; then \
	( cd CPV ; $(MAKE) TLDEPS= all || exit 1) ; fi

ph : phlibs
	if test -d PHonon; then \
	( cd PHonon; $(MAKE) TLDEPS= all || exit 1) ; fi

hp : hplibs
	if test -d HP; then \
	( cd HP; $(MAKE) TLDEPS= all || exit 1) ; fi

neb : pwlibs
	if test -d NEB; then \
	( cd NEB; $(MAKE) TLDEPS= all || exit 1) ; fi

tddfpt : lrmods
	if test -d TDDFPT; then \
	( cd TDDFPT; $(MAKE) TLDEPS= all || exit 1) ; fi

pp : pwlibs
	if test -d PP ; then \
	( cd PP ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

pwcond : pwlibs
	if test -d PWCOND ; then \
	( cd PWCOND ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

acfdt : phlibs
	if test -d ACFDT ; then \
	( cd ACFDT ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

gwl : phlibs
	if test -d GWW ; then \
	( cd GWW ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

gipaw : pwlibs
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

d3q : phlibs
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

ld1 : bindir mods
	if test -d atomic ; then \
	( cd atomic ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

xspectra : pwlibs
	if test -d XSpectra ; then \
	( cd XSpectra ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

couple : pw cp
	if test -d COUPLE ; then \
	( cd COUPLE ; $(MAKE) TLDEPS= all || exit 1 ) ; fi

epw: phlibs
	if test -d EPW ; then \
	( cd EPW ; $(MAKE) all || exit 1; \
		cd ../bin; ln -fs ../EPW/bin/epw.x . ); fi

all_currents:
	if test -d QEHeat ; then \
	( cd QEHeat ; $(MAKE) all || exit 1; ) ; fi

travis : pwall epw
	if test -d test-suite ; then \
	( cd test-suite ; make run-travis || exit 1 ) ; fi

gui :
	@if test -d GUI/PWgui ; then \
	    cd GUI/PWgui ; \
	    $(MAKE) TLDEPS= init; \
	    echo ; \
	    echo "  PWgui has been built in ./GUI/PWgui/. You may try it either as:  "; \
	    echo "         ./GUI/PWgui/pwgui" ; \
	    echo "     or"; \
	    echo "         cd ./GUI/PWgui";\
	    echo "         ./pwgui" ; \
	    echo ; \
	else \
	    echo ; \
	    echo "  Sorry, gui works only for git sources !!!" ; \
	    echo ; \
	fi

pwall : pw neb ph pp pwcond acfdt

all   : pwall cp ld1 tddfpt hp xspectra gwl 

###########################################################
# Auxiliary targets used by main targets:
# compile modules, libraries, directory for binaries, etc
###########################################################

pwlibs: bindir mods libks_solvers dftd3
	if test -d PW ; then \
	( cd PW ; $(MAKE) pw-lib || exit 1) ; fi

phlibs: pwlibs lrmods
	if test -d PHonon; then \
	( cd PHonon; $(MAKE) ph-lib || exit 1) ; fi

hplibs: pwlibs lrmods
	if test -d HP; then \
	( cd HP; $(MAKE) hp-lib || exit 1) ; fi

gwwlib : phlibs
	if test -d GWW ; then \
	( cd GWW ; $(MAKE) gwwa || exit 1 ) ; fi

pw4gwwlib : phlibs
	if test -d GWW ; then \
	( cd GWW ; $(MAKE) pw4gwwa || exit 1 ) ; fi

mods : libfox libutil libla libfft libupf libmbd librxc
	( cd Modules ; $(MAKE) TLDEPS= all || exit 1 )

libks_solvers : libutil libla
	( cd KS_Solvers ; $(MAKE) TLDEPS= all || exit 1 )

libla : liblapack libutil libcuda
	( cd LAXlib ; $(MAKE) TLDEPS= all || exit 1 )

libfft : 
	( cd FFTXlib ; $(MAKE) TLDEPS= all || exit 1 )

librxc : 
	( cd XClib ; $(MAKE) TLDEPS= all || exit 1 )

libutil : 
	( cd UtilXlib ; $(MAKE) TLDEPS= all || exit 1 )

libupf : libutil libcuda
	( cd upflib ; $(MAKE) TLDEPS= all || exit 1 )

lrmods : mods pwlibs
	( cd LR_Modules ; $(MAKE) TLDEPS= all || exit 1 )

dftd3 : mods
	( cd dft-d3 ; $(MAKE) TLDEPS= all || exit 1 )

bindir :
	test -d bin || mkdir bin

#############################################################
# Targets for external libraries
############################################################

liblapack: 
	cd install ; $(MAKE) -f extlibs_makefile $@

libfox: 
	cd install ; $(MAKE) -f extlibs_makefile $@

libcuda: 
	cd install ; $(MAKE) -f extlibs_makefile $@

libmbd:
	cd install ; $(MAKE) -f extlibs_makefile $@

#########################################################
# plugins
#########################################################

w90: bindir liblapack
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

want: liblapack
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

yambo: liblapack
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

#############################################################
# 'make install' works with "configure --prefix=PREFIX"
# - If the PREFIX/bin directory does not exists it creates it
#############################################################

install : 
	mkdir -p $(PREFIX)/bin ; \
	for x in `find * ! -path "test-suite/*" -name *.x -type f` ; do \
		cp -v $$x $(PREFIX)/bin/ ; done
	@echo -e '\nQuantum ESPRESSO binaries are installed in $(PREFIX)/bin\n'

#########################################################
# Other targets: clean up
#########################################################

# remove object files and executables
clean : 
	touch make.inc 
	for dir in \
		CPV LAXlib FFTXlib XClib UtilXlib upflib Modules PP PW EPW KS_Solvers \
		NEB ACFDT COUPLE GWW XSpectra PWCOND dft-d3 \
		atomic LR_Modules upflib \
		dev-tools extlibs Environ TDDFPT PHonon HP GWW Doc GUI \
		QEHeat \
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
	- (cd install ; rm -rf config.log configure.msg config.status \
		configure.h make_wannier90.inc autom4te.cache )
	- (cd include; rm -rf configure.in qe_cdefs.h )
	- rm -f espresso.tar.gz
	- rm -rf make.inc
	- rm -rf FoX
	- rm -rf MBD 
# remove everything not in the original distribution
distclean : veryclean
	- cd pseudo; ./clean_ps ; cd -
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
	    $(MAKE) TLDEPS= clean init pwgui-source; \
	    mv PWgui-*.tgz ../.. ; \
	else \
	    echo ; \
	    echo "  Sorry, tar-gui works only for git sources !!!" ; \
	    echo ; \
	fi

tar-qe-modes :
	@if test -d GUI/QE-modes ; then \
	    cd GUI/QE-modes ; \
	    $(MAKE) TLDEPS= veryclean tar; \
	    mv QE-modes-*.tar.gz ../.. ; \
	else \
	    echo ; \
	    echo "  Sorry, tar-qe-modes works only for git sources !!!" ; \
	    echo ; \
	fi

# NOTICE about "make doc": in order to build the .html and .txt
# documentation in Doc, "tcl", "tcllib", "xsltproc" are needed;
# in order to build the .pdf files in Doc, "pdflatex" is needed;
# in order to build html files for the user guide,
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

depend:
	@echo 'Checking dependencies...'
	- ( if test -x install/makedeps.sh ; then install/makedeps.sh ; fi)
