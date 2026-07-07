# Copyright (C) 2001-2025 Quantum ESPRESSO Foundation
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

-include make.inc

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
	@echo '               (compiles w90 as well)'
	@echo '  kcw          KCW code: implementation of Koopmans functionals in primitive cell'
	@echo '  pioud        Path Integral Molecular Dynamics with PIOUD algorithm'
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

# The syntax "( cd PW ; $(MAKE) all || exit 1)" below
# guarantees that error code 1 is returned in case of error and make stops
# If "|| exit 1" is not present, the error code from make in subdirectories
# is not returned and make goes on even if compilation has failed

pw : pwlibs
	if test -d PW ; then \
	( cd PW ; $(MAKE) all || exit 1) ; fi

cp : bindir mods
	if test -d CPV ; then \
	( cd CPV ; $(MAKE) all || exit 1) ; fi

ph : phlibs
	if test -d PHonon; then \
	( cd PHonon; $(MAKE) all || exit 1) ; fi

hp : hplibs
	if test -d HP; then \
	( cd HP; $(MAKE) all || exit 1) ; fi

neb : pwlibs
	if test -d NEB; then \
	( cd NEB; $(MAKE) all || exit 1) ; fi

tddfpt : lrmods
	if test -d TDDFPT; then \
	( cd TDDFPT; $(MAKE) all || exit 1) ; fi

pp : pwlibs
	if test -d PP ; then \
	( cd PP ; $(MAKE) all || exit 1 ) ; fi

pwcond : pwlibs
	if test -d PWCOND ; then \
	( cd PWCOND ; $(MAKE) all || exit 1 ) ; fi

acfdt : phlibs
	if test -d ACFDT ; then \
	( cd ACFDT ; $(MAKE) all || exit 1 ) ; fi

gwl : phlibs
	if test -d GWW ; then \
	( cd GWW ; $(MAKE) all || exit 1 ) ; fi

gipaw : pwlibs
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

d3q : phlibs
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

ld1 : bindir mods
	if test -d atomic ; then \
	( cd atomic ; $(MAKE) all || exit 1 ) ; fi

xspectra : pwlibs
	if test -d XSpectra ; then \
	( cd XSpectra ; $(MAKE) all || exit 1 ) ; fi

couple : pw cp
	if test -d COUPLE ; then \
	( cd COUPLE ; $(MAKE) all || exit 1 ) ; fi

epw: pw ph pp ld1 libw90 
	if test -d EPW ; then \
	( cd EPW ; $(MAKE) all || exit 1; \
		cd ../bin; ln -fs ../EPW/bin/epw.x . ); fi

all_currents:
	if test -d QEHeat ; then \
	( cd QEHeat ; $(MAKE) all || exit 1; ) ; fi

travis : pwall epw
	if test ! $(TOPDIR) -ef $(BUILDDIR) ; then \
	   cp -r $(TOPDIR)/test-suite $(BUILDDIR) ; \
	   cp -r $(TOPDIR)/pseudo $(BUILDDIR) ; fi
	if test -d test-suite ; then \
	( cd test-suite ; make run-travis || exit 1 ) ; fi

kcw : pwlibs lrmods kcwlib
	if test -d KCW ; then \
	( cd KCW ; $(MAKE) all || exit 1 ) ; fi

pioud : pw pwlibs 
	if test -d PIOUD ; then \
	( cd PIOUD ; $(MAKE) all || exit 1 ) ; fi

gui : bindir
	@if test ! $(TOPDIR) -ef $(BUILDDIR) ; then \
	   echo "make $@ not supported in out-of-source builds" ; \
	else \
	   if test -d GUI/PWgui ; then \
	       cd GUI/PWgui ; \
	       $(MAKE) init; \
	       echo ; \
	       echo "  ------------------------------------------------------------"; \
	       echo "  PWgui was built in ./GUI/PWgui/ and a link was made in bin/."; \
	       echo "  ------------------------------------------------------------"; \
	       echo "  Try it either as:  "; \
	       echo "         ./GUI/PWgui/pwgui" ; \
	       echo "     or"; \
	       echo "         ./bin/pwgui";\
	       echo ; \
	   else \
	       echo ; \
	       echo "  Sorry, GUI/PWgui directory does not exist !" ; \
	       echo ; \
	   fi ; \
	fi

pwall : pw neb ph pp pwcond acfdt

all   : pwall cp ld1 tddfpt hp xspectra gwl kcw pioud

###########################################################
# Auxiliary targets used by main targets:
# compile modules, libraries, directory for binaries, etc
###########################################################

pwlibs: bindir mods libks_solvers dftd3
	if test -d PW ; then \
	( cd PW ; $(MAKE) pwlibs || exit 1) ; fi

phlibs: pwlibs lrmods
	if test -d PHonon; then \
	( cd PHonon; $(MAKE) phlibs || exit 1) ; fi

hplibs: pwlibs lrmods
	if test -d HP; then \
	( cd HP; $(MAKE) hplibs || exit 1) ; fi

gwwlib : phlibs
	if test -d GWW ; then \
	( cd GWW ; $(MAKE) gwwa || exit 1 ) ; fi

kcwlib : pwlibs lrmods
	if test -d KCW ; then \
	( cd KCW ; $(MAKE) kcwlib || exit 1 ) ; fi

pw4gwwlib : phlibs
	if test -d GWW ; then \
	( cd GWW ; $(MAKE) pw4gwwa || exit 1 ) ; fi

mods : $(FOX) libutil libla libfft libupf libmbd librxc
	( cd Modules ; $(MAKE) all || exit 1 )

libks_solvers : libutil libla
	( cd KS_Solvers ; $(MAKE) all || exit 1 )

libla : $(LAPACK) libutil libdevx
	( cd LAXlib ; $(MAKE) all || exit 1 )

libfft : 
	( cd FFTXlib ; $(MAKE) all || exit 1 )

librxc : 
	( cd XClib ; $(MAKE) all || exit 1 )

libutil : 
	( cd UtilXlib ; $(MAKE) all || exit 1 )

libupf : libutil libdevx
	( cd upflib ; $(MAKE) all || exit 1 )

lrmods : mods pwlibs
	( cd LR_Modules ; $(MAKE) all || exit 1 )

dftd3 : mods
	( cd dft-d3 ; $(MAKE) all || exit 1 )

bindir :
	test -d bin || mkdir bin

#############################################################
# Targets for external libraries
############################################################

libdevx:
	( cd install ; $(MAKE) -f extlibs_makefile $@ || exit 1 )

libmbd:
	( cd install ; $(MAKE) -f extlibs_makefile $@ || exit 1 )

libw90:
	( cd install ; $(MAKE) -f extlibs_makefile $@ || exit 1 )

# next two targets are obsolescent if not obsolete
liblapack: 
	cd install ; $(MAKE) -f oldlibs_makefile $@

libfox: 
	cd install ; $(MAKE) -f oldlibs_makefile $@

#########################################################
# plugins
#########################################################

w90: libw90

want: $(LAPACK)
	( cd install ; $(MAKE) -f plugins_makefile $@ || exit 1 )

yambo: $(LAPACK)
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
		LAXlib FFTXlib XClib UtilXlib upflib Modules KS_Solvers \
		dft-d3 LR_Modules PW CPV PP PHonon HP EPW NEB TDDFPT GWW \
		XSpectra PWCOND atomic QEHeat KCW PIOUD COUPLE Doc GUI \
		dev-tools ACFDT Environ \
	; do \
	    if test -d $$dir ; then \
		( cd $$dir ; \
		$(MAKE) clean ) \
	    fi \
	done
	- @(cd install ; $(MAKE) -f plugins_makefile clean)
	- @(cd install ; $(MAKE) -f extlibs_makefile clean)
	- /bin/rm -rf bin/*.x tempdir

# remove files produced by "configure" as well
veryclean : clean
	-@if test ! $(TOPDIR) -ef $(BUILDDIR) ; then \
	   echo "make $@ not supported in out-of-source builds" ; \
	   echo "just re-create $(BUILDDIR) and re-run configure" ; \
	else \
	- @(cd install ; $(MAKE) -f plugins_makefile veryclean) ; \
	- (cd install ; rm -rf config.log configure.msg config.status \
		make_wannier90.inc autom4te.cache ) ; \
	- rm -f espresso.tar.gz ; \
	- rm -rf make.inc ; \
	- rm -rf MBD wannier90 devxlib ;\
	- rm -rf FoX lapack ; \
	fi
# remove everything not in the original distribution
# place deinit at the very end such that makefiles clean up as much as possible.
distclean : veryclean
	-@if test ! $(TOPDIR) -ef $(BUILDDIR) ; then \
	   echo "make $@ not supported in out-of-source builds" ; \
	else \
		cd pseudo; ./clean_ps ; cd - ;\
		(cd install ; $(MAKE) -f extlibs_makefile $@) ;\
		(cd install ; $(MAKE) -f plugins_makefile $@) ;\
		git submodule deinit --all --force  ;\
	fi

# find line: do not include unneeded stuff  
tar :
	@if test ! $(TOPDIR) -ef $(BUILDDIR) ; then \
		echo "make $@ not supported in out-of-source builds" ; \
	else \
		if test -f espresso.tar.gz ; then /bin/rm espresso.tar.gz ; fi ;\
		find ./ -type f | grep -v -e /.svn/ -e'/\.' -e'\.o$$' -e'\.mod$$'\
			-e /.git/ -e'\.a$$' -e'\.d$$' -e'\.i$$' -e'_tmp\.f90$$' -e'\.x$$' \
			-e'~$$' -e'\./GUI' -e '\./tempdir' | xargs tar rvf espresso.tar ;\
		gzip espresso.tar ;\
	fi

#########################################################
# Tools for the developers
#########################################################
tar-gui :
	@if test ! $(TOPDIR) -ef $(BUILDDIR) ; then \
	   echo "make $@ not supported in out-of-source builds" ; \
	else \
		if test -d GUI/PWgui ; then \
		    cd GUI/PWgui ; \
		    $(MAKE) clean init pwgui-source; \
		    mv PWgui-*.tgz ../.. ; \
		else \
		    echo ; \
		    echo "  Sorry, tar-gui works only for git sources !!!" ; \
		    echo ; \
		fi ;\
	fi

tar-qe-modes :
	@if test ! $(TOPDIR) -ef $(BUILDDIR) ; then \
	   echo "make $@ not supported in out-of-source builds" ; \
	else \
		if test -d GUI/QE-modes ; then \
		    cd GUI/QE-modes ; \
		    $(MAKE) veryclean tar; \
		    mv QE-modes-*.tar.gz ../.. ; \
		else \
		    echo ; \
		    echo "  Sorry, tar-qe-modes works only for git sources !!!" ; \
		    echo ; \
		fi ;\
	fi

# NOTICE about "make doc": in order to build the .html and .txt
# documentation in Doc, "tcl", "tcllib", "xsltproc" are needed;
# in order to build the .pdf files in Doc, "pdflatex" is needed;
# in order to build html files for the user guide,
# "latex2html" and "convert" (from Image-Magick) are needed.
doc : 
	if test ! $(TOPDIR) -ef $(BUILDDIR) ; then \
	   echo "make $@ not supported in out-of-source builds" ; \
	else \
	   if test -d Doc ; then \
	   ( cd Doc ; $(MAKE) all ) ; fi ;\
	   for dir in */Doc; do \
	   ( if test -f $$dir/Makefile ; then \
	   ( cd $$dir; $(MAKE) all ) ; fi ) ;  done ; \
	fi

doc_clean :
	if test -d Doc ; then \
	( cd Doc ; $(MAKE) clean ) ; fi
	for dir in */Doc; do \
	( if test -f $$dir/Makefile ; then \
	( cd $$dir; $(MAKE) clean ) ; fi ) ;  done

depend:
	echo 'Checking dependencies...'
	-@if test ! $(TOPDIR) -ef $(BUILDDIR) ; then \
	    $(TOPDIR)/install/makedeps.sh $(BUILDDIR) ; \
	else \
	    install/makedeps.sh ; \
	fi
