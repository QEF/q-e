include make.sys

default :
	@echo 'to install, type at the shell prompt:'
	@echo '  ./configure'
	@echo '  make target'
	@echo 'where target is one of the following:'
	@echo '  pw           basic code for scf, structure optimization, MD'
	@echo '  cp           CP code: CP MD with ultrasoft pseudopotentials'
	@echo '  ph           phonon code'
	@echo '  neb          code for Nudged Elastic Band method'
	@echo '  tddfpt       time dependent dft code'
	@echo '  pp           postprocessing programs'
	@echo '  gamma        Gamma-only version of phonon code'
	@echo '  pwcond       ballistic conductance'
	@echo '  d3           third-order derivatives'
	@echo '  vdw          vdW calculation'
	@echo '  gipaw        magnetic response (NMR, EPR, ...)'
	@echo '  w90          Maximally localised Wannier Functions'
	@echo '  want         Quantum Transport with Wannier functions'
	@echo '  plumed       Patch for calculating free-energy paths with pw or cp'
	@echo '  gww          GW with Wannier Functions'
	@echo '  yambo        electronic excitations with plane waves'
	@echo '  tools        misc tools for data analysis'
	@echo '  ld1          utilities for pseudopotential generation'
	@echo '  upf          utilities for pseudopotential conversion'
	@echo '  xspectra     X-ray core-hole spectroscopy calculations '
	@echo '  pwall        same as "make pw ph pp gamma pwcond d3 tools"'
	@echo '  all          same as "make pwall cp ld1 upf tddfpt"'
	@echo '  clean        remove executables and objects'
	@echo '  veryclean    revert distribution to the original status'
	@echo '  tar          create a tarball of the source tree'
	@echo '  tar-gui      create a tarball of the GUI sources'
	@echo '  doc          build documentation'
	@echo '  log          create ChangeLog and ChangeLog.html files'
	@echo '  links        create links to all executables in bin/'

###########################################################
# Main targets
###########################################################
pw : bindir mods liblapack libblas libs libiotk
	if test -d PW ; then \
	( cd PW ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

cp : bindir mods liblapack libblas libs libiotk
	if test -d CPV ; then \
	( cd CPV ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= cp ; \
	else $(MAKE) $(MFLAGS) TLDEPS= cp ; fi ) ; fi

ph : bindir mods libs pw
	if test -d PH ; then \
	( cd PH ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

neb : bindir mods libs pw
	if test -d NEB ; then \
	( cd NEB ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

tddfpt : bindir mods libs pw ph
	if test -d TDDFPT ; then \
	( cd TDDFPT ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

pp : bindir mods libs pw
	if test -d PP ; then \
	( cd PP ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

gamma : bindir mods libs pw
	if test -d Gamma ; then \
	( cd Gamma ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

pwcond : bindir mods libs pw pp
	if test -d PWCOND ; then \
	( cd PWCOND ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

d3 : bindir mods libs pw ph
	if test -d D3 ; then \
	( cd D3 ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

vdw : bindir mods libs pw ph pp
	if test -d VdW ; then \
	( cd VdW ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

acfdt : bindir mods libs pw ph
	if test -d ACFDT ; then \
	( cd ACFDT ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

gipaw : bindir mods libs pw
	if test -d GIPAW ; then \
	( cd GIPAW ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

gww   : bindir pw ph
	if test -d GWW ; then \
	( cd GWW ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

tools : bindir mods libs pw
	if test -d pwtools ; then \
	( cd pwtools ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

ld1 : bindir mods libs
	if test -d atomic ; then \
	( cd atomic ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

upf : mods libs
	if test -d upftools ; then \
	( cd upftools ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

pw_export : libiotk bindir mods libs pw
	if test -d PP ; then \
	( cd PP ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= pw_export.x ; \
	else $(MAKE) $(MFLAGS) TLDEPS= pw_export.x ; fi ) ; fi

xspectra : bindir mods libs pw gipaw
	if test -d XSpectra ; then \
	( cd XSpectra ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

pwall : pw neb ph pp gamma pwcond d3 vdw tools acfdt
all   : pwall cp ld1 upf gww tddfpt

###########################################################
# Auxiliary targets used by main targets:
# compile modules, libraries, directory for binaries, etc
###########################################################
mods : libiotk
	( cd Modules ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi )
libs : mods
	( cd clib ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi )
	( cd flib ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= $(FLIB_TARGETS) ; \
	else $(MAKE) $(MFLAGS) TLDEPS= $(FLIB_TARGETS) ; fi )

eelib : mods
	( cd EE ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all  ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi )
bindir :
	test -d bin || mkdir bin

#############################################################
# Targets for external libraries
############################################################
libblas:
	if test -e extlibs/archive/blas-1.tar.gz ; then \
	( cd extlibs ; $(MAKE) $(MFLAGS) $@) ; fi

liblapack:
	if test -e extlibs/archive/lapack-3.2.tar.gz ; then \
	( cd extlibs ; $(MAKE) $(MFLAGS) $@) ; fi

mglib:
	if test -e extlibs/archive/multigrid.tar ; then \
	( cd extlibs ; $(MAKE) $(MFLAGS) $@) ; fi

libiotk:
	if test -e extlibs/archive/iotk-1.1.beta.tar.gz ; then \
	( cd extlibs ; $(MAKE) $(MFLAGS) $@) ; fi
# In case of trouble with iotk and compilers, replace with the following
#	if test -e extlibs/archive/iotk-1.1.beta.tar.gz ; then \
#	( cd extlibs ; $(MAKE) $(MFLAGS) FFLAGS="$(FFLAGS_NOOPT)" $@) ; fi

#########################################################
# plugins
#########################################################

w90: bindir libblas liblapack
	cd install ; $(MAKE) $(MFLAGS) -f plugins_makefile $@

want : touch-dummy
	cd install ; $(MAKE) $(MFLAGS) -f plugins_makefile $@

yambo: touch-dummy
	cd install ; $(MAKE) $(MFLAGS) -f plugins_makefile $@

plumed: touch-dummy
	cd install ; $(MAKE) $(MFLAGS) -f plugins_makefile $@

touch-dummy :
	$(dummy-variable)

#########################################################
# Links and copies. "make links" is likely obsolete.
# "make inst INSTALLDIR=/some/place" will copy all 
# available executables to /some/place/ (must exist and
# be writable), prepending "qe_" to all executables (e.g.:
# /some/place/qe_pw.x). This allows installation of QE
# into system directories with no danger of name conflicts
#########################################################
inst : 
	( for exe in */*.x ../GWW/*.x ; do \
	   file=`basename $$exe`; if test "$(INSTALLDIR)" != ""; then \
		cp $(PWD)/$$exe $(INSTALLDIR)/qe_$$file ; fi ; \
	done )

links : bindir
	( cd bin/ ; \
	rm -f *.x ; \
	for exe in ../*/*.x ; do \
	    if test ! -L $$exe ; then ln -fs $$exe . ; fi \
	done \
	)


#########################################################
# Other targets: clean up
#########################################################

# remove object files and executables
clean :
	touch make.sys 
	for dir in \
		CPV D3 Gamma Modules PH PP PW PWCOND VdW ACFDT EE \
		atomic clib flib pwtools upftools iotk GIPAW XSpectra \
		dev-tools GWW extlibs TDDFPT \
	; do \
	    if test -d $$dir ; then \
		( cd $$dir ; \
		if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= clean ; \
		else $(MAKE) $(MFLAGS) TLDEPS= clean ; fi ) \
	    fi \
	done
	- cd install ; $(MAKE) $(MFLAGS) --file=plugins_makefile clean
	- /bin/rm -rf bin/*.x tmp
	- cd tests; /bin/rm -rf CRASH *.out *.out2 

# remove configuration files too
distclean veryclean : clean
	- test -d extlibs && ( cd extlibs ; $(MAKE) veryclean)
	- rm -rf make.sys
	- rm -rf install/patch-plumed
	- cd install ; rm -f config.log configure.msg config.status autom4te.cache \
	CPV/version.h ChangeLog* intel.pcl */intel.pcl
	- rm -f espresso.tar.gz
	- cd examples ; ./make_clean
	- cd atomic_doc ; ./make_clean
	- for dir in Doc doc-def; do \
	    test -d $$dir && ( cd $$dir ; $(MAKE) $(MFLAGS) TLDEPS= clean ) \
	done
	- test -d GUI && ( cd GUI ;  $(MAKE) $(MFLAGS) TLDEPS= veryclean )

tar :
	@if test -f espresso.tar.gz ; then /bin/rm espresso.tar.gz ; fi
	# do not include unneeded stuff 
	find ./ -type f | grep -v -e /CVS/ -e /results/ -e'/\.' -e'\.o$$' \
             -e'\.mod$$' -e'\.a$$' -e'\.d$$' -e'\.i$$' -e'\.F90$$' -e'\.x$$' \
	     -e'~$$' -e'\./GUI' | xargs tar rvf espresso.tar
	gzip espresso.tar

#########################################################
# Tools for the developers
#########################################################
# TAR-GUI works only if we have CVS-sources !!!
tar-gui :
	@if test -d GUI/PWgui ; then \
	    cd GUI/PWgui ; \
	    if test "$(MAKE)" = "" ; then \
		make $(MFLAGS) TLDEPS= clean cvsinit pwgui-source; \
	    else $(MAKE) $(MFLAGS) TLDEPS= clean cvsinit pwgui-source; fi; \
	    mv PWgui-*.tgz ../.. ; \
	else \
	    echo ; \
	    echo "  Sorry, tar-gui works only for CVS-sources !!!" ; \
	    echo ; \
	fi

# NOTICE about "make doc": in order to build the .html and .txt
# documentation in Doc, "tcl", "tcllib", "xsltproc" are needed;
# in order to build the .pdf files in Doc, "pdflatex" is needed;
# in order to build html files for user guide and developer manual,
# "latex2html" and "convert" (from Image-Magick) are needed.
doc :
	if test -d Doc ; then \
	( cd Doc ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi
	if test -d doc-def ; then \
	( cd doc-def ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi


# produce a log file of all changes - will work only if you have read
# access to the cvs server. For usage by developers
log :
	-perl ./cvs2cl.pl
	-perl ./cvs2cl.pl --xml --header /dev/null --stdout \
		| perl ./cl2html.pl --entries 0 > ChangeLog.html

depend:
	@echo 'Checking dependencies...'
	- ( if test -x install/makedeps.sh ; then install/makedeps.sh ; fi)

