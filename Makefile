default :
	@echo 'to install, type at the shell prompt:'
	@echo '  ./configure'
	@echo '  make target'
	@echo 'where target is one of the following:'
	@echo '  pw           basic code for scf, structure optimization, MD'
	@echo '  cp           CP code: CP MD with ultrasoft pseudopotentials'
	@echo '  ph           phonon code'
	@echo '  pp           postprocessing programs'
	@echo '  gamma        Gamma-only version of phonon code'
	@echo '  pwcond       ballistic conductance'
	@echo '  d3           third-order derivatives'
	@echo '  vdw          vdW calculation'
	@echo '  gipaw        magnetic response (NMR, EPR, ...)'
	@echo '  tools        misc tools for data analysis'
	@echo '  ld1          utilities for pseudopotential generation'
	@echo '  upf          utilities for pseudopotential conversion'
	@echo '  pwall        same as "make pw ph pp gamma pwcond d3 tools"'
	@echo '  all          same as "make pwall cp ld1 upf"'
	@echo '  clean        remove executables and objects'
	@echo '  veryclean    revert distribution to the original status'
	@echo '  tar          create a tarball of the source tree'
	@echo '  tar-gui      create a tarball of the GUI sources'
	@echo '  log          create ChangeLog and ChangeLog.html files'

pw : bindir mods libs libiotk
	if test -d PW ; then \
	( cd PW ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

cp : bindir mods libs libiotk
	if test -d CPV ; then \
	( cd CPV ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= cp ; \
	else $(MAKE) $(MFLAGS) TLDEPS= cp ; fi ) ; fi

ph : bindir mods libs pw
	if test -d PH ; then \
	( cd PH ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
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

gipaw : bindir mods libs pw
	if test -d GIPAW ; then \
	( cd GIPAW ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

tools : bindir mods libs pw
	if test -d pwtools ; then \
	( cd pwtools ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

ld1 : bindir mods libs pw
	if test -d atomic ; then \
	( cd atomic ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

upf : mods libs
	if test -d upftools ; then \
	( cd upftools ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi ) ; fi

libiotk :
	if test -d iotk ; then \
	( cd iotk ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= libiotk.a ; \
	else $(MAKE) $(MFLAGS) TLDEPS= libiotk.a ; fi ) ; fi

pw_export : libiotk bindir mods libs pw
	if test -d PP ; then \
	( cd PP ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= pw_export.x ; \
	else $(MAKE) $(MFLAGS) TLDEPS= pw_export.x ; fi ) ; fi

pwall : pw ph pp gamma pwcond d3 vdw tools
all   : pwall cp ld1 upf 

mods : libiotk
	( cd Modules ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi )
libs : mods
	( cd clib ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi )
	( cd flib ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= all ; \
	else $(MAKE) $(MFLAGS) TLDEPS= all ; fi )
bindir :
	test -d bin || mkdir bin

# remove object files and executables
clean :
	touch make.sys 
	for dir in \
		CPV D3 Gamma Modules PH PP PW PWCOND VdW\
		atomic clib flib pwtools upftools iotk GIPAW \
	; do \
	    if test -d $$dir ; then \
		( cd $$dir ; \
		if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= clean ; \
		else $(MAKE) $(MFLAGS) TLDEPS= clean ; fi ) \
	    fi \
	done
	- /bin/rm -rf bin/*.x

# remove configuration files too
distclean veryclean : clean
	- /bin/rm -rf make.sys \
		      config.log configure.msg config.status autom4te.cache \
		      espresso.tar.gz CPV/version.h ChangeLog* \
		      intel.pcl */intel.pcl
	- cd examples ; ./make_clean
	- cd atomic_doc ; ./make_clean
	- if test -d GUI ; then \
	( cd GUI ; if test "$(MAKE)" = "" ; then make $(MFLAGS) TLDEPS= veryclean ; \
		else $(MAKE) $(MFLAGS) TLDEPS= veryclean ; fi ) \
	  fi

tar :
	tar cvf espresso.tar \
	    License README* */README* Makefile */Makefile */makefile.* */make.depend \
	    configure configure.ac config.guess config.sub configure.msg.in \
            install-sh make.sys.in \
	    makedeps.sh moduledep.sh includedep.sh ifcmods.sh \
	    */*.f90 */*.c */*.f clib/*.h include/*.h* upftools/UPF \
	    pwtools/*.awk pwtools/*.sh
	# remove unneeded stuff from iotk
	find iotk -type f | grep -v -e /CVS/ -e'\.o$$' -e'\.mod$$' -e'\.a$$' \
	    -e'\.d$$' -e'\.i$$' -e'\.F90$$' | xargs tar rvf espresso.tar
	# archive a few entire directories, but without CVS subdirs
	find install Doc atomic_doc examples pseudo -type f \
		| grep -v -e /CVS/ -e /results | xargs tar rvf espresso.tar
	gzip espresso.tar

# TAR-GUI works only if we have CVS-sources !!!
tar-gui :
	@if test -d GUI/PWgui ; then \
	    cd GUI/PWgui ; \
	    if test "$(MAKE)" = "" ; then \
		make $(MFLAGS) TLDEPS= clean cvsinit pwgui-source-notcl; \
	    else $(MAKE) $(MFLAGS) TLDEPS= clean cvsinit pwgui-source-notcl; fi; \
	    mv PWgui-*.tgz ../.. ; \
	else \
	    echo ; \
	    echo "  Sorry, tar-gui works only for CVS-sources !!!" ; \
	    echo ; \
	fi

log :
	-perl ./cvs2cl.pl
	-perl ./cvs2cl.pl --xml --header /dev/null --stdout \
		| perl ./cl2html.pl --entries 0 > ChangeLog.html

links : bindir
	( cd bin/ ; \
	for exe in \
	    ../CPV/cp.x \
	    ../D3/d3.x \
	    ../CPV/cppp.x \
	    ../Gamma/phcg.x \
	    ../PH/ph.x ../PH/dynmat.x ../PH/matdyn.x ../PH/q2r.x \
	    ../PP/average.x ../PP/bands.x ../PP/dos.x \
	      ../PP/efg.x ../PP/plotband.x ../PP/plotrho.x ../PP/pmw.x \
	      ../PP/pp.x ../PP/projwfc.x ../PP/pw2casino.x ../PP/pw2wan.x \
	      ../PP/voronoy.x ../PP/pw_export.x \
	    ../PW/pw.x \
	    ../PWCOND/pwcond.x \
	    ../atomic/ld1.x \
	    ../pwtools/band_plot.x ../pwtools/dist.x \
	      ../pwtools/ev.x ../pwtools/fqha.x ../pwtools/kpoints.x \
	      ../pwtools/path_int.x ../pwtools/pwi2xsf.x \
	; do \
	      if test -f $$exe ; then ln -fs $$exe . ; fi \
	done \
	)

depend:
	@echo 'Checking dependencies...'
	- ( if test -x ./makedeps.sh ; then ./makedeps.sh ; fi)
