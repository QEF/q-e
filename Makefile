default :
	@echo 'to install, type at the shell prompt:'
	@echo '  ./configure'
	@echo '  make target'
	@echo 'where target is one of the following:'
	@echo '  pw           basic code for scf, structure optimization, MD'
	@echo '  fpmd         FPMD code for Car-Parrinello MD'
	@echo '  cp           CP code: CP MD with ultrasoft pseudopotentials'
	@echo '  ph           phonon code'
	@echo '  pp           postprocessing programs'
	@echo '  gamma        Gamma-only version of phonon code'
	@echo '  nc           non collinear magnetic version of pw code'
	@echo '  pwcond       ballistic conductance'
	@echo '  d3           third-order derivatives'
	@echo '  tools        misc tools for data analysis'
	@echo '  upf          utilities for pseudopotential conversion'
	@echo '  pwall        same as "make pw ph pp gamma nc pwcond d3 tools"'
	@echo '  all          same as "make pwall fpmd cp upf"'
	@echo '  links        creates links to executables in bin/'
	@echo '  clean        remove executables and objects'
	@echo '  veryclean    revert distribution to the original status'
	@echo '  tar          create a tarball of the source tree'
	@echo '  tar-gui      create a tarball of the GUI sources'

pw : modules libs
	if test -d PW   ; then ( cd PW   ; make all ) ; fi
fpmd : modules libs
	if test -d FPMD ; then ( cd FPMD ; make all ) ; fi
cp : modules libs
	if test -d CPV  ; then ( cd CPV  ; make all ) ; fi

ph : modules libs pw
	if test -d PH ; then ( cd PH ; make all ) ; fi
pp : modules libs pw
	if test -d PP ; then ( cd PP ; make all ) ; fi
gamma : modules libs pw
	if test -d Gamma  ; then ( cd Gamma  ; make all ) ; fi
nc : modules libs pw
	if test -d PWNC   ; then ( cd PWNC   ; make all ) ; fi
pwcond : modules libs pw pp
	if test -d PWCOND ; then ( cd PWCOND ; make all ) ; fi
d3 : modules libs pw ph
	if test -d D3 ; then ( cd D3 ; make all ) ; fi

tools : modules libs pw
	if test -d pwtools  ; then ( cd pwtools  ; make all ) ; fi
upf : modules libs
	if test -d upftools ; then ( cd upftools ; make all ) ; fi

pwall : pw ph pp gamma nc pwcond d3 tools
all   : pwall fpmd cp upf 

modules :
	( cd Modules; make all )
libs : modules
	( cd clib ; make all )
	( cd flib ; make all )

# create link only if file exists
links :
	test -d bin || mkdir bin
	( cd bin/ ; \
	  for exe in \
	      ../PW/pw.x ../PW/memory.x ../PH/ph.x ../D3/d3.x \
	      ../Gamma/phcg.x ../CPV/cp.x ../FPMD/par2.x \
	      ../PP/average.x ../PP/bands.x ../PP/chdens.x ../PP/dos.x \
	      ../PP/plotrho.x ../PP/pp.x ../PP/projwfc.x ../PP/voronoy.x \
	      ../PP/plotband.x ../PWNC/pwnc.x ../PWCOND/pwcond.x \
	      ../pwtools/band_plot.x ../pwtools/dynmat.x ../pwtools/fqha.x \
	      ../pwtools/matdyn.x ../pwtools/q2r.x ../pwtools/dist.x \
	      ../pwtools/ev.x ../pwtools/kpoints.x ../pwtools/path_int.x \
	  ; do \
	      if test -f $$exe ; then ln -fs $$exe . ; fi \
	  done \
	)

# remove object files and executables
clean :
	touch make.rules make.sys # make complains if they aren't there
	#                         # same with .dependencies below
	for dir in PW PWNC PH PP D3 PWCOND Gamma pwtools upftools \
		   Modules install clib flib FPMD CPV ; do \
	    if test -d $$dir ; then \
		( cd $$dir ; touch .dependencies ; make clean ) \
	    fi \
	done

# remove configuration files too
veryclean : clean
	- /bin/rm -rf make.rules make.sys */.dependencies \
		      config.log config.status bin/*.x \
		      autom4te.cache pw.tar.gz FPMD/version.h \
		      intel.pcl */intel.pcl
	- if test -d GUI ; then ( cd GUI; make veryclean ) ; fi

tar :
	tar cvf pw.tar \
	    License README */README README.cvs README.configure \
            INSTALL Makefile */Makefile \
	    configure configure.ac config.guess config.sub install-sh \
	    makedeps.sh moduledep.sh make.rules.in make.sys.in configure.old \
	    */*.f90 */*.c */*.f clib/*.h include/*.h* upftools/UPF \
	    pwtools/*.awk pwtools/*.sh
	# archive a few entire directories, but without CVS subdirs
	find install *docs *_examples pseudo -type f \
		| grep -v -e /CVS/ -e /results/ | xargs tar rvf pw.tar
	gzip pw.tar


# TAR-GUI works only if we have CVS-sources !!!
tar-gui :
	@if test -d GUI/PWgui ; then \
	    cd GUI/PWgui ; \
	    make clean cvsinit pwgui-source-notcl ; \
	    mv PWgui-*.tgz ../.. ; \
	else \
	    echo ; \
	    echo "  Sorry, tar-gui works only for CVS-sources !!!" ; \
	    echo ; \
	fi
