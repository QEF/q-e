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
	@echo '  raman        raman tensor'
	@echo '  tools        misc tools for data analysis'
	@echo '  ld1          utilities for pseudopotential generation'
	@echo '  upf          utilities for pseudopotential conversion'
	@echo '  pwall        same as "make pw ph pp gamma nc pwcond d3 raman tools"'
	@echo '  all          same as "make pwall fpmd cp ld1 upf"'
	@echo '  clean        remove executables and objects'
	@echo '  veryclean    revert distribution to the original status'
	@echo '  tar          create a tarball of the source tree'
	@echo '  tar-gui      create a tarball of the GUI sources'

pw : mods libs
	if test -d PW   ; then ( cd PW   ; make all ) ; fi

fpmd : mods libs
	if test -d FPMD ; then ( cd FPMD ; make all ) ; fi

cp : mods libs
	if test -d CPV  ; then ( cd CPV  ; make all ) ; fi

ph : mods libs pw
	if test -d PH ; then ( cd PH ; make all ) ; fi

pp : mods libs pw
	if test -d PP ; then ( cd PP ; make all ) ; fi

gamma : mods libs pw
	if test -d Gamma  ; then ( cd Gamma  ; make all ) ; fi

nc : mods libs pw
	if test -d PWNC   ; then ( cd PWNC   ; make all ) ; fi

pwcond : mods libs pw pp
	if test -d PWCOND ; then ( cd PWCOND ; make all ) ; fi

d3 : mods libs pw ph
	if test -d D3 ; then ( cd D3 ; make all ) ; fi

raman : mods libs pw ph
	if test -d Raman ; then ( cd Raman ; make all ) ; fi

tools : mods libs pw
	if test -d pwtools  ; then ( cd pwtools  ; make all ) ; fi

ld1 : mods libs pw
	if test -d atomic ; then ( cd atomic ; make all ) ; fi

upf : mods libs
	if test -d upftools ; then ( cd upftools ; make all ) ; fi

pwall : pw ph pp gamma nc pwcond d3 raman tools
all   : pwall fpmd cp ld1 upf 

mods :
	( cd Modules; make all )
libs : mods
	( cd clib ; make all )
	( cd flib ; make all )

# remove object files and executables
clean :
	touch make.rules make.sys 
	# make complains if they aren't there; same with .dependencies below
	for dir in PW PWNC PH PP D3 PWCOND Gamma pwtools upftools atomic \
		   Modules install clib flib FPMD CPV Raman ; do \
	    if test -d $$dir ; then \
		( cd $$dir ; touch .dependencies ; make clean ) \
	    fi \
	done
	- /bin/rm -rf bin/*.x

# remove configuration files too
veryclean : clean
	- /bin/rm -rf make.rules make.sys */.dependencies \
		      config.log config.status autom4te.cache \
		      pw.tar.gz FPMD/version.h \
		      intel.pcl */intel.pcl
	- cd examples ; ./make_clean
	- if test -d GUI ; then ( cd GUI; make veryclean ) ; fi

tar :
	tar cvf pw.tar \
	    License README* */README* Makefile */Makefile \
	    configure configure.ac config.guess config.sub install-sh \
	    makedeps.sh moduledep.sh make.rules.in make.sys.in configure.old \
	    */*.f90 */*.c */*.f clib/*.h include/*.h* upftools/UPF \
	    pwtools/*.awk pwtools/*.sh bin
	# archive a few entire directories, but without CVS subdirs
	find install Doc *docs examples pseudo -type f \
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
