default:
	@echo 'to install, type at the shell prompt:'
	@echo '  ./configure.new'
	@echo '  make target'
	@echo 'where target is one of the following:'
	@echo '  pw           basic code for scf, structure optimization, MD'
	@echo '  fpmd         FPMD code for Car-Parrinello MD'
	@echo '  cp           CP code: CP MD with ultrasoft pseudopotentials'
	@echo '  ph           phonon code'
	@echo '  pp           postprocessing programs'
	@echo '  gamma        Gamma-only version of pw and ph'
	@echo '  nc           non collinear magnetic version of pw code'
	@echo '  pwneb        version of pw with NEB'
	@echo '  pwcond       ballistic conductance'
	@echo '  d3           third-order derivatives'
	@echo '  tools        misc tools for data analysis'
	@echo '  upf          utilities for pseudopotential conversion'
	@echo '  all          same as "make pw ph pp gamma pwcond d3 tools"'
	@echo '  links        creates links to executables in bin/'
	@echo '  clean        remove executables and objects'
	@echo '  veryclean    revert distribution to the original status'
	@echo '  tar          create a tarball of the source tree'

pw: modules libs
	( cd PW; make all )
fpmd: modules libs
	( cd FPMD; make all )
cp: modules libs
	( cd CPV; make all )

ph: pw
	( cd PH; make all )
pp: pw
	( cd PP; make all )
gamma: pw
	( cd Gamma; make all )
nc: pw
	( cd PWNC;  make all )
pwneb: pw
	( cd NEB; make all )

pwcond: pw pp
	( cd PWCOND; make all )
d3: pw ph
	( cd D3; make all )

tools: libs
	( cd pwtools ; make all )
upf: libs
	( cd upftools ; make all )

all: pw ph pp gamma pwcond d3 tools
all_: all fpmd cp nc pwneb upf

modules:
	( cd Modules; make all )
libs: modules
	( cd clib; make all )
	( cd flib; make all )

# create link only if file exists
links:
	test -d bin || mkdir bin
	( cd bin/ ; \
	  for exe in ../PW/pw.x ../PW/memory.x ../NEB/pwneb.x ../PH/ph.x \
	   	     ../D3/d3.x ../Gamma/pwg.x ../Gamma/phcg.x ../CPV/cp.x \
		../FPMD/par2.x ../PP/average.x ../PP/bands.x ../PP/chdens.x \
		../PP/dos.x ../PP/plotrho.x ../PP/pp.x ../PP/projwfc.x \
		../PP/voronoy.x ../PP/plotband.x ../PWCOND/pwcond.x \
		../pwtools/band_plot.x ../pwtools/dynmat.x ../pwtools/fqha.x \
		../pwtools/matdyn.x ../pwtools/q2r.x ../pwtools/dist.x \
		../pwtools/ev.x ../pwtools/kpoints.x ; do \
	    if test -f $$exe ; then ln -fs $$exe . ; fi \
	  done \
	)

# remove object files and executables
clean:
	touch make.rules make.sys # make complains if they aren't there
				  # same with .dependencies below
	for dir in PW PWNC NEB PH PP D3 PWCOND Gamma pwtools upftools \
		   Modules install clib flib FPMD CPV ; do \
	  if test -d $$dir ; then \
	    ( cd $$dir ; touch .dependencies ; make clean_ ) \
	  fi \
	done

# remove configuration files too
veryclean: clean
	- /bin/rm -rf make.rules make.sys */.dependencies \
		      config.log config.status */dum1 */dum2 bin/*.x \
		      autom4te.cache pw.tar.gz FPMD/version.h

tar:
	tar cvf pw.tar \
	    License README */README README.cvs INSTALL Makefile */Makefile \
	    configure.new configure.ac config.guess config.sub install-sh \
	    makedeps.sh moduledep.sh make.rules.in make.sys.in configure \
	    */*.f90 */*.c */*.f clib/*.h include/*.h* upftools/UPF \
	    pwtools/*.awk pwtools/*.sh
	# archive a few entire directories, but without CVS subdirs
	find install *docs *_examples pseudo -type f | \
		grep -v /CVS/ | xargs tar rvf pw.tar
	gzip pw.tar

