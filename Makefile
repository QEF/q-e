what:
	@echo 'First step:'
	@echo './configure "your-system"'
	@echo './configure with no arguments gives a list of supported systems'
	@echo 'Then: Edit the file make.sys , if necessary'
	@echo 'Then: make "task", where task is one of the following:'
	@echo '   pw        (basic code for scf, struct. optimization, MD)'
	@echo '   pp        (postprocessing programs)'
	@echo '   ph        (phonon code)'
	@echo '   d3        (third-order derivatives)'
	@echo '   gamma     (Gamma-only version of pw and ph)'
	@echo '   upf       (utilities for pseudopotential conversion)'
	@echo '   tools     (misc tools for data analysis)'
	@echo '   tar       (create a tar file containing the distribution)'
	@echo '   clean     (remove executables and object)'
	@echo '   veryclean (revert distribution to the original status)'
	@echo '   fpmd      (FPMD code for Car-Parrinello MD)'
	@echo '   cpv       (CPV code: CP MD with ultrasoft pseudopotentials)'
	@echo '   links     (creates links to executables in bin/)'

tools:
	( cd pwtools ; make all )

upf:
	( cd upftools ; make all )

gamma: modules pw
	( cd Gamma; make all )

all: pp ph d3

d3: modules pw ph
	( cd D3; make all )

ph: modules pw
	( cd PH; make all )

pp: modules pw
	( cd PP; make all )

pw: modules
	( cd PW; make all )

modules:
	( cd Modules; make all )

libs: modules
	( cd clib; make all ); ( cd flib; make all );

fpmd: modules libs
	( cd FPMD; make all )

cpv: modules libs
	( cd CPV; make all )

links:
	( cd bin/ ; ln -fs ../PW/pw.x ../PH/ph.x ../D3/d3.x ../Gamma/pwg.x ../Gamma/phcg.x ../CPV/cp.x ../FPMD/par2.x ../PP/average.x ../PP/bands.x ../PP/chdens.x ../PP/dos.x ../PP/plotrho.x ../PP/pp.x ../PP/projwfc.x ../PP/voronoy.x ../pwtools/band_plot.x ../pwtools/dynmat.x ../pwtools/fqha.x ../pwtools/matdyn.x ../pwtools/q2r.x . )

clean :
	( cd PW ; make clean ) ; \
	( cd PH ; make clean ) ; \
	( cd PP ; make clean ) ; \
	( cd D3 ; make clean ) ; \
	( cd Gamma ; make clean ) ; \
	( cd pwtools ; make clean ) ; \
	( cd upftools ; make clean ) ; \
        ( cd Modules ; make clean ) ; \
        ( cd install ; make clean ) ; \
	( cd clib ; make clean ) ; \
	( cd flib ; make clean ) ; \
	( cd FPMD ; make clean ) ; \
	( cd CPV ; make clean ) ;

veryclean: clean
	/bin/rm make.rules make.sys */.dependencies */dum1 */dum2 bin/*

tar:
	tar -cf - License Makefile README include/machine.h* \
                  install/* \
                  PW/Makefile PW/*.f90 \
                  PH/Makefile PH/*.f90 \
                  PP/Makefile PP/*.f90 \
                  D3/Makefile D3/*.f90 \
                  Modules/Makefile Modules/*.f90 \
                  pwtools/Makefile pwtools/*.f90 \
                  upftools/Makefile upftools/*.f90 upftools/UPF \
	          pwdocs pw_examples/ pseudo/ | gzip > pw.tar.gz
