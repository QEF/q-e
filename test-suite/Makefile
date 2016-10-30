# Copyright (C) 2001-2016 Quantum ESPRESSO group
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
#
# Maintainer: Filippo Spiga (filippo.spiga@quantum-espresso.org)

include ./ENVIRONMENT

default :
	@echo 'To run Quantum ESPRESSO test-suite, type at the shell prompt:'
	@echo ' '
	@echo '  make target'
	@echo ' '
	@echo 'where <target> identifies an action'
	@echo ' run-tests                   : run-tests-serial (default)'
	@echo ' run-tests-serial            : run all tests in SERIAL' 
	@echo ' run-tests-XX-serial         : run tests for XX={pw,cp} only, as above' 
	@echo ' run-tests-parallel          : run all tests in PARALLEL (4 MPI)'
	@echo ' run-tests-XX-parallel       : run tests for XX={pw,cp} only, as above' 
	@echo ' run-custom-test testdir=DIR : run test in DIR only (serial)'
	@echo ' compare                     : compare last output with reference'
	@echo ' clean                       : clean stdout/sderr of all tests'
	@echo ' '
	@echo 'For additional advanced commands and settings please manually inspect'
	@echo 'ENVIRONMENT and Makefile files'

prolog :
	@sed "s|XXXXXX|$(ESPRESSO_ROOT)|g" < userconfig.tmp > userconfig

pseudo :
	@./check_pseudo.sh pw_
	@./check_pseudo.sh cp_


run-tests : run-tests-serial

run-tests-serial : run-tests-cp-serial run-tests-pw-serial

run-tests-parallel : run-tests-cp-parallel run-tests-pw-parallel


run-tests-cp-serial : prolog pseudo
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=cp_all

run-tests-cp-parallel : prolog pseudo
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=cp_all

run-tests-pw-serial : prolog pseudo
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=pw_all

run-tests-pw-parallel : prolog pseudo
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=pw_all

run-tests-epw-serial : prolog pseudo
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=epw_all

run-tests-epw-parallel : prolog pseudo
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=epw_all


run-custom-test : run-custom-test-serial

run-custom-test-serial : prolog pseudo
	@if test -d $(testdir); then \
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=$(testdir) ; fi

run-custom-test-parallel : prolog pseudo
	@if test -d $(testdir); then \
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=$(testdir) ; fi


run-tests-xsd-pw-serial : prolog pseudo
	cd xsd_pw ; for file in ./*.in ; do \
  	echo "Running using $$file" ; \
		env QE_USE_MPI=0 ../run-pw.sh -in $$file &> $$file.out ; \
		python ../validate_xsd_pw.py $$file; \
    done ; cd ..

create-reference-cp : pseudo prolog
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --category=cp_all make-benchmarks

create-reference-pw : pseudo prolog
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --category=pw_all make-benchmarks

create-reference-epw : pseudo prolog
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --category=epw_all make-benchmarks


compare : compare-cp compare-pw compare-epw

compare-cp : prolog
	${TESTCODE_DIR}/bin/testcode.py --category=cp_all  --verbose compare 2>&1 | tee out.CP.`date +%Y%m%d_%H%M%S`

compare-pw : prolog
	${TESTCODE_DIR}/bin/testcode.py --category=pw_all  --verbose compare 2>&1 | tee out.PW.`date +%Y%m%d_%H%M%S`

compare-epw : prolog
	${TESTCODE_DIR}/bin/testcode.py --category=epw_all  --verbose compare 2>&1 | tee out.EPW.`date +%Y%m%d_%H%M%S`

clean:
	@for x in `find pw_* -name "test*"`; do rm -rf $$x; done
	@for x in `find cp_* -name "test*"`; do rm -rf $$x; done
	@for x in `find epw_* -name "test*"`; do rm -rf $$x; done
# Special case for 'pw_vc-relax' test-case
	@rm -f pw_vc-relax/ave pw_vc-relax/e pw_vc-relax/tv pw_vc-relax/p pw_vc-relax/eal pw_vc-relax/avec
# Special case for 'pw_langevin' test-case
	@rm -f pw_langevin/trajectory-pwscf.xyz
# Special cases for CP	
	@rm -f cp_si/qe_input.xml cp_sio2/qe_input.xml cp_o2/qe_input.xml cp_h2o/qe_input.xml
	@rm -f userconfig

# Calling 'purge' means remove all reference outputs... be careful!
purge: clean
	@for x in `find pw_* -name "benchmark*"`; do rm -rf $$x; done
	@for x in `find cp_* -name "benchmark*"`; do rm -rf $$x; done
