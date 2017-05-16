#!/bin/bash                                                                                                           #
# Copyright (C) 2001 Quantum ESPRESSO                                                                                 #
# This program is free software; you can redistribute it and/or                                                       # modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2                                                      # of the License. See the file `License' in the root directory
# of the present distribution.

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
	@echo ' run-custom-test-parallel testdir=DIR : run test in DIR only (parallel 4 MPI)'
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

run-tests-ph-serial : prolog pseudo
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=ph_base

run-tests-ph-parallel : prolog pseudo
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=ph_base

run-tests-epw-serial : prolog pseudo
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=epw_all

run-tests-epw-parallel : prolog pseudo
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=epw_all

run-travis : prolog pseudo
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=epw_base

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
	@for x in `find ph_* -name "test*"`; do rm -rf $$x; done
	@for x in `find cp_* -name "test*"`; do rm -rf $$x; done
	@for x in `find epw_* -name "test*"`; do rm -rf $$x; done
	@for x in `find pw_* -name "CRASH"`; do rm -rf $$x; done
	@for x in `find ph_* -name "CRASH*"`; do rm -rf $$x; done
	@for x in `find cp_* -name "CRASH"`; do rm -rf $$x; done
	@for x in `find epw_* -name "CRASH"`; do rm -rf $$x; done
# Special case for 'pw_vc-relax' test-case
	@rm -f pw_vc-relax/ave pw_vc-relax/e pw_vc-relax/tv pw_vc-relax/p pw_vc-relax/eal pw_vc-relax/avec
# Special case for 'noexe_langevin' test-case
	@rm -f noexe_langevin/trajectory-pwscf.xyz
# Special case for PH
	@for x in `find ph_* -name "*dyn*"`; do rm -rf $$x; done
# Cleaning EPW final files
	@rm -rf epw_*/_ph0 epw_*/*.save epw_*/*.wfc* epw_*/*.dyn* \
		epw_*/*.qdos_25.00 epw_*/*.P epw_*/*.freq epw_*/*.wout \
		epw_*/*.fmt epw_*/*.epmat_wanep epw_*/*.kgmap epw_*/*.egnv epw_*/*.ikmap \
		epw_*/*.nnkp epw_*/*.imag_iso_25.00 epw_*/*.pade_iso_25.00 epw_*/*.epmatwp* \
		epw_*/*.acon_iso_25.00 epw_*/*.crystal.fmt epw_*/*.a2f_iso epw_*/*.ephmat* \
		epw_*/*.kmap epw_*/*.lambda_k_pairs epw_*/*.ukk epw_*/*.imag_aniso_25.00 \
		epw_*/*.pade_aniso_gap0_25.00 epw_*/*.win epw_*/*.imag_aniso_gap0_25.00 \
		epw_*/*.chk epw_*/*.phdos_proj epw_*/*..a2f epw_*/*.imag_aniso_gap_FS_25.00 \
		epw_*/*.epmatwp epw_*/*.phdos epw_*/*.epwane epw_*/*.ephmat \
		epw_*/*.fmt epw_*/*.H epw_*/*.epb* epw_*/*self epw_*/*.a2f epw_*/*.01 \
		epw_*/*.kpt epw_*/*.gnu epw_*/*.eig epw_*/*.lambda_FS
# Special cases for EPW
	@rm -rf epw_base/save epw_super/save

# Special case EPW
	@rm -rf epw_trev/save	
	@rm -f userconfig

# Calling 'purge' means remove all reference outputs... be careful!
purge: clean
	@for x in `find pw_* -name "benchmark*"`; do rm -rf $$x; done
	@for x in `find cp_* -name "benchmark*"`; do rm -rf $$x; done
