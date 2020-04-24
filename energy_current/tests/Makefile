#!/bin/bash
#
# Copyright (C) 2001 Quantum ESPRESSO
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
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
	@echo ' pseudo                      : download needed PPs into ESPRESSO_PSEUDO'
	@echo ' clean                       : clean stdout/sderr of all tests'
	@echo ' '
	@echo 'For additional advanced commands and settings please manually inspect'
	@echo 'ENVIRONMENT and Makefile files'

prolog :
	@sed "s|XXXXXX|$(ESPRESSO_ROOT)|g" < userconfig.tmp > userconfig


run-tests : run-tests-serial

run-tests-serial : run-tests-all_currents-serial

run-tests-parallel : run-tests-all_currents-parallel


run-tests-all_currents-serial : prolog
	@if ! test -d $(ESPRESSO_TMPDIR); then mkdir $(ESPRESSO_TMPDIR); fi 
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=all_currents_all

run-tests-all_currents-parallel : prolog
	@if ! test -d $(ESPRESSO_TMPDIR); then mkdir $(ESPRESSO_TMPDIR); fi 
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=all_currents_all


run-custom-test : run-custom-test-serial

run-custom-test-serial : prolog
	@if test -d $(testdir); then \
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=$(testdir) ; fi

run-custom-test-parallel : prolog
	@if test -d $(testdir); then \
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=$(testdir) ; fi



create-reference-all_currents : prolog
	@if ! test -d $(ESPRESSO_TMPDIR); then mkdir $(ESPRESSO_TMPDIR); fi 
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --category=all_currents_all make-benchmarks



compare : compare-all_currents

compare-all_currents : prolog
	${TESTCODE_DIR}/bin/testcode.py --category=all_currents_all  --verbose compare 2>&1 | tee out.ALL_CURRENTS.`date +%Y%m%d_%H%M%S`


clean:
# output files
	@for x in `find energy_current_* -name "test*"`; do rm -rf $$x; done
# error files
	@for x in `find energy_current_* -name "CRASH"`; do rm -rf $$x; done
	@rm -f userconfig

# Calling 'purge' means remove all reference outputs... be careful!
purge: clean
	@for x in `find energy_current_* -name "benchmark*"`; do rm -rf $$x; done
