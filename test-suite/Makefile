# Copyright (C) 2001-2016 Quantum ESPRESSO group
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
#
# Maintainer: Filippo Spiga (filippo.spiga@quantum-espresso.org)
#
# Makefile to run regression numerical tests. Simple is good

include ./ENVIRONMENT

default : run-tests

prolog :
	@sed "s|XXXXXX|$(ESPRESSO_ROOT)|g" < userconfig.tmp > userconfig

pseudo :
	@./check_pseudo.sh pw_
	@./check_pseudo.sh cp_

run-tests : prolog pseudo
	${TESTCODE_DIR}/bin/testcode.py --processors=0 --category=_default_ 

run-tests-cp : prolog pseudo
	${TESTCODE_DIR}/bin/testcode.py --processors=0 --category=cp_all

run-tests-pw : prolog pseudo
	${TESTCODE_DIR}/bin/testcode.py --processors=0 --category=pw_all

run-custom-test : prolog pseudo
	@if test -d $(testdir); then \
	${TESTCODE_DIR}/bin/testcode.py --processors=0 --category=$(testdir) ; fi

create-reference : pseudo prolog
	${TESTCODE_DIR}/bin/testcode.py --processors=0 --category=_default_ make-benchmarks

compare :
	${TESTCODE_DIR}/bin/testcode.py --processors=0 --category=_default_  --verbose compare 2>&1 | tee out.`date +%Y%m%d_%H%M%S`

clean: prolog
	@for x in `find pw_* -name "test*"`; do rm -rf $$x; done
	@for x in `find cp_* -name "test*"`; do rm -rf $$x; done
# Special case for 'pw_vc-relax' test-case
	@rm -f pw_vc-relax/ave pw_vc-relax/e pw_vc-relax/tv pw_vc-relax/p pw_vc-relax/eal pw_vc-relax/avec
# Special case for 'pw_langevin' test-case
	@rm -f pw_langevin/trajectory-pwscf.xyz
# Soecial cases for CP	
	@rm -f cp_si/qe_input.xml cp_sio2/qe_input.xml cp_o2/qe_input.xml cp_h2o/qe_input.xml
	@rm userconfig

# Calling 'purge' means remove all reference outputs... be careful!
purge: clean
	@for x in `find pw_* -name "benchmark*"`; do rm -rf $$x; done
	@for x in `find cp_* -name "benchmark*"`; do rm -rf $$x; done
