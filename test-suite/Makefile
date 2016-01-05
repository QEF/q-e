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

# --test-id={v<release-version>, r<revision-number>}
create-reference : pseudo
	${TESTCODE_DIR}/bin/testcode.py --processors=0 --category=_default_  --test-id=$(REFERENCE_VERSION) make-benchmarks

compare :
	${TESTCODE_DIR}/bin/testcode.py --processors=0 --category=_default_  --verbose compare 2>&1 | tee out.`date +%Y%m%d_%H%M%S`

clean: prolog
	@for x in `find pw_* -name "test*"`; do rm -rf $$x; done
	@for x in `find cp_* -name "test*"`; do rm -rf $$x; done
# Special case for 'pw_vc-relax' test-case
	@rm -f rm pw_vc-relax/e pw_vc-relax/tv pw_vc-relax/p pw_vc-relax/eal pw_vc-relax/avec
	@rm userconfig

purge: clean
	@for x in `find pw_* -name "benchmark*"`; do rm -rf $$x; done
	@for x in `find cp_* -name "benchmark*"`; do rm -rf $$x; done
