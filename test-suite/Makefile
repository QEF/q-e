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
	@echo ' run-tests-XX-serial         : run tests for XX={pw,cp,ph,epw,hp} only, as above' 
	@echo ' run-tests-parallel          : run all tests in PARALLEL (4 MPI)'
	@echo ' run-tests-XX-parallel       : run tests for XX={pw,cp,ph,epw,hp} only, as above' 
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

pseudo :
	@./check_pseudo.sh pw_
	@./check_pseudo.sh cp_
	@./check_pseudo.sh epw_
	@./check_pseudo.sh tddfpt_
	@./check_pseudo.sh hp_
	@./check_pseudo.sh ph_

run-tests : run-tests-serial

run-tests-serial : run-tests-pw-serial run-tests-cp-serial run-tests-ph-serial run-tests-epw-serial run-tests-hp-serial  

run-tests-parallel : run-tests-pw-parallel run-tests-cp-parallel run-tests-ph-parallel run-tests-epw-parallel run-tests-hp-parallel


run-tests-cp-serial : prolog pseudo
	@if ! test -d $(ESPRESSO_TMPDIR); then mkdir $(ESPRESSO_TMPDIR); fi 
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=cp_all

run-tests-cp-parallel : prolog pseudo
	@if ! test -d $(ESPRESSO_TMPDIR); then mkdir $(ESPRESSO_TMPDIR); fi 
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=cp_all

run-tests-pw-serial : prolog pseudo
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=pw_all

run-tests-pw-parallel : prolog pseudo
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=pw_all

run-tests-ph-serial : prolog pseudo
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=ph_all

run-tests-ph-parallel : prolog pseudo
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=ph_all

run-tests-epw-serial : prolog pseudo
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=epw_partial

run-tests-epw-parallel : prolog pseudo
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=epw_all

run-tests-zg-serial : prolog
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=zg_all

run-tests-zg-parallel : prolog
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=zg_all

run-tests-hp-serial : prolog pseudo
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=hp_all

run-tests-hp-parallel : prolog pseudo
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=hp_all

run-tests-pp-serial : prolog pseudo
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=pp_all

run-tests-pp-parallel : prolog pseudo
	env QE_USE_MPI=1 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=pp_all

run-travis : clean prolog pseudo
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
	@if ! test -d $(ESPRESSO_TMPDIR); then mkdir $(ESPRESSO_TMPDIR); fi 
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --category=cp_all make-benchmarks

create-reference-pw : pseudo prolog
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --category=pw_all make-benchmarks

create-reference-epw : pseudo prolog
	env QE_USE_MPI=0 ${TESTCODE_DIR}/bin/testcode.py --category=epw_all make-benchmarks


compare : compare-cp compare-pw compare-ph compare-epw

compare-cp : prolog
	${TESTCODE_DIR}/bin/testcode.py --category=cp_all  --verbose compare 2>&1 | tee out.CP.`date +%Y%m%d_%H%M%S`

compare-pw : prolog
	${TESTCODE_DIR}/bin/testcode.py --category=pw_all  --verbose compare 2>&1 | tee out.PW.`date +%Y%m%d_%H%M%S`

compare-ph : prolog
	${TESTCODE_DIR}/bin/testcode.py --category=ph_all  --verbose compare 2>&1 | tee out.PH.`date +%Y%m%d_%H%M%S`

compare-epw : prolog
	${TESTCODE_DIR}/bin/testcode.py --category=epw_all  --verbose compare 2>&1 | tee out.EPW.`date +%Y%m%d_%H%M%S`

clean:
# output files
	@rm -rf pw_*/test* ph_*/test* cp_*/test* hp_*/test* pp_*/test* epw_*/test* zg_*/test*
# error files
	@rm -rf pw_*/CRASH ph_*/CRASH* cp_*/CRASH hp_*/CRASH pp_*/CRASH epw_*/CRASH zg_*/CRASH
# pwscf data files
	@rm -rf pw_*/pwscf.* ph_*/pwscf.* pp_*/pwscf.* epw_*/pwscf.*
# cp data files
	@rm -rf cp_si/si_*.* cp_o2/o2_*.* cp_h2o*/h2o_*.* cp_sio2/sio2_*.*
	@rm -rf cp_si/si.??? cp_o2/o2.??? cp_h2o*/h2o.??? cp_sio2/sio2.???
# Special case for 'pw_vc-relax' test-case
	@rm -f pw_vc-relax/ave pw_vc-relax/e pw_vc-relax/tv pw_vc-relax/p pw_vc-relax/eal pw_vc-relax/avec
# Special case for 'noexe_langevin' test-case
	@rm -f noexe_langevin/trajectory-pwscf.xyz
# Special case for PH
	@for x in `find ph_* -name "*dyn"`; do rm -rf $$x; done
# Cleaning PH final files
	@rm -rf ph_*/_ph0 ph_*/*.save ph_*/*.xml ph_*/a2F.* ph_*/*.dyn? ph_*/*.dyn?? ph_*/alpha2F.dat \
		ph_*/dyna2F ph_*/elph.gamma* ph_*/gam.lines ph_*/lambda ph_*/lambda.dat \
		ph_*/elph_dir ph_*/Al444.* ph_*/matdyn.modes ph_*/phonon.dos
	@rm -rf ph_metal/aluminum.* ph_base/carbon.* ph_base/nickel.* ph_base/nickelpaw.* ph_base/silicon.* \
	        ph_U_metal_us/iron.* ph_U_insulator_us/bn.* ph_U_metal_paw/nickel.* ph_U_insulator_paw/bn.*
# Special case for PH dvscf_q2r
	@rm -rf ph_interpol*/w_pot/ ph_interpol*/*.dyn*.elph.* ph_interpol*/*.wfc* ph_interpol*/*.dyn.freq
# Special case for PH ahc
	@rm -rf ph_ahc*/wpot*/ ph_ahc*/selfen_*.dat ph_ahc*/diam.modes? ph_ahc*/*.ifc ph_ahc*/ahc_dir* \
	  ph_ahc*/*.wfc*
# Cleaning HP final files
	@rm -rf hp_*/HP hp_*/*hub* hp_*/*.xml hp_*/*.wfc* hp_*/*.save hp_*/*Hubbard_parameters* \
		hp_*/*HubbardV* hp_*/*mix* hp_*/*parameters.out
# Cleaning EPW final files
	@rm -rf epw_*/_ph0 epw_*/*.save epw_*/*.wfc* epw_*/*.dyn* \
		epw_*/*.qdos_* epw_*/*.P epw_*/*.wout epw_*/*.fmt epw_*/*.pade* \
		epw_*/*.epmat_wanep epw_*/*.kgmap epw_*/*.nnkp epw_*/*.epmatwp* \
		epw_*/*.acon_iso_* epw_*/*.a2f* epw_*/*.phdos* epw_*/*.lambda \
		epw_*/*.kmap epw_*/*.lambda_k_pairs epw_*/*.ukk epw_*/*.imag_aniso* \
		epw_*/*.win epw_*/*.imag_iso* epw_*/*.chk epw_*/*.freq  \
		epw_*/*.epmatwp epw_*/*.epwane epw_*/*.ephmat epw_*/epsilon2* \
		epw_*/*.H epw_*/*.epb* epw_*/*self epw_*/*.res* epw_*/*.01 \
		epw_*/*.kpt epw_*/*.gnu epw_*/*.eig epw_*/*.lambda_FS epw_*/*.sigma_restart* \
		epw_*/specfun* epw_*/linewidth.*  epw_*/*.tau_restart* epw_*/lambda.phself* \
		epw_*/*_elcond_* epw_*/fort.* epw_*/decay.* epw_*/scattering_rate* \
		epw_*/*.fc.* epw_*/*_band.dat epw_*/*F_restart* epw_*/*F_restart_CB* \
		epw_*/*.mmn epw_*/*.epmatkq* epw_*/*.epmatkqcb* epw_*/sparse* \
		epw_*/*.Fin_restart* epw_*/*.Fin_restartcb* epw_*/*.acon_iso_* \
		epw_*/*.bvec epw_*/*.fc epw_*/*.xml epw_*/*.labelinfo.dat \
		epw_*/Fepmatkq1 epw_*/Fepmatkqcb1 epw_*/Fsparse epw_*/Fsparsecb \
		epw_*/inv_taucb_tmp epw_*/inv_tau_tmp epw_*/*.cube epw_*/EPW.bib
 
# Special cases for EPW
	@rm -rf epw_base/save epw_super/save epw_pl/save/ epw_polar/save/ \
                epw_trev/save epw_trev_uspp/save epw_trev_paw/save \
		epw_mob/save epw_mob_ibte/save epw_mob_ibte_sym/save epw_mob_polar/save

# Cleaning ZG final files
	@rm -f zg_conf/equil_pos.txt zg_conf/ZG-*_0.020.dat

	@rm -f userconfig

# Calling 'purge' means remove all reference outputs... be careful!
purge: clean
	@for x in `find pw_* -name "benchmark*"`; do rm -rf $$x; done
	@for x in `find cp_* -name "benchmark*"`; do rm -rf $$x; done
