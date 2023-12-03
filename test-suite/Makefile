
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
	@echo ' run-tests-XX [NPROCS=N]        : run tests [optionally, on N processors]'
	@echo ' for XX=pw, cp, ph, epw, hp, tddfpt, kcw, all_currents, pp, zg, xsd-pw'
	@echo ' run-tests [NPROCS=N]           : run all tests above (except the last four)'
	@echo ' run-custom-test testdir=DIR [NPROCS=N]: run test in DIR only'
	@echo ' compare                     : compare last output with reference'
	@echo ' pseudo                      : download needed PPs into ESPRESSO_PSEUDO'
	@echo ' clean                       : clean stdout/sderr of all tests'
	@echo ' '
	@echo 'For additional advanced commands and settings please manually inspect'
	@echo 'ENVIRONMENT and Makefile files'

prolog :
	@sed "s|XXXXXX|$(ESPRESSO_ROOT)|g" < userconfig.tmp > userconfig
	@if ! test -d $(ESPRESSO_TMPDIR); then mkdir $(ESPRESSO_TMPDIR); fi

pseudo :
	@./check_pseudo.sh pw_
	@./check_pseudo.sh cp_
	@./check_pseudo.sh epw_
	@./check_pseudo.sh tddfpt_
	@./check_pseudo.sh hp_
	@./check_pseudo.sh ph_
	@./check_pseudo.sh QEHeat_
	@./check_pseudo.sh kcw_
	@./check_pseudo.sh oscdft_

run-tests : run-tests-pw run-tests-cp run-tests-ph run-tests-epw run-tests-hp run-tests-tddfpt run-tests-kcw

run-tests-cp : prolog pseudo
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --verbose --category=cp_all

run-tests-pw : prolog pseudo
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --verbose --category=pw_all

run-tests-ph : prolog pseudo
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --verbose --category=ph_all

run-tests-epw : prolog pseudo
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --verbose --category=epw_all

run-tests-zg : prolog
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --verbose --category=zg_all

run-tests-hp : prolog pseudo
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --verbose --category=hp_all

run-tests-pp : prolog pseudo
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --verbose --category=pp_all

run-tests-tddfpt : prolog pseudo
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --verbose --category=tddfpt_all

run-tests-image : prolog pseudo
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --verbose --category=image_all

run-tests-kcw : prolog pseudo
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --verbose  --category=kcw_all

run-tests-all_currents : prolog pseudo
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --verbose --category=all_currents_all

run-tests-oscdft : prolog pseudo
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --verbose --category=oscdft_all

run-travis : clean prolog pseudo
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --verbose --category=epw_base

run-custom-test : prolog pseudo
	@if test -d $(testdir); then \
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --verbose --category=$(testdir) ; fi

run-tests-xsd-pw : prolog pseudo
	cd xsd_pw ; for file in ./*.in ; do \
	echo "Running using $$file" ; \
		env QE_USE_MPI=$(NPROCS) ../run-pw.sh -in $$file &> $$file.out ; \
		python ../validate_xsd_pw.py $$file; \
    done ; cd ..

create-reference-cp : pseudo prolog
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --category=cp_all make-benchmarks

create-reference-pw : pseudo prolog
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --category=pw_all make-benchmarks

create-reference-epw : pseudo prolog
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --category=epw_all make-benchmarks

create-reference-all_currents : pseudo prolog
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --category=all_currents_all make-benchmarks

create-reference-oscdft : pseudo prolog
	env QE_USE_MPI=$(NPROCS) ${TESTCODE_DIR}/bin/testcode.py --category=oscdft_all make-benchmarks

compare : compare-cp compare-pw compare-ph compare-epw compare-tddfpt compare-hp compare-all_currents

compare-cp : prolog
	${TESTCODE_DIR}/bin/testcode.py --category=cp_all  --verbose compare 2>&1 | tee out.CP.`date +%Y%m%d_%H%M%S`

compare-pw : prolog
	${TESTCODE_DIR}/bin/testcode.py --category=pw_all  --verbose compare 2>&1 | tee out.PW.`date +%Y%m%d_%H%M%S`

compare-ph : prolog
	${TESTCODE_DIR}/bin/testcode.py --category=ph_all  --verbose compare 2>&1 | tee out.PH.`date +%Y%m%d_%H%M%S`

compare-epw : prolog
	${TESTCODE_DIR}/bin/testcode.py --category=epw_all  --verbose compare 2>&1 | tee out.EPW.`date +%Y%m%d_%H%M%S`

compare-tddfpt : prolog
	${TESTCODE_DIR}/bin/testcode.py --category=tddfpt_all  --verbose compare 2>&1 | tee out.TDDFPT.`date +%Y%m%d_%H%M%S`

compare-hp : prolog
	${TESTCODE_DIR}/bin/testcode.py --category=hp_all  --verbose compare 2>&1 | tee out.HP.`date +%Y%m%d_%H%M%S`

compare-all_currents : prolog pseudo
	${TESTCODE_DIR}/bin/testcode.py --category=all_currents_all  --verbose compare 2>&1 | tee out.ALL_CURRENTS.`date +%Y%m%d_%H%M%S`


clean:
# output files
	@rm -rf pw_*/test* ph_*/test* cp_*/test* hp_*/test* pp_*/test* epw_*/test* zg_*/test* tddfpt_*/test* kcw_*/test*
# error files
	@rm -rf pw_*/CRASH ph_*/CRASH* cp_*/CRASH hp_*/CRASH pp_*/CRASH epw_*/CRASH zg_*/CRASH tddfpt_*/CRASH kcw_*/CRASH
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
	@rm -rf ph_*/*.wfc* ph_*/*.hub*
# Special case for PH dvscf_q2r
	@rm -rf ph_interpol*/w_pot/ ph_interpol*/*.dyn*.elph.* ph_interpol*/*.dyn.freq
# Special case for PH ahc
	@rm -rf ph_ahc*/wpot*/ ph_ahc*/selfen_*.dat ph_ahc*/diam.modes? ph_ahc*/*.ifc ph_ahc*/ahc_dir*
# Cleaning HP final files
	@rm -rf hp_*/HP hp_*/*hub* hp_*/*.xml hp_*/*.wfc* hp_*/*.save hp_*/*Hubbard_parameters* \
		hp_*/*HubbardV* hp_*/*mix* hp_*/*parameters.out
# Cleaning KCW final files
	@rm -rf	kcw_*/out kcw_*/*.nnkp kcw_*/*.wout kcw_*/*.amn kcw_*/*.mmn kcw_*/*.eig kcw_*/*band* kcw_*/*chk \
                kcw_*/*.dat kcw_*/*.xyz kcw_*/*.mat 
# Cleaning EPW final files
	@rm -rf epw_*/_ph0 epw_*/*.save epw_*/*.wfc* epw_*/*.dyn* \
		epw_*/*.qdos_* epw_*/*.P epw_*/*.wout epw_*/crystal.fmt epw_*/*.pade* \
		epw_*/*.epmat_wanep epw_*/*.kgmap epw_*/*.nnkp epw_*/*.epmatwp* \
		epw_*/*.acon_iso_* epw_*/*.a2f* epw_*/*.phdos* epw_*/*.lambda \
		epw_*/*.kmap epw_*/*.lambda_k_pairs epw_*/*.ukk epw_*/*.imag_aniso* \
		epw_*/*.win epw_*/*.imag_iso* epw_*/*.chk epw_*/*.freq  \
		epw_*/*.epmatwp epw_*/*.epwane epw_*/*.ephmat epw_*/epsilon2* \
		epw_*/*.H epw_*/*.epb* epw_*/*self epw_*/*.res* epw_*/*.01 \
		epw_*/*.kpt epw_*/*.gnu epw_*/*.lambda_FS epw_*/*.sigma_restart* \
		epw_*/specfun* epw_*/linewidth.*  epw_*/*.tau_restart* epw_*/lambda.phself* \
		epw_*/*_elcond_* epw_*/fort.* epw_*/decay.* epw_*/scattering_rate* \
		epw_*/*.fc.* epw_*/*_band.dat epw_*/*F_restart* epw_*/*F_restart_CB* \
		epw_*/*.mmn epw_*/*.epmatkq* epw_*/*.epmatkqcb* epw_*/sparse* \
		epw_*/*.Fin_restart* epw_*/*.Fin_restartcb* epw_*/*.acon_iso_* \
		epw_*/*.bvec epw_*/*.fc epw_*/*.xml epw_*/*.labelinfo.dat \
		epw_*/Fepmatkq1 epw_*/Fepmatkqcb1 epw_*/Fsparse epw_*/Fsparsecb \
		epw_*/inv_taucb_tmp epw_*/inv_tau_tmp epw_*/*.cube epw_*/EPW.bib \
		epw_*/*_hr.dat epw_*/*_wsvec.dat epw_*/epwdata.fmt epw_*/IBTEvel_sup.fmt \
		epw_*/inv_taucb.fmt epw_*/inv_tau.fmt epw_*/restart.fmt epw_*/vmedata.fmt \
		epw_*/selecq.fmt epw_*/dmedata.fmt epw_*/band.eig epw_*/ksdata.fmt \
        epw_*/inv_tau_freq.fmt epw_*/inv_taucb_freq.fmt epw_*/inv_taucb_mode.fmt \
        epw_*/mobility_nk.fmt epw_*/mobility_nuq.fmt epw_*/inv_tau_mode.fmt \
        epw_*/*.dos epw_*/*.dmat* epw_*/*.qmap epw_*/*.dwmatwe* epw_*/*.symk \
        epw_*/*.sthmatwe* epw_*/*.xqc* epw_*/*.dgmatwe* epw_*/*.cpmew*

# Special cases for EPW
	@rm -rf epw_base/save epw_super/save epw_pl/save/ epw_polar/save/ \
		epw_trev/save epw_trev_uspp/save epw_trev_paw/save \
		epw_mob/save epw_mob_ibte/save epw_mob_ibte_sym/save epw_mob_polar/save \
		epw_qdpt/save epw_2D/save epw_hall/save epw_qdpt/*.dat epw_wfpt/save \
		epw_wfpt/matdyn.modes epw_plrn/*.plrn epw_plrn/plrn_tmp epw_plrn/save \
		epw_plrn/dtau.plrn.xsf

# Cleaning IMAGE final files
	@rm -rf image_*/_ph0 image_*/*.save image_*/*.xml image_*/*.dyn? image_*/*.dyn??  \
                image_*/*.wfc*

# Cleaning ZG final files
	@rm -f zg_conf/equil_pos.txt zg_conf/ZG-*_0.020.dat

# Cleaning TDDFPT final files
	@rm -rf tddfpt_*/*.plot_chi.dat tddfpt_*/*.plot_eps.dat tddfpt_*/out/ tddfpt_*/turbo_spectrum.out
# Cleaning OS-CDFT final files
	@rm -rf oscdft_*/*.out oscdft_*/*.save oscdft_*/*.xml oscdft_*/test* oscdft_*/*.md oscdft_*/*.bfgs oscdft_*/*.mix* oscdft_*/*.wfc* oscdft_*/oscdft.in oscdft_*/*.tmp
# all_currents files
	@for x in `find QEHeat_* -name "test*"`; do rm -rf $$x; done
	@for x in `find QEHeat_* -name "CRASH"`; do rm -rf $$x; done
	@rm -rf QEHeat_*/save \
                QEHeat_*/current_hz*
	@rm -f userconfig


# Calling 'purge' means remove all reference outputs... be careful!
purge: clean
	@for x in `find pw_* -name "benchmark*"`; do rm -rf $$x; done
	@for x in `find cp_* -name "benchmark*"`; do rm -rf $$x; done
	@for x in `find QEHeat_* -name "benchmark*"`; do rm -rf $$x; done


run-tests-serial : run-tests-pw-serial run-tests-cp-serial run-tests-ph-serial run-tests-epw-serial run-tests-hp-serial run-tests-tddfpt-serial run-tests-kcw-serial
	@echo '$@ obsolete, use "make run-tests" instead'
run-tests-pw-serial : run-tests-pw
	@echo '$@ obsolete, use "make run-tests-pw" instead'
run-tests-cp-serial : run-tests-cp
	@echo '$@ obsolete, use "make run-tests-cp" instead'
run-tests-ph-serial : run-tests-ph
	@echo '$@ obsolete, use "make run-tests-ph" instead'
run-tests-epw-serial : run-tests-epw
	@echo '$@ obsolete, use "make run-tests-epw" instead'
run-tests-hp-serial : run-tests-hp
	@echo '$@ obsolete, use "make run-tests-hp" instead'
run-tests-tddfpt-serial : run-tests-tddfpt
	@echo '$@ obsolete, use "make run-tests-tddfpt" instead'
run-tests-kcw-serial : run-tests-kcw
	@echo '$@ obsolete, use "make run-tests-kcw" instead'
run-tests-pp-serial : run-tests-pp
	@echo '$@ obsolete, use "make run-tests-pp" instead'
run-tests-zg-serial : run-tests-zg
	@echo '$@ obsolete, use "make run-tests-zg" instead'

run-tests-parallel : run-tests-pw-parallel run-tests-cp-parallel run-tests-ph-parallel run-tests-epw-parallel run-tests-hp-parallel run-tests-tddfpt-parallel run-tests-kcw-parallel
	@echo '$@ obsolete, use "make run-tests NPROCS=4" instead'
run-tests-pw-parallel : prolog pseudo
	@echo '$@ obsolete, use "make run-tests-pw NPROCS=4" instead'
	env QE_USE_MPI=4 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=pw_all
run-tests-cp-parallel : prolog pseudo
	@echo '$@ obsolete, use "make run-tests-cp NPROCS=4" instead'
	env QE_USE_MPI=4 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=cp_all
run-tests-ph-parallel : prolog pseudo
	@echo '$@ obsolete, use "make run-tests-ph NPROCS=4" instead'
	env QE_USE_MPI=4 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=ph_all
run-tests-epw-parallel : prolog pseudo
	@echo '$@ obsolete, use "make run-tests-epw NPROCS=4" instead'
	env QE_USE_MPI=4 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=epw_all
run-tests-hp-parallel : prolog pseudo
	@echo '$@ obsolete, use "make run-tests-hp NPROCS=4" instead'
	env QE_USE_MPI=4 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=hp_all
run-tests-tddfpt-parallel : prolog pseudo
	@echo '$@ obsolete, use "make run-tests-tddfpt NPROCS=4" instead'
	env QE_USE_MPI=4 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=tddfpt_all
run-tests-kcw-parallel : prolog pseudo
	@echo '$@ obsolete, use "make run-tests-kcw NPROCS=4" instead'
	env QE_USE_MPI=4 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=kcw_all
run-tests-pp-parallel : prolog pseudo
	@echo '$@ obsolete, use "make run-tests-pp NPROCS=4" instead'
	env QE_USE_MPI=4 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=pp_all
run-tests-zg-parallel : prolog pseudo
	@echo '$@ obsolete, use "make run-tests-zg NPROCS=4" instead'
	env QE_USE_MPI=4 ${TESTCODE_DIR}/bin/testcode.py --verbose --category=zg_all
