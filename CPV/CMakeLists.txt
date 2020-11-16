set(sources
    src/berryion.f90
    src/bforceion.f90
    src/cell_nose.f90
    src/cg.f90
    src/cg_sub.f90
    src/cglib.f90
    src/chargedensity.f90
    src/cp_autopilot.f90
    src/cp_emass.f90
    src/cp_interfaces.f90
    src/cp_restart_new.f90
    src/cplib.f90
    src/cpr_mod.f90
    src/cpr.f90
    src/dealloc.f90
    src/dforceb.f90
    src/efermi.f90
    src/efield.f90
    src/eigs0.f90
    src/electrons.f90
    src/electrons_nose.f90
    src/energies.f90
    src/ensemble_dft.f90
    src/exch_corr.f90
    src/exx_cg.f90
    src/exx_es.f90
    src/exx_gs.f90
    src/exx_module.f90
    src/exx_pair.f90
    src/exx_psi.f90
    src/exx_vofr.f90
    src/forces.f90
    src/fromscra.f90
    src/gram.f90
    src/gtable.f90
    src/init.f90
    src/init_run.f90
    src/inner_loop_cold.f90
    src/input.f90
    src/ions_nose.f90
    src/ions_positions.f90
    src/ksstates.f90
    src/ldaU_cp.f90
    src/ldaU.f90
    src/ldaUpen.f90
    src/mainvar.f90
    src/cpr_loop.f90
    src/cplib_meta.f90
    src/metaxc.f90
    src/modules.f90
    src/move_electrons.f90
    src/newd.f90
    src/nl_base.f90
    src/nlcc.f90
    src/ortho_base.f90
    src/ortho.f90
    src/phasefactor.f90
    src/plugin_add_potential.f90
    src/plugin_clean.f90
    src/plugin_clock.f90
    src/plugin_energy.f90
    src/plugin_ext_forces.f90
    src/plugin_get_potential.f90
    src/plugin_init_base.f90
    src/plugin_init_cell.f90
    src/plugin_initialization.f90
    src/plugin_init_ions.f90
    src/plugin_int_forces.f90
    src/plugin_print_energies.f90
    src/plugin_print_info.f90
    src/plugin_read_input.f90
    src/plugin_utilities.f90
    src/potentials.f90
    src/pres_ai_mod.f90
    src/print_out.f90
    src/printout_base.f90
    src/problem_size.f90
    src/pseudo_base.f90
    src/pseudopot.f90
    src/pseudopot_sub.f90
    src/qmatrixd.f90
    src/qqberry.f90
    src/restart.f90
    src/restart_sub.f90
    src/runcp.f90
    src/sic.f90
    src/smallbox.f90
    src/smallbox_grid.f90
    src/smallbox_gvec.f90
    src/smallbox_subs.f90
    src/smallbox_lib.f90
    src/spline.f90
    src/spinsq.f90
    src/stop_run.f90
    src/stress.f90
    src/vol_clu.f90
    src/vofrho.f90
    src/wannier_base.f90
    src/wannier.f90
    src/wave.f90
    src/wave_types.f90
    src/wf.f90
    src/makov_payne.f90
    src/entropy.f90)

qe_add_library(qe_cpv ${sources})
add_library(QE::CPV ALIAS qe_cpv)
target_link_libraries(qe_cpv
    PRIVATE
        QE::LAX
        QE::UTILX
        QE::Modules
        QE::FOX
        QE::FFTX
        QE::UPF
        QE::OpenMP_Fortran
        QE::MPI_Fortran)

###########################################################
# cp.x
###########################################################
set(sources src/cprstart.f90)
qe_add_executable(qe_cpv_exe ${sources})
set_target_properties(qe_cpv_exe PROPERTIES OUTPUT_NAME cp.x)
target_link_libraries(qe_cpv_exe
    PRIVATE
        QE::Modules
        QE::LAX
        QE::CPV)

###########################################################
# manycp.x
###########################################################
set(sources src/manycp.f90)
qe_add_executable(qe_cpv_manycp_exe ${sources})
set_target_properties(qe_cpv_manycp_exe PROPERTIES OUTPUT_NAME manycp.x)
target_link_libraries(qe_cpv_manycp_exe
    PRIVATE
        QE::LAX
        QE::Modules
        QE::UTILX
        QE::CPV)

###########################################################
# cppp.x
###########################################################
set(sources src/cppp.f90)
qe_add_executable(qe_cpv_cppp_exe ${sources})
set_target_properties(qe_cpv_cppp_exe PROPERTIES OUTPUT_NAME cppp.x)
target_link_libraries(qe_cpv_cppp_exe
    PRIVATE
        QE::Modules
        QE::CPV)

###########################################################
# wfdd.x
###########################################################
set(sources src/wfdd.f90)
qe_add_executable(qe_cpv_wfdd_exe ${sources})
set_target_properties(qe_cpv_wfdd_exe PROPERTIES OUTPUT_NAME wfdd.x)
target_link_libraries(qe_cpv_wfdd_exe
    PRIVATE
        QE::LAPACK
        QE::Modules
        QE::CPV)

###########################################################

qe_install_targets(
    # Libraries
    qe_cpv
    # Exacutables
    qe_cpv_exe
    qe_cpv_manycp_exe
    qe_cpv_cppp_exe
    qe_cpv_wfdd_exe
)
