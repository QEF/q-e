set(src_cpv
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
    src/wave_base.f90
    src/wave_types.f90
    src/wf.f90
    src/makov_payne.f90
    src/entropy.f90)
qe_enable_cuda_fortran("${src_cpv}")

qe_add_library(qe_cpv ${src_cpv})
target_link_libraries(qe_cpv
    PRIVATE
        qe_lax
        qe_modules
        qe_fftx
        qe_upflib
        qe_openmp_fortran
        qe_mpi_fortran
        qe_xclib)
if(QE_ENABLE_CUDA)
    target_link_libraries(qe_cpv
        PRIVATE
            qe_devxlib)
endif()

###########################################################
# cp.x
###########################################################
set(src_cp_x src/cprstart.f90)
qe_enable_cuda_fortran("${src_cp_x}")
qe_add_executable(qe_cpv_exe ${src_cp_x})
set_target_properties(qe_cpv_exe 
    PROPERTIES 
        OUTPUT_NAME cp.x)
target_link_libraries(qe_cpv_exe
    PRIVATE
        qe_modules
        qe_lax
        qe_cpv
        qe_lapack
        qe_xclib)

###########################################################
# manycp.x
###########################################################
set(src_manycp_x src/manycp.f90)
qe_enable_cuda_fortran("${src_manycp_x}")
qe_add_executable(qe_cpv_manycp_exe ${src_manycp_x})
set_target_properties(qe_cpv_manycp_exe 
    PROPERTIES 
        OUTPUT_NAME manycp.x)
target_link_libraries(qe_cpv_manycp_exe
    PRIVATE
        qe_lax
        qe_modules
        qe_cpv
        qe_lapack
        qe_xclib)

###########################################################
# cppp.x
###########################################################
set(src_cppp_x src/cppp.f90)
qe_add_executable(qe_cpv_cppp_exe ${src_cppp_x})
set_target_properties(qe_cpv_cppp_exe 
    PROPERTIES 
        OUTPUT_NAME cppp.x)
target_link_libraries(qe_cpv_cppp_exe
    PRIVATE
        qe_modules
        qe_cpv
        qe_clib
        qe_lapack
        qe_xclib)

###########################################################
# wfdd.x
###########################################################
set(src_wfdd_x src/wfdd.f90)
qe_add_executable(qe_cpv_wfdd_exe ${src_wfdd_x})
set_target_properties(qe_cpv_wfdd_exe 
    PROPERTIES 
        OUTPUT_NAME wfdd.x)
target_link_libraries(qe_cpv_wfdd_exe
    PRIVATE
        qe_lapack
        qe_modules
        qe_cpv
        qe_xclib)

###########################################################

qe_install_targets(
    # Libraries
    qe_cpv
    # Exacutables
    qe_cpv_exe
    qe_cpv_manycp_exe
    qe_cpv_cppp_exe
    qe_cpv_wfdd_exe)
