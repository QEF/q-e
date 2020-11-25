set(sources
    src/hp_allocate_q.f90
    src/hp_bcast_input.f90
    src/hpcom.f90
    src/hp_clean_q.f90
    src/hp_close_q.f90
    src/hp_calc_chi.f90
    src/hp_check_pert.f90
    src/hp_check_type.f90
    src/hp_dnstot_sum_q.f90
    src/hp_dnsq.f90
    src/hp_dealloc_q.f90
    src/hp_dealloc_1.f90
    src/hp_dealloc_2.f90
    src/hp_dvpsi_pert.f90
    src/hp_efermi_shift.f90
    src/hp_find_inequiv_sites.f90
    src/hp_generate_grids.f90
    src/hp_init.f90
    src/hp_init_q.f90
    src/hp_load_q.f90
    src/hp_ns_trace.f90
    src/hp_openfil_q.f90
    src/hp_psymdvscf.f90
    src/hp_postproc.f90
    src/hp_prepare_q.f90
    src/hp_print_clock.f90
    src/hp_q_points.f90
    src/hp_R_points.f90
    src/hp_rotate_dnsq.f90
    src/hp_readin.f90
    src/hp_read_chi.f90
    src/hp_read_dnsq.f90
    src/hp_run_nscf.f90
    src/hp_solve_linear_system.f90
    src/hp_symdnsq.f90
    src/hp_symdvscf.f90
    src/hp_setup_q.f90
    src/hp_summary.f90
    src/hp_summary_q.f90
    src/hp_stop_smoothly.f90
    src/hp_write_chi.f90
    src/hp_write_chi_full.f90
    src/hp_write_dnsq.f90)

qe_add_library(qe_hp ${sources})
add_library(QE::HP ALIAS qe_hp)
target_link_libraries(qe_hp
    PRIVATE
        QE::PW
        QE::LR_Modules
        QE::Modules
        QE::FFTX
        QE::UTILX
        QE::UPF)

###########################################################
# hp.x
###########################################################
set(sources src/hp_main.f90)
qe_add_executable(qe_hp_exe ${sources})
set_target_properties(qe_hp_exe PROPERTIES OUTPUT_NAME hp.x)
target_link_libraries(qe_hp_exe
    PRIVATE
        QE::PW
        QE::LAX
        QE::Modules
        QE::HP
        QE::UTILX)

###########################################################

qe_install_targets(
    qe_hp
    qe_hp_exe)