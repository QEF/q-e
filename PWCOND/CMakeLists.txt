###########################################################
# pwcond.x
###########################################################
set(sources
    src/allocate_cond.f90
    src/bessj.f90
    src/compbs.f90
    src/compbs_2.f90
    src/condcom.f90
    src/condmain.f90
    src/cond_out.f90
    src/cond_restart.f90
    src/do_cond.f90
    src/eigenchnl.f90
    src/form_zk.f90
    src/four.f90
    src/free_mem.f90
    src/gep_x.f90
    src/gramsh.f90
    src/hev_ab.f90
    src/init_cond.f90
    src/init_orbitals.f90
    src/init_gper.f90
    src/integrals.f90
    src/jbloch.f90
    src/kbloch.f90
    src/local.f90
    src/local_set.f90
    src/openfil_cond.f90
    src/plus_u_setup.f90
    src/poten.f90
    src/print_clock_pwcond.f90
    src/realus_scatt.f90
    src/rotproc.f90
    src/save_cond.f90
    src/scat_states_plot.f90
    src/scatter_forw.f90
    src/summary_band.f90
    src/summary_tran.f90
    src/sunitary.f90
    src/transmit.f90
)

qe_add_executable(qe_pwcond_exe ${sources})
set_target_properties(qe_pwcond_exe PROPERTIES OUTPUT_NAME pwcond.x)
target_link_libraries(qe_pwcond_exe
    PRIVATE
        QE::LAPACK
        QE::PW
        QE::Modules)

###########################################################

qe_install_targets(qe_pwcond_exe)
