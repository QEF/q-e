set(sources
    src/addusdens1d.f90
    src/add_shift_cc.f90
    src/add_shift_lc.f90
    src/add_shift_us.f90
    src/atomic_wfc_nc_proj.f90
    src/chdens_module.f90
    src/chdens_bspline.f90
    src/compute_ppsi.f90
    src/compute_sigma_avg.f90
    src/cube.f90
    src/dosg.f90
    src/do_initial_state.f90
    src/do_shift_ew.f90
    src/d_matrix_nc.f90
    src/d_matrix_so.f90
    src/elf.f90
    src/fermisurfer_common.f90
    src/ggen1d.f90
    src/local_dos.f90
    src/local_dos_mag.f90
    src/local_dos1d.f90
    src/openfil_pp.f90
    src/partialdos.f90
    src/paw_postproc.f90
    src/projections_mod.f90
    src/projwfc_box.f90
    src/punch_plot.f90
    src/read_proj.f90
    src/smallgk.f90
    src/stm.f90
    src/stop_pp.f90
    src/sum_band_kin.f90
    src/sym_band.f90
    src/vasp_read_chgcar_mod.f90
    src/vasp_xml_module.f90
    src/wannier_proj.f90
    src/wannier_enrg.f90
    src/work_function.f90
    src/write_p_avg.f90
    src/write_proj.f90
    src/write_io_header.f90
    src/write_hamiltonians.f90
    src/xc_vdW_scale_mod.f90
    src/projwfc.f90)

qe_add_library(qe_pp ${sources})
add_library(QE::PP ALIAS qe_pp)
target_compile_definitions(qe_pp PRIVATE "__WANLIB")
target_link_libraries(qe_pp
    PRIVATE
        QE::OpenMP_Fortran
        QE::FOX
        QE::PW
        QE::Modules)

###########################################################
# pp.x
###########################################################
set(sources src/postproc.f90)
qe_add_executable(qe_pp_exe ${sources})
set_target_properties(qe_pp_exe PROPERTIES OUTPUT_NAME pp.x)
target_link_libraries(qe_pp_exe
    PRIVATE
        QE::PW
        QE::Modules
        QE::PP)

###########################################################
# open_grid.x
###########################################################
set(sources src/open_grid.f90)
qe_add_executable(qe_pp_opengrid_exe ${sources})
set_target_properties(qe_pp_opengrid_exe PROPERTIES OUTPUT_NAME open_grid.x)
target_link_libraries(qe_pp_opengrid_exe
    PRIVATE
        QE::PW
        QE::UTILX
        QE::Modules
        QE::PP)

###########################################################
# average.x
###########################################################
set(sources src/average.f90)
qe_add_executable(qe_pp_average_exe ${sources})
set_target_properties(qe_pp_average_exe PROPERTIES OUTPUT_NAME average.x)
target_link_libraries(qe_pp_average_exe
    PRIVATE
        QE::PW
        QE::FFTX
        QE::Modules
        QE::PP)

###########################################################
# bands.x
###########################################################
set(sources src/bands.f90)
qe_add_executable(qe_pp_bands_exe ${sources})
set_target_properties(qe_pp_bands_exe PROPERTIES OUTPUT_NAME bands.x)
target_link_libraries(qe_pp_bands_exe
    PRIVATE
        QE::PW
        QE::Modules
        QE::PP)

###########################################################
# dos.x
###########################################################
set(sources src/dos.f90)
qe_add_executable(qe_pp_dos_exe ${sources})
set_target_properties(qe_pp_dos_exe PROPERTIES OUTPUT_NAME dos.x)
target_link_libraries(qe_pp_dos_exe
    PRIVATE
        QE::PW
        QE::PP)

###########################################################
# pawplot.x
###########################################################
set(sources src/pawplot.f90)
qe_add_executable(qe_pp_pawplot_exe ${sources})
set_target_properties(qe_pp_pawplot_exe PROPERTIES OUTPUT_NAME pawplot.x)
target_link_libraries(qe_pp_pawplot_exe
    PRIVATE
        QE::PW
        QE::Modules
        QE::PP)

###########################################################
# plan_avg.x
###########################################################
set(sources src/plan_avg.f90)
qe_add_executable(qe_pp_planavg_exe ${sources})
set_target_properties(qe_pp_planavg_exe PROPERTIES OUTPUT_NAME plan_avg.x)
target_link_libraries(qe_pp_planavg_exe
    PRIVATE
        QE::PW
        QE::Modules
        QE::PP)

###########################################################
# plotband.x
###########################################################
set(sources src/plotband.f90)
qe_add_executable(qe_pp_plotband_exe ${sources})
set_target_properties(qe_pp_plotband_exe PROPERTIES OUTPUT_NAME plotband.x)
target_link_libraries(qe_pp_plotband_exe
    PRIVATE
        QE::UTILX
        QE::Modules
        QE::PP)

###########################################################
# plotproj.x
###########################################################
set(sources src/plotproj.f90)
qe_add_executable(qe_pp_plotproj_exe ${sources})
set_target_properties(qe_pp_plotproj_exe PROPERTIES OUTPUT_NAME plotproj.x)
target_link_libraries(qe_pp_plotproj_exe
    PRIVATE
        QE::Modules
        QE::PP)

###########################################################
# plotrho.x
###########################################################
set(sources src/plotrho.f90)
qe_add_executable(qe_pp_plotrho_exe ${sources})
set_target_properties(qe_pp_plotrho_exe PROPERTIES OUTPUT_NAME plotrho.x)
target_link_libraries(qe_pp_plotrho_exe
    PRIVATE
        QE::Modules
        QE::PP)

###########################################################
# pmw.x
###########################################################
set(sources src/poormanwannier.f90)
qe_add_executable(qe_pp_pmw_exe ${sources})
set_target_properties(qe_pp_pmw_exe PROPERTIES OUTPUT_NAME pmw.x)
target_link_libraries(qe_pp_pmw_exe
    PRIVATE
        QE::LAPACK
        QE::PW
        QE::UTILX
        QE::Modules
        QE::PP)

###########################################################
# xctest.x
###########################################################
set(sources src/xctest.f90)
qe_add_executable(qe_pp_xctest_exe ${sources})
set_target_properties(qe_pp_xctest_exe PROPERTIES OUTPUT_NAME xctest.x)
target_link_libraries(qe_pp_xctest_exe
    PRIVATE
        QE::UTILX
        QE::Modules
        QE::CPV
        QE::PP)

###########################################################
# projwfc.x
###########################################################
set(sources src/projwfc.f90)
qe_add_executable(qe_pp_projwfc_exe ${sources})
set_target_properties(qe_pp_projwfc_exe PROPERTIES OUTPUT_NAME projwfc.x)
target_link_libraries(qe_pp_projwfc_exe
    PRIVATE
        QE::LAPACK
        QE::PW
        QE::Modules
        QE::PP)

###########################################################
# pw2wannier90.x
###########################################################
set(sources src/pw2wannier90.f90)
qe_add_executable(qe_pp_pw2wannier90_exe ${sources})
set_target_properties(qe_pp_pw2wannier90_exe PROPERTIES OUTPUT_NAME pw2wannier90.x)
target_link_libraries(qe_pp_pw2wannier90_exe
    PRIVATE
        QE::LAPACK
        QE::PW
        QE::UTILX
        QE::Modules
        QE::PP)

###########################################################
# pw2critic.x
###########################################################
set(sources src/pw2critic.f90)
qe_add_executable(qe_pp_pw2critic_exe ${sources})
set_target_properties(qe_pp_pw2critic_exe PROPERTIES OUTPUT_NAME pw2critic.x)
target_link_libraries(qe_pp_pw2critic_exe
    PRIVATE
        QE::PW
        QE::Modules
        QE::PP)

###########################################################
# wfck2r.x
###########################################################
set(sources src/wfck2r.f90)
qe_add_executable(qe_pp_wfck2r_exe ${sources})
set_target_properties(qe_pp_wfck2r_exe PROPERTIES OUTPUT_NAME wfck2r.x)
target_link_libraries(qe_pp_wfck2r_exe
    PRIVATE
        QE::PW
        QE::Modules
        QE::PP)

###########################################################
# initial_state.x
###########################################################
set(sources src/initial_state.f90)
qe_add_executable(qe_pp_initial_state_exe ${sources})
set_target_properties(qe_pp_initial_state_exe PROPERTIES OUTPUT_NAME initial_state.x)
target_link_libraries(qe_pp_initial_state_exe
    PRIVATE
        QE::PW
        QE::Modules
        QE::PP)

###########################################################
# pw2gw.x
###########################################################
set(sources src/pw2gw.f90)
qe_add_executable(qe_pp_pw2gw_exe ${sources})
set_target_properties(qe_pp_pw2gw_exe PROPERTIES OUTPUT_NAME pw2gw.x)
target_link_libraries(qe_pp_pw2gw_exe
    PRIVATE
        QE::PW
        QE::UTILX
        QE::Modules
        QE::PP)

###########################################################
# sumpdos.x
###########################################################
set(sources src/sumpdos.f90)
qe_add_executable(qe_pp_sumpdos_exe ${sources})
set_target_properties(qe_pp_sumpdos_exe PROPERTIES OUTPUT_NAME sumpdos.x)
target_link_libraries(qe_pp_sumpdos_exe
    PRIVATE
        QE::PP)

###########################################################
# epsilon.x
###########################################################
set(sources src/epsilon.f90)
qe_add_executable(qe_pp_epsilon_exe ${sources})
set_target_properties(qe_pp_epsilon_exe PROPERTIES OUTPUT_NAME epsilon.x)
target_link_libraries(qe_pp_epsilon_exe
    PRIVATE
        QE::PW
        QE::Modules
        QE::PP)

###########################################################
# wannier_ham.x
###########################################################
set(sources src/wannier_ham.f90)
qe_add_executable(qe_pp_wannierham_exe ${sources})
set_target_properties(qe_pp_wannierham_exe PROPERTIES OUTPUT_NAME wannier_ham.x)
target_link_libraries(qe_pp_wannierham_exe
    PRIVATE
        QE::CLIB
        QE::PW
        QE::UTILX
        QE::Modules
        QE::PP)

###########################################################
# wannier_plot.x
###########################################################
set(sources src/wannier_plot.f90)
qe_add_executable(qe_pp_wannierplot_exe ${sources})
set_target_properties(qe_pp_wannierplot_exe PROPERTIES OUTPUT_NAME wannier_plot.x)
target_link_libraries(qe_pp_wannierplot_exe
    PRIVATE
        QE::PW
        QE::UTILX
        QE::Modules
        QE::PP)

###########################################################
# molecularpdos.x
###########################################################
set(sources src/molecularpdos.f90)
qe_add_executable(qe_pp_molecularpdos_exe ${sources})
set_target_properties(qe_pp_molecularpdos_exe PROPERTIES OUTPUT_NAME molecularpdos.x)
target_link_libraries(qe_pp_molecularpdos_exe
    PRIVATE
        QE::UTILX
        QE::Modules
        QE::PP)

###########################################################
# pw2bgw.x
###########################################################
set(sources src/pw2bgw.f90)
qe_add_executable(qe_pp_pw2bgw_exe ${sources})
set_target_properties(qe_pp_pw2bgw_exe PROPERTIES OUTPUT_NAME pw2bgw.x)
target_link_libraries(qe_pp_pw2bgw_exe
    PRIVATE
        QE::PW
        QE::Modules
        QE::PP)

###########################################################
# fs.x
###########################################################
set(sources src/fermisurface.f90)
qe_add_executable(qe_pp_fermisurface_exe ${sources})
set_target_properties(qe_pp_fermisurface_exe PROPERTIES OUTPUT_NAME fs.x)
target_link_libraries(qe_pp_fermisurface_exe
    PRIVATE
        QE::PW
        QE::UTILX
        QE::Modules
        QE::PP)

###########################################################
# fermi_velocity.x
###########################################################
set(sources src/fermi_velocity.f90)
qe_add_executable(qe_pp_fermivelocity_exe ${sources})
set_target_properties(qe_pp_fermivelocity_exe PROPERTIES OUTPUT_NAME fermi_velocity.x)
target_link_libraries(qe_pp_fermivelocity_exe
    PRIVATE
        QE::PW
        QE::UTILX
        QE::Modules
        QE::PP)

###########################################################
# fermi_proj.x
###########################################################
set(sources src/fermi_proj.f90)
qe_add_executable(qe_pp_fermiproj_exe ${sources})
set_target_properties(qe_pp_fermiproj_exe PROPERTIES OUTPUT_NAME fermi_proj.x)
target_link_libraries(qe_pp_fermiproj_exe
    PRIVATE
        QE::PW
        QE::Modules
        QE::PP)

###########################################################
# ppacf.x
###########################################################
set(sources src/ppacf.f90)
qe_add_executable(qe_pp_ppacf_exe ${sources})
set_target_properties(qe_pp_ppacf_exe PROPERTIES OUTPUT_NAME ppacf.x)
target_link_libraries(qe_pp_ppacf_exe
    PRIVATE
        QE::LAPACK
        QE::PW
        QE::PP)

###########################################################
# benchmark_libxc.x
###########################################################
set(sources src/xctest_qe_libxc.f90)
qe_add_executable(qe_pp_benchmarklibxc_exe ${sources})
set_target_properties(qe_pp_benchmarklibxc_exe PROPERTIES OUTPUT_NAME xctest_qe_libxc.x)
target_link_libraries(qe_pp_benchmarklibxc_exe
    PRIVATE
        QE::PP)

###########################################################
# ef.x
###########################################################
set(sources 
    simple_transport/src/cryst_to_car.f90
    simple_transport/src/lint.f90
    simple_transport/src/ef.f90
)
qe_add_executable(qe_pp_st_ef_exe ${sources})
set_target_properties(qe_pp_st_ef_exe PROPERTIES OUTPUT_NAME ef.x)
target_link_libraries(qe_pp_st_ef_exe
    PRIVATE
        QE::OpenMP_Fortran
        QE::Modules)

###########################################################
# dos_sp.x
###########################################################
set(sources 
    simple_transport/src/smearing_mod.f90 
    simple_transport/src/cryst_to_car.f90 
    simple_transport/src/lint.f90
    simple_transport/src/vband_ibz.f90 
    simple_transport/src/dos.f90
)
qe_add_executable(qe_pp_st_dos_exe ${sources})
set_target_properties(qe_pp_st_dos_exe PROPERTIES OUTPUT_NAME dos_sp.x)
target_link_libraries(qe_pp_st_dos_exe
    PRIVATE
        QE::OpenMP_Fortran
        QE::Modules)

###########################################################
# fermi_int_0.x
###########################################################
set(sources 
    simple_transport/src/cryst_to_car.f90 
    simple_transport/src/lint.f90 
    simple_transport/src/vband_ibz.f90 
    simple_transport/src/fermi_int_0.f90
)
qe_add_executable(qe_pp_st_fermiint0_exe ${sources})
set_target_properties(qe_pp_st_fermiint0_exe PROPERTIES OUTPUT_NAME fermi_int_0.x)
target_link_libraries(qe_pp_st_ef_exe
    PRIVATE
        QE::OpenMP_Fortran)

###########################################################
# fermi_int_1.x
###########################################################
set(sources 
    simple_transport/src/cryst_to_car.f90 
    simple_transport/src/reducegrid.f90 
    simple_transport/src/lint.f90 
    simple_transport/src/vband_ibz.f90 
    simple_transport/src/fermi_int_1.f90
)
qe_add_executable(qe_pp_st_fermiint1_exe ${sources})
set_target_properties(qe_pp_st_fermiint1_exe PROPERTIES OUTPUT_NAME fermi_int_1.x)
target_link_libraries(qe_pp_st_ef_exe
    PRIVATE
        QE::OpenMP_Fortran)

###########################################################

qe_install_targets(
    # Libraries
    qe_pp
    # Executables
    qe_pp_exe
    qe_pp_opengrid_exe
    qe_pp_average_exe
    qe_pp_bands_exe
    qe_pp_dos_exe
    qe_pp_pawplot_exe
    qe_pp_planavg_exe
    qe_pp_plotband_exe
    qe_pp_plotproj_exe
    qe_pp_plotrho_exe
    qe_pp_pmw_exe
    qe_pp_xctest_exe
    qe_pp_projwfc_exe
    qe_pp_pw2wannier90_exe
    qe_pp_pw2critic_exe
    qe_pp_wfck2r_exe
    qe_pp_initial_state_exe
    qe_pp_pw2gw_exe
    qe_pp_sumpdos_exe
    qe_pp_epsilon_exe
    qe_pp_wannierham_exe
    qe_pp_wannierplot_exe
    qe_pp_molecularpdos_exe
    qe_pp_pw2bgw_exe
    qe_pp_fermivelocity_exe
    qe_pp_fermisurface_exe
    qe_pp_fermiproj_exe
    qe_pp_ppacf_exe
    qe_pp_benchmarklibxc_exe
    # Simple Transport
    qe_pp_st_ef_exe
    qe_pp_st_dos_exe
    qe_pp_st_fermiint0_exe
    qe_pp_st_fermiint1_exe
)
