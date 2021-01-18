set(src_pp
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
qe_enable_cuda_fortran("${src_pp}")

qe_add_library(qe_pp ${src_pp})
target_compile_definitions(qe_pp PRIVATE "__WANLIB")
target_link_libraries(qe_pp
    PRIVATE
        qe_openmp_fortran
        qe_pw
        qe_modules
        qe_fftx
        qe_upflib
        qe_lax
        qe_xclib)

###########################################################
# pp.x
###########################################################
set(src_pp_x src/postproc.f90)
qe_add_executable(qe_pp_exe ${src_pp_x})
set_target_properties(qe_pp_exe PROPERTIES OUTPUT_NAME pp.x)
target_link_libraries(qe_pp_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_pp
        qe_fftx
        qe_xclib)

###########################################################
# open_grid.x
###########################################################
set(src_open_grid_x src/open_grid.f90)
qe_add_executable(qe_pp_opengrid_exe ${src_open_grid_x})
set_target_properties(qe_pp_opengrid_exe PROPERTIES OUTPUT_NAME open_grid.x)
target_link_libraries(qe_pp_opengrid_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_pp
        qe_fftx
        qe_upflib
        qe_xclib)

###########################################################
# average.x
###########################################################
set(src_average_x src/average.f90)
qe_add_executable(qe_pp_average_exe ${src_average_x})
set_target_properties(qe_pp_average_exe PROPERTIES OUTPUT_NAME average.x)
target_link_libraries(qe_pp_average_exe
    PRIVATE
        qe_pw
        qe_fftx
        qe_modules
        qe_pp
        qe_upflib
        qe_xclib)

###########################################################
# bands.x
###########################################################
set(src_bands_x src/bands.f90)
qe_add_executable(qe_pp_bands_exe ${src_bands_x})
set_target_properties(qe_pp_bands_exe PROPERTIES OUTPUT_NAME bands.x)
target_link_libraries(qe_pp_bands_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_pp
        qe_upflib)

###########################################################
# dos.x
###########################################################
set(src_dos_x src/dos.f90)
qe_add_executable(qe_pp_dos_exe ${src_dos_x})
set_target_properties(qe_pp_dos_exe PROPERTIES OUTPUT_NAME dos.x)
target_link_libraries(qe_pp_dos_exe
    PRIVATE
        qe_pw
        qe_pp
        qe_modules)

###########################################################
# pawplot.x
###########################################################
set(src_pawplot_x src/pawplot.f90)
qe_add_executable(qe_pp_pawplot_exe ${src_pawplot_x})
set_target_properties(qe_pp_pawplot_exe PROPERTIES OUTPUT_NAME pawplot.x)
target_link_libraries(qe_pp_pawplot_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_pp
        qe_upflib
        qe_fftx
        qe_xclib)

###########################################################
# plan_avg.x
###########################################################
set(src_plan_avg_x src/plan_avg.f90)
qe_add_executable(qe_pp_planavg_exe ${src_plan_avg_x})
set_target_properties(qe_pp_planavg_exe PROPERTIES OUTPUT_NAME plan_avg.x)
target_link_libraries(qe_pp_planavg_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_pp
        qe_upflib
        qe_fftx)

###########################################################
# plotband.x
###########################################################
set(src_plotband_x src/plotband.f90)
qe_add_executable(qe_pp_plotband_exe ${src_plotband_x})
set_target_properties(qe_pp_plotband_exe PROPERTIES OUTPUT_NAME plotband.x)
target_link_libraries(qe_pp_plotband_exe
    PRIVATE
        qe_modules
        qe_pp)

###########################################################
# plotproj.x
###########################################################
set(src_plotproj_x src/plotproj.f90)
qe_add_executable(qe_pp_plotproj_exe ${src_plotproj_x})
set_target_properties(qe_pp_plotproj_exe PROPERTIES OUTPUT_NAME plotproj.x)
target_link_libraries(qe_pp_plotproj_exe
    PRIVATE
        qe_modules
        qe_pp)

###########################################################
# plotrho.x
###########################################################
set(src_plotrho_x src/plotrho.f90)
qe_add_executable(qe_pp_plotrho_exe ${src_plotrho_x})
set_target_properties(qe_pp_plotrho_exe PROPERTIES OUTPUT_NAME plotrho.x)
target_link_libraries(qe_pp_plotrho_exe
    PRIVATE
        qe_modules
        qe_pp)

###########################################################
# pmw.x
###########################################################
set(src_pmw_x src/poormanwannier.f90)
qe_add_executable(qe_pp_pmw_exe ${src_pmw_x})
set_target_properties(qe_pp_pmw_exe PROPERTIES OUTPUT_NAME pmw.x)
target_link_libraries(qe_pp_pmw_exe
    PRIVATE
        qe_lapack
        qe_pw
        qe_modules
        qe_pp
        qe_upflib)

###########################################################
# xctest.x
###########################################################
set(src_xctest_x src/xctest.f90)
qe_add_executable(qe_pp_xctest_exe ${src_xctest_x})
set_target_properties(qe_pp_xctest_exe PROPERTIES OUTPUT_NAME xctest.x)
target_link_libraries(qe_pp_xctest_exe
    PRIVATE
        qe_modules
        qe_cpv
        qe_pp
        qe_xclib)

###########################################################
# projwfc.x
###########################################################
set(src_projwfc_x src/projwfc.f90)
qe_add_executable(qe_pp_projwfc_exe ${src_projwfc_x})
set_target_properties(qe_pp_projwfc_exe PROPERTIES OUTPUT_NAME projwfc.x)
target_link_libraries(qe_pp_projwfc_exe
    PRIVATE
        qe_lapack
        qe_pw
        qe_modules
        qe_pp
        qe_lax
        qe_upflib)

###########################################################
# pw2wannier90.x
###########################################################
set(src_pw2wannier90_x src/pw2wannier90.f90)
qe_add_executable(qe_pp_pw2wannier90_exe ${src_pw2wannier90_x})
set_target_properties(qe_pp_pw2wannier90_exe PROPERTIES OUTPUT_NAME pw2wannier90.x)
target_link_libraries(qe_pp_pw2wannier90_exe
    PRIVATE
        qe_lapack
        qe_pw
        qe_modules
        qe_pp
        qe_fftx
        qe_upflib
        qe_xclib)

###########################################################
# pw2critic.x
###########################################################
set(src_pw2critic_x src/pw2critic.f90)
qe_add_executable(qe_pp_pw2critic_exe ${src_pw2critic_x})
set_target_properties(qe_pp_pw2critic_exe PROPERTIES OUTPUT_NAME pw2critic.x)
target_link_libraries(qe_pp_pw2critic_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_pp
        qe_fftx)

###########################################################
# wfck2r.x
###########################################################
set(src_wfck2r_x src/wfck2r.f90)
qe_add_executable(qe_pp_wfck2r_exe ${src_wfck2r_x})
set_target_properties(qe_pp_wfck2r_exe PROPERTIES OUTPUT_NAME wfck2r.x)
target_link_libraries(qe_pp_wfck2r_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_pp
        qe_fftx)

###########################################################
# initial_state.x
###########################################################
set(src_initial_state_x src/initial_state.f90)
qe_add_executable(qe_pp_initial_state_exe ${src_initial_state_x})
set_target_properties(qe_pp_initial_state_exe PROPERTIES OUTPUT_NAME initial_state.x)
target_link_libraries(qe_pp_initial_state_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_pp
        qe_upflib)

###########################################################
# pw2gw.x
###########################################################
set(src_pw2gw_x src/pw2gw.f90)
qe_add_executable(qe_pp_pw2gw_exe ${src_pw2gw_x})
set_target_properties(qe_pp_pw2gw_exe PROPERTIES OUTPUT_NAME pw2gw.x)
target_link_libraries(qe_pp_pw2gw_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_pp
        qe_fftx
        qe_upflib
        qe_mpi_fortran
        qe_xclib)

###########################################################
# sumpdos.x
###########################################################
set(src_sumpdos_x src/sumpdos.f90)
qe_add_executable(qe_pp_sumpdos_exe ${src_sumpdos_x})
set_target_properties(qe_pp_sumpdos_exe PROPERTIES OUTPUT_NAME sumpdos.x)
target_link_libraries(qe_pp_sumpdos_exe
    PRIVATE
        qe_pp)

###########################################################
# epsilon.x
###########################################################
set(src_epsilon_x src/epsilon.f90)
qe_add_executable(qe_pp_epsilon_exe ${src_epsilon_x})
set_target_properties(qe_pp_epsilon_exe PROPERTIES OUTPUT_NAME epsilon.x)
target_link_libraries(qe_pp_epsilon_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_pp
        qe_upflib)

###########################################################
# wannier_ham.x
###########################################################
set(src_wannier_ham_x src/wannier_ham.f90)
qe_add_executable(qe_pp_wannierham_exe ${src_wannier_ham_x})
set_target_properties(qe_pp_wannierham_exe PROPERTIES OUTPUT_NAME wannier_ham.x)
target_link_libraries(qe_pp_wannierham_exe
    PRIVATE
        qe_clib
        qe_pw
        qe_modules
        qe_pp)

###########################################################
# wannier_plot.x
###########################################################
set(src_wannier_plot_x src/wannier_plot.f90)
qe_add_executable(qe_pp_wannierplot_exe ${src_wannier_plot_x})
set_target_properties(qe_pp_wannierplot_exe PROPERTIES OUTPUT_NAME wannier_plot.x)
target_link_libraries(qe_pp_wannierplot_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_pp
        qe_fftx)

###########################################################
# molecularpdos.x
###########################################################
set(src_molecularpdos_x src/molecularpdos.f90)
qe_add_executable(qe_pp_molecularpdos_exe ${src_molecularpdos_x})
set_target_properties(qe_pp_molecularpdos_exe PROPERTIES OUTPUT_NAME molecularpdos.x)
target_link_libraries(qe_pp_molecularpdos_exe
    PRIVATE
        qe_modules
        qe_pp)

###########################################################
# pw2bgw.x
###########################################################
set(src_pw2bgw_x src/pw2bgw.f90)
qe_add_executable(qe_pp_pw2bgw_exe ${src_pw2bgw_x})
set_target_properties(qe_pp_pw2bgw_exe PROPERTIES OUTPUT_NAME pw2bgw.x)
target_link_libraries(qe_pp_pw2bgw_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_pp
        qe_upflib
        qe_fftx
        qe_mpi_fortran
        qe_xclib)

###########################################################
# fs.x
###########################################################
set(src_fs_x src/fermisurface.f90)
qe_add_executable(qe_pp_fermisurface_exe ${src_fs_x})
set_target_properties(qe_pp_fermisurface_exe PROPERTIES OUTPUT_NAME fs.x)
target_link_libraries(qe_pp_fermisurface_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_pp)

###########################################################
# fermi_velocity.x
###########################################################
set(src_fermi_velocity_x src/fermi_velocity.f90)
qe_add_executable(qe_pp_fermivelocity_exe ${src_fermi_velocity_x})
set_target_properties(qe_pp_fermivelocity_exe PROPERTIES OUTPUT_NAME fermi_velocity.x)
target_link_libraries(qe_pp_fermivelocity_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_pp)

###########################################################
# fermi_proj.x
###########################################################
set(src_fermi_proj_x src/fermi_proj.f90)
qe_add_executable(qe_pp_fermiproj_exe ${src_fermi_proj_x})
set_target_properties(qe_pp_fermiproj_exe PROPERTIES OUTPUT_NAME fermi_proj.x)
target_link_libraries(qe_pp_fermiproj_exe
    PRIVATE
        qe_pw
        qe_modules
        qe_pp
        qe_upflib)

###########################################################
# ppacf.x
###########################################################
set(src_ppacf_x src/ppacf.f90)
qe_add_executable(qe_pp_ppacf_exe ${src_ppacf_x})
set_target_properties(qe_pp_ppacf_exe PROPERTIES OUTPUT_NAME ppacf.x)
target_link_libraries(qe_pp_ppacf_exe
    PRIVATE
        qe_lapack
        qe_pw
        qe_pp
        qe_modules
        qe_fftx
        qe_upflib
        qe_xclib)

###########################################################
# ef.x
###########################################################
set(src_ef_x 
    simple_transport/src/cryst_to_car.f90
    simple_transport/src/lint.f90
    simple_transport/src/ef.f90)
qe_add_executable(qe_pp_st_ef_exe ${src_ef_x})
set_target_properties(qe_pp_st_ef_exe PROPERTIES OUTPUT_NAME ef.x)
target_link_libraries(qe_pp_st_ef_exe
    PRIVATE
        qe_openmp_fortran
        qe_modules)

###########################################################
# dos_sp.x
###########################################################
set(src_dos_sp_x
    simple_transport/src/smearing_mod.f90 
    simple_transport/src/cryst_to_car.f90 
    simple_transport/src/lint.f90
    simple_transport/src/vband_ibz.f90 
    simple_transport/src/dos.f90)
qe_add_executable(qe_pp_st_dos_exe ${src_dos_sp_x})
set_target_properties(qe_pp_st_dos_exe PROPERTIES OUTPUT_NAME dos_sp.x)
target_link_libraries(qe_pp_st_dos_exe
    PRIVATE
        qe_openmp_fortran
        qe_modules)

###########################################################
# fermi_int_0.x
###########################################################
set(src_fermi_int_0_x
    simple_transport/src/cryst_to_car.f90 
    simple_transport/src/lint.f90 
    simple_transport/src/vband_ibz.f90 
    simple_transport/src/fermi_int_0.f90)
qe_add_executable(qe_pp_st_fermiint0_exe ${src_fermi_int_0_x})
set_target_properties(qe_pp_st_fermiint0_exe PROPERTIES OUTPUT_NAME fermi_int_0.x)
target_link_libraries(qe_pp_st_ef_exe
    PRIVATE
        qe_openmp_fortran)

###########################################################
# fermi_int_1.x
###########################################################
set(src_fermi_int_1_x
    simple_transport/src/cryst_to_car.f90 
    simple_transport/src/reducegrid.f90 
    simple_transport/src/lint.f90 
    simple_transport/src/vband_ibz.f90 
    simple_transport/src/fermi_int_1.f90)
qe_add_executable(qe_pp_st_fermiint1_exe ${src_fermi_int_1_x})
set_target_properties(qe_pp_st_fermiint1_exe PROPERTIES OUTPUT_NAME fermi_int_1.x)
target_link_libraries(qe_pp_st_ef_exe
    PRIVATE
        qe_openmp_fortran)

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
    # Simple Transport
    qe_pp_st_ef_exe
    qe_pp_st_dos_exe
    qe_pp_st_fermiint0_exe
    qe_pp_st_fermiint1_exe)
