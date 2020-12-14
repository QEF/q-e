set(sources_gww
    gww/basic_structures.f90
    gww/compact_product.f90
    gww/contour.f90
    gww/create_hf.f90
    gww/create_quasi_particle.f90
    gww/create_quasi_particle_off.f90
    gww/do_contour.f90
    gww/do_polarization_lanczos.f90
    gww/do_self_lanczos.f90
    gww/do_self_lanczos_time.f90
    gww/do_self_lanczos_full.f90
    gww/energies_gww.f90
    gww/expansion.f90
    gww/fft_gw.f90
    gww/fit_multipole.f90
    gww/fit_polynomial.f90
    gww/go_dressed_w.f90
    gww/go_exchange.f90
    gww/go_fft.f90
    gww/go_green.f90
    gww/go_polarization.f90
    gww/green_function.f90
    gww/gv_time.f90
    gww/input_gw.f90
    gww/lanczos_polarization.f90
    gww/para_gww.f90
    gww/polarization.f90
    gww/read_data_pw.f90
    gww/self_energy.f90
    gww/self_energy_storage.f90
    gww/start_end.f90
    gww/times_gw.f90
    gww/vcprim.f90)

set(sources_pw4gww
    pw4gww/allocate_wannier.f90
    pw4gww/produce_wannier_gamma.f90
    pw4gww/start_pw4gww.f90
    pw4gww/stop_pp.f90
    pw4gww/openfil_pw4gww.f90
    pw4gww/dft_exchange.f90
    pw4gww/wfc_real.f90
    pw4gww/full.f90
    pw4gww/energies_xc.f90
    pw4gww/wannier.f90
    pw4gww/write_wannier_matrix.f90
    pw4gww/rotate_wannier.f90
    pw4gww/self_lanczos.f90
    pw4gww/pola_lanczos.f90
    pw4gww/matrix_wannier_gamma.f90
    pw4gww/calculate_wing.f90
    pw4gww/fake_conduction.f90
    pw4gww/o_1psi.f90
    pw4gww/o_bands.f90
    pw4gww/optimal.f90
    pw4gww/pola_partial.f90
    pw4gww/semicore.f90
    pw4gww/semicore_read.f90
    pw4gww/v_basis.f90
    pw4gww/wannier_uterms.f90
    pw4gww/write_vpot_matrix.f90
    pw4gww/diago_cg.f90
    pw4gww/o_rinitcgg.f90
    pw4gww/diago_cg_g.f90
    pw4gww/contour_terms.f90
    pw4gww/fft_custom.f90
    pw4gww/exchange_custom.f90
    pw4gww/mp_wave_parallel.f90
    pw4gww/wannier_bse.f90
    pw4gww/hpsi_pw4gww.f90
    pw4gww/cgsolve_all_gamma.f90
    pw4gww/realus.f90)

set(sources_bse
    bse/openfil_bse.f90
    bse/start_bse.f90
    bse/bse_wannier.f90
    bse/transitions.f90
    bse/read_export.f90
    bse/bse_basic_structure.f90
    bse/exciton.f90
    bse/diago_exc.f90
    bse/exchange_exc.f90
    bse/rotate_wannier.f90
    bse/write_wannier_matrix.f90
    bse/direct_v_exc.f90
    bse/direct_w_exc.f90
    bse/direct_www.f90
    bse/exc_h_a.f90
    bse/sdescent.f90
    bse/check_basis.f90
    bse/conjgrad.f90
    bse/conj_grad_stuff.f90
    bse/find_eig.f90
    bse/plot_excwfn.f90
    bse/tspace.f90
    bse/absorption.f90
    bse/spectrum.f90
    bse/absorption.f90
    bse/cgsolve.f90
    bse/h_h.f90
    bse/dvpsi_bse.f90
    bse/lanczos.f90
    bse/print_spectrum.f90
    bse/print_bse.f90
    bse/contract_w_exc.f90
    bse/qpe_exc.f90
    bse/qpcorrections.f90)

set(sources_head
    head/bcast_ph_input.f90
    head/close_phq.f90
    head/lanczos_k.f90
    head/openfilq.f90
    head/phq_readin.f90
    head/solve_head.f90)

set(sources_minpack
    minpack/dpmpar.f90
    minpack/fdjac2.f90
    minpack/lmder1.f90
    minpack/lmdif1.f90
    minpack/mlagzo.f90
    minpack/qrfac.f90
    minpack/test.f90
    minpack/enorm.f90
    minpack/lmder.f90
    minpack/lmdif.f90
    minpack/lmpar.f90
    minpack/mlegzo.f90
    minpack/qrsolv.f90)

set(sources_simple
    simple/stop_pp.f90
    simple/read_export.f90
    simple/openfile_simple.f90
    simple/input_simple.f90
    simple/wfc_basis.f90
    simple/product_basis.f90
    simple/v_product.f90
    simple/epe.f90
    simple/gk_sort_limit.f90
    simple/khamiltonian.f90
    simple/init_us_2_max.f90
    simple/commutator.f90)

set(sources_simple_bse
    simple_bse/simple_object.f90
    simple_bse/start_end.f90
    simple_bse/input_simple_exc.f90
    simple_bse/simple_eigen.f90
    simple_bse/diago_exc_sd.f90
    simple_bse/hamiltonian.f90
    simple_bse/derived_object.f90
    simple_bse/diago_exc_cg.f90
    simple_bse/spectrum.f90
    simple_bse/lanczos.f90
    simple_bse/build_eemat.f90)

set(sources_simple_ip
    simple_ip/simple_ip_objects.f90
    simple_ip/start_end.f90
    simple_ip/input_simple_ip.f90
    simple_ip/tetra_mod1.f90
    simple_ip/interpolation.f90
    simple_ip/diagonalization.f90
    simple_ip/dielectric.f90)

qe_add_library(qe_gww ${sources_gww})
add_library(QE::GWW ALIAS qe_gww)
target_link_libraries(qe_gww
    PRIVATE
        QE::FFTX
    PUBLIC
        QE::UTILX
	QE::Modules)

qe_add_library(qe_gww_pw4gww ${sources_pw4gww})
add_library(QE::GWW::pw4gww ALIAS qe_gww_pw4gww)
target_link_libraries(qe_gww_pw4gww
    PRIVATE
        QE::SCALAPACK
        QE::PW
        QE::GWW
        QE::FFTX
        QE::UPF
        QE::LAX)

qe_add_library(qe_gww_bse ${sources_bse})
add_library(QE::GWW::bse ALIAS qe_gww_bse)
target_link_libraries(qe_gww_bse
    PUBLIC
        QE::GWW
        QE::PW
        QE::GWW::pw4gww
        QE::UPF)

qe_add_library(qe_gww_head ${sources_head})
add_library(QE::GWW::head ALIAS qe_gww_head)
target_link_libraries(qe_gww_head
    PRIVATE
        QE::PW
        QE::LR_Modules
        QE::Modules
        QE::PHonon::PH
        QE::UTILX
        QE::UPF
        QE::FFTX)

qe_add_library(qe_gww_minpack ${sources_minpack})
add_library(QE::GWW::minpack ALIAS qe_gww_minpack)
target_link_libraries(qe_gww_minpack
    PRIVATE
        QE::Modules)

qe_add_library(qe_gww_simple ${sources_simple})
add_library(QE::GWW::simple ALIAS qe_gww_simple)
target_link_libraries(qe_gww_simple
    PRIVATE
        QE::PW
        QE::GWW
        QE::UPF
        QE::FFTX)

qe_add_library(qe_gww_simplebse ${sources_simple_bse})
add_library(QE::GWW::simple_bse ALIAS qe_gww_simplebse)
target_link_libraries(qe_gww_simplebse
    PRIVATE
        QE::UTILX
        QE::Modules)

qe_add_library(qe_gww_simpleip ${sources_simple_ip})
add_library(QE::GWW::simple_ip ALIAS qe_gww_simpleip)
target_link_libraries(qe_gww_simpleip
    PRIVATE
        QE::UTILX
        QE::Modules
        QE::PW
)

###########################################################
# graph.x
###########################################################
set(sources util/graph.f90)
qe_add_executable(qe_gww_util_grap_exe ${sources})
set_target_properties(qe_gww_util_grap_exe PROPERTIES OUTPUT_NAME graph.x)

###########################################################
# abcoeff_to_eps.x
###########################################################
set(sources util/abcoeff_to_eps.f90)
qe_add_executable(qe_gww_util_abcoefftoeps_exe ${sources})
set_target_properties(qe_gww_util_abcoefftoeps_exe PROPERTIES OUTPUT_NAME abcoeff_to_eps.x)

###########################################################
# memory_pw4gww.x
###########################################################
set(sources util/memory_pw4gww.f90)
qe_add_executable(qe_gww_util_memorypw4gww_exe ${sources})
set_target_properties(qe_gww_util_memorypw4gww_exe PROPERTIES OUTPUT_NAME memory_pw4gww.x)

#############################################################
#bse_main.x
############################################################
set(sources bse/bse_main.f90)
qe_add_executable(qe_gww_bse_bse_main_exe ${sources})
set_target_properties(qe_gww_bse_bse_main_exe PROPERTIES OUTPUT_NAME bse_main.x)
target_link_libraries(qe_gww_bse_bse_main_exe
    PRIVATE
       QE::GWW::bse )
################################################################
#gww.x
################################################################
set (sources gww/gww.f90)
qe_add_executable(qe_gww_gww_exe ${sources})
set_target_properties(qe_gww_gww_exe PROPERTIES OUTPUT_NAME gww.x)
target_link_libraries(qe_gww_gww_exe
       PRIVATE
       QE::GWW
       QE::GWW::minpack
       QE::PW
       )

################################################################
#gww_fit.x 
################################################################
set (sources gww/gww_fit.f90)
qe_add_executable(qe_gww_gww_fit_exe  ${sources})
set_target_properties(qe_gww_gww_fit_exe PROPERTIES OUTPUT_NAME gww.x)
target_link_libraries(qe_gww_gww_fit_exe
       PRIVATE 
       QE::GWW
       QE::GWW::minpack
       )

################################################################
#head.x 
################################################################
set (sources head/head.f90)
qe_add_executable(qe_gww_head_exe ${sources})
set_target_properties(qe_gww_head_exe PROPERTIES OUTPUT_NAME head.x)
target_link_libraries(qe_gww_head_exe
       PRIVATE 
       QE::GWW::head
       QE::GWW::minpack
       QE::PHonon::PH
       QE::Modules
       )

###########################################################
#simple_bse.x
###########################################################
set (sources simple_bse/simple_bse.f90)
qe_add_executable(qe_gww_simple_bse_exe ${sources})
set_target_properties(qe_gww_simple_bse_exe PROPERTIES OUTPUT_NAME simple_bse.x)
target_link_libraries(qe_gww_simple_bse_exe
    PRIVATE
        QE::GWW::simple_bse
)

###############################################################
#simple_ip.x
################################################################
set (sources simple_ip/simple_ip.f90)
qe_add_executable(qe_gww_simple_ip_exe ${sources}) 
set_target_properties(qe_gww_simple_ip_exe PROPERTIES OUTPUT_NAME simple_ip.x)
target_link_libraries(qe_gww_simple_ip_exe
    PRIVATE
        QE::GWW::simple_ip
)

###########################################################

qe_install_targets(
    # Libraries
    qe_gww
    qe_gww_pw4gww
    qe_gww_bse
    qe_gww_head
    qe_gww_minpack
    qe_gww_simple
    qe_gww_simplebse
    qe_gww_simpleip
    # Executables
    qe_gww_util_grap_exe
    qe_gww_util_abcoefftoeps_exe
    qe_gww_util_memorypw4gww_exe)
