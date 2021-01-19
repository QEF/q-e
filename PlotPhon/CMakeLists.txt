###########################################################
# k_for_bands.x
###########################################################
set(src_k_for_bands_x SRC/K_for_bands.f90)
qe_add_executable(qe_plotphon_kforbands_exe ${src_k_for_bands_x})
set_target_properties(qe_plotphon_kforbands_exe 
    PROPERTIES 
        OUTPUT_NAME k_for_bands.x)

###########################################################
# k_for_bands.x
###########################################################
set(src_k_for_bands_x SRC/Bands_to_gnuplot.f90)
qe_add_executable(qe_plotphon_bandstognuplot_exe ${src_k_for_bands_x})
set_target_properties(qe_plotphon_bandstognuplot_exe 
    PROPERTIES 
        OUTPUT_NAME bands_to_gnuplot.x)

###########################################################
# E_min_max.x
###########################################################
set(src_E_min_max_x SRC/E_min_max.f90)
qe_add_executable(qe_plotphon_eminmax_exe ${src_E_min_max_x})
set_target_properties(qe_plotphon_eminmax_exe 
    PROPERTIES 
        OUTPUT_NAME E_min_max.x)

###########################################################

qe_install_targets(
    # Executables
    qe_plotphon_kforbands_exe
    qe_plotphon_bandstognuplot_exe
    qe_plotphon_eminmax_exe)
