set(sources
    pbecor.c 
    beefun.c
)

qe_add_library(qe_libbeef ${sources})
add_library(QE::libbeef ALIAS qe_libbeef)

###########################################################

qe_install_targets(qe_libbeef)
