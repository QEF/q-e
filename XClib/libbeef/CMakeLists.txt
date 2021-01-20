set(sources_libbeef
    pbecor.c
    beefun.c
)

qe_add_library(qe_libbeef ${sources_libbeef})
target_link_libraries(qe_libbeef PRIVATE qe_lapack)

qe_install_targets(qe_libbeef)
