set(sources
    pbecor.c 
    beefun.c
)

qe_add_library(qe_libbeef ${sources})
add_library(QE::libbeef ALIAS qe_libbeef)
if(CMAKE_C_COMPILER_ID STREQUAL "PGI")
	target_compile_options(qe_libbeef PUBLIC "-c11") 
    
    endif()


###########################################################

qe_install_targets(qe_libbeef)
