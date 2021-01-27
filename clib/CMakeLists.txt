set(sources
    customize_signals.c
    stack.c
    c_mkdir.c
    copy.c
    cptimer.c
    eval_infix.c
    fletcher32.c
    md5.c
    md5_from_file.c
    memstat.c
    ptrace.c
    sockets.c
    qmmm_aux.c
    qsort.c)

qe_add_library(qe_clib ${sources})

###########################################################

qe_install_targets(qe_clib)