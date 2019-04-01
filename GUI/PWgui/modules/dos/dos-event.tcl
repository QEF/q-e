tracevar bz_sum w {
    switch -glob -- [varvalue bz_sum] {
        *tetrahedra* {
            widget ngauss  disable
            widget degauss disable
        }
        default {
            widget ngauss  enable
            widget degauss enable
        }
    }
}
varset bz_sum -value {}
