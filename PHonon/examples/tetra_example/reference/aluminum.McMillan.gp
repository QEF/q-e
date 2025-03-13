unset key
set xlabel "mu*"
set ylabel "T_c [K]"
set xrange [0.0:    0.31966E+00]
plot     0.25429E+03*exp(   -0.14546E+01/(    0.39867E+00-x*    0.12472E+01))
pause -1
