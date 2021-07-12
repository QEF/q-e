set xrange [0.4:0.55]
# set arrow
set key left
set arrow from 0.45,0.042 to 0.5,0.042 heads
set label "ZPR = 50 meV" at 0.47,0.044
# ZPR stands for zero-point-renormalization
plot "JDOS_Gaus_333_equil.dat" u 1:(sqrt($2/19562791)) t "Equil" w l lw 4,\
     "JDOS_Gaus_333_0K.dat" u 1:(sqrt($2/22841581)) t "0K" w l lw 4
