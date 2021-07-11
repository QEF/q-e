set xrange [0.4:0.55]
# set arrow
set key left
set arrow from 0.442,0.042 to 0.497,0.042 heads
set label "ZPR = 55 meV" at 0.47,0.044
# ZPR stands for zero-point-renormalization
plot "JDOS_Gaus_555_equil.dat" u 1:(sqrt($2/16473383)) t "Equil" w l lw 4,\
     "JDOS_Gaus_555_0K.dat" u 1:(sqrt($2/19203399)) t "0K" w l lw 4

# Normalize by Dividing with the number transitions involved in that particular enegy range
