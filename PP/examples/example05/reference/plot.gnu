#!/usr/bin/gnuplot

set xrange [0.0:3.1868]
unset xtics
set yrange [-21.0:15.0]
set ylabel 'E - E_F  (eV)'

set border lw 2
set style arrow 1 nohead front lw 2 lc rgb 'black'
set label 'G' at graph 0,graph -0.03 center
set arrow from 0.8292,graph 0 to 0.8292,graph 1 as 1
set label 'L' at 0.8292, graph -0.03 center
set arrow from 1.1223,graph 0 to 1.1223,graph 1 as 1
set label "K'" at 1.1223,graph -0.03 center
set arrow from 1.8878,graph 0 to 1.8878,graph 1 as 1
set label 'T' at 1.8878, graph -0.03 center
set arrow from 2.3208,graph 0 to 2.3208,graph 1 as 1
set label 'G' at 2.3208, graph -0.03 center
set label 'X' at graph 1,graph -0.03 center

plot 'feo_af.bands.all' u 1:($2 - 10.9920):3 w l palette lw 2 notitle, \
     0.0 lt 1 lw 2 lc rgb 'grey50' notitle

