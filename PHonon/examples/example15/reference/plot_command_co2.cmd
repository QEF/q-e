 set lmargin 8
 set rmargin 3
 
 set multiplot
 set key left top
 set origin 0.0,0.5
 set size 1,0.5
 set yrange [0:] 
 set format x ""
 set tmargin 1
 plot 'plotdata_co2.dat' u ($2):($3) title ' CO2-RAMAN' w i lw 2

 set key left bottom
 set origin 0.0,0.0
 set size 1,0.587
 set yrange [0:] reverse
 set format x
 set xlabel "Wavenumber [cm-1]"
 set bmargin 3
 set ylabel "Intensity" offset 0,5
 plot 'plotdata_co2.dat' u ($2):($4) title 'CO2-IR' w i lw 2 lc 2 
 set nomultiplot
