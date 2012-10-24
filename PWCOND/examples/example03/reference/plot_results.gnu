##
# Script to visualize the results of the DFT+U PWcond example
##

set style line 11 lt 1 lc rgbcolor 'blue' lw 2 pt 6
set style line 12 lt 2 lc rgbcolor 'cyan' lw 1 pt 6
set style line 21 lt 1 lc rgbcolor 'red' lw 2 pt 2
set style line 22 lt 2 lc rgbcolor 'magenta' lw 1 pt 2
set style arrow 1 nohead lt 3 lc rgbcolor 'dark-green' lw 1.5
set style arrow 2 nohead lt 1 lc rgbcolor 'black' lw 1

#1. compare CBS of Au chain within LDA and LDA+U
set key center right
set xlabel 'Re(k_z)'
set label 'Im(k_z)=0' at 0.02, -2.5 left
set label 'Im(k_z)' at first -0.25, screen 0.03 center
set label 'Re(k_z)=0' at -0.02, -2.5 right
set label 'Re(k_z)+Im(k_z)' at 0.75, screen 0.03 center
set label 'Re(k_z)=0.5' at 0.52, -2.5 left
set ylabel 'E - E_F  (eV)'
set arrow from graph 0,first 0 to graph 1, first 0 as 1
set arrow from 0,graph 0 to 0,graph 1 as 2
set arrow from 0.5,graph 0 to 0.5,graph 1 as 2
set xrange [-0.5:1.0]
plot '< awk "{if(\$1>0.0) print}" bands.Auwire.re'  w p ls 11 title 'U=0',\
     'bands.Auwire.im'  w p ls 11 notitle,\
     '< awk "{if(\$1>0.0) print}" bandsU.Auwire.re' w p ls 21 title 'U=3',\
     'bandsU.Auwire.im' w p ls 21 notitle

unset arrow
unset label
pause -1  "Hit return to continue"

## extract the number of channels
! echo "# channels" > nch.tmp
! grep Nchannels COatAuwire.cond.out  | cut -d\= -f2 >> nch.tmp
! echo "# channels" > nchU.tmp
! grep Nchannels COatAuwireU.cond.out | cut -d\= -f2 >> nchU.tmp

#2. compare the ballistic transmission for CO@Au chain
set xlabel 'E - E_F  (eV)'
set ylabel 'Transmittance'
set arrow from 0,graph 0 to 0,graph 1 as 1
set xrange [-1.0:1.0]
plot 'trans.AuwireCO'  u 1:(0.5*$2) w lp ls 11 title 'T(U=0)', \
     '< paste trans.AuwireCO nch.tmp'  u 1:3 w lp ls 12 title 'N(U=0)', \
     'transU.AuwireCO' u 1:(0.5*$2) w lp ls 21 title 'T(U=3)',\
     '< paste transU.AuwireCO nchU.tmp'  u 1:3 w lp ls 22 title 'N(U=3)'


