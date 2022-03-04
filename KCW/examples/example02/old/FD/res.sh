rm -fr ene_vs_occ.dat
for charge in 0.0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1; do 
occ=`echo 1.0-$charge | bc -l`
ene=`grep ! h2o.scf_charge$charge.out | awk '{print $5}'`
echo $occ $ene >> ene_vs_occ.dat
done 

homo=`grep high h2o.scf_charge0.0.out | awk '{print $5}'`
e0=`grep ! h2o.scf_charge0.0.out | awk '{printf "%12.8f", $5*13.6057}'`
d2=`grep relaxed ../h2o.kc_screen.out | head -4 | tail -1 | awk '{print $6}'`
echo $homo
echo $e0

cat > plot.gnu << EOF
pl 'ene_vs_occ.dat' u 1:(\$2*13.6057)-($e0) w lp tit 'E(N)-E0', $homo*(x-1)+0.5*$d2*13.6057*(x-1)**2 w l tit '2nd order approx'
pause -1
pl 'ene_vs_occ.dat' u 1:(\$2*13.6057)-($e0 + $homo*(\$1-1)) w lp tit 'E-E0-dE/dN (x-1)', 0.5*$d2*13.6057*(x-1)**2 w l tit '2nd order approx'
pause -1
EOF

gnuplot plot.gnu



