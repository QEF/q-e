set tit "with spin [nspin=2 OR (nspin=4 and domag)]" 
pl 'nspin2/3_hamiltonian/si.kcw_bands.dat' w lp pt 7 lc rgb 'black' tit "nspin=2" ,'nspin4_noSOC_MAG/3_hamiltonian/si.kcw_bands.dat' w l lc rgb 'red' tit 'nspin=4 noSOC MAG'
pause -1 

set tit "w/0 spin [nspin=1 OR (nspin=4 and NOT domag)]" 
pl 'nspin1/3_hamiltonian/si.kcw_bands.dat' w lp pt 7 lc rgb 'black' tit "nspin=2" ,'nspin4_noSOC_noMAG/3_hamiltonian/si.kcw_bands.dat' w l lc rgb 'red' tit 'nspin=4 noSOC MAG'
pause -1 
