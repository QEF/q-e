home=`pwd`

cd example01/
./run_example
./compare.sh
echo "DONE"
rm -fr home/nicola/Scrivania/CODES/git/koopmans/quantum_espresso/qe_koopmans/tempdir
cd $home

echo "Running example02"
cd example02/
./run_example
./compare.sh
echo "DONE"
rm -fr home/nicola/Scrivania/CODES/git/koopmans/quantum_espresso/qe_koopmans/tempdir
cd $home 

cd example01/old/kpoints/nspin1/
echo "Running example01/old/kpoints/nspin1"
sh run.sh
./compare.sh
echo "DONE"
rm -fr out
cd $home

echo "Running example01/old/kpoints/nspin2/spin_up"
cd example01/old/kpoints/nspin2/spin_up
sh run.sh
./compare.sh
echo "DONE"
rm -fr out
cd $home

echo "Running example01/old/kpoints/nspin2/spin_dw"
cd example01/old/kpoints/nspin2/spin_dw
sh run.sh
./compare.sh
echo "DONE"
rm -fr out
cd $home

