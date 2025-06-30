#!/bin/bash
#SBATCH -J lif                   # Job name
#SBATCH -p development           # Queue (partition) name
#SBATCH -N 1                     # Total # of nodes
#SBATCH -n 56                   # Total # of mpi tasks
#SBATCH -t 00:30:00              # Run time (hh:mm:ss)
#SBATCH --mail-type=all          # Send email at begin and end of job
#SBATCH -A DMR21002              # Project/Allocation name (req'd if you have more than 1)

module load hdf5

. ../setup.sh

###Move the KS wave functions for gauge consistency###
rm -rf lif.save 
ln -sf ../espresso/05-Wfn_fi/lif.save lif.save 

###Move all the BSE eigenvectors to the current folder###
if [ -d "eigv" ]; then
  echo "BSE eigenvector exists"
else
  echo "Start copying BSE eigenvectors"

  mkdir eigv
  for i in {1..64}
  do
      j=$((i-1))
      cd eigv
      mkdir q_$j
      cd q_$j
      ln -sf ../../../finiteQ_grid/Q$i/14-absorption/eigenvectors.h5 .
      cd ../../
  done

  echo "Done copy BSE eigenvectors"
fi
##########################################################

$MPIRUN  $QEPATH/epw.x -npool 56 -inp epw1.in > epw1.out
$MPIRUN  $QEPATH/epw.x -npool 56 -inp epw2.in > epw2.out
$MPIRUN  $QEPATH/epw.x -npool 56 -inp epw3.in > epw3.out
