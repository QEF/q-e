#!/bin/bash
#SBATCH -J lif                   # Job name
#SBATCH -p development          # Queue (partition) name
#SBATCH -N 1                     # Total # of nodes
#SBATCH -n 56                   # Total # of mpi tasks
#SBATCH -t 00:30:00              # Run time (hh:mm:ss)
#SBATCH --mail-type=all          # Send email at begin and end of job
#SBATCH -A DMR21002              # Project/Allocation name (req'd if you have more than 1)

module load hdf5

. ../setup.sh

export OMP_NUM_THREADS=1

cp ../pseudos/*.upf .

$MPIRUN  $QEPATH/pw.x -nk 14 -inp scf.in > scf.out
$MPIRUN  $QEPATH/ph.x -nk 14 -inp ph.in > ph.out

echo "Next, please run pp.py"
echo "For exciton polaron, please do not run another nscf calculations."
echo "Copy all the files in prefix.save in ../espresso/05-Wfn_fi to current_folder/prefix.save."
echo "The run.epw.sh script should help accomplish already."
