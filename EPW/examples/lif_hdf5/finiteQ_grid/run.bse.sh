#!/bin/bash
#SBATCH -J lif             # Job name
#SBATCH -p development           # Queue (partition) name
#SBATCH -N 1                  # Total # of nodes
#SBATCH -n 56              # Total # of mpi tasks
#SBATCH -t 01:00:00           # Run time (hh:mm:ss)
#SBATCH --mail-type=all       # Send email at begin and end of job
#SBATCH -A DMR21002    # Project/Allocation name (req'd if you have more than 1)
##SBATCH --mail-user=zhenbang.dai@austin.utexas.edu     # This's now in comment mode with one more # tag

# Load all modules needed

module load phdf5
module load mvapich2-x

for i in {1..64}
do
    echo "Starting BSE at the ${i}-th(/64) Q point..."

    cd Q$i
    bash link.sh

    cd 13-kernel
    bash run.sh
    cd ..

    cd 14-absorption
    bash run.sh
    cd ../13-kernel
    rm bsemat.h5
    cd ../../
done

