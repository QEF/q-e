#!/bin/bash
#SBATCH -J lif                 # Job name
#SBATCH -p development         # Queue (partition) name
#SBATCH -N 1                   # Total # of nodes
#SBATCH -n 56                  # Total # of mpi tasks
#SBATCH -t 00:30:00            # Run time (hh:mm:ss)
#SBATCH --mail-type=all        # Send email at begin and end of job
#SBATCH -A DMR21002            # Project/Allocation name (req'd if you have more than 1)
##SBATCH --mail-user=zhenbang.dai@austin.utexas.edu     # This's now in comment mode with one more # tag

# Load modules required for QE 
module load hdf5

bash espresso/link.sh
cd espresso/01-Density
bash run.sh
cd ../..

cd espresso/02-Wfn
mkdir -p lif.save
bash run.sh
cd ../..

cd espresso/03-Wfnq
mkdir -p lif.save
bash run.sh
cd ../..

cd espresso/04-Wfn_co
mkdir -p lif.save
bash run.sh
cd ../..

cd espresso/05-Wfn_fi
mkdir -p lif.save
bash run.sh
cd ../..

cd espresso/06-Wfnq_fi
mkdir -p lif.save
bash run.sh
cd ../..

# Load modules required for BerkeleyGW 
module load phdf5
module load mvapich2-x

cd 11-epsilon
bash run.sh
cd ..

cd 12-sigma
bash run.sh
cd ..

cd 13-kernel
bash run.sh
cd ..

cd 14-absorption
bash run.sh
cd ..


