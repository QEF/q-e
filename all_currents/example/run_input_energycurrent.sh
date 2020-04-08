#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=182000
#SBATCH --time=00:05:00
#SBATCH --account=Sis20_baroni
#SBATCH --partition=skl_usr_dbg
#SBATCH --job-name=example

SCRIPTS_DIR="./scripts"
QEBIN_DIR="/marconi_work/Sis19_baroni/qe-thermalcurrents2/bin"
NPROCS=8

module purge
module load autoload
module load intel/pe-xe-2018--binary
module load intelmpi/2018--binary
module load mkl/2018--binary
module load scalapack/2.0.2--intelmpi--2018--binary
module load python/2.7.12
module load numpy/1.14.0--python--2.7.12

#######################################
# run ONLY all_currents
mpirun -np ${NPROCS} ${QEBIN_DIR}/all_currents.x -in input_energycurrent > input_energycurrent_out


# compute and write total current
#python ${SCRIPTS_DIR}/read_hartree_currents.py current_hz >> current_hz

echo "Current computed:  current_hz2"
cat current_hz2

