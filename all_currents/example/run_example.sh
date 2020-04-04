#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=182000
#SBATCH --time=00:30:00
#SBATCH --account=Sis19_baroni
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

######################################
# generate input files
${SCRIPTS_DIR}/create_pw_input2 < tabella
${SCRIPTS_DIR}/create_hz_input < tabella

######################################
# run PW 1 & 2
mpirun -np ${NPROCS} ${QEBIN_DIR}/pw.x -in input_pw_1 > out_pw_1
mpirun -np ${NPROCS} ${QEBIN_DIR}/pw.x -in input_pw_2 > out_pw_2

#######################################
# run HARTREE-init
mpirun -np ${NPROCS} ${QEBIN_DIR}/all_currents.x -in input_hz_init > out_hz_init

# run HARTREE-run
mpirun -np ${NPROCS} ${QEBIN_DIR}/all_currents.x -in input_hz_run > out_hz_run

# compute and write total current
python ${SCRIPTS_DIR}/read_hartree_currents.py current_hz >> current_hz

echo "Current computed:  current_hz"
cat current_hz

