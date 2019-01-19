#!/bin/bash
#SBATCH -A nstaff
#SBATCH -t 02:00:00
#SBATCH -C gpu
#SBATCH --exclusive
#SBATCH --gres=gpu:8

#load modules
module load gcc
module load pgi
#. ../buildscripts/env_gori.sh

#set up env
execdir=../buildscripts/install/6.3/knl/pgi/bin
mpidir=/global/homes/t/tkurth/src/openmpi_ucx/install_pgi/ompi

#some parameters
rankspernode=1
totalranks=$(( ${rankspernode} * ${SLURM_NNODES} ))

#openmp stuff
export OMP_NUM_THREADS=1 #$(( 40 / ${rankspernode} ))
#export OMP_PLACES=threads
#export OMP_PROC_BIND=true

#debug
export PGI_TERM='trace'
export USEGPU=no

#run
MPI_RUN="srun --mpi=pmi2 -N ${SLURM_NNODES} -n ${totalranks} -c $(( 80 / ${rankspernode} )) --cpu_bind=cores"
#cmd="${execdir}/pw.x -ndiag ${totalranks} -nbgrp ${totalranks} -in nano.in"
cmd="${execdir}/pw.x -in nano.in"
echo ${cmd}
#${MPI_RUN} ${execdir}/pw.x -in nano.in 
${cmd}
