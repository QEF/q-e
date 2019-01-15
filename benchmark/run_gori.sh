#!/bin/bash
#SBATCH -A nstaff
#SBATCH -t 02:00:00
#SBATCH -C gpu
#SBATCH --exclusive
#SBATCH --gres=gpu:8

. ../buildscripts/env_gori.sh

#set up env
execdir=../buildscripts/install/6.3/v100/bin
mpidir=/global/homes/t/tkurth/src/openmpi_ucx/install_pgi/ompi

#some parameters
rankspernode=1
totalranks=$(( ${rankspernode} * ${SLURM_NNODES} ))

#openmp stuff
export OMP_NUM_THREADS=$(( 40 / ${rankspernode} ))
export OMP_PLACES=threads
#export OMP_PROC_BIND=true
export MKL_FAST_MEMORY_LIMIT=0

#run
MPI_RUN="srun --mpi=pmi2 -N ${SLURM_NNODES} -n ${totalranks} -c $(( 80 / ${rankspernode} )) --cpu_bind=cores"
${MPI_RUN} ${execdir}/pw.x -ndiag ${totalranks} -nbgrp ${totalranks} -in small.in 
