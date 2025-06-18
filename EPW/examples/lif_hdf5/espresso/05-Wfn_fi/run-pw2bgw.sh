#!/bin/bash


. ../../setup.sh

srun -n 1 $PW2BGW -in wfn.pp.in &> wfn.pp.out

