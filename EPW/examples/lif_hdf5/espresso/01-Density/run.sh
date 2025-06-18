#!/bin/bash


. ../../setup.sh

$MPIRUN $PW $PWFLAGS -in scf.in &> scf.out

