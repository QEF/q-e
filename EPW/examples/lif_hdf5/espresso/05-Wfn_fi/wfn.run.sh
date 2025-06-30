#!/bin/bash


. ../../setup.sh

cp -f ../01-Density/lif.save/data-file-schema.xml lif.save/data-file-schema.xml
$MPIRUN $PW $PWFLAGS -in wfn.in &> wfn.out

