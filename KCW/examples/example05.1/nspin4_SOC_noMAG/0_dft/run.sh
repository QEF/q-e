#!/bin/bash

printf "  PWSCF   ... "
pw.x -in scf.pwi > scf.pwo
echo " done"
printf "  PWbands ... "
pw.x -in bands.pwi > bands.pwo
echo " done"
printf "  BANDS   ... "
bands.x -in bands.in > bands.out
echo " done"

