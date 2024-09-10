#!/bin/bash

printf "  DFT scf       ... "
pw.x -in scf.pwi > scf.pwo
echo "  done"
printf "  DFT nscf      ... "
pw.x -in nscf.pwi > nscf.pwo
echo "  done"
printf "  Wann occ      ... "
cd occ
wannier90.x -pp wann
pw2wannier90.x -in pw2wan.p2wi > pw2wan.p2wo
wannier90.x wann
echo "  done"
cd ../
printf "  Wann emp      ... "
cd emp
wannier90.x -pp wann
pw2wannier90.x -in pw2wan.p2wi > pw2wan.p2wo
wannier90.x wann
echo "  done"
cd ../

./link_wann.sh

printf "  KCW interface ... "
kcw.x -in kc.w2ki > kc.w2ko
echo "  done"
