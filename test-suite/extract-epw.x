# Copyright (C) 2001 Quantum ESPRESSO
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
#
# SP: This can lead to issue if you reach the OS pipe buffer
#     You can increase the buffer in /proc/sys/fs/pipe-max-size
  
fname=$1
args=$(echo $fname | awk -F= '{print $NF}')

##echo $fname > /home/sponce/program/espresso/test-suite/tmp.txt 
###    echo $args >> /home/sponce/program/espresso/test-suite/tmp.txt

### if [[ "$args" == "1" ]]
##  then
# SCF
e1=`grep ! $fname | tail -1 | awk '{printf "%12.6f\n", $5}'`
n1=`grep 'convergence has' $fname | tail -1 | awk '{print $6}'`
f1=`grep "Total force" $fname | head -1 | awk '{printf "%8.4f\n", $4}'`
p1=`grep "P= " $fname | tail -1 | awk '{print $6}'`
### fi


# NSCF
ef1=`grep "the Fermi energy is" $fname | awk '{print $5}'`
eh1=`grep "highest occupied" $fname | awk '{print $7}'`
el1=`grep "highest occupied" $fname | awk '{print $8}'`
tf1=`grep " P = " $fname | head -1 | awk '{printf "%7.5f", $3}'`

# PH
diel=`grep -A 4 '  Dielectric constant in cartesian' $fname | grep -v '  Dielectric constant' | awk '{print $2; print $3; print $4 }'`
born=`grep "     E[x-z]  ( " $fname | awk '{print $3; print $4; print $5}'`
phfreq=`grep "     freq (.*THz" $fname | awk '{print $5; print $8}'`

# Q2R
qdir=`grep " q= " $fname | awk '{print $2; print $3; print $4}'`

# EPW
q1=`grep "   q(" $fname | awk '{print $6; print $7; print $8}'`
dos1=`grep "DOS =" $fname | awk '{print $3}'`
e2=`grep " E(" $fname | awk '{print $4}'`
rsig=`grep "Re\[Sigma\]=" $fname | awk '{print $7}'` 
isig=`grep "Im\[Sigma\]=" $fname | awk '{print $10}'` 
z1=`grep " Z=" $fname | awk '{print $13}'`
lam=`grep "lam= " $fname | awk '{print $15}'`
lambda=`grep "     lambda___(" $fname | awk '{print $4}'`
lambda_tr=`grep "  lambda_tr(" $fname | awk '{print $4}'`
gamma=`grep " gamma___=" $fname | awk '{print $6}'`
omega=`grep " omega=" $fname | awk '{print $9}'`
lam_tot=`grep " lambda :" $fname | awk '{print $3}'`
lam_tr=`grep " lambda_tr :" $fname | awk '{print $3}'`
logavg=`grep " logavg =" $fname | awk '{print $3}'`
l_a2F=`grep "l_a2F =" $fname | awk '{print $6}'`
efm=`grep "at Ef=" $fname | awk '{print $8}'`
lam_max=`grep "lambda_max =" $fname | awk '{print $3}'`
lam_kmax=`grep "lambda_k_max =" $fname | awk '{print $6}'`
elph=`grep "Electron-phonon coupling strength =" $fname | awk '{print $5}'`
allDyn=`grep "Estimated Allen-Dynes Tc =" $fname | awk '{print $5}'`
bcsgap=`grep "Estimated BCS superconducting gap =" $fname | awk '{print $6}'`
pi=`grep "Re[Pi]=" $fname | awk '{print $4; print $7; print $10}'`
mobvb=`grep "Mobility VB Fermi level" $fname | awk '{print $5}'`
mobcb=`grep "Mobility CB Fermi level" $fname | awk '{print $5}'`
density=`grep " x-axis" $fname | awk '{print $1; print $2; print $3}'`
mobx=`grep " x-axis" $fname | awk '{print $4}'`
mobav=`grep "   avg" $fname | awk '{print $1}'`
mobxZ=`grep " x-axis [Z]" $fname | awk '{print $1; print $2; print $3; print $4}'`
indabs=`grep "  (cm-1)" $fname | awk '{print $1; print $2; print $3; print $4}'` 

if test "$efm" != ""; then
        echo efm
        for x in $efm; do echo $x; done
fi

if test "$lam_max" != ""; then
        echo lam_max
        echo $lam_max
fi

if test "$lam_kmax" != ""; then
        echo lam_kmax
        echo $lam_kmax
fi


if test "$elph" != ""; then
        echo elph
        echo $elph
fi

if test "$allDyn" != ""; then
        echo allDyn
        echo $allDyn
fi

if test "$bcsgap" != ""; then
        echo bcsgap
        echo $bcsgap
fi

if test "$mobvb" != ""; then
        echo mobvb
        for x in $mobvb; do echo $x; done
fi

if test "$mobcb" != ""; then
        echo mobcb
        for x in $mobcb; do echo $x; done
fi

if test "$mobx" != ""; then
        echo mobx
        for x in $mobx; do echo $x; done
fi

if test "$mobav" != ""; then
        echo mobav
        for x in $mobav; do echo $x; done
fi

if test "$density" != ""; then
        echo density
        for x in $density; do echo $x; done
fi

if test "$mobxZ" != ""; then
        echo mobxZ
        for x in $mobxZ; do echo $x; done
fi

if test "$indabs" != ""; then
        echo indabs
        for x in $indabs; do echo $x; done
fi

if test "$qdir" != ""; then
        echo qdir
        for x in $qdir; do echo $x; done
fi

if test "$q1" != ""; then
        echo q1
        for x in $q1; do echo $x; done
fi

if test "$dos1" != ""; then
        echo dos1
        for x in $dos1; do echo $x; done
fi

if test "$e2" != ""; then
        echo e2
        for x in $e2; do echo $x; done
fi

if test "$rsig" != ""; then
        echo rsig
        for x in $rsig; do echo $x; done
fi

if test "$isig" != ""; then
        echo isig
        for x in $isig; do echo $x; done
fi

if test "$z1" != ""; then
        echo z1
        for x in $z1; do echo $x; done
fi

if test "$lam" != ""; then
        echo lam
        for x in $lam; do echo $x; done
fi

if test "$lamda" != ""; then
        echo lamda
        for x in $lamda; do echo $x; done
fi

if test "$lambda_tr" != ""; then
        echo lambda_tr
        for x in $lambda_tr; do echo $x; done
fi

if test "$gamma" != ""; then
        echo gamma
        for x in $gamma; do echo $x; done
fi

if test "$omega" != ""; then
        echo omega
        for x in $omega; do echo $x; done
fi


if test "$lam_tot" != ""; then
        echo lam_tot
        echo $lam_tot
fi

if test "$lam_tr" != ""; then
        echo lam_tr
        echo $lam_tr
fi

if test "$logavg" != ""; then
        echo logavg
        echo $logavg
fi

if test "$l_a2F" != ""; then
        echo l_a2F
        echo $l_a2F
fi

if test "$e1" != ""; then
	echo e1
	echo $e1
fi

if test "$n1" != ""; then
	echo n1
	echo $n1
fi

if test "$f1" != ""; then
	echo f1
	echo $f1
fi

if test "$p1" != ""; then
	echo p1
	echo $p1
fi

if test "$ef1" != ""; then
	echo ef1
        for x in $ef1; do echo $x; done
fi

if test "$eh1" != ""; then
        echo eh1
        for x in $eh1; do echo $x; done
fi

if test "$el1" != ""; then
        echo el1
        for x in $el1; do echo $x; done
fi

if test "$tf1" != ""; then
        echo tf1
        for x in $tf1; do echo $x; done
fi

if test "$diel" != ""; then
        echo diel
        for x in $diel; do echo $x; done
fi

if test "$born" != ""; then
        echo born
        for x in $born; do echo $x; done
fi

if test "$phfreq" != ""; then
        echo phfreq
        for x in $phfreq; do echo $x; done
fi

if test "$pi" != ""; then
        echo pi
        for x in $pi; do echo $x; done
fi


