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

scf=$(echo $fname | awk '/pw/{print 1}' )
turbodav=$(echo $fname | awk '/-dav.in/{print 1}' )
turbolancz=$(echo $fname | awk '/-lan.in/{print 1}' )
eels=$(echo $fname | awk '/eels.in/{print 1}' )
magnons=$(echo $fname | awk '/magnons.in/{print 1}' )
turbospec=$(echo $fname | awk '/pp.in/{print 1}' )
turbospeceels=$(echo $fname | awk '/pp_eels.in/{print 1}' )
turbospecmagnons=$(echo $fname | awk '/pp_magnons.in/{print 1}' )

# SCF
if [ "$scf" = "1" ]; then
        e1=`grep ! $fname | tail -1 | awk '{printf "%12.6f\n", $5}'`
        n1=`grep 'convergence has' $fname | tail -1 | awk '{print $6}'`
        f1=`grep "Total force" $fname | head -1 | awk '{printf "%8.4f\n", $4}'`
        p1=`grep "P= " $fname | tail -1 | awk '{print $6}'`
fi

# turbo_davidson.x
if [ "$turbodav" = "1" ]; then
        e3=`grep "The transition energy is" $fname | awk '{print $(NF-3)}' | sed 's/(//g' `
        n3=` grep      "Davidson diagonalization has finished" $fname |              awk '{print $(NF-1)}' `
        nb1=`grep -A 1 "Davidson diagonalization has finished" $fname | tail -n 1  | awk '{print $NF}' `
        nb2=`grep -A 2 "Davidson diagonalization has finished" $fname | tail -n 1  | awk '{print $NF}' `
        chi=`grep --no-group-separator  -A  2  "Chi_i_i"       $fname | grep -v Chi_i_i   | tr '\n' ' ' `
fi

# turbo_lanczos.x
if [ "$turbolancz" = "1" ]; then
        alpha=`grep "alpha(" $fname | awk '{print $2}'`
        beta=`grep "beta( " $fname | awk '{print $3}'`
        gamma=`grep "gamma(" $fname | awk '{print $2}'`
fi

#turbo_eels.x
if [ "$eels" = "1" ]; then
        nblanczos=`grep "Number of Lanczos iterations" $fname | awk '{print $6}'`
        epsiloninv=`grep "epsilon^-1" $fname | awk '{print $3; print $4}'`
        chirr=`grep "chirr" $fname | awk '{print $3; print $4}'`
fi

# turbo_magnons.x
if [ "$magnons" = "1" ]; then
        nblanczos=`grep "Number of Lanczos iterations" $fname | awk '{print $6}'`
fi

# turbo_spectrum.x
if [ "$turbospec" = "1" ]; then
        rechi=`grep "chi_1_1=" $fname | awk '{print $3}'`
        imchi=`grep "chi_1_1=" $fname | awk '{print $4}'`
fi

if [ "$turbospeceels" = "1" ]; then
        freq=`tail -n +2 $fname | awk '{print $1}'`
        reepsm1=`tail -n +2 $fname | awk '{print $2}'`
        imepsm1=`tail -n +2 $fname | awk '{print $3}'`
        reeps=`tail -n +2 $fname | awk '{print $4}'`
        imeps=`tail -n +2 $fname | awk '{print $5}'`
fi

if [ "$turbospecmagnons" = "1" ]; then
       imchi=`grep "chi_2_2=" $fname | awk '{print $4}'`
fi

if test "$nblanczos" != ""; then
        echo nblanczos
        for x in $nblanczos; do echo $x; done
fi

if test "$epsiloninv" != ""; then
        echo epsiloninv
        for x in $epsiloninv; do echo $x; done
fi

if test "$chirr" != ""; then
        echo chirr
        for x in $chirr; do echo $x; done
fi

if test "$rechi" != ""; then
        echo rechi
        for x in $rechi; do echo $x; done
fi

if test "$imchi" != ""; then
        echo imchi
        for x in $imchi; do echo $x; done
fi

if test "$freq" != ""; then
        echo freq
        for x in $freq; do echo $x; done
fi

if test "$reepsm1" != ""; then
        echo reepsm1
        for x in $reepsm1; do echo $x; done
fi

if test "$imepsm1" != ""; then
        echo imepsm1
        for x in $imepsm1; do echo $x; done
fi

if test "$reeps" != ""; then
        echo reeps
        for x in $reeps; do echo $x; done
fi

if test "$imeps" != ""; then
        echo imeps
        for x in $imeps; do echo $x; done
fi

if test "$average" != ""; then
        echo average
        for x in $average; do echo $x; done
fi

if test "$averageosc" != ""; then
        echo averageosc
        for x in $averageosc; do echo $x; done
fi

if test "$alpha" != ""; then
        echo alpha
        for x in $alpha; do echo $x; done
fi

if test "$beta" != ""; then
        echo beta
        for x in $beta; do echo $x; done
fi

if test "$gamma" != ""; then
        echo gamma
        for x in $gamma; do echo $x; done
fi

if test "$e1" != ""; then
	echo e1
	echo $e1
fi

if test "$e3" != ""; then
	echo e3
	for x in $e3; do echo $x; done 
fi

if test "$n3" != ""; then
	echo n3
	echo $n3
fi

if test "$nb1" != ""; then
	echo nb1
	echo $nb1
fi

if test "$nb2" != ""; then
	echo nb2
	echo $nb2
fi

if test "$chi" != ""; then
	echo chi
	for x in $chi; do echo $x; done 
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
