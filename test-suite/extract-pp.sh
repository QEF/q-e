# Copyright (C) 2020 Quantum ESPRESSO Foundation
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
 
fname=$1

# SCF
e1=`grep ! $fname | tail -1 | awk '{printf "%12.6f\n", $5}'`

# PPACF
efock=`grep 'Fock energy' $fname | awk '{print $4}'`
exlda=`grep 'LDA Exchange' $fname | awk '{print $3}'`
eclda=`grep 'LDA Correlation' $fname | awk '{print $3}'`
exc=`grep 'Exchange + Correlation' $fname | awk '{print $4}'`
etcl=`grep 'T_c^LDA' $fname | awk '{print $2}'`
etnl=`grep 'T_c^nl' $fname | awk '{print $2}'`
ekc=`grep 'Kinetic-correlation Energy' $fname | awk '{print $3}'`
enl=`grep 'Non-local energy' $fname | awk '{print $4}'`

if test "$e1" != ""; then
        echo e1
        echo $e1
fi
if test "$efock" != ""; then
        echo efock
        echo $efock
fi
if test "$exlda" != ""; then
        echo exlda
        echo $exlda
fi
if test "$eclda" != ""; then
        echo eclda
        echo $eclda
fi
if test "$exc" != ""; then
        echo exc
        echo $exc
fi
if test "$etcl" != ""; then
        echo etcl
        echo $etcl
fi
if test "$etnl" != ""; then
        echo etnl
        echo $etnl
fi
if test "$ekc" != ""; then
        echo ekc
        echo $ekc
fi
if test "$enl" != ""; then
        echo enl
        echo $enl
fi
