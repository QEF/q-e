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

# PP
efock=`grep 'Fock energy' $fname | awk '{print $4}'`

if test "$e1" != ""; then
        echo e1
        echo $e1
fi
if test "$efock" != ""; then
        echo efock
        echo $efock
fi
