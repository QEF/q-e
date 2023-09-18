# Copyright (C) 2001 Quantum ESPRESSO
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
  
fname=$1

# zg_conf
conf=`grep "ZG_conf:" $fname | awk '{print $3; print $4; print $5}'`

if test "$conf" != ""; then
        echo conf
        for x in $conf; do echo $x; done
fi
