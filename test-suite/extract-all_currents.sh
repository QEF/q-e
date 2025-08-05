#!/bin/bash
# Copyright (C) 2001 Quantum ESPRESSO
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

paste -d " "  <(awk 'BEGIN{h=1}$1=="total" && $3=="current:"{if(h==1) {print "jx jy jz ";h=0;} print $4,$5,$6}' $1 )  <(awk 'BEGIN{h=1} $1=="center" && $4=="velocity" && $7=="1:" {if(h==1) {print "jx1 jy1 jz1 ";h=0;}print $8,$9,$10}' $1 ) <(awk 'BEGIN{h=1}$1=="center" && $4=="velocity" && $7=="2:" {if(h==1) {print "jx2 jy2 jz2 ";h=0;} print $8,$9,$10}' $1 )  <(awk 'BEGIN{h=1}$1=="center" && $4=="velocity" && $7=="3:" {if(h==1) {print "jx3 jy3 jz3 ";h=0;}print $8,$9,$10}' $1 )    <(awk 'BEGIN{h=1}$1=="center" && $4=="velocity" && $7=="4:" {if(h==1) {print "jx4 jy4 jz4 ";h=0;}print $8,$9,$10}' $1 )  <(awk 'BEGIN{h=1}$1=="center" && $4=="velocity" && $7=="5:" {if(h==1) {print "jx5 jy5 jz5 ";h=0;}print $8,$9,$10}' $1 )
exit 0
