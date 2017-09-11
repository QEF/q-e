# Copyright (C) 2001 Quantum ESPRESSO
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
 
fname=$1

# SCF
e1=`grep ! $fname | tail -1 | awk '{printf "%12.6f\n", $5}'`
n1=`grep 'convergence has' $fname | tail -1 | awk '{print $6}'`
f1=`grep "Total force" $fname | head -1 | awk '{printf "%8.4f\n", $4}'`
p1=`grep "P= " $fname | tail -1 | awk '{print $6}'`

# PH
diel=`grep -A 4 '  Dielectric constant in cartesian' $fname | grep -v '  Dielectric constant' | awk '{print $2; print $3; print $4 }'`
born=`grep "     E[x-z]  ( " $fname | awk '{print $3; print $4; print $5}'`
phfreq=`grep "     freq (.*THz" $fname | awk '{print $5; print $8}'`
dos=`grep "DOS =" $fname | awk '{print $3; print $8}'`
lambda=$(awk '
# reset counters
/Diagonalizing the dynamical matrix/{
  num_freq = 0
  ifreq = 0
}
# determine degenerate frequencies
# last frequency of a degenerate subspace has list number of degeneracies
/freq.*THz/{
  freq[num_freq] = $8
  degen[num_freq] = 1
  for (i = num_freq - 1; i >= 0; i--) {
    if (degen[i] >= 1) {
      if (freq[i] == freq[num_freq]) {
        degen[num_freq] = degen[i] + 1
        degen[i] = 0
      }
      break
    }
  }
  num_freq++
}
# store lambda in array and sort over degenerate subspace
/lambda/{
  # read lambda and gamma
  lambda[ifreq] = $3
  gamma[ifreq]  = $5
  # we use the maximum to deactivate found entries
  max_gamma = (max_gamma < gamma[ifreq] ? gamma[ifreq] : max_gamma)
  # sort over degenerate subspace
  for (i = 0; i < degen[ifreq]; i++) {
    # find minimum gamma
    j = 0
    min_pos   = ifreq - j
    min_gamma = gamma[min_pos]
    for (j = 1; j < degen[ifreq]; j++) {
      if (min_gamma > gamma[ifreq - j]) {
        min_pos   = ifreq - j
        min_gamma = gamma[min_pos]
      }
    }
    # print sorted array
    print lambda[min_pos]"\n"min_gamma
    # remove lambda from list
    gamma[min_pos] = max_gamma + 1.0
  }
  ifreq = (ifreq + 1) % num_freq
}' $fname)

# Q2R
qpt=`grep "q= " $fname | awk '{print $2; print $3; print $4}'` 

# LAMBDA
lambda2=`grep "lambda =" $fname | awk '{print $3; print $5; print $9 ;print $12; print $15}'`


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

if test "$dos" != ""; then
        echo dos
        for x in $dos; do echo $x; done
fi

if test "$lambda" != ""; then
        echo lambda
        for x in $lambda; do echo $x; done
fi

if test "$qpt" != ""; then
        echo qpt
        for x in $qpt; do echo $x; done
fi

if test "$lambda2" != ""; then
        echo lambda2
        for x in $lambda2; do echo $x; done
fi
