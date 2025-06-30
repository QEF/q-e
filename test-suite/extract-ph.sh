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
band=`sed -n "/bands (ev)/{n;n;p}" $fname | awk '{print $1; print $2; print $3; print $4; print $5 }' | head -$num_band`

# NSCF
ef1=`grep Fermi $fname | head -$max_iter | awk '{print $5}'`
eh1=`grep "highest occupied" $fname | tail -1 | awk '{print $5}'`
ehl1=`grep "highest occupied, lowest unoccupied" $fname | tail -1 | awk '{print $7; print $8}'`
tf1=`grep " P = " $fname | head -1 | awk '{printf "%7.5f", $3}'`

# PH
qpoint=`grep '  Calculation of q ' $fname | awk '{print $5; print $6; print $7 }'`
diel=`grep -A 4 '  Dielectric constant in cartesian' $fname | grep -v '  Dielectric constant' | awk '{print $2; print $3; print $4 }'`
born=`grep "     E[x-z]  ( " $fname | awk '{print $3; print $4; print $5}'`
# phfreq=`grep "     freq (.*THz" $fname | awk '{print $5; print $8}'`
# in the version below, phfreq only contains frequencies in cm^-1, not in THz
phfreq=`grep "     freq (.*THz" $fname | awk '{print $8}'`
dos=`grep "DOS =" $fname | awk '{print $3; print $8}'`
quad=`grep -A 2 'quadrupole.fmt' $fname | tail -1 | awk '{print $6}'`
epsil=`grep -A 6 'epsilon.fmt' $fname | tail -1 | awk '{print $5}'`
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
dnsscf_e=`grep -A 5 "DNS_SCF SYMMETRIZED IN ELECTRIC FIELD IN CARTESIAN COORDINATES" $fname  | tail -3`
dnsscf_ph=`grep -A 5 "DNS_SCF SYMMETRIZED IN CARTESIAN COORDINATES" $fname  | tail -3`

# Q2R
qpt=`grep "q= " $fname | awk '{print $2; print $3; print $4}'`

# MATDYN
born_diff=`grep "Norm of the difference between" $fname | awk '{print $NF}'`
# DYNMAT 
dynmat_freqs=`awk -e '/# mode/ {found++}; found && NF {print $0}' $fname | grep ^[[:space:]]*[1-9] | awk -e '{print $2}'`
# LAMBDA
lambda2=`grep "lambda =" $fname | awk '{print $3; print $5; print $9 ;print $12; print $15}'`

#DVSCF_Q2R
rlatt_cart=`grep "rlatt_cart" $fname | awk '{print $2; print $3; print $4}'`
rlatt_crys=`grep "rlatt_crys" $fname | awk '{print $2; print $3; print $4}'`
sum_w_pot=`grep "sum_w_pot" $fname | awk '{print $2; print $3; print $4}'`

#POSTAHC
postahc_selfen=`sed -n '/Begin postahc output/,/End postahc output/p' $fname \
    | head -n -3 | tail -n +6 | awk '{print $3; print $4; print $5 ; print $6 ; print $7}'`

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

if test "$band" != ""; then
        echo band
        for x in $band; do echo $x; done
fi

if test "$ef1" != ""; then
  echo ef1
  for x in $ef1; do echo $x; done
fi

if test "$eh1" != ""; then
        echo eh1
        for x in $eh1; do echo $x; done
fi

if test "$ehl1" != ""; then
        echo ehl1
        for x in $ehl1; do echo $x; done
fi

if test "$tf1" != ""; then
        echo tf1
        for x in $tf1; do echo $x; done
fi

if test "$qpoint" != ""; then
        echo qpoint
        for x in $qpoint; do echo $x; done
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

if test "$dynmat_freqs" != ""; then
        echo dynmat_freqs 
        for x in $dynmat_freqs; do echo $x; done
fi

if test "$dos" != ""; then
        echo dos
        for x in $dos; do echo $x; done
fi

if test "$quad" != ""; then
        echo quad 
        for x in $quad; do echo $quad; done
fi

if test "$epsil" != ""; then
        echo epsil 
        for x in $epsil; do echo $epsil; done
fi

if test "$lambda" != ""; then
        echo lambda
        for x in $lambda; do echo $x; done
fi

if test "$dnsscf_e" != ""; then
        echo dnsscf_e
        for x in $dnsscf_e; do echo $x; done
fi

if test "$dnsscf_ph" != ""; then
        echo dnsscf_ph
        for x in $dnsscf_ph; do echo $x; done
fi

if test "$qpt" != ""; then
        echo qpt
        for x in $qpt; do echo $x; done
fi

if test "$born_diff" != ""; then
        echo born_diff
        for x in $born_diff; do echo $x; done
fi

if test "$lambda2" != ""; then
        echo lambda2
        for x in $lambda2; do echo $x; done
fi

if test "$rlatt_cart" != ""; then
        echo rlatt_cart
        for x in $rlatt_cart; do echo $x; done
fi
if test "$rlatt_crys" != ""; then
        echo rlatt_crys
        for x in $rlatt_crys; do echo $x; done
fi
if test "$sum_w_pot" != ""; then
        echo sum_w_pot
        for x in $sum_w_pot; do echo $x; done
fi
if test "$postahc_selfen" != ""; then
        echo postahc_selfen
        for x in $postahc_selfen; do echo $x; done
fi
