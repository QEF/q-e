#!/bin/sh

# Checks are presently implemented only for pw.x and for the following cases:
# calculation='scf', 'relax', 'md'
# The following quantites are verified against reference output :
#    the converged total energy
#    the total force ( sqrt(\sum_i f_i^2)) if calculated;
#    the pressure P if calculated
# Presently: input data *.in, reference results *.res, output *.out

if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

ESPRESSO_ROOT=$HOME/espresso/
ESPRESSO_TMPDIR=./tmp/
ESPRESSO_PSEUDO=$ESPRESSO_ROOT/pseudo/

# no need to specify outdir and pseudo_dir in all *.in files
export ESPRESSO_TMPDIR ESPRESSO_PSEUDO

if test ! -d $ESPRESSO_TMPDIR ; then
   mkdir $ESPRESSO_TMPDIR
fi

# usage : ./check.j [input data files]
# With no arguments, checks all *.in files
if test $# = 0
then
    files=`/bin/ls *.in`
else
    files=$*
fi

for file in $files
do
  name=`basename $file .in`
  $ECHO "Checking $name...\c"
  ###
  $ESPRESSO_ROOT/bin/pw.x < $name.in > $name.out
  ###
  if test $? != 0; then
     $ECHO "FAILED with error condition!"
     $ECHO "Input: $name.in, Output: $name.out, Reference: $name.ref"
     $ECHO "Aborting"
     exit 1
  fi
  #
  if test -f $name.ref ; then
     # reference file exists
     # get reference total energy (cut to 6 significant digits)
     e0=`grep ! $name.ref | tail -1 | awk '{printf "%12.6f\n", $5}'`
     # get reference force (cut to 4 significant digits)
     f0=`grep "Total force = " $name.ref | tail -1 | awk '{printf "%8.4f\n", $4}'`
     # get reference pressure
     p0=`grep "P= " $name.ref | tail -1 | awk '{print $6}'`
     #
     e1=`grep ! $name.out | tail -1 | awk '{printf "%12.6f\n", $5}'`
     f1=`grep "Total force = " $name.out | tail -1 | awk '{printf "%8.4f\n", $4}'`
     p1=`grep "P= " $name.out | tail -1 | awk '{print $6}'`
     #
     if test "$e1" = "$e0"; then
        if test "$f1" = "$f0"; then
           if test "$p1" = "$p0"; then
              $ECHO  "passed"
           fi
        fi
     fi
     if test "$e1" != "$e0"; then
        $ECHO "discrepancy in total energy detected"
        $ECHO "Reference: $e0, You got: $e1"
     fi
     if test "$f1" != "$f0"; then
        $ECHO "discrepancy in force detected"
        $ECHO "Reference: $f0, You got: $f1"
     fi
     if test "$p1" != "$p0"; then
        $ECHO "discrepancy in pressure detected"
        $ECHO "Reference: $p0, You got: $p1"
     fi
  else
     # reference does not exist
     $ECHO  "not checked, reference file not available "
  fi
done
