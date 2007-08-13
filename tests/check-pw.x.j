#!/bin/sh

# Checks for pw.x are presently implemented only for the following
# calculations: 'scf', 'relax', 'md', 'nscf' (see below for the latter)
# The following quantites are verified against reference output :
#    the converged total energy
#    the module of the force ( sqrt(\sum_i f_i^2)) if calculated;
#    the pressure P if calculated
# Input data: *.in, reference results: *.res, output: *.out
# ./check-pw.x.j checks all *.in files
# ./check-pw.x.j "some file(s)" checks the specified files
# Example: 
# ./check-pw.x.j atom*.in lsda*
#
# For 'nscf' case, the data is in file $name.in2, where $name.in is the
# data for the scf calculation that must be executed before the nscf one.
# Output is written to $name.out2 and checked vs reference data $name.res2
# The following quantities are compared: Fermi energy, or HOMO and LUMO

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
# With an argument, checks files (ending with .in) matching the argument

if test $# = 0
then
    files=`/bin/ls *.in`
else
    files=`/bin/ls $*| grep "\.in$"`
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
     # get reference initial force (cut to 4 significant digits)
     f0=`grep "Total force = " $name.ref | head -1 | awk '{printf "%8.4f\n", $4}'`
     # get reference pressure
     p0=`grep "P= " $name.ref | tail -1 | awk '{print $6}'`
     #
     e1=`grep ! $name.out | tail -1 | awk '{printf "%12.6f\n", $5}'`
     f1=`grep "Total force = " $name.out | head -1 | awk '{printf "%8.4f\n", $4}'`
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
     $ECHO  "not checked, reference file not available "
  fi
  #
  # now check subsequent non-scf step if required
  # look for $name.in2
  n=2
  if test -f $name.in$n; then
     $ECHO "Checking $name, step $n ...\c"
     ###
     $ESPRESSO_ROOT/bin/pw.x < $name.in$n > $name.out$n
     ###
     if test $? != 0; then
        $ECHO "FAILED with error condition!"
        $ECHO "Input: $name.in$n, Output: $name.out$n, Reference: $name.ref$n"
        $ECHO "Aborting"
        exit 1
     fi
     #
     if test -f $name.ref$n ; then
        # reference file exists
        # get reference Fermi energy
        ef0=`grep Fermi $name.ref$n | awk '{print $5}'`
        # get reference HOMO and LUMO
        eh0=`grep "highest occupied" $name.ref$n | awk '{print $7}'`
        el0=`grep "highest occupied" $name.ref$n | awk '{print $8}'`
        #
        ef1=`grep Fermi $name.out$n | awk '{print $5}'`
        eh1=`grep "highest occupied" $name.out$n | awk '{print $7}'`
        el1=`grep "highest occupied" $name.out$n | awk '{print $8}'`
        #
        if test "$ef1" = "$ef0"; then
           if test "$eh1" = "$eh0"; then
              if test "$el1" = "$el0"; then
                 $ECHO  "passed"
              fi
           fi
        fi
        if test "$ef1" != "$ef0"; then
           $ECHO "discrepancy in Fermi energy detected"
           $ECHO "Reference: $ef0, You got: $ef1"
        fi
        if test "$eh1" != "$eh0"; then
           $ECHO "discrepancy in HOMO detected"
           $ECHO "Reference: $eh0, You got: $eh1"
        fi
        if test "$el1" != "$el0"; then
           $ECHO "discrepancy in LUMO detected"
           $ECHO "Reference: $el0, You got: $el1"
        fi
     else
        $ECHO  "not checked, reference file not available "
     fi
  fi
done
