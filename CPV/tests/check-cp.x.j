#!/bin/sh

# Automated checks for cp.x - CC 2009
#
. ../../environment_variables
#
# You shouldn't need to modify anything below this line.
#
# Some specific quantities are checked against a reference output
# Checks are implemented for the following calculations: 
#   ' fill in here '
# (see below for the three latter)
#
# Input data: *.in, reference results: *.res, output: *.out
# ./check-cp.x.j checks all *.in files
# ./check-cp.x.j "some file(s)" checks the specified files
# Example:
#    ./check-cp.x.j h2o*.in lsda*
# If you want to save a copy in file "logfile":
#    ./check-cp.x.j h2o*.in lsda* | tee logfile
#
# The quantites that are verified are:
#    the last value of total energy, forces and stress

# taken from examples - not sure it is really needed
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi
#
ESPRESSO_ROOT=`cd ../.. ; pwd`
#ESPRESSO_TMPDIR=$ESPRESSO_ROOT/CPV/tmp/
ESPRESSO_TMPDIR=$TMP_DIR
ESPRESSO_PSEUDO=$ESPRESSO_ROOT/pseudo/

# no need to specify outdir and pseudo_dir in all *.in files
export ESPRESSO_TMPDIR ESPRESSO_PSEUDO

if test ! -d $ESPRESSO_TMPDIR
then
   mkdir $ESPRESSO_TMPDIR
fi

# this is the current directory, where the test is executed
TESTDIR=`pwd`

# With no arguments, checks all *.in files
# With an argument, checks files (ending with .in1, .in2, ecc...) matching the argument

if test $# = 0
then
    files=`/bin/ls *.in1`
else
    files=`/bin/ls $*| grep "\.in[1-9]"`
fi

########################################################################
# function to get pseudopotentials from the web if missing
########################################################################
get_pp () {
    ppfiles=`grep UPF $1 | awk '{print $3}'`
    for ppfile in $ppfiles
    do
        if ! test -f $ESPRESSO_PSEUDO/$ppfile ; then
            $ECHO "Downloading $ppfile to $ESPRESSO_PSEUDO...\c"
            $WGET $ESPRESSO_PSEUDO/$ppfile $NETWORK_PSEUDO/$ppfile 2> /dev/null
            if test $? != 0; then
                $ECHO "failed!"
                $ECHO "test $1 will not be executed"
                # status=1
            else
                $ECHO "success"
                # status=0
            fi
        fi
    done
}
########################################################################
# function to test NEB calculations - usage: check_neb "file prefix"
# obsolete - will be moved to NEB-specific tests
########################################################################
check_neb () {
  # get reference number of neb iterations
  n0=`grep 'neb: convergence' $1.ref | awk '{print $1}'`
  # get reference activation energy (truncated to 4 significant digits)
  e0=`grep 'activation energy' $1.ref | tail -1 | awk '{printf "%8.4f\n", $5}'`
  #
  n1=`grep 'neb: convergence' $1.out | awk '{print $1}'`
  e1=`grep 'activation energy' $1.out | tail -1 | awk '{printf "%8.4f\n", $5}'`
  if test "$e1" = "$e0"
  then
    if test "$n1" = "$n0"
    then
      $ECHO  "passed"
    fi
  fi
  if test "$e1" != "$e0"
  then
    $ECHO "discrepancy in activation energy detected"
    $ECHO "Reference: $e0, You got: $e1"
  fi
  if test "$n1" != "$n0"
  then
    $ECHO "discrepancy in number of neb iterations detected"
    $ECHO "Reference: $n0, You got: $n1"
  fi
}
########################################################################
# function to test scf calculations - usage: check_scf "file prefix"
########################################################################
check_cp () {
  fname=$1.ref$2
  # get reference total energy (cut to 6 significant digits)
  e0=`grep "total energy =" $fname | tail -1 | awk '{printf "%18.6f\n", $4}'`
  # get reference number for stress matrix
  s0=`grep -A 3 "Total stress" $fname | tail -3 | tr '\n' ' ' | awk '{ printf "%-18.8f", $1+$2+$3+$4+$5+$6+$7+$8+$9 }'`
  # get reference eigenvalues
  v0u=`grep -A 2 "Eigenvalues (eV).*spin.*1" $fname | tail -1 | awk '{ for(i=1;i<=NF;i++) { v=v+$i; } print v }'` 
  v0d=`grep -A 2 "Eigenvalues (eV).*spin.*2" $fname | tail -1 | awk '{ for(i=1;i<=NF;i++) { v=v+$i; } print v }'` 
  # get average temperature over the step of the current execution
  t0=`grep -A 6 "Averaged Physical Quantities"  $fname | tail -1 | awk '{ print $4 }'`
  # note that only the final energy, pressure, number of iterations, 
  # and only the initial force are tested - hopefully this should 
  # cover the various MD and optimization cases as well as simple scf
  #
  fname=$1.out$2
  e1=`grep "total energy =" $fname | tail -1 | awk '{printf "%18.6f\n", $4}'`
  s1=`grep -A 3 "Total stress" $fname | tail -3 | tr '\n' ' ' | awk '{ printf "%-18.8f", $1+$2+$3+$4+$5+$6+$7+$8+$9 }'`
  v1u=`grep -A 2 "Eigenvalues (eV).*spin.*1" $fname | tail -1 | awk '{ for(i=1;i<=NF;i++) { v=v+$i; } print v }'` 
  v1d=`grep -A 2 "Eigenvalues (eV).*spin.*2" $fname | tail -1 | awk '{ for(i=1;i<=NF;i++) { v=v+$i; } print v }'` 
  t1=`grep -A 6 "Averaged Physical Quantities"  $fname | tail -1 | awk '{ print $4 }'`
  #
  #echo $e1
  #echo $s1
  #echo $v1
  #echo $t1
  #
  if test "$e1" = "$e0"
  then
    if test "$s1" = "$s0"
    then
      if test "$v1u" = "$v0u"
      then
        if test "$v1u" = "$v0u"
        then
          if test "$t1" = "$t0"
          then
            $ECHO  " $2 passed"
          fi
        fi
      fi
    fi
  fi
  if test "$e1" != "$e0"
  then
    $ECHO "discrepancy in total energy detected"
    $ECHO "Reference: $e0, You got: $e1"
  fi
  if test "$s1" != "$s0"
  then
    $ECHO "discrepancy in stress detected"
    $ECHO "Reference: $s0, You got: $s1"
  fi
  if test "$v1u" != "$v0u"
  then
    $ECHO "discrepancy in eigenvalues detected"
    $ECHO "Reference: $v0u, You got: $v1u"
  fi
  if test "$v1d" != "$v0d"
  then
    $ECHO "discrepancy in eigenvalues detected"
    $ECHO "Reference: $v0d, You got: $v1d"
  fi
  if test "$t1" != "$t0"
  then
    $ECHO "discrepancy in average temperature"
    $ECHO "Reference: $t0, You got: $t1"
  fi
}

########################################################################
# function to get wall times - usage: get_times "file prefix"
########################################################################
get_times () {
  # convert from "1h23m45.6s" to seconds
  # the following line prevents cases such as "2m 7.5s"
  grep 'WALL$' $1.ref$2 | sed 's/m /m0/' > $1.tmp
  # in order to get cpu instead of wall time, replace $3 to $5
  tref=`awk '{ str = $5; h = m = s = 0;
                  if (split(str, x, "h") == 2) { h = x[1]; str = x[2]; }
                  if (split(str, x, "m") == 2) { m = x[1]; str = x[2]; }
                  if (split(str, x, "s") == 2) { s = x[1]; str = x[2]; }
                  t += h * 3600 + m * 60 + s; }
                END { printf("%.2f\n", t); }' \
               $1.tmp`
  # as above for file *.out
  grep 'WALL$' $1.out$2 | sed 's/m /m0/' > $1.tmp 
  tout=`awk '{ str = $5; h = m = s = 0;
                  if (split(str, x, "h") == 2) { h = x[1]; str = x[2]; }
                  if (split(str, x, "m") == 2) { m = x[1]; str = x[2]; }
                  if (split(str, x, "s") == 2) { s = x[1]; str = x[2]; }
                  t += h * 3600 + m * 60 + s; }
                END { printf("%.2f\n", t); }' \
               $1.tmp`
  /bin/rm $1.tmp
  # accumulate data
  totref=`echo $totref $tref | awk '{print $1+$2}'`
  totout=`echo $totout $tout | awk '{print $1+$2}'`
}

for file in $files
do
  name=`basename $file .in1`
  $ECHO "Checking $name...\c"
  ###
  # run the code in the scratch directory
  #
  cd $ESPRESSO_TMPDIR
  #
  steps=""
  #
  for i in 1 2 3 4 5 6 7 8 9
  do
    if test -f $TESTDIR/$name.in$i ; then
      get_pp $TESTDIR/$name.in$i 
      $ECHO ".$i.\c"
      steps=`echo $steps $i`
      $PARA_PREFIX $ESPRESSO_ROOT/bin/cp.x $PARA_POSTFIX \
                               -i $TESTDIR/$name.in$i > $TESTDIR/$name.out$i
      if test $? != 0; then
        $ECHO "FAILED with error condition!"
        $ECHO "Input: $name.in$i, Output: $name.out$i, Reference: $name.ref$i"
        $ECHO "Aborting"
        exit 1
      fi
    fi
  done 
  #
  cd $TESTDIR
  #
  echo
  #
  for i in $steps
  do
    if test -f $name.ref$i ; then
      # reference file exists
      if grep 'neb: convergence achieved' $name.ref$i > /dev/null; then
         #
         # Specific test for NEB
         #
 	check_neb $name
         #
      else
         #
         # Test for scf/relax/md/vc-relax
         #
  	check_cp $name $i
         #echo check
         #
      fi
      #
      # extract wall time statistics
      #
      get_times $name $i
      #
    else
      $ECHO  "not checked, reference file not available "
    fi
  done
  #

done

$ECHO  "Total wall time (s) spent in this run: " $totout
$ECHO  "Reference                            : " $totref

