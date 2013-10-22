#!/bin/sh

# Automated checks for coupling to QE codes 
#
. ../../environment_variables
#
# You shouldn't need to modify anything below this line.
#
#
# this takes some existing test inputs and primarily checks
# that the wrappers for fortran and c work correctly with
# different sets of processor counts and a subcommunicator. 
#
# taken from examples - not sure it is really needed
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

ESPRESSO_ROOT=`cd ../../ ; pwd`
ESPRESSO_TMPDIR=$ESPRESSO_ROOT/tmp/
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
# With an argument, checks files (ending with .in) matching the argument

if test $# = 0
then
    pwfiles=`/bin/ls *.pw.in`
    cpfiles=`/bin/ls *.cp.in`
else
    pwfiles=`/bin/ls $*| grep "\.pw.in$"`
    cpfiles=`/bin/ls $*| grep "\.cp.in$"`
fi

########################################################################
# function to get pseudopotentials from the web if missing
########################################################################
get_pp () {
    ppfiles=`grep UPF $1.in | awk '{print $3}'`
    for ppfile in $ppfiles
    do
	if ! test -f $ESPRESSO_PSEUDO/$ppfile ; then
	    $ECHO "Downloading $ppfile to $ESPRESSO_PSEUDO...\c"
	    $WGET  $ESPRESSO_PSEUDO/$ppfile $NETWORK_PSEUDO/$ppfile 2> /dev/null
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
# function to test scf calculations - usage: check_scf "file prefix"
########################################################################
check_scf () {
  # get reference total energy (cut to 6 significant digits)
  e0=`grep ! $1.ref | tail -1 | awk '{printf "%12.6f\n", $5}'`
  # get reference number of scf iterations
  n0=`grep 'convergence has' $1.ref | tail -1 | awk '{print $6}'`
  # get reference initial force (cut to 4 significant digits)
  f0=`grep "Total force" $1.ref | head -1 | awk '{printf "%8.4f\n", $4}'`
  # get reference pressure
  p0=`grep "P= " $1.ref | tail -1 | awk '{print $6}'`
  #
  # note that only the final energy, pressure, number of iterations, 
  # and only the initial force are tested - hopefully this should 
  # cover the various MD and optimization cases as well as simple scf
  #
  e1=`grep ! $1.out | tail -1 | awk '{printf "%12.6f\n", $5}'`
  n1=`grep 'convergence has' $1.out | tail -1 | awk '{print $6}'`
  f1=`grep "Total force" $1.out | head -1 | awk '{printf "%8.4f\n", $4}'`
  p1=`grep "P= " $1.out | tail -1 | awk '{print $6}'`
  #
  if test "$e1" = "$e0"
  then
    if test "$n1" = "$n0"
    then
      if test "$f1" = "$f0"
      then
        if test "$p1" = "$p0"
        then
          $ECHO  "passed"
        fi
      fi
    fi
  fi
  if test "$e1" != "$e0"
  then
    $ECHO "discrepancy in total energy detected"
    $ECHO "Reference: $e0, You got: $e1"
  fi
  if test "$n1" != "$n0"
  then
    $ECHO "discrepancy in number of scf iterations detected"
    $ECHO "Reference: $n0, You got: $n1"
  fi
    if test "$f1" != "$f0"
    then
    $ECHO "discrepancy in force detected"
    $ECHO "Reference: $f0, You got: $f1"
  fi
  if test "$p1" != "$p0"
  then
    $ECHO "discrepancy in pressure detected"
    $ECHO "Reference: $p0, You got: $p1"
  fi
}

########################################################################
# function to test nscf calculations - usage: check_nscf "file prefix" "number"
########################################################################
check_nscf () {
  # get reference Fermi energy
  ef0=`grep Fermi $1.ref$2 | awk '{print $5}'`
  # get reference HOMO and LUMO
  eh0=`grep "highest occupied" $1.ref$2 | awk '{print $7}'`
  el0=`grep "highest occupied" $1.ref$2 | awk '{print $8}'`
  # get total polarization (for Berry's phase calculation)
  tf0=`grep " P = " $1.ref$2 | head -1 | awk '{printf "%7.5f", $3}'`
  #
  ef1=`grep Fermi $name.out$n | awk '{print $5}'`
  eh1=`grep "highest occupied" $1.out$2 | awk '{print $7}'`
  el1=`grep "highest occupied" $1.out$2 | awk '{print $8}'`
  tf1=`grep " P = " $1.out$2 | head -1 | awk '{printf "%7.5f", $3}'`
  #
  if test "$ef1" = "$ef0"
  then
    if test "$eh1" = "$eh0"
    then
      if test "$el1" = "$el0"
      then
        if test "$tf1" = "$tf0"
        then
          $ECHO  "passed"
        fi
      fi
    fi
  fi
  if test "$ef1" != "$ef0"
  then
    $ECHO "discrepancy in Fermi energy detected"
    $ECHO "Reference: $ef0, You got: $ef1"
  fi
  if test "$eh1" != "$eh0"
  then
    $ECHO "discrepancy in HOMO detected"
    $ECHO "Reference: $eh0, You got: $eh1"
  fi
  if test "$el1" != "$el0"
  then
    $ECHO "discrepancy in LUMO detected"
    $ECHO "Reference: $el0, You got: $el1"
  fi
  if test "$tf1" != "$tf0"
  then
    $ECHO "discrepancy in polarization detected"
    $ECHO "Reference: $tf0, You got: $tf1"
  fi
}
########################################################################
# function to test cp calculations - usage: check_cp "file prefix"
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
  grep 'WALL$' $1.ref | sed 's/m /m0/' > $1.tmp 
  # in order to get cpu instead of wall time, replace $3 to $5
  tref=`awk '{ str = $5; h = m = s = 0;
                  if (split(str, x, "h") == 2) { h = x[1]; str = x[2]; }
                  if (split(str, x, "m") == 2) { m = x[1]; str = x[2]; }
                  if (split(str, x, "s") == 2) { s = x[1]; str = x[2]; }
                  t += h * 3600 + m * 60 + s; }
                END { printf("%.2f\n", t); }' \
               $1.tmp`
  # as above for file *.out
  grep 'WALL$' $1.out | sed 's/m /m0/' > $1.tmp 
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

########################################################################
# Perform here required checks on pw.x based tests
########################################################################
for file in $pwfiles
do
  name=`basename $file .in`
  get_pp $name
  $ECHO "Checking $name..."
  ###
  # run the code in the scratch directory
  #
  for res in 0 2 4 16
  do
    $ECHO "f2pw.x with $res reserved MPI ranks...\c"
    cd $ESPRESSO_TMPDIR
    $PARA_PREFIX $ESPRESSO_ROOT/COUPLE/examples/f2pw.x $PARA_POSTFIX \
        -i $TESTDIR/$name.in -nk 2 -ndiag 4 -nres $res > $TESTDIR/$name.f-$res.out
    if test $? != 0; then
      $ECHO "FAILED with error condition!"
      $ECHO "Input: $name.in, Output: $name.f-$res.out, Reference: $name.f-$res.ref"
      $ECHO "Aborting"
      exit 1
    fi
    cd $TESTDIR
    ###
    if test -f $name.f-$res.ref ; then
      # reference file exists
      # Test for scf/relax/md/vc-relax
      #
      check_scf $name.f-$res
      #
      # extract wall time statistics
      #
      get_times $name.f-$res
      #
    else
      $ECHO  "not checked, reference file not available "
    fi
    $ECHO "c2pw.x with $res reserved MPI ranks...\c"
    cd $ESPRESSO_TMPDIR
    $PARA_PREFIX $ESPRESSO_ROOT/COUPLE/examples/c2pw.x $PARA_POSTFIX \
        -i $TESTDIR/$name.in -nk 2 -nd 4 -nres $res > $TESTDIR/$name.c-$res.out
    if test $? != 0; then
      $ECHO "FAILED with error condition!"
      $ECHO "Input: $name.in, Output: $name.c-$res.out, Reference: $name.c-$res.ref"
      $ECHO "Aborting"
      exit 1
    fi
    cd $TESTDIR
    ###
    if test -f $name.c-$res.ref ; then
      # reference file exists
      # Test for scf/relax/md/vc-relax
      #
      check_scf $name.c-$res
      #
      # extract wall time statistics
      #
      get_times $name.c-$res
      #
    else
      $ECHO  "not checked, reference file not available "
    fi
  done
  #
done

########################################################################
# Perform here required checks on cp.x based tests
########################################################################
for file in $cpfiles
do
  name=`basename $file .in`
  get_pp $name
  $ECHO "Checking $name..."
  ###
  # run the code in the scratch directory
  #
  for res in 0 2 4 16
  do
    $ECHO "f2cp.x with $res reserved MPI ranks...\c"
    cd $ESPRESSO_TMPDIR
    $PARA_PREFIX $ESPRESSO_ROOT/COUPLE/examples/f2cp.x $PARA_POSTFIX \
        -i $TESTDIR/$name.in -ndiag 4 -nres $res > $TESTDIR/$name.f-$res.out
    if test $? != 0; then
      $ECHO "FAILED with error condition!"
      $ECHO "Input: $name.in, Output: $name.f-$res.out, Reference: $name.f-$res.ref"
      $ECHO "Aborting"
      exit 1
    fi
    cd $TESTDIR
    ###
    if test -f $name.f-$res.ref ; then
      # reference file exists
      # Test for scf/relax/md/vc-relax
      #
      check_cp $name.f-$res
      #
      # extract wall time statistics
      #
      get_times $name.f-$res
      #
    else
      $ECHO  "not checked, reference file not available "
    fi
    $ECHO "c2cp.x with $res reserved MPI ranks...\c"
    cd $ESPRESSO_TMPDIR
    $PARA_PREFIX $ESPRESSO_ROOT/COUPLE/examples/c2cp.x $PARA_POSTFIX \
        -i $TESTDIR/$name.in -nd 4 -nres $res > $TESTDIR/$name.c-$res.out
    if test $? != 0; then
      $ECHO "FAILED with error condition!"
      $ECHO "Input: $name.in, Output: $name.c-$res.out, Reference: $name.c-$res.ref"
      $ECHO "Aborting"
      exit 1
    fi
    cd $TESTDIR
    ###
    if test -f $name.c-$res.ref ; then
      # reference file exists
      # Test for scf/relax/md/vc-relax
      #
      check_cp $name.c-$res
      #
      # extract wall time statistics
      #
      get_times $name.c-$res
      #
    else
      $ECHO  "not checked, reference file not available "
    fi
  done
  #
done

$ECHO  "Total wall time (s) spent in this run: " $totout
$ECHO  "Reference                            : " $totref

