#!/bin/sh

# Automated checks for neb.x - PG 2011
# Same logic of "check-pw.x.j"
# You shouldn't need to modify anything below this line.

# taken from examples - not sure it is really needed
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

ESPRESSO_ROOT=`cd ../.. ; pwd`
. $ESPRESSO_ROOT/environment_variables
ESPRESSO_TMPDIR=$ESPRESSO_ROOT/tempdir/
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
    files=`/bin/ls *.in`
else
    files=`/bin/ls $*| grep "\.in$"`
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
	    $WGET  $ESPRESSO_PSEUDO/$ppfile \
                http://www.quantum-espresso.org/pseudo/1.3/UPF/$ppfile \
		2> /dev/null
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
# Perform here required checks
########################################################################
for file in $files
do
  name=`basename $file .in`
  get_pp $name
  $ECHO "Checking $name...\c"
  ###
  # run the code in the scratch directory
  #
  cd $ESPRESSO_TMPDIR
  $PARA_PREFIX $ESPRESSO_ROOT/bin/neb.x $PARA_POSTFIX \
        -inp $TESTDIR/$name.in > $TESTDIR/$name.out 2> /dev/null
  if test $? != 0; then
     $ECHO "FAILED with error condition!"
     $ECHO "Input: $name.in, Output: $name.out, Reference: $name.ref"
     $ECHO "Aborting"
     exit 1
  fi
  #
  cd $TESTDIR
  ###
  if test -f $name.ref ; then
     # reference file exists
     if grep 'neb: convergence achieved' $name.ref > /dev/null; then
        # Specific test for NEB
	check_neb $name
        #
     fi
     #
     # extract wall time statistics
     #
     get_times $name
     #
  else
     $ECHO  "not checked, reference file not available "
  fi
  #
done

$ECHO  "Total wall time (s) spent in this run: " $totout
$ECHO  "Reference                            : " $totref

