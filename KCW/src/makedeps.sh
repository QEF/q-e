#!/bin/sh
# compute dependencies for the PWscf directory tree

# make sure there is no locale setting creating unneeded differences.
LC_ALL=C
export LC_ALL

# run from directory where this script is
#cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname

dirs="./"


# the following command removes a trailing slash
DIR=`echo ${dir%/}`

# the following would also work
#DIR=`echo $dir | sed "s,/$,,"`

# set inter-directory dependencies - only directories containing
# modules that are used, or files that are included, by routines
# in directory DIR should be listed in DEPENDS

LEVEL1=..
LEVEL2=../..
# default
DEPENDS="$LEVEL1/include" 
# for convenience, used later
DEPEND1="$LEVEL1/include $LEVEL1/iotk/src $LEVEL1/FFTXlib/src $LEVEL1/LAXlib $LEVEL1/UtilXlib"
DEPEND2="$LEVEL2/include $LEVEL2/iotk/src $LEVEL2/FFTXlib/src $LEVEL2/LAXlib $LEVEL2/UtilXlib $LEVEL2/Modules $LEVEL2/upflib $LEVEL2/XClib"

DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL2/LR_Modules"

# generate dependencies file (only for directories that are present)
   
../../install/moduledep.sh $DEPENDS > make.depend
../../install/includedep.sh $DEPENDS >> make.depend

# handle special cases: hardware-specific monitoring tools
sed '/@\/cineca\/prod\/hpm\/include\/f_hpm.h@/d;/@ifcore@/d' make.depend > make.depend.tmp
# handle special cases: modules for C-fortran binding, system utilities
sed '/@iso_c_binding@/d;/@f90_unix_env@/d;s/fft_scalar.*.o/fft_scalar.o/' make.depend.tmp > make.depend
rm -f make.depend.tmp

# check for missing dependencies 
if grep @ make.depend
then
   notfound=1
   echo WARNING: dependencies not found in directory $DIR
else
   echo directory $DIR : ok
fi

if test "$notfound" = ""
then
    echo all dependencies updated successfully
fi
