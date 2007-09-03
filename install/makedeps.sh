#!/bin/sh
# compute dependencies for the PWscf directory tree

# make sure there is no locale setting creating unneeded differences.
LC_ALL=C
export LC_ALL

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
TOPDIR=`pwd`

if test $# = 0
then
    dirs=" Modules clib PW CPV flib pwtools upftools PP PWCOND \
           Gamma PH D3 atomic GIPAW VdW"
else
    dirs=$*
fi

for DIR_ in $dirs
do
    # set inter-directory dependencies
    DIR=`echo $DIR_ | sed 's?/??' `
    DEPENDS="../include ../flib ../iotk/src"
    case $DIR in 
	PW | CPV | flib | pwtools | upftools | atomic )
		  DEPENDS="$DEPENDS ../Modules"             ;;
	PP | PWCOND | Gamma | PH | GIPAW )
		  DEPENDS="$DEPENDS ../Modules ../PW"       ;;
	D3 | VdW ) 
                  DEPENDS="$DEPENDS ../Modules ../PW ../PH" ;;
    esac

    # generate dependencies file
    if test -d $TOPDIR/$DIR
    then
	cd $TOPDIR/$DIR
	$TOPDIR/moduledep.sh $DEPENDS > make.depend
	$TOPDIR/includedep.sh $DEPENDS >> make.depend
    fi

    # handle special cases
    mv make.depend make.depend.tmp
    sed '/@\/cineca\/prod\/hpm\/include\/f_hpm.h@/d' \
        make.depend.tmp > make.depend

    if test "$DIR" = "clib"
    then
        mv make.depend make.depend.tmp
        sed 's/@fftw.c@/fftw.c/' make.depend.tmp > make.depend
    fi

    if test "$DIR" = "Modules"
    then
	cat make.depend | sed 's/fft_scalar.o : @mkl_dfti.f90@//' > make.depend.tmp 
	cat make.depend.tmp | sed 's/fft_scalar.o : @mkl_dft_type@//' > make.depend
   fi

    rm -f make.depend.tmp

    # check for missing dependencies
    if grep @ make.depend
    then
	notfound=1
	echo WARNING: dependencies not found in directory $DIR
    else
        echo directory $DIR : ok
    fi
done
if test "$notfound" = ""
then
    echo all dependencies updated successfully
fi
