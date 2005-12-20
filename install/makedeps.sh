#!/bin/sh
# compute dependencies for the PWscf directory tree

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
TOPDIR=`pwd`

for DIR in Modules clib PW CPV flib pwtools upftools PP PWCOND \
           Gamma PH D3 atomic Nmr CPVIB
do
    # set inter-directory dependencies
    case $DIR in
	Modules | clib )
		  DEPENDS="../include ../iotk/src"                        ;;
	PW | CPV | flib | pwtools | upftools | atomic )
		  DEPENDS="../include ../Modules ../iotk/src"             ;;
	PP | PWCOND | Gamma | PH )
		  DEPENDS="../include ../Modules ../PW ../iotk/src"       ;;
	D3 | Nmr) DEPENDS="../include ../Modules ../PW ../PH ../iotk/src" ;;
        CPVIB )   DEPENDS="../include ../Modules ../PW ../iotk/src ../CPV"
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
