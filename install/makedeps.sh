#!/bin/sh
# compute dependencies for the PWscf directory tree

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
TOPDIR=`pwd`

for DIR in Modules PW CPV flib pwtools upftools PP PWCOND PWNC \
           Gamma PH D3 Raman atomic Nmr
do
    # set inter-directory dependencies
    case $DIR in
	Modules )                         DEPENDS=""                       ;;
	PW | CPV | flib | pwtools | upftools | atomic )
	                                  DEPENDS="../Modules"             ;;
	PP | PWCOND | PWNC | Gamma | PH ) DEPENDS="../Modules ../PW"       ;;
	D3 | Raman | Nmr)                 DEPENDS="../Modules ../PW ../PH" ;;
    esac

    # generate dependencies file
    if test -d $TOPDIR/$DIR
    then
	cd $TOPDIR/$DIR
	$TOPDIR/moduledep.sh $DEPENDS > make.depend
    fi

    # check for missing dependencies
    if grep @ make.depend
    then
	echo WARNING: modules not found in directory $DIR
    fi
done
