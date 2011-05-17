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
           Gamma PH D3 atomic VdW XSpectra \
	   GWW/gww GWW/pw4gww GWW/head ACFDT NEB Solvent" 
          
else
    dirs=$*
fi

for DIR_ in $dirs
do
    DIR=`echo $DIR_ | sed 's?/??' `
    # set inter-directory dependencies - only directories containing
    # modules that are used, or files that are included, by routines
    # in directory DIR should be listed in DEPENDS
    DEPENDS="../include ../iotk/src"
    case $DIR in 
        EE | flib | upftools | atomic | CPV | Solvent | PW )
                  DEPENDS="$DEPENDS ../Modules "      ;;
	PP | PWCOND | Gamma | PH | pwtools )
		  DEPENDS="$DEPENDS ../Modules ../PW" ;;
	D3 | VdW | ACFDT ) 
                  DEPENDS="$DEPENDS ../Modules ../PW ../PH" ;;
	XSpectra  )
		  DEPENDS="$DEPENDS ../Modules ../PW"  ;;
        GWW/pw4gww )
                  DEPENDS="../../include ../../iotk/src ../../Modules \
		  ../../PW " ;;
	GWW/gww )
                  DEPENDS="../../include ../../iotk/src ../../Modules " ;;
	GWW/head )
                  DEPENDS="../../include ../../iotk/src ../../Modules \
		  ../../PW ../../PH ../pw4gww " ;;
	NEB )
		  DEPENDS="$DEPENDS ../Modules ../PW" ;;

    esac

    # generate dependencies file (only for directories that are present)
    if test -d $TOPDIR/../$DIR
    then
	cd $TOPDIR/../$DIR
       
	$TOPDIR/moduledep.sh $DEPENDS > make.depend
	$TOPDIR/includedep.sh $DEPENDS >> make.depend

        # handle special cases
        sed '/@\/cineca\/prod\/hpm\/include\/f_hpm.h@/d' \
            make.depend > make.depend.tmp
        sed '/@iso_c_binding@/d' make.depend.tmp > make.depend

        if test "$DIR" = "clib"
        then
            mv make.depend make.depend.tmp
            sed 's/@fftw.c@/fftw.c/' make.depend.tmp > make.depend
        fi

        if test "$DIR" = "PW"
        then
            mv make.depend make.depend.tmp
            sed '/@solvent_base@/d' make.depend.tmp > make.depend
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
    fi
done
if test "$notfound" = ""
then
    echo all dependencies updated successfully
fi
