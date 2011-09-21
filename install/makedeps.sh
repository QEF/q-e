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
    dirs=" Modules clib PW CPV/src flib pwtools upftools PP PWCOND/src \
           PHonon/Gamma PHonon/PH PHonon/D3 atomic/src VdW XSpectra/src \
	   GWW/gww GWW/pw4gww GWW/head ACFDT NEB/src Solvent" 
          
else
    dirs=$*
fi

for dir in $dirs; do

    # the following command removes a trailing slash
    DIR=`echo ${dir%/}`

    # the following would also work
    #DIR=`echo $dir | sed "s,/$,,"`

    # set inter-directory dependencies - only directories containing
    # modules that are used, or files that are included, by routines
    # in directory DIR should be listed in DEPENDS
    LEVEL1=..
    LEVEL2=../..
    DEPENDS="$LEVEL1/include $LEVEL1/iotk/src"
    case $DIR in 
        EE | flib | upftools | Solvent | PW )
             DEPENDS="$LEVEL1/include $LEVEL1/iotk/src $LEVEL1/Modules" ;;
	PP | PWCOND | pwtools )
             DEPENDS="$LEVEL1/include $LEVEL1/iotk/src $LEVEL1/Modules \
                      $LEVEL1/PW" ;;
	VdW | ACFDT ) 
             DEPENDS="$LEVEL1/include $LEVEL1/iotk/src $LEVEL1/Modules \
                      $LEVEL1/PW $LEVEL1/PHonon/PH" ;;
	CPV/src | atomic/src | GWW/gww )
             DEPENDS="$LEVEL2/include $LEVEL2/iotk/src $LEVEL2/Modules" ;;
	PHonon/PH | PHonon/Gamma | XSpectra/src  | PWCOND/src | GWW/pw4gww | NEB/src )
             DEPENDS="$LEVEL2/include $LEVEL2/iotk/src $LEVEL2/Modules \
                      $LEVEL2/PW" ;;
	PHonon/D3 )
	     DEPENDS="$LEVEL2/include $LEVEL2/iotk/src $LEVEL2/Modules \
	              $LEVEL2/PW $LEVEL2/PHonon/PH" ;;	
	GWW/head )
             DEPENDS="$LEVEL2/include $LEVEL2/iotk/src $LEVEL2/Modules \
                      $LEVEL2/PW $LEVEL2/PHonon/PH $LEVEL1/pw4gww " ;;
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
