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
    dirs=" Modules clib PW/src CPV/src flib PW/tools upftools PP/src PWCOND/src\
           PHonon/Gamma PHonon/PH PHonon/D3 PHonon/FD atomic/src XSpectra/src \
           ACDFT NEB/src TDDFPT/src GIPAW/src GWW/pw4gww GWW/gww GWW/head" 
          
elif
    test $1 = "-addson" 
then
    echo "The script for adding new dependencies is running"
    echo "Usage: $0 -addson DIR DEPENDENCY_DIRS"
    echo "$0 assumes that the new dependencies are in $TOPDIR/../"
#    ninput=$#
#    echo "number of input arguments: $ninput"
    dirs=$2
    shift
    shift
    add_deps=$*
    echo "dependencies in $add_deps will be searched for $dirs"
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
        flib | upftools )
             DEPENDS="$LEVEL1/include $LEVEL1/iotk/src $LEVEL1/Modules" ;;
	PP/src  )
             DEPENDS="$LEVEL2/include $LEVEL2/iotk/src $LEVEL2/Modules \
                      $LEVEL2/PW/src" ;;
	ACFDT ) 
             DEPENDS="$LEVEL1/include $LEVEL1/iotk/src $LEVEL1/Modules \
                      $LEVEL1/PW/src $LEVEL1/PHonon/PH" ;;
	PW/src )
	     DEPENDS="$LEVEL2/include $LEVEL2/iotk/src $LEVEL2/Modules" ;;
	PW/tools | PWCOND/src | PHonon/FD )
	     DEPENDS="$LEVEL2/include $LEVEL2/PW/src $LEVEL2/iotk/src $LEVEL2/Modules" ;;
	CPV/src | atomic/src | GWW/gww )
             DEPENDS="$LEVEL2/include $LEVEL2/iotk/src $LEVEL2/Modules" ;;
	PHonon/PH | PHonon/Gamma | XSpectra/src  | PWCOND/src | GWW/pw4gww | NEB/src | GIPAW/src )
             DEPENDS="$LEVEL2/include $LEVEL2/iotk/src $LEVEL2/Modules \
                      $LEVEL2/PW/src" ;;
	PHonon/D3 )
	     DEPENDS="$LEVEL2/include $LEVEL2/iotk/src $LEVEL2/Modules \
	              $LEVEL2/PW/src $LEVEL2/PHonon/PH" ;;	
        GWW/pw4gww )
            DEPENDS="$LEVEL2/include $LEVEL2/iotk/src $LEVEL2/Modules \
                       $LEVEL2/PW/src  " ;;
	GWW/gww )
            DEPENDS="$LEVEL2/include $LEVEL2/iotk/src $LEVEL2/Modules " ;;
        GWW/head )
             DEPENDS="$LEVEL2/include $LEVEL2/iotk/src $LEVEL2/Modules \
                      $LEVEL2/PW/src $LEVEL2/PHonon/PH " ;;
	TDDFPT/src )
             DEPENDS="$LEVEL2/include $LEVEL2/iotk/src $LEVEL2/Modules \
                      $LEVEL2/PW/src $LEVEL2/PHonon/PH" ;;
    *)
# if addson needs a make.depend file
	DEPENDS="$DEPENDS $add_deps"

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
        sed '/@iso_c_binding@/d;/@ifcore@/d' make.depend.tmp > make.depend

        if test "$DIR" = "Modules"
        then
            sed '/@mpi@/d;/@elpa1@/d' make.depend > make.depend.tmp
            sed '/@mkl_dfti/d;/@fftw3.f03@/d' make.depend.tmp > make.depend
        fi

        if test "$DIR" = "clib"
        then
            mv make.depend make.depend.tmp
            sed 's/@fftw.c@/fftw.c/' make.depend.tmp > make.depend
        fi

        if test "$DIR" = "PW/src" || test "$DIR" = "TDDFPT/src"
        then
            sed '/@environ_/d;/@solvent_tddfpt@/d' make.depend > make.depend.tmp
            sed '/fft_defs.h@/d' make.depend.tmp > make.depend
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
    else
       echo directory $DIR : not present in $TOPDIR 
    fi
done
if test "$notfound" = ""
then
    echo all dependencies updated successfully
fi
