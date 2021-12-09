#!/bin/sh
# compute dependencies for the PWscf directory tree

# make sure there is no locale setting creating unneeded differences.
LC_ALL=C
export LC_ALL

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
TOPDIR=$(cd ../../install/; pwd)

dirs=.

for dir in $dirs; do

    # the following command removes a trailing slash
    DIR=`echo ${dir%/}`

    # the following would also work
    #DIR=`echo $dir | sed "s,/$,,"`

    # set inter-directory dependencies - only directories containing
    # modules that are used, or files that are included, by routines
    # in directory DIR should be listed in DEPENDS
    DEPENDS=". ../../XClib ../../Modules ../../PW/src  ../../PHonon/PH ../../iotk/src ../../FFTXlib ../../LAXlib ../../LR_Modules ../../UtilXlib ../../upflib ../../FoX/finclude/"
    # generate dependencies file (only for directories that are present)
    if test -d $DIR
    then
	cd $DIR
       
	$TOPDIR/moduledep.sh $DEPENDS > make.depend
	$TOPDIR/includedep.sh $DEPENDS >> make.depend

        # handle special cases
        sed '/@\/cineca\/prod\/hpm\/include\/f_hpm.h@/d' \
            make.depend > make.depend.tmp
        sed '/@iso_c_binding@/d;/@ifcore@/d;/@iso_fortran_env@/d' make.depend.tmp > make.depend
        sed -i '/@fox_dom@/d;/@fox_wxml@/d;/@m_common_io@/d;/@fox_sax@/d;/@fox_common@/d'  make.depend

        if test "$DIR" = "Modules"
        then
            sed '/@mpi@/d' make.depend > make.depend.tmp
            sed '/@elpa1@/d' make.depend.tmp > make.depend
        fi

        if test "$DIR" = "clib"
        then
            mv make.depend make.depend.tmp
            sed 's/@fftw.c@/fftw.c/' make.depend.tmp > make.depend
        fi

        if test "$DIR" = "PW/src"
        then
            mv make.depend make.depend.tmp
            sed '/@environ_base@/d' make.depend.tmp > make.depend
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
