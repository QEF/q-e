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
# this is the list of all directories for which we want to find dependencies
# upon include files *.h or *.fh or modules. Note that libraries that are 
# externally maintained should not go into this list

    dirs=" LAXlib FFTXlib UtilXlib clib \
           KS_Solvers/Davidson KS_Solvers/Davidson_RCI KS_Solvers/CG \
	   KS_Solvers/PPCG KS_Solvers/ParO  KS_Solvers/DENSE  \
           upflib Modules LR_Modules PW/src CPV/src PW/tools PP/src PWCOND/src \
           PHonon/Gamma PHonon/PH PHonon/FD HP/src atomic/src \
           EPW/src XSpectra/src ACFDT/src NEB/src TDDFPT/src \
           GWW/pw4gww GWW/gww GWW/head GWW/bse GWW/simple \
	   GWW/simple_bse GWW/simple_ip" 
          
elif
    test $1 = "-addson" 
then
    echo "The script for adding new dependencies is running"
    echo "Usage: $0 -addson DIR DEPENDENCY_DIRS"
    echo "$0 assumes that the new dependencies are in $TOPDIR/../"
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
    # (directory DIR itself should not be listed in DEPENDS)
    LEVEL1=..
    LEVEL2=../..
    # default
    DEPENDS="$LEVEL1/include" 
    # for convenience, used later
    DEPEND1="$LEVEL1/include $LEVEL1/FFTXlib $LEVEL1/LAXlib $LEVEL1/UtilXlib \
	     $LEVEL1/upflib"
    DEPEND3="$LEVEL2/include $LEVEL2/FFTXlib $LEVEL2/LAXlib $LEVEL2/UtilXlib"
    DEPEND2="$DEPEND3 $LEVEL2/upflib $LEVEL2/Modules"
    case $DIR in 
        Modules )
             DEPENDS="$DEPEND1" ;;
        LR_Modules )
             DEPENDS="$DEPEND1 $LEVEL1/Modules $LEVEL1/PW/src" ;;
	ACFDT/src ) 
             DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL2/PHonon/PH $LEVEL2/LR_Modules" ;;
	atomic/src | GWW/gww )
	     DEPENDS="$DEPEND2" ;;
	PW/src | CPV/src )
	     DEPENDS="$DEPEND2 ../../KS_Solvers/Davidson ../../KS_Solvers/CG ../../KS_Solvers/PPCG ../../KS_Solvers/ParO ../../KS_Solvers/DENSE ../../dft-d3" ;;
	KS_Solvers/Davidson | KS_Solvers/Davidson_RCI | KS_Solvers/CG | KS_Solvers/PPCG | KS_Solvers/ParO | KS_Solvers/DENSE )
	     DEPENDS="$DEPEND3" ;;
	PW/tools | PP/src | PWCOND/src | GWW/pw4gww | NEB/src )
	     DEPENDS="$DEPEND2 $LEVEL2/PW/src" ;;
	PHonon/FD | PHonon/PH | PHonon/Gamma | HP/src | TDDFPT/src | XSpectra/src  | GIPAW/src )
	     DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL2/LR_Modules" ;;
        EPW/src )
             DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL2/LR_Modules $LEVEL2/PHonon/PH $LEVEL2/Modules" ;; 
	GWW/head )
	     DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL2/PHonon/PH $LEVEL2/LR_Modules" ;;	
	GWW/bse )
	     DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL2/PHonon/PH $LEVEL2/LR_Modules $LEVEL2/GWW/pw4gww $LEVEL2/GWW/gww" ;;
	GWW/simple )
	     DEPENDS="$DEPEND2 $LEVEL2/PW/src $LEVEL2/GWW/pw4gww $LEVEL2/GWW/gww" ;;
	GWW/simple_bse )
	     DEPENDS="$DEPEND2 $LEVEL2/GWW/gww" ;;
	GWW/simple_ip)
	     DEPENDS="$DEPEND2" ;;
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

        # handle special cases: modules for C-fortran binding,
        #   	                hdf5, MPI, FoX, libxc
        sed '/@iso_c_binding@/d' make.depend > tmp; mv tmp make.depend
        sed '/@hdf5@/d' make.depend > tmp; mv tmp make.depend
        sed '/@mpi@/d'  make.depend > tmp; mv tmp make.depend
        sed '/@fox_dom@/d;/@fox_wxml@/d;/@m_common_io@/d' make.depend > tmp; mv tmp make.depend
        sed '/@xc_version.h@/d;/@xc_f03_lib_m@/d' make.depend > tmp; mv tmp make.depend

        if test "$DIR" = "FFTXlib"
        then
            # more special cases: modules for FFTs, GPU, OpenMP
            sed '/@omp_lib@/d' make.depend > tmp; mv tmp make.depend
            sed '/@mkl_dfti/d' make.depend > tmp; mv tmp make.depend
            sed '/@fftw3.f/d;s/@fftw.c@/fftw.c/' make.depend > tmp; mv tmp make.depend
            sed '/@cudafor@/d;/@cufft@/d;/@flops_tracker@/d' make.depend > tmp; mv tmp make.depend
        fi

        if test "$DIR" = "LAXlib"
        then
            # more special cases: modules for ELPA, GPUs
            sed '/@elpa1@/d;/@elpa@/d' make.depend > tmp; mv tmp make.depend
            sed '/@cudafor@/d;/@cusolverdn@/d;/@gbuffers@/d' make.depend > tmp; mv tmp make.depend
            sed '/@zhegvdx_gpu@/d;/@dsyevd_gpu@/d;/@dsygvdx_gpu@/d' make.depend > tmp; mv tmp make.depend
            sed '/@cublas@/d;/@eigsolve_vars@/d;/@nvtx_inters@/d' make.depend > tmp ; mv tmp make.depend
            sed '/@device_fbuff_m@/d' make.depend > tmp ; mv tmp make.depend
        fi

        if test "$DIR" = "UtilXlib"
        then
            sed '/@ifcore@/d' make.depend > tmp; mv tmp make.depend
            sed '/@cudafor@/d' make.depend> tmp; mv tmp make.depend
        fi

        if test "$DIR" = "Modules"
        then
            sed '/@mbd@/d' make.depend > tmp; mv tmp make.depend
        fi

        if test "$DIR" = "PW/src" || test "$DIR" = "TDDFPT/src"
        then
            sed '/@environ_/d;/@solvent_tddfpt@/d' make.depend > tmp; mv tmp make.depend
        fi

        if test "$DIR" = "CPV/src"
        then
            sed '/@f90_unix_proc@/d' make.depend > tmp; mv tmp make.depend
        fi

        if test "$DIR" = "EPW/src"
        then
            sed '/@f90_unix_io@/d' make.depend > tmp; mv tmp make.depend
            sed '/@f90_unix_env@/d' make.depend> tmp; mv tmp make.depend
            sed '/@w90_io@/d' make.depend > tmp; mv tmp make.depend
            sed '/@ifport@/d' make.depend > tmp; mv tmp make.depend
        fi

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
