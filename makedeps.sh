#!/bin/sh
# compute dependencies

# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
TOPDIR=`pwd`

# directory that only depends on itself
DIR="Modules"
DEPENDS=""
if test -d $TOPDIR/$DIR ; then
  cd $TOPDIR/$DIR
  $TOPDIR/moduledep.sh $DEPENDS > .dependencies
fi

# directories that depend on Modules/
DIRS="PW FPMD CPV flib"
DEPENDS="../Modules"
for DIR in $DIRS ; do
  if test -d $TOPDIR/$DIR ; then
    cd $TOPDIR/$DIR
    $TOPDIR/moduledep.sh $DEPENDS > .dependencies
  fi
done

# directories that depend on Modules/, PW/
DIRS="PP PWCOND PWNC Gamma PH NEB"
DEPENDS="../Modules ../PW"
for DIR in $DIRS ; do
  if test -d $TOPDIR/$DIR ; then
    cd $TOPDIR/$DIR
    $TOPDIR/moduledep.sh $DEPENDS > .dependencies
  fi
done

# directories that depend on Modules/, PW/, PH/
DIRS="D3"
DEPENDS="../Modules ../PW ../PH"
for DIR in $DIRS ; do
  if test -d $TOPDIR/$DIR ; then
    cd $TOPDIR/$DIR
    $TOPDIR/moduledep.sh $DEPENDS > .dependencies
  fi
done
