#!/bin/sh
# compute dependencies

TOPDIR=`pwd`

# directory that only depends on itself
cd $TOPDIR/Modules
$TOPDIR/moduledep.sh > .dependencies

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
DIRS="PP PWCOND PWNC Gamma PH"
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
