#!/bin/bash
# PATCH SCRIPT FOR QE Makefiles with addson
#

# this script needs to be launched from the root directory of the host code
# two arguments are neede where the addson source code is and
# where the source code has to be linked in order to be
# compiled by QE

# This script has been adapted from an original patch script
# of plumed (www.plumed-code.org)

destination="$PWD"
echo "root directory of host package: $destination"

LINKED_FILES="$1/*.f90"
WHERE_LINKS="$2/"

echo "LINKED_FILES are: $LINKED_FILES"
echo "WHERE_LINKS are: $WHERE_LINKS"

function to_do_before_patch () {
  echo > /dev/null
}

function to_do_after_patch () {
  {
    echo -n "ADDSON_OBJECTS="
    for file in $destination/$LINKED_FILES
      do f=${file##*/}
      echo " \\"
      echo -n "	${f%.f90}.o"
    done
    echo
    echo -n "ADDSON_SRC="
    for file in $destination/$LINKED_FILES
      do f=${file##*/}
      echo " \\"
      echo -n "	${f%.f90}.f90"
    done
    echo
    echo
  } > $destination/$WHERE_LINKS/ADDSON.inc
  cp $destination/make.sys $destination/make.sys.original
}

function to_do_before_revert () {
  rm $destination/$WHERE_LINKS/addson.inc
}

function to_do_after_revert () {
  echo > /dev/null
  mv $destination/make.sys.original $destination/make.sys
}

#########

NAME="$0"
echo "NAME $NAME "
if test -e $destination/install/addsontool.sh ; then
  source $destination/install/addsontool.sh
else
  echo "missing file addsontool.sh in install directory"
  EXIT
fi
