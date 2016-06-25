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
#echo "root directory of host package: $destination"

# bisogna prendere il nome del plugin in input
ADDSON_NAME="$1"

LINKED_FILES="$2/*.f90"
WHERE_LINKS="$3/"

echo "The NAME of the addson is: $ADDSON_NAME"
echo "LINKED_FILES are: $LINKED_FILES"
echo "WHERE_LINKS are: $WHERE_LINKS"

function to_do_before_patch () {
  echo > /dev/null
  cp $destination/make.inc $destination/make.inc.pre$ADDSON_NAME
  cp $destination/$WHERE_LINKS/Makefile $destination/$WHERE_LINKS/Makefile.pre$ADDSON_NAME
  if test -e $destination/$WHERE_LINKS/make.depend ; then 
    cp $destination/$WHERE_LINKS/make.depend $destination/$WHERE_LINKS/make.depend.pre$ADDSON_NAME
  fi
}

function to_do_after_patch () {
  {
    echo -n "${ADDSON_NAME}_OBJECTS=" 
    for file in $destination/$LINKED_FILES
      do f=${file##*/}
      echo " \\"
      echo -n "	${f%.f90}.o"
    done
    echo
    echo -n "${ADDSON_NAME}_SRC="
    for file in $destination/$LINKED_FILES
      do f=${file##*/}
      echo " \\"
      echo -n "	${f%.f90}.f90"
    done
    echo
    echo
  } >> $destination/$WHERE_LINKS/$ADDSON_NAME.inc
}

function to_do_before_revert () {
  rm $destination/$WHERE_LINKS/$ADDSON_NAME.inc
 echo > /dev/null
}

function to_do_after_revert () {
  echo > /dev/null
  mv $destination/make.inc.pre$ADDSON_NAME $destination/make.inc
  mv $destination/$WHERE_LINKS/Makefile.pre$ADDSON_NAME $destination/$WHERE_LINKS/Makefile
  if test -e $destination/$WHERE_LINKS/make.depend.pre$ADDSON_NAME ; then \
  mv $destination/$WHERE_LINKS/make.depend.pre$ADDSON_NAME $destination/$WHERE_LINKS/make.depend ; fi
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
