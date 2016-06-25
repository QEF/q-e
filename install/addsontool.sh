#!/bin/bash

# everything is performed in the destination directory
# this script has to be run by addsonpatch.sh
# do not execut it manually

cd "$destination"

if [ "$#" -eq 0 ];
then
 echo "[ USAGE :"
 echo "./install/addsonpatch.sh  ADDSON_NAME WHERE_SOURCE WHERE_LINKS (-patch) (-revert)   "
 echo " addsonpatch.sh has to be run from the Quantum ESPRESSO root directory"
 echo "WHERE_SOURCE is the relative path to the sources of the Addson code "
 echo "WHERE_LINKS is the relative path to the QE directory where the addson sources have to be linked"
 echo "at the moment it only allows for pure f90 routines to be linked in flib"
 echo "or pure f90 modules to be linked in Modules"
 echo " -patch  : apply patch to Makefiles " 
 echo " -revert : revert Makefiles to original "
 echo " ]"
 exit
fi

case "$4" in
(-patch)

  echo "* I will try to patch needed files for integrated compilation ..."

  if test -e "${ADDSON_NAME}_PATCH" ; then
    echo "-- File $destination/${ADDSON_NAME}_PATCH exists"
    echo "-- I guess you have already patched $ADDSON_NAME"
    echo "-- Please unpatch it first, or start from a clean source tree"
    echo "-- See you later..."
    echo "* ABORT"
    exit
  fi
  echo "#Please do not remove or modify this file"                    >  ${ADDSON_NAME}_PATCH
  echo "#It is keeps track of the steps for patching $ADDSON package" >> ${ADDSON_NAME}_PATCH
  
#-------------------
  echo "-- Executing pre script"

  command -v patch &>/dev/null || { echo "I require patch command but it's not installed. Aborting." >&2; exit 1;  }

#------------------- check if GNU patch works
  cat > test_patch1 << \EOF
  alfa
  beta
EOF

  cat > test_patch2 << \EOF
  alfa
  gamma
EOF

  cat > test_patch3 << \EOF_EOF
  patch -c -l -b -F 3 --suffix=.pre "./test_patch1" << \EOF
EOF_EOF

  diff -c test_patch1 test_patch2 >> test_patch3

  echo EOF >> test_patch3

  bash test_patch3 &> test_patch4

  status=$?
  if [ $status -ne 0 ]
  then
    echo "patch does not work! Error message:"
    echo "**********"
    cat test_patch4
    echo "**********"
    echo "Please install a recent version of the GNU patch utility and try again."
    exit
  fi

  rm test_patch1 test_patch2 test_patch3 test_patch4
  if [ -e test_patch1.pre ]
  then
    rm test_patch1.pre
  fi
#-------------------------------------------

  command -v sed &>/dev/null || { echo "I require sed command but it's not installed. Aborting." >&2; exit 1;  }

#------------------- check if GNU sed works
  cat > test_sed1 << \EOF
  alfa
  beta
EOF

  cat > test_sed2 << \EOF
  alfa
  gamma
  beta
EOF

  sed '/alfa/ a\
  gamma' test_sed1 > tmp.1

  mv tmp.1 test_sed1

  diff -c test_sed1 test_sed2 >> test_sed3

#  echo EOF >> test_sed3

  bash test_sed3 &> test_sed4

  status=$?
  if [ $status -ne 0 ]
  then
    echo "sed does not work! Error message:"
    echo "**********"
    cat test_sed4
    echo "**********"
    echo "Please install a recent version of the GNU sed utility and try again."
    exit
  fi

  rm test_sed1 test_sed2 test_sed3 test_sed4
# -----------------------------------------
# -----------------------------------------
  to_do_before_patch

  echo "-- Setting up symlinks"
  for file in $destination/$LINKED_FILES ; do
    base="${file##*/}"
    if test -e $destination/$WHERE_LINKS/$base ; then
      echo "PATCH ERROR: file $base is already in $WHERE_LINKS"
      exit 1
    fi
#    echo "$destination/$WHERE_LINKS/$base"
    ln -s $file $destination/$WHERE_LINKS/$base
  done
  
  tmp_var=\$\(${ADDSON_NAME}_OBJECTS\)
  
  echo "-- modifying $WHERE_LINKS/Makefile"
  sed < $destination/$WHERE_LINKS/Makefile.pre$ADDSON_NAME > $destination/$WHERE_LINKS/tmp.1 '/make.inc/ a\
  include '"${ADDSON_NAME}"'.inc \
   '
  sed < $destination/$WHERE_LINKS/tmp.1 > $destination/$WHERE_LINKS/Makefile '/= \\/ a\
  '"${tmp_var}"' \\'
  
  rm $destination/$WHERE_LINKS/tmp.1

  echo "-- Executing post script"
  to_do_after_patch

  echo "- DONE!"
;;
(-revert)
  echo "* I will try to revert ..."
  echo "-- Executing pre script"

  to_do_before_revert


  echo "-- Removing symlinks"
  for file in $destination/$LINKED_FILES ; do
    base="${file##*/}"
    if test -e $destination/$WHERE_LINKS/$base ; then \
#      echo "$destination/$WHERE_LINKS/$base" ; \
      rm $destination/$WHERE_LINKS/$base ; \
    else
      echo "where_links base: $destination/$WHERE_LINKS/$base"
      echo "PATCH WARNING: file $base is not in $destination/$WHERE_LINKS"
    fi
  done



  echo "-- Restoring .pre$ADDSON_NAME files"
  PREADDSON=$(find . -name "*.pre*")
  if ! test "$PREADDSON" ; then
    echo "-- I cannot find any .pre$ADDSON_NAME file"
    echo "* ABORT"
    exit
  fi
  
  rm ${ADDSON_NAME}_PATCH

  echo "-- Executing post script"
  to_do_after_revert

  echo "* DONE!"
;;
(*)
  echo "Missing input argument"
esac
