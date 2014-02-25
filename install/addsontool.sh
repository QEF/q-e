#!/bin/bash

# everything is performed in the destination directory
cd "$destination"

if [ "$#" -eq 0 ];
then
 echo "USAGE :"
 echo "$NAME  ADDSON_NAME WHERE_SOURCE WHERE_LINKS (-patch) (-revert)   "
 echo "WHERE_SOURCE is the relative path to the sources of the Addson code "
 echo "WHERE_LINKS is the relative path to the QE directory where the addson sources have to be linked"
 echo " -patch  : apply patch to Makefiles " 
 echo " -revert : revert Makefiles to original "
 echo " )"
 exit
fi

case "$4" in
(-patch)
  echo "* I will try to patch needed files for integrated compilation ..."
  
#-------------------
  echo "-- Executing pre script"
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
  sed < $destination/$WHERE_LINKS/Makefile.pre$ADDSON_NAME > $destination/$WHERE_LINKS/tmp.1 '/make.sys/ a\
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


  echo "-- Executing post script"
  to_do_after_revert

  echo "* DONE!"
;;
(*)
  echo "Missing input argument"
esac
