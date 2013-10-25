#!/bin/bash

# everything is performed in the destination directory
cd "$destination"

if [ "$#" -eq 0 ];
then
 echo "USAGE :"
 echo "$NAME  WHERE_SOURCE WHERE_LINKS (-patch) (-revert)   "
 echo "WHERE_SOURCE is the relative path to the sources of the Addson code "
 echo "WHERE_LINKS is the relative path to the QE directory where the addson sources have to be linked"
 echo " -patch  : apply patch to Makefiles " 
 echo " -revert : revert Makefiles to original "
 echo " )"
 exit
fi

case "$3" in
(-patch)
  echo "* I will try to patch needed files for integrated compilation ..."
  
#-------------------

  to_do_before_patch

  echo "-- Setting up symlinks"
  for file in $destination/$LINKED_FILES ; do
    base="${file##*/}"
    if test -e $destination/$WHERE_LINKS/$base ; then
      echo "PATCH ERROR: file $base is already in $WHERE_LINKS"
      exit 1
    fi
    ln -s $file $destination/$WHERE_LINKS/$base
  done

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
    if test -e $destination/$WHERE_LINKS/$base ; then
      rm $destination/$WHERE_LINKS/$base 
    else
      echo "where_links base: destination$WHERE_LINKS/$base"
      echo "PATCH WARNING: file $base is not in $destination/$WHERE_LINKS"
    fi
  done



  echo "-- Restoring .original files"
  PREADDSON=$(find . -name "*.original")
  if ! test "$PREADDSON" ; then
    echo "-- I cannot find any .original file"
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
