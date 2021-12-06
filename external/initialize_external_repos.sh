#!/bin/bash

# this script use git to initialize all the external package repositories.

commit_hash_file=./submodule_commit_hash_records

if [ ! -f ${commit_hash_file} ]; then
  echo "Failed to find $commit_hash_file, are you in the QE TOPDIR/external subdirectory"
  exit
fi

if ! which git; then
  echo "Failed to find git!"
  exit
fi

function init_component()
{
  SUBMODULE_NAME=$1
  echo "initializing $SUBMODULE_NAME ..."

  if [ -d $SUBMODULE_NAME/.git ]; then
    echo "  $SUBMODULE_NAME subdirectory exists. Skipping."
    echo "  If you would like to force a reinitializaton, remove $SUBMODULE_NAME subdirectory."
    return
  fi

  RECORD_HASH=`grep $SUBMODULE_NAME $commit_hash_file | awk '{print $1}'`
  echo "    hash $RECORD_HASH"

  SUBMODULE_URL=`git config --file ../.gitmodules --get submodule.external/$SUBMODULE_NAME.URL`
  echo "  Cloning ${SUBMODULE_URL} into ${SUBMODULE_NAME}."

  git init ${SUBMODULE_NAME}
  cd ${SUBMODULE_NAME}
    git remote add origin ${SUBMODULE_URL}
    git fetch --depth 1 origin ${RECORD_HASH}
    git checkout -b recorded_HEAD FETCH_HEAD
  cd ..
}

for component in `awk '{print $2}' $commit_hash_file`
do
  init_component $component
done
