#!/bin/sh 

# Purpose: install the qe-modes package

prefix=${prefix:-$HOME/.emacs.d/}
prefix=${prefix%/}

fulldir=$prefix/qe-modes


have_install=1
check_install=`type install`
if test $? -gt 0; then have_install=0; fi

if test ! -d $prefix; then
    mkdir -p $prefix;
    if test $? -gt 0; then  exit 1; fi
fi


if test $have_install -eq 1; then

    echo "
installing qe-modes in directory :   $fulldir
"
    install -d $fulldir
    install -m 644 qe-modes/* $fulldir/
    
else
    echo "
copying qe-modes to directory :   $fulldir
"
    if test ! -d $fulldir; then mkdir $fulldir; fi
    cp -r qe-modes/* $fulldir/
fi
