#!/bin/sh 

# Purpose: install the QE-modes package

prefix=${prefix:-$HOME/.emacs.d/}
prefix=${prefix%/}

fulldir=$prefix/qe-modes


have_install=1
check_install=`type install`
if test $? -gt 0; then have_install=0; fi

# check if $fulldir exists 

if test -d $fulldir; then
    echo -n "
Directory $fulldir to which QE-modes will be installed exists.
Do you want to overwrite (y|n): "
    read ans    
    if test "x$ans" != xy; then
        echo "Aborting."
        exit
    else
        echo "Removing the old installation *.elc files ..."
        rm -f $fulldir/*.elc
    fi
fi

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

# byte-compile the *.el files in $fulldir

echo -n "Byte compiling QE-mode Elisp files in $fulldir ... "
(
    cd $fulldir
    emacs -Q --batch --eval "(add-to-list 'load-path \".\")" -f batch-byte-compile *.el
    if test $? -eq 0; then
        echo [OK]
    fi
)
